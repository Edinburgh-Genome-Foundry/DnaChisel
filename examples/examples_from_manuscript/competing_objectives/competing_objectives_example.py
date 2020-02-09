from copy import deepcopy
from collections import OrderedDict
import dnachisel as dc
import pandas
import seaborn as sns

sns.set()

# DEFINE HOW OPTIMIZATION PROBLEMS ARE CREATED

specifications = {
    "~keep": dc.AvoidChanges(),
    "~no(CG)": dc.AvoidPattern("CG"),
    "~codon_optimize": dc.CodonOptimize(species="e_coli"),
    "~unique_kmers": dc.UniquifyAllKmers(20),
    "~gc(39%)": dc.EnforceGCContent(target=0.39, window=200),
}
class_to_label = {
    spec.__class__: label for label, spec in specifications.items()
}
sequence = dc.load_record("record.gb")


def create_problem(boost_profile):
    location = dc.Location(1000, 9247)
    objectives = []
    for spec_name, boost in boost_profile.items():
        spec = specifications[spec_name]
        spec = spec.copy_with_changes(boost=boost, location=location)
        objectives.append(spec)
    return dc.DnaOptimizationProblem(
        sequence,
        constraints=[dc.EnforceTranslation(location=location)],
        objectives=objectives,
    )


def get_scores_from_problem(problem):
    """Extract the score of each objective post-optimization."""
    return {
        class_to_label[e.specification.__class__]: e.score
        for e in problem.objectives_evaluations().evaluations
    }


# DEFINE THE DIFFERENT OPTIMIZATION PROFILES USED IN THIS EXPERIMENT

profiles = OrderedDict()
for spec in specifications:
    profiles[spec + " alone"] = {
        s: (1 if s == spec else 0) for s in specifications
    }
profiles["equal weights"] = {spec: 1 for spec in specifications}
for spec in specifications:
    profiles[spec + " 5x"] = {
        s: (5 if s == spec else 1) for s in specifications
    }

# PERFORM OPTIMIZATIONS FOR EACH PROFILE

results = OrderedDict()

for name, profile in profiles.items():
    print(name)
    problem = create_problem(profile)
    problem.resolve_constraints()
    problem.optimize()
    results[name] = get_scores_from_problem(problem)
problem = create_problem(profiles["equal weights"])
results["initial sequence"] = get_scores_from_problem(problem)

# PLOT THE RESULTS

df = -pandas.DataFrame.from_dict(results).T[::-1]
df.to_excel("optimization_data.xlsx")
ax = df.plot.barh(legend=False)
ax.legend(bbox_to_anchor=(1.05, 1), frameon=False)
ax.figure.set_size_inches((5, 12))
ax.set_xscale("log")
ax.figure.savefig("competing_objectives_optimization.svg", bbox_inches="tight")

