# Example with competing objectives


**Note:** install seaborn (pip install seaborn) for to run this example.


This example explores different variations of a problem where a gene is optimized using synonymous mutations in order to increase the sum of 5 competing objective scores:

- ``unique_kmers(20)`` ensures that no k-mer of size 20 or more have homologies elsewhere in the sequence (the negative score is the number of non-unique k-mers)
- ``~gc(39%/100bp)`` is an objective targeting 39% GC on every 100bp window (the negative score is the sum of all absolute deviations from 39% on all 100bp windows)
- ``~no(CG)`` seeks to minimize the occurences of "CG" in the sequence (the negative score is the number of CG patterns in the sequence).
- ``~codon_optimize`` seeks to maximize the CAI gene's CAI for *E. coli* (the negative score is the sum of each codon's deviation from its optimal relative frequency).
- ``~keep`` is an objective fostering parsimony in sequence changes by penalizing differences (the negative score is the number of changes in the sequence).

The script runs variations of this problem where:
- All objectives are considered at once, with equal weight ("equal weights")
- Only a single objective is considered ("~keep only", "~gc(39%) only", etc.)
- One objective is given 5 times more weight than the others ("~keep 5x", etc.).

The results are plotted in ``competing_objectives_optimization.svg``

Some remarks:

- When optimized separately as single objectives, all objectives can be fully optimized (e.g. no CG at all, or no non-unique k-mers), with the exception of ``~gc(39%)`` as it is not possible to ensure exactly 39% GC on every 100bp window doing only synonymous mutations (but the score is still improved by a factor of 7.5, from 309 to 41)
- Applying a 5x boost to an objective has the expected effect: this objective becomes much more optimized, at the expense of others.
- The ``unique_kmers`` seems to be a preferred target of the algorithm due to the easiness with which k-mers can be removed (a single mutation can remove several non-unique kmers at once).

![](competing_objectives_optimization.svg)