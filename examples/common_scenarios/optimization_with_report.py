"""Example of optimization from a Genbank record, with report

In this example, the Genbank record data/example_sequence.gbk
gets optimized and a report is created at reports/optimization_with_report/
"""
import os
from dnachisel import DnaOptimizationProblem

genbank_path = os.path.join("data", "example_sequence.gbk")
report_path = "optimization_with_report.zip"  # could also be a folder

problem = DnaOptimizationProblem.from_record(genbank_path)
success, message, _ = problem.optimize_with_report(target=report_path)

print(message + " A report was generated in " + report_path)
