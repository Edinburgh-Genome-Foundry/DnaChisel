"""Basic demo of the high-level method optimize_with_report."""
import os
from dnachisel import DnaOptimizationProblem

genbank_path = os.path.join("data", "example_sequence.gbk")
report_folder = os.path.join("reports", "optimization_with_report")

problem = DnaOptimizationProblem.from_record(genbank_path)
success, message, _ = problem.optimize_with_report(target=report_folder)

print(message + " A report was generated in " + report_folder)
