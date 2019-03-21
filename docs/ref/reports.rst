Reports
~~~~~~~~~~~~~~~~~~~~

.. mermaid::
   :align: center

   graph TD;
      D[DnaOptimizationProblem] --> o[optimize_with_report] 
      o -->|constraints clash| N[write_no_solution_error]
      o -->|no constraints clash| W[write_optimization_report]

.. automodule:: dnachisel.reports.optimization_reports
   :members:

.. automodule:: dnachisel.reports.constraints_reports
   :members:
