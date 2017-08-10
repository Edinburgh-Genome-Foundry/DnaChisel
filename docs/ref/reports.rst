Reports
~~~~~~~~~~~~~~~~~~~~

.. mermaid::

   graph TD;
      o[optimization_with_report] <--> D[DnaOptimizationProblem]
      o -->|constraints clash| N[write_no_solution_error]
      o -->|no constraints clash| W[write_optimization_report]

.. automodule:: dnachisel.reports.optimization_reports
   :members:
.. automodule:: dnachisel.reports.constraints_reports
  :members:
