class NoSolutionError(Exception):
    """Exception returned when a DnaOptimizationProblem aborts.
    This means that the constraints are found to be unsatisfiable.
    """

    def __init__(self, message, problem, constraint=None, location=None):
        """Initialize."""
        Exception.__init__(self, message)
        self.message = message
        self.problem = problem
        self.constraint = constraint
        self.location = location

    def __str__(self):
        return self.message
