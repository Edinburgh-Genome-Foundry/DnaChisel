class NoSolutionError(Exception):
    """Exception returned when a DnaOptimizationProblem aborts.
    This means that the constraints are found to be unsatisfiable.
    """

    def __init__(self, message, problem, constraint=None, location=None):
        """Initialize."""
        # Passing all of our args into the superclass constructor allows
        # this exception to roundtrip through pickle;
        # https://stackoverflow.com/a/41809333
        Exception.__init__(self, message, problem, constraint, location)
        self.message = message
        self.problem = problem
        self.constraint = constraint
        self.location = location

    def __str__(self):
        return self.message
