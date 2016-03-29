

class Objective:
    pass


class CodonOptimizationObjective(Objective):

    def __init__(self, window, organism):
        self.window = window
        self.organism = organism

    def evaluate(self, sequence):
        pass
