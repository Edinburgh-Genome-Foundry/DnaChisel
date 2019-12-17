class SpecificationSet:
    """Generic class for writing Specs which are actually made of more specs.

    Behaves as a Specification when it comes to instanciation, reading it
    from annotated records, etc. but the initialization actually creates a
    dictionnary of standard Specifications in the DNAOptimizationProblem
    """

    def register_specifications(self, specifications):
        for name, spec in specifications.items():
            spec.parent_specification = self
            spec.name_in_parent = name
        self.specifications = specifications
