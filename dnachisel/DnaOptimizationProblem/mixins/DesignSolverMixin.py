class DesignSolverMixin:
    def _init_disign(self):
        self.scores = {}
        self.levels = {}
        self._set_scores()
        self._set_level()
        
    def _set_scores(self):
        for i in range(len(self.design_space.feature_label)):
            self.scores[self.design_space.feature_label[
                i]] = self.design_space.feature_specification[i].evaluate(
                    self).score

    def _set_level(self):
        for i in range(len(self.design_space.feature_label)):
            feature_label = self.design_space.feature_label[i]
            self.levels[
                feature_label +
                "_Level"] = self.design_space.range_set_list[i].get_level(
                    self.scores[feature_label])

    def is_match_design(self, desired_design):
        levels = [
            str(self.levels[feature + '_Level'])
            for feature in self.design_space.feature_label
        ]

        return '.'.join(levels) == desired_design
