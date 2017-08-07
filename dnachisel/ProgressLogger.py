from tqdm import tqdm
from copy import deepcopy

class ProgressLogger:

    def __init__(self, init_state=None):

        self.state = {}
        if init_state is not None:
            self.state.update(init_state)
    def callback(self, **kw):
        pass

    def __call__(self, **kw):
        self.state.update(kw)
        self.callback(**kw)

class ProgressBarLogger(ProgressLogger):

    def __init__(self, bars_descriptions):
        bars_descriptions = [
            bar if isinstance(bar, (tuple, list)) else
            (bar, bar.lower() + '_index', bar.lower() + '_total')
            for bar in bars_descriptions
        ]
        self.bars_descriptions = bars_descriptions
        self.bars = [
            None
            for i in range(len(bars_descriptions))
        ]
        self.key_to_bar = {
            key: (i, role)
            for i, (label, index, total) in enumerate(bars_descriptions)
            for key, role in zip((index, total), ('index', 'total'))
        }
        ProgressLogger.__init__(self, {
            key: 0
            for (label, index, total) in bars_descriptions
            for key in (index, total)
        })
        self.previous_state = deepcopy(self.state)

    def new_bar(self, bar_id):
        desc, _, total = self.bars_descriptions[bar_id]
        total = self.state[total]
        self.bars[bar_id] = tqdm(total=total, desc=desc)

    def callback(self, **kw):

        for label, value in kw.items():

            if label not in self.key_to_bar:
                continue

            old_value = self.previous_state[label]
            self.previous_state[label] = value

            bar_id, role = self.key_to_bar[label]
            if role == 'total':
                self.new_bar(bar_id)
            else:
                if self.bars[bar_id] is None:
                    self.new_bar(bar_id)
                self.bars[bar_id].update(max(0, value - old_value))
                _, _, total = self.bars_descriptions[bar_id]
                if value >= self.state[total]:
                    self.bars[bar_id].close()
                    self.bars[bar_id] = None

class DnaOptimizationProgressBar(ProgressBarLogger):

    def __init__(self, specs=True, locations=True, mutations=False):
        bars = []
        if specs:
            bars += ['Objective', 'Constraint']
        if locations:
            bars.append('Location')
        if mutations:
            bars.append('Mutations')
        ProgressBarLogger.__init__(self, bars_descriptions=bars)


if __name__ == '__main__':
    import time
    logger = ProgressBarLogger(bars_descriptions=[('I', 'i', 'total_i'),
                                                  ('J', 'j', 'total_j'),
                                                  ('K', 'k', 'total_k')])
    logger(total_i=9, total_j=10, total_k=11)
    for i in range(10):
        logger(i=i)
        for j in range(11):
            logger(j=j)
            for k in range(12):
                logger(k=k)
                time.sleep(0.05)
