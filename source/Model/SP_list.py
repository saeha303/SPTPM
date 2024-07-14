import copy

class SP_list:
    def __init__(self, period, seq_eventInsIds):
        self.period = period
        self.seq_eventInsIds = copy.deepcopy(seq_eventInsIds)