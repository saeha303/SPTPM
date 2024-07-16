from .Relation import *
from .Utils import *

class Per_ds:
    def __init__(self, early_start, latest_start, list_instance):
        self.early_start=early_start
        self.latest_start=latest_start
        self.list_instance = list_instance

    def getTimes(self):
        return [self.early_start, self.latest_start]

    def __str__(self):
        return (f"Per_ds(early_start={self.early_start}, "
                f"latest_start={self.latest_start})")

    def __repr__(self):
        return self.__str__()