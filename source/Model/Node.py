from .Relation import *
from .Utils import *
from .EventInstance import *
from .Event import*
from .Pattern import*


class Node:
    def __init__(self, ids_event, bitmap, patterns=None,perPatterns=None):
        self.ids_event = ids_event
        self.bitmap = bitmap
        self.patterns = patterns
        # ==================================
        self.perPatterns = perPatterns

    def get_ids_events(self):
        return tuple(self.ids_event)

    def get_bitmap(self):
        return self.bitmap

    def get_support(self):
        return len(self.bitmap)

    def get_patterns(self):
        return self.patterns

    def to_dict(self, event_table, maxper, minPS, event_instance_table=None):
        result = {}
        event_labels = self.get_ids_events()
        pattern_result = []
        # there can be only one pattern in self.patterns for single event
        for pattern in self.patterns:
            if pattern.isSeasonal:
                temp_obj = pattern.to_dict(event_labels, maxper, minPS,
                                      event_instance_table)
                obj=temp_obj['periodic_intervals']
                pattern_result.append(obj)

        # =================================================
        event_name = ','.join(event_labels)

        if pattern_result:
            result['name_node'] = event_name
            result['patterns'] = pattern_result

        # return result

        return pattern_result

    def __str__(self):
        return ','.join(self.ids_event)
