from distutils.command.config import config
from Model import *


class MIFTPMining:
    def __init__(self, database, satisfiedAtts, maxper, minSR, minPS, epsilon, minOverlap, mindist, maxdist, maxPatternSize):
        self.satisfiedAtts = satisfiedAtts
        self.epsilon = epsilon
        self.minoverlap = minOverlap
        self.maxPatternSize = maxPatternSize
        self.EventInstanceTable = {}
        self.EventTable = {}
        self.Nodes = {}
        self.num_patterns = 0
        self.maxPer = maxper
        self.minSR = minSR
        self.minPS = minPS
        self.mindist = mindist
        self.maxdist = maxdist

        for idx, row in enumerate(database):
            start, end, id_event, sID = row
            sID = int(sID)
            start = int(start)
            end = int(end)
            eventIns = EventInstance(idx, start, end)
            self.EventInstanceTable[idx] = eventIns

            if id_event not in self.EventTable:
                event = Event(id_event, {sID: [idx]})
                self.EventTable[id_event] = event
            elif sID not in self.EventTable[id_event].seq_eventInsIds:
                self.EventTable[id_event].seq_eventInsIds[sID] = [idx]
            else:
                self.EventTable[id_event].seq_eventInsIds[sID].append(idx)

    def Find1Event(self):
        self.Nodes['level1'] = {}
        for id_event in self.EventTable:
            bitmap = self.EventTable[id_event].get_bitmap()
            SR, ESR = Utils.calculateSR_ESR(tuple(bitmap), self.maxPer, self.minPS, self.mindist, self.maxdist)

            list_pattern_candidates = []

            if ESR >= self.minSR:
                pattern = Pattern([], self.EventTable[id_event].seq_eventInsIds)
                pattern.SR = SR
                list_pattern_candidates.append(pattern)

                if SR >= self.minSR:
                    list_pattern_candidates[-1].isSeasonal = True
                    self.num_patterns += 1

            if list_pattern_candidates:
                node = Node(tuple([id_event]), bitmap, list_pattern_candidates)
                self.Nodes['level1'][(id_event)] = node

    def Find2FrequentPatterns(self):
        self.Nodes['level2'] = {}
        for node1_id in self.Nodes['level1']:
            att1 = node1_id[:node1_id.index('-')]
            #att1 = node1_id[:node1_id.index('*')]
            if att1 not in self.satisfiedAtts:
                continue
            node1 = self.Nodes['level1'][node1_id]
            node1_bitmap = set(node1.get_bitmap())

            for node2_id in self.Nodes['level1']:
                att2 = node2_id[:node2_id.index('-')]
                #att2 = node2_id[:node2_id.index('*')]
                if att2 not in self.satisfiedAtts[att1]:
                    continue
                node2 = self.Nodes['level1'][node2_id]
                node2_bitmap = set(node2.get_bitmap())

                bitmap = node1_bitmap.intersection(node2_bitmap)

                _, ESR = Utils.calculateSR_ESR(tuple(bitmap), self.maxPer, self.minPS, self.mindist, self.maxdist)
                if ESR >= self.minSR:

                    id_event1 = node1.ids_event[0]
                    id_event2 = node2.ids_event[0]
                    list_patterns = self.Find2Patterns(id_event1, id_event2, bitmap)

                    list_pattern_candidates = []

                    for pattern in list_patterns:
                        SR, ESR = Utils.calculateSR_ESR(tuple(pattern.get_bitmap()), self.maxPer, self.minPS, self.mindist, self.maxdist)
                        if ESR >= self.minSR:
                            list_pattern_candidates.append(pattern)
                            pattern.SR = SR
                            if SR >= self.minSR:
                                list_pattern_candidates[-1].isSeasonal = True
                                self.num_patterns += 1

                    if list_pattern_candidates:
                        node = Node((id_event1, id_event2), bitmap, list_pattern_candidates)
                        self.Nodes['level2'][(id_event1, id_event2)] = node

    def Find2Patterns(self, id_event1, id_event2, bitmap):
        e1_f_e2_instances = {}  # follow
        e1_c_e2_instances = {}  # contain
        e1_o_e2_instances = {}  # overlap
        for sID in bitmap:
            list_event_instances1 = self.EventTable[id_event1].get_list_instance_at_sequence_id(sID)
            list_event_instances2 = self.EventTable[id_event2].get_list_instance_at_sequence_id(sID)
            for id_event_instances1 in list_event_instances1:
                time_1 = self.EventInstanceTable[id_event_instances1]
                for id_event_instances2 in list_event_instances2:
                    time_2 = self.EventInstanceTable[id_event_instances2]
                    relation = Utils.check_relation(time_1, time_2, id_event1, id_event2, self.epsilon, self.minoverlap)

                    if relation is Relation.Follows:
                        if sID in e1_f_e2_instances:
                            e1_f_e2_instances[sID].append((id_event_instances1, id_event_instances2))
                        else:
                            e1_f_e2_instances[sID] = [(id_event_instances1, id_event_instances2)]

                    elif relation is Relation.Contains:
                        if sID in e1_c_e2_instances:
                            e1_c_e2_instances[sID].append((id_event_instances1, id_event_instances2))
                        else:
                            e1_c_e2_instances[sID] = [(id_event_instances1, id_event_instances2)]

                    elif relation is Relation.Overlaps:
                        if sID in e1_o_e2_instances:
                            e1_o_e2_instances[sID].append((id_event_instances1, id_event_instances2))
                        else:
                            e1_o_e2_instances[sID] = [(id_event_instances1, id_event_instances2)]

        list_patterns = []

        follow_pattern = Pattern([Relation.Follows], e1_f_e2_instances)
        list_patterns.append(follow_pattern)

        overlap_pattern = Pattern([Relation.Overlaps], e1_o_e2_instances)
        list_patterns.append(overlap_pattern)


        contain_pattern = Pattern([Relation.Contains], e1_c_e2_instances)
        list_patterns.append(contain_pattern)

        return list_patterns

    def FindKFrequentPatterns(self):
        level = 3
        while self.maxPatternSize == -1 or level <= self.maxPatternSize:
            level_name = 'level{}'.format(level)
            self.Nodes[level_name] = {}
            k_1_Freq = self.Nodes['level{}'.format(level-1)]
            if len(k_1_Freq) == 0:
                break
            for k_1_id_events in k_1_Freq:
                k_1_node = k_1_Freq[k_1_id_events]
                k_1_bitmap = set(k_1_node.get_bitmap())

                for _1_id_event in self.Nodes['level1']:
                    _1_node = self.Nodes['level1'][_1_id_event]
                    _1_bitmap = set(_1_node.get_bitmap())

                    bitmap = k_1_bitmap.intersection(_1_bitmap)

                    _, ESR = Utils.calculateSR_ESR(tuple(bitmap), self.maxPer, self.minPS, self.mindist, self.maxdist)

                    if ESR >= self.minSR:

                        list_patterns = self.FindKPatterns(k_1_id_events, _1_id_event, level) 
                        list_pattern_candidates = [] 

                        for pattern in list_patterns: 
                            SR, ESR = Utils.calculateSR_ESR(tuple(pattern.get_bitmap()), self.maxPer, self.minPS, self.mindist, self.maxdist)
                            if ESR  >= self.minSR:
                                list_pattern_candidates.append(pattern)
                                pattern.SR = SR  
                                if SR >= self.minSR: 
                                    list_pattern_candidates[-1].isSeasonal = True 
                                    self.num_patterns += 1 

                        if list_pattern_candidates:
                            temp = list(k_1_id_events)
                            temp.append(_1_id_event)
                            k_id_events = tuple(temp)
                            
                            node = Node(k_id_events, bitmap, list_pattern_candidates)
                            self.Nodes[level_name][k_id_events] = node
            
            level += 1

    def FindKPatterns(self, k_1_id_events, _1_id_event, level):
        for prev_id in k_1_id_events:
            pair_id = (prev_id, _1_id_event)
            if pair_id not in self.Nodes['level2']:
                return []

        k_1_patterns = self.Nodes['level{}'.format(level-1)][k_1_id_events].get_patterns()
        # ABC <- D
        # check CD -> AD -> BD

        # phase 1: check C in ABC and CD
        last_2_id_events = (k_1_id_events[-1], _1_id_event)  # CD
        patterns_last_2_events = self.Nodes['level2'][last_2_id_events].get_patterns()

        temp_pattern_candidates = {}

        for a_pattern in k_1_patterns:
            a_pattern_bitmap = set(a_pattern.get_bitmap())

            _1_bitmap = set(self.Nodes['level1'][_1_id_event].get_bitmap())
            temp_bitmap = a_pattern_bitmap.intersection(_1_bitmap)

            # find relation between a Pattern of ABC and Patterns of CD
            for last_pattern in patterns_last_2_events:
                last_pattern_bitmap = set(last_pattern.get_bitmap())

                temp_bitmap = tuple(a_pattern_bitmap.intersection(last_pattern_bitmap))
                _, ESR = Utils.calculateSR_ESR(temp_bitmap, self.maxPer, self.minPS, self.mindist, self.maxdist)

                if ESR < self.minSR:
                    continue

                pattern_sID_instances = {}
                for sequenceID in temp_bitmap:
                    current_list_instances = a_pattern.get_instance_at_sequence_id(sequenceID)
                    last_list_instances = last_pattern.get_instance_at_sequence_id(sequenceID)

                    for current_instance in current_list_instances:
                        for last_instance in last_list_instances:
                            if current_instance[-1] == last_instance[0]:
                                if sequenceID in pattern_sID_instances:
                                    pattern_sID_instances[sequenceID].append(
                                        tuple(list(current_instance) + [last_instance[-1]]))
                                else:
                                    pattern_sID_instances[sequenceID] = [
                                        tuple(list(current_instance) + [last_instance[-1]])]

                pattern_name = tuple(a_pattern.get_list_relation() + last_pattern.get_list_relation())
                temp_pattern_candidates[pattern_name] = pattern_sID_instances

        # Phase 2:
        if len(temp_pattern_candidates) == 0:
            return []

        for index, id_event in enumerate(k_1_id_events[:-1]):
            two_events = (id_event, _1_id_event)  # AD, BD
            patterns_two_events = self.Nodes['level2'][two_events].get_patterns()

            new_pattern_candidates_update = {}

            for two_pattern in patterns_two_events:  # loop parttern in AD and check event instances in AD and current pattern candidates
                two_pattern_bitmap = set(two_pattern.get_bitmap())

                for pattern_candidate in temp_pattern_candidates:
                    current_pattern_instance = temp_pattern_candidates[pattern_candidate]
                    current_pattern_bitmap = set(current_pattern_instance.keys())

                    temp_bitmap = tuple(current_pattern_bitmap.intersection(two_pattern_bitmap))
                    
                    _, ESR = Utils.calculateSR_ESR(temp_bitmap, self.maxPer, self.minPS, self.mindist, self.maxdist)

                    if ESR < self.minSR:
                        continue

                    pattern_sID_instances = {}

                    for sequenceID in temp_bitmap:
                        current_list_instances = current_pattern_instance[sequenceID]
                        two_list_instances = two_pattern.get_instance_at_sequence_id(sequenceID)

                        for current_instance in current_list_instances:
                            for two_instance in two_list_instances:
                                if current_instance[index] == two_instance[0] and current_instance[-1] == two_instance[-1]:
                                    if sequenceID in pattern_sID_instances:
                                        pattern_sID_instances[sequenceID].append(current_instance)
                                    else:
                                        pattern_sID_instances[sequenceID] = [current_instance]

                    pattern_name = tuple(list(pattern_candidate) + two_pattern.get_list_relation())
                    new_pattern_candidates_update[pattern_name] = pattern_sID_instances

            temp_pattern_candidates = new_pattern_candidates_update  # update temp pattern for next check

        list_pattern = []

        for pattern_name, pattern_sID_instances in temp_pattern_candidates.items():
            tmp = Pattern(pattern_name, pattern_sID_instances)
            list_pattern.append(tmp)

        patterns = list_pattern

        return patterns
 
