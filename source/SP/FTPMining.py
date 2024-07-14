from Model import *
from itertools import chain
import math
class FTPMining:
    def __init__(self, database, maxper, minSR, minPS, epsilon, minOverlap, mindist, maxdist, maxPatternSize,minper,minconf,minoccur,granularity_of_G,fineness):
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
        # =====================================
        self.minper=minper
        self.minconf=minconf
        self.minoccur=minoccur
        self.fineness=fineness
        self.lowestG=0
        self.G_granularity=granularity_of_G

        # calculating SUP set for all the events
        for idx, row in enumerate(database):
            start, end, id_event, sID = row
            sID = int(sID)
            start = int(start)
            end = int(end)
            if idx==0:
                self.lowestG=start
            # eventIns doesn't have sID, cool
            eventIns = EventInstance(idx, start, end)
            # EventInstanceTable has the same data as the database, just instead of rows it has event instances as element
            self.EventInstanceTable[idx] = eventIns

            if id_event not in self.EventTable:
                # if it's a completely new event not seen before
                event = Event(id_event, {sID: [idx]})
                # EventTable stores the events in the format 'time series-value' (C:0) and their instances as in 'corresponding H granule (key):position of instance in the EventInstanceTable(iddx)'
                self.EventTable[id_event] = event
            elif sID not in self.EventTable[id_event].seq_eventInsIds:
                # the event is present but not the granule
                self.EventTable[id_event].seq_eventInsIds[sID] = [idx]
            else:
                # that 1 0 1 case, you have already seen the first 1 and the sequence is (C:1, [G1,G1]), for the second one you need to add a new entry (C:1, [G3,G3])
                self.EventTable[id_event].seq_eventInsIds[sID].append(idx)

    def Find1Event(self):
        # HLH1 in case of 'level1'
        self.Nodes['level1'] = {}
        for id_event in self.EventTable:
            # get all the H granules for the corresponding id_event
            bitmap = self.EventTable[id_event].get_bitmap()
            # SR,ESR=number of seasons and maxSeason respectively, now compare them with minSeason (minSR)
            SR, ESR = Utils.calculateSR_ESR(tuple(bitmap), self.maxPer, self.minPS, self.mindist, self.maxdist)

            list_pattern_candidates = []
            # if maxSeason>=minSeason
            if ESR >= self.minSR:
                # pattern has the whole dictionary of granule:idx entries
                # at most one pattern for every id_event
                pattern = Pattern([], self.EventTable[id_event].seq_eventInsIds)
                pattern.SR = SR
                # candidate patterns
                list_pattern_candidates.append(pattern)

                if SR >= self.minSR:
                    list_pattern_candidates[-1].isSeasonal = True
                    self.num_patterns += 1

            if list_pattern_candidates:
                # candidate single-event patterns, not only the frequent ones, frequent ones will have their isSeasonal field equal to true
                # why is id_event a tuple though
                # at most one node for every id_event
                node = Node(tuple([id_event]), bitmap, list_pattern_candidates)
                # 'level1' means HLH1
                self.Nodes['level1'][(id_event)] = node

    def Find2FrequentPatterns(self):
        # 'level2' means HLH2
        self.Nodes['level2'] = {}
        for node1_id in self.Nodes['level1']:
            node1 = self.Nodes['level1'][node1_id]
            # getting the H granules for node1
            node1_bitmap = set(node1.get_bitmap())

            for node2_id in self.Nodes['level1']:
                node2 = self.Nodes['level1'][node2_id]
                node2_bitmap = set(node2.get_bitmap())
                # the intersection part of two candidate single events
                bitmap = node1_bitmap.intersection(node2_bitmap)

                _, ESR = Utils.calculateSR_ESR(tuple(bitmap), self.maxPer, self.minPS, self.mindist, self.maxdist)
                # check is exceeds maxSeason
                if ESR >= self.minSR:

                    id_event1 = node1.ids_event[0]
                    id_event2 = node2.ids_event[0]
                    # in single-event, we didn't have to think about relation between two events, here we do, that's why a new list
                    # why aren't we checking if they are the same event?
                    list_patterns = self.Find2Patterns(id_event1, id_event2, bitmap) 

                    list_pattern_candidates = [] 
                    # just 3 patterns in list_patterns, one for follows, one for contains and the other for overlaps - for node1 and node2
                    # that means we are gonna have one of each of the follows, contains and overlaps pattern involving two events in the worst case
                    for pattern in list_patterns:  
                        SR, ESR = Utils.calculateSR_ESR(tuple(pattern.get_bitmap()), self.maxPer, self.minPS, self.mindist, self.maxdist)
                        if ESR >= self.minSR:
                            list_pattern_candidates.append(pattern)
                            pattern.SR = SR  
                            if SR >= self.minSR: 
                                list_pattern_candidates[-1].isSeasonal = True 
                                self.num_patterns += 1 

                    if list_pattern_candidates:
                        # list_pattern_candidates can at most have 3 entries for 3 types of relations
                        node = Node((id_event1, id_event2), bitmap, list_pattern_candidates)
                        self.Nodes['level2'][(id_event1, id_event2)] = node 

    def Find2Patterns(self, id_event1, id_event2, bitmap):
        # 2d matrices, row indicates H granules that we got by intersecting id_event1 and id_event2, each has a list of G granule pairs
        # that have some relation between them
        e1_f_e2_instances = {}  # follow
        e1_c_e2_instances = {}  # contain
        e1_o_e2_instances = {}  # overlap
        for sID in bitmap:
            # assume sID is H1, then get all the positions of the G granules in EventInstanceTable for both events
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
        # maximum 3-event pattern, GOOD FOR US
        level = 3
        # self.maxPatternSize == -1 is for ASTPM, for ESTPM the default value is 4
        while self.maxPatternSize == -1 or level <= self.maxPatternSize:
            level_name = 'level{}'.format(level)
            self.Nodes[level_name] = {}
            # level2 nodes for every eligible event-pair
            k_1_Freq = self.Nodes['level{}'.format(level-1)]
            # if there is no (level-1)-pattern, break
            if len(k_1_Freq) == 0:
                break
            # k_1_id_events is an event-pair
            for k_1_id_events in k_1_Freq:
                # only 1 node
                k_1_node = k_1_Freq[k_1_id_events]
                k_1_bitmap = set(k_1_node.get_bitmap())

                for _1_id_event in self.Nodes['level1']:
                    _1_node = self.Nodes['level1'][_1_id_event]
                    _1_bitmap = set(_1_node.get_bitmap())
                    # works with outer H
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
            # if (r(k-n)(k),E(k-n),E(k)) is not in HLH2 for all n=[1,k-1], then return
            if pair_id not in self.Nodes['level2']:
                return []
        # P(k-1)
        k_1_patterns = self.Nodes['level{}'.format(level-1)][k_1_id_events].get_patterns()
        # ABC <- D
        # check CD -> AD -> BD

        # phase 1: check C in ABC and CD
        last_2_id_events = (k_1_id_events[-1], _1_id_event)  # CD,(E(k-1),E(k))
        patterns_last_2_events = self.Nodes['level2'][last_2_id_events].get_patterns()

        temp_pattern_candidates = {}
        # follows, contains, overlaps one by one
        for a_pattern in k_1_patterns:
            a_pattern_bitmap = set(a_pattern.get_bitmap())

            _1_bitmap = set(self.Nodes['level1'][_1_id_event].get_bitmap())
            # works with inner H granules for each pattern
            temp_bitmap = a_pattern_bitmap.intersection(_1_bitmap)

            # find relation between a Pattern of ABC and Patterns of CD
            for last_pattern in patterns_last_2_events:
                last_pattern_bitmap = set(last_pattern.get_bitmap())
                # works with inner H granules for each pattern along with this last_pattern
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
                            # if index matches
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
                                # cde and ce, so c and c, e and e
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
    
    # =========================================================================== Our code
    def binary_search(arr, target):
        low = 0
        high = len(arr) - 1

        while low <= high:
            mid = (low + high) // 2
            if arr[mid] == target:
                return mid
            elif arr[mid] < target:
                low = mid + 1
            else:
                high = mid - 1

        return -1
    def Find1PeriodicPatterns(self):
        
        for id_event in self.Nodes['level1']:
            # season means node here, supposed to be only one node for each event (chosen ones, ehm)
            for node in self.Nodes['level1'][id_event]:
                candidates={}
                # eligible_seasons is a 2d matrix with each row having the H granules
                temp_seasons=node.to_dict(self.event_table,self.maxPer,self.minPS)
                eligible_seasons=temp_seasons['patterns']
                # flattened_list=list(chain.from_iterable(eligible_seasons))
                
                for season_granules in eligible_seasons:
                    temp_start_end_list=[]
                    for sID in season_granules:
                        current_instances=self.EventTable[id_event].seq_eventInsIds[sID]
                        for instance in current_instances:
                            if instance.start==instance.end:
                                temp_start_end_list.append(instance.start)
                            else:
                                temp_start_end_list.append(instance.start)
                                temp_start_end_list.append(instance.end)
                    H_start=season_granules[0]
                    H_end=season_granules[-1]

                    total_elements = len(season_granules)
                    last=math.floor((total_elements*self.fineness)/self.minoccur)
                    period_start=self.minper*self.G_granularity
                    period_end=(last+1)*self.G_granularity

                    inner_start=self.lowestG+H_start*self.fineness*self.G_granularity
                    end=self.lowestG+(H_end+1)*self.fineness*self.G_granularity

                    g_start_idx=self.EventTable[id_event].seq_eventInsIds[H_start]
                    g_start_instance=self.EventInstanceTable[g_start_idx]
                    p=g_start_instance.start
                    for period in range(period_start,period_end,self.G_granularity):
                        inner_end=inner_start+period
                        for start_pos in range(p,p+period,self.G_granularity):
                            sID_index=0
                            perfect_periodicity=0
                            count=0
                            for i in (start_pos,end,period):
                                current_H=season_granules[sID_index]
                                H_granule_of_i=((i-self.lowestG)/self.G_granularity)//(self.fineness-1)
                                if H_granule_of_i!=current_H:
                                    sID_index+=1
                                else:
                                    perfect_periodicity+=1
                                    result=self.binary_search(temp_start_end_list,i)
                                    if result!=-1:
                                        count+=1
                            conf=count/perfect_periodicity
                            if conf>=self.minconf:
                                if id_event not in candidates:
                                    candidates[id_event] = {}
                                if period not in candidates[id_event]:
                                    candidates[id_event][period] = []
                                candidates[id_event][period].append(start_pos)
                node.perPatterns=candidates
                
        
    def FindKPeriodicPatterns(self):
        getKsp(self.maxPer, self.EventInstanceTable, self.fineness, minsup, self.minoccur, self.granularity_of_G, pattern, self.maxPatternSize, self.minper)

 