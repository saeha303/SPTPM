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
                # print(node)

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
    def binary_search(self,arr, target):
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
            # print(id_event)
            # season means node here, supposed to be only one node for each event (chosen ones, ehm)
            # for node in self.Nodes['level1'][id_event]:
            node=self.Nodes['level1'][id_event]
            candidates={}
            # eligible_seasons is a 2d matrix with each row having the H granules
            temp_seasons=node.to_dict(self.EventTable,self.maxPer,self.minPS)
            temp_obj=temp_seasons['patterns']
            eligible_seasons=[]
            if temp_obj:
                eligible_seasons=temp_obj[0]['periodic_intervals']
            # print(eligible_seasons) 
            # flattened_list=list(chain.from_iterable(eligible_seasons))
            if eligible_seasons:
                # print("am def here")
                for season_granules in eligible_seasons:
                    # print("am here")
                    temp_start_end_list=[]
                    for sID in season_granules:
                        current_instances=self.EventTable[id_event].seq_eventInsIds[sID]
                        # print(current_instances)
                        for instance_idx in current_instances:
                            instance=self.EventInstanceTable[instance_idx]
                            if instance.start==instance.end:
                                temp_start_end_list.append(instance.start)
                            else:
                                temp_start_end_list.append(instance.start)
                                temp_start_end_list.append(instance.end)
                    H_start=season_granules[0]
                    H_end=season_granules[-1]
                    # print(temp_start_end_list)
                    total_elements = len(season_granules)
                    last=math.floor((total_elements*self.fineness)/self.minoccur)
                    # print(last)
                    period_start=self.minper*self.G_granularity
                    period_end=(last+1)*self.G_granularity
                    # print(period_start,period_end)
                    inner_start=self.lowestG+H_start*self.fineness*self.G_granularity
                    end=self.lowestG+(H_end+1)*self.fineness*self.G_granularity
                    # print(end)
                    g_start_idx=self.EventTable[id_event].seq_eventInsIds[H_start]
                    g_start_instance=self.EventInstanceTable[g_start_idx[0]]
                    p=g_start_instance.start
                    for period in range(period_start,period_end,self.G_granularity):
                        # print("period:")
                        # print(period)
                        inner_end=inner_start+period
                        # print("start_pos:")
                        for start_pos in range(p,p+period,self.G_granularity):
                            # print(start_pos, end)
                            sID_index=0
                            perfect_periodicity=0
                            count=0
                            # print("i:")
                            for i in range(start_pos,end,period):
                                # print(i)
                                # current_H=season_granules[sID_index]
                                # H_granule_of_i=int(((i-self.lowestG)/self.G_granularity)//(self.fineness))
                                # print("H")
                                # print(current_H,H_granule_of_i)
                                # if H_granule_of_i!=current_H:
                                #     # from here
                                #     sID_index+=1
                                # else:
                                perfect_periodicity+=1
                                result=self.binary_search(temp_start_end_list,i)
                                if result!=-1:
                                    count+=1
                            # print("end of i")
                            print(count,perfect_periodicity)
                            conf=count/perfect_periodicity
                            if conf>=self.minconf:
                                if id_event not in candidates:
                                    candidates[id_event] = {}
                                if period not in candidates[id_event]:
                                    candidates[id_event][period] = []
                                candidates[id_event][period].append(start_pos)
                    pattern=node.get_patterns()    
                    pattern[0].per_patterns=candidates[id_event]
            print(candidates)
            node.perPatterns=candidates
                
        
    def FindKPeriodicPatterns(self):
        level = 2
        while self.maxPatternSize == -1 or level <= self.maxPatternSize:
            level_name = 'level{}'.format(level)
            nodes = self.Nodes[level_name]
            for id_event in nodes:
                node=nodes[id_event]
                patterns = node.get_patterns()
                for pattern in patterns:
                    if pattern.isSeasonal == True:
                        self.getKsp(self.maxPer, self.EventInstanceTable, self.fineness, self.minconf, self.minoccur, self.G_granularity, pattern, level, self.minper)
            level += 1

    def sid_to_season(self, sid_list, maxper):
        season={}
        count=0
        first=sid_list[0]
        for sid in sid_list[1:]:
            per = sid - first

            if per <= maxper:
                if count in season:
                    season[count].append(sid)
                else:
                    season[count] = [sid]
            else:
                count+=1
                season[count] = [sid]

            first = sid
        return season

    def get_multiple(self, diff_list, value):
        idx = -1
        if len(diff_list) == 0:
            return idx
        for index, diff in enumerate(diff_list):
            if value % diff[0] == 0:
                idx = index
                break
        return idx

    def get_gcd(self, diff_list, value):
        idx = -1
        if len(diff_list) == 0:
            return idx
        for index, diff in enumerate(diff_list):
            if math.gcd(diff[0], value) != 1:
                idx = index
                break
        return idx

    #for k-pattern timestamp iteration
    def generate_pairs(self, n):
        # 2d matrix
        pairs = []
        for a in range(2, n + 1):
            pairs.append((a - 1, a))
            for b in range(1, a - 1):
                pairs.append((b, a))

                # so for 123, pairs will be 2,1//3,2//3,1
        return pairs

    #adjust for k pattern
    def getKsp(self, maxper, event_instance_table, fineness, minconf, minoccur, granularity_of_G, pattern, level, minper):
        start_times = {} # contains season number as key and Per_ds object as value
        sid_list = list(pattern.get_bitmap())
        seasons = self.sid_to_season(sid_list, maxper)
        total_support_ps = {} #for each season, keep the count of events

        pairs = self.generate_pairs(level)

        for count, season in seasons.items():
            ts_count = 0 #count number of events in a season
            for sid in season:
                relation_symbol = pattern.get_list_relation()
                list_instances = pattern.get_instance_at_sequence_id(sid)
            
                if len(list_instances) == 1:
                    # not for the case 101: len(list_instances) == 1
                    pair_count = 0
                    size = 1
                    start = -1
                    end = -1
                    time3 = -1
                    while(relation_symbol):
                        split_relation = relation_symbol[:size]
                        relation_symbol = relation_symbol[size:]
                        relation = split_relation[0]
                        obj1 = event_instance_table[list_instances[0][pairs[pair_count][0]-1]]
                        obj2 = event_instance_table[list_instances[0][pairs[pair_count][1]-1]]
                        #check the timing
                        if start == -1:
                            start = obj1.start
                            time3 = start
                            end = obj2.end

                        if obj2.end < end:
                            end = obj2.end

                        #find latest starting time
                        if relation is Relation.Follows:
                            temp_time3 = obj1.end - granularity_of_G
                        elif relation is Relation.Contains:
                            temp_time3 = obj2.end - granularity_of_G
                        elif relation is Relation.Overlaps:
                            temp_time3 = obj2.start - granularity_of_G
                        
                        if temp_time3 >= start and temp_time3 < end and temp_time3 <= time3: #latest starting time
                            time3 = temp_time3

                        ts_count += 1
                        pair_count += 1
                        for i in range(0, size-1):
                            relation = split_relation[i+1]
                            obj1 = event_instance_table[list_instances[0][pairs[pair_count][0]-1]]
                            obj2 = event_instance_table[list_instances[0][pairs[pair_count][1]-1]]
                            #check the timing
                            if obj2.end < end:
                                end = obj2.end

                            if relation is Relation.Follows:
                                temp_time3 = obj1.end - granularity_of_G
                            elif relation is Relation.Contains:
                                temp_time3 = obj2.end - granularity_of_G
                            elif relation is Relation.Overlaps:
                                temp_time3 = obj2.start - granularity_of_G
                            
                            if temp_time3 >= start and temp_time3 < end and temp_time3 <= time3:
                                time3 = temp_time3
                            
                            ts_count += 1
                            pair_count += 1
                        size += 1
                    times = Per_ds(start, time3, list_instances)

                    if count in start_times:
                        start_times[count].append(times)
                    else:
                        start_times[count] = [times]

                else: 
                    # case 101
                    for instances in list_instances:
                        pair_count = 0
                        size = 1
                        start = -1
                        end = -1
                        time3 = -1
                        while(relation_symbol):
                            split_relation = relation_symbol[:size]
                            relation_symbol = relation_symbol[size:]
                            relation = split_relation[0]
                            obj1 = event_instance_table[instances[pairs[pair_count][0]-1]]
                            obj2 = event_instance_table[instances[pairs[pair_count][1]-1]]
                            #check the timing
                            if start == -1:
                                start = obj1.start
                                time3 = start
                                end = obj2.end

                            if obj2.end < end:
                                end = obj2.end

                            if relation is Relation.Follows:
                                temp_time3 = obj1.end - granularity_of_G
                            elif relation is Relation.Contains:
                                temp_time3 = obj2.end - granularity_of_G
                            elif relation is Relation.Overlaps:
                                temp_time3 = obj2.start - granularity_of_G
                            
                            if temp_time3 >= start and temp_time3 < end and temp_time3 <= time3:
                                time3 = temp_time3

                            ts_count += 1
                            pair_count += 1
                            for i in range(0, size-1):
                                relation = split_relation[i+1]
                                obj1 = event_instance_table[instances[pairs[pair_count][0]-1]]
                                obj2 = event_instance_table[instances[pairs[pair_count][1]-1]]
                                #check the timing
                                if obj2.end < end:
                                    end = obj2.end

                                if relation is Relation.Follows:
                                    temp_time3 = obj1.end - granularity_of_G
                                elif relation is Relation.Contains:
                                    temp_time3 = obj2.end - granularity_of_G
                                elif relation is Relation.Overlaps:
                                    temp_time3 = obj2.start - granularity_of_G
                                
                                if temp_time3 >= start and temp_time3 < end and temp_time3 <= time3:
                                    time3 = temp_time3
                                ts_count += 1
                                pair_count += 1
                            size += 1
                        times = Per_ds(start, time3, list_instances)
                        if count in start_times:
                            start_times[count].append(times)
                        else:
                            start_times[count] = [times]
            total_support_ps[count] = ts_count
            #print(start_times)
        #print(total_support_ps)

        #check if periodic, need to adjust according to the seasons ds
        all_sp_list = {} #store seasonal periodic patterns for each season
        sp_count = 0 #counts how many seasonal periodic patterns have been found

        for count, season in start_times.items():
            max_count = int(len(season) - len(season) * minconf) #minconf is a conf given in decimal, less than 1
            #print(max_count)
            cursor = 1
            per_list = {} #for a particular season, contains cursor as key and diff_list as value
            for events in season[:max_count+1]:
                diff_list = {}
                for next_events in season[cursor:]:
                    #print('here?')
                    for event1time in range(events.early_start, events.latest_start+granularity_of_G, granularity_of_G):
                        diff_list[event1time] = []
                        for event2time in range(next_events.early_start, next_events.latest_start+granularity_of_G, granularity_of_G):
                            diff = event2time - event1time
                            if diff < minper:
                                continue

                            idx_m = self.get_multiple(diff_list[event1time], diff)
                            idx_g = self.get_gcd(diff_list[event1time], diff)

                            if idx_m == -1 and idx_g == -1: #if no multiplier or gcd exists
                                if diff > maxper*fineness*granularity_of_G: #diff is outside of periodicity
                                    continue
                                else:
                                    #diff_ds = SP_list(diff, event2time)
                                    #diff_list[event1time] = [diff_ds]
                                    diff_list[event1time] = [diff]

                            elif idx_m != -1: #multiple is in priority
                                last = diff_list[idx_m][-1]
                                if (diff - last) > maxper*fineness*granularity_of_G:
                                    continue
                                else:
                                    #diff_ds = SP_list(diff, event2time)
                                    #diff_list[event1time][idx_m].append(diff_ds)
                                    diff_list[event1time][idx_m].append(diff)

                            elif idx_g != -1: #check gcd
                                last = diff_list[idx_g][-1]
                                if (diff - last) > maxper*fineness*granularity_of_G:
                                    continue
                                else:
                                    #diff_ds = SP_list(diff, event2time)
                                    #diff_list[event1time][idx_g].append(diff_ds)
                                    diff_list[event1time][idx_g].append(diff)

                per_list[cursor] = diff_list  
                cursor += 1 
            all_sp_list[count] = per_list    
        print(all_sp_list)
        #find support in each season. if threshold is reached, increase sp_count
        #keep a ds to store the periodic timestamps
        candidates = {}
        for count, pers in all_sp_list.items():
            is_sp = False
            for all_diffs in pers.values():
                #print(all_diffs)
                for time, diffs in all_diffs.items():
                    conf = (len(diffs)+1)/total_support_ps[count]
                    #print(conf)
                    if conf >= minconf:
                        
                        if len(diffs) == 0:
                            continue
                        elif len(diffs)==1:
                            period = diffs[0]
                        else:
                            period = math.gcd(diffs[0], diffs[1])

                        is_sp = True

                        if period not in candidates:
                            candidates[period] = []
                        candidates[period].append(time)
            if is_sp == True:
                sp_count += 1

        print(minoccur)
        #if sp_count >= minseason, then the pattern is seasonal peridic
        if sp_count >= minoccur:
            print('herehere')
            pattern.per_patterns = candidates
 