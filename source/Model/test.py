from .Relation import *
from .Per_ds import *
from .SP_list import *
import math


'''
separate the pattern into seasons -
put the times according to the season -
iterate through the season and times
new function for checking multiple and gcd, store the index of list having the multiple or gcd
'''

def sid_to_season(sid_list, maxper):
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

def get_multiple(diff_list, value):
    idx = -1
    if len(diff_list) == 0:
        return idx
    for index, diff in enumerate(diff_list):
        if value % diff[0] == 0:
            idx = index
            break
    return idx

def get_gcd(diff_list, value):
    idx = -1
    if len(diff_list) == 0:
        return idx
    for index, diff in enumerate(diff_list):
        if math.gcd(diff[0], value) != 1:
            idx = index
            break
    return idx

#for k-pattern timestamp iteration
def generate_pairs(n):
    # 2d matrix
    pairs = []
    for a in range(2, n + 1):
        pairs.append((a - 1, a))
        for b in range(1, a - 1):
            pairs.append((b, a))

            # so for 123, pairs will be 2,1//3,2//3,1
    return pairs

#adjust for k pattern
def getKsp(maxper, event_instance_table, fineness, minconf, minoccur, granularity_of_G, pattern, maxPatternSize, minper = None):
    start_times = {} # contains season number as key and Per_ds object as value
    sid_list = pattern.get_bitmap()
    seasons = sid_to_season(sid_list, maxper)
    total_support_ps = {} #for each season, keep the count of events

    pairs = generate_pairs(maxPatternSize)

    for count, season in seasons:
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
                    obj1 = event_instance_table[list_instances[0][pairs[pair_count][0]]]
                    obj2 = event_instance_table[list_instances[0][pairs[pair_count][1]]]
                    #check the timing
                    if start == -1:
                        start = obj1.start
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
                    
                    if time3 == -1: #latest starting time
                        time3 = temp_time3
                    elif temp_time3 < end and temp_time3 > time3: #if latest starting time is greater than previous
                        time3 = temp_time3

                    ts_count += 1
                    pair_count += 1
                    for i in range(0, size-1):
                        relation = split_relation[i+1]
                        obj1 = event_instance_table[list_instances[0][pairs[pair_count][0]]]
                        obj2 = event_instance_table[list_instances[0][pairs[pair_count][1]]]
                        #check the timing
                        if obj2.end < end:
                            end = obj2.end

                        if relation is Relation.Follows:
                            temp_time3 = obj1.end - granularity_of_G
                        elif relation is Relation.Contains:
                            temp_time3 = obj2.end - granularity_of_G
                        elif relation is Relation.Overlaps:
                            temp_time3 = obj2.start - granularity_of_G
                        
                        if time3 == -1:
                            time3 = temp_time3
                        elif temp_time3 < end and temp_time3 > time3:
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
                        obj1 = event_instance_table[instances[pairs[pair_count][0]]]
                        obj2 = event_instance_table[instances[pairs[pair_count][1]]]
                        #check the timing
                        if start == -1:
                            start = obj1.start
                            end = obj2.end

                        if obj2.end < end:
                            end = obj2.end

                        if relation is Relation.Follows:
                            temp_time3 = obj1.end - granularity_of_G
                        elif relation is Relation.Contains:
                            temp_time3 = obj2.end - granularity_of_G
                        elif relation is Relation.Overlaps:
                            temp_time3 = obj2.start - granularity_of_G
                        
                        if time3 == -1:
                            time3 = temp_time3
                        elif temp_time3 < end and temp_time3 < time3:
                            time3 = temp_time3

                        ts_count += 1
                        pair_count += 1
                        for i in range(0, size-1):
                            relation = split_relation[i+1]
                            obj1 = event_instance_table[instances[pairs[pair_count][0]]]
                            obj2 = event_instance_table[instances[pairs[pair_count][1]]]
                            #check the timing
                            if obj2.end < end:
                                end = obj2.end

                            if relation is Relation.Follows:
                                temp_time3 = obj1.end - granularity_of_G
                            elif relation is Relation.Contains:
                                temp_time3 = obj2.end - granularity_of_G
                            elif relation is Relation.Overlaps:
                                temp_time3 = obj2.start - granularity_of_G
                            
                            if time3 == -1:
                                time3 = temp_time3
                            elif temp_time3 < end and temp_time3 < time3:
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

    #check if periodic, need to adjust according to the seasons ds
    all_sp_list = {} #store seasonal periodic patterns for each season
    sp_count = 0 #counts how many seasonal periodic patterns have been found

    for count, season in start_times:
        max_count = len(season) - len(season) * minconf #minconf is a conf given in decimal, less than 1
        cursor = 1
        per_list = {} #for a particular season, contains cursor as key and diff_list as value
        for events in season[:max_count+1]:
            diff_list = {}
            for next_events in season[cursor:]:
                for event1time in range(events.early_start, events.latest_start, granularity_of_G):
                    diff_list[event1time] = []
                    for event2time in range(next_events.early_start, next_events.latest_start, granularity_of_G):
                        diff = event2time - event1time
                        if minper and diff < minper:
                            continue

                        idx_m = get_multiple(diff_list[event1time], diff)
                        idx_g = get_gcd(diff_list[event1time], diff)

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

    #find support in each season. if threshold is reached, increase sp_count
    #keep a ds to store the periodic timestamps
    candidates = {}
    for count, pers in all_sp_list:
        is_sp = False
        for all_diffs in pers:
            for time, diff_lists in all_diffs:
                for diffs in diff_lists
                    conf = (len(diffs)+1)/total_support_ps[count]
                    if conf >= minconf:
                        is_sp = True
                        if len(diffs)==1:
                            period = diffs[0]
                        else:
                            period = math.gcd(diffs[0], diffs[1])

                        if period not in candidates:
                            candidates[period] = []
                        candidates[period].append(time)
        if is_sp == True:
            sp_count += 1

    #if sp_count >= minseason, then the pattern is seasonal peridic
    if sp_count >= minoccur:
        pattern.per_patterns = candidates
