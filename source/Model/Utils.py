from .Relation import *
from datetime import *
import math 

class Utils:
    def check_relation(time_a, time_b, label_a, label_b, epsilon, min_overlap):
        relation = None
        if time_a.start == time_b.start:
            if time_a.end == time_b.end:
                # same relation check with same span is done here, why > though
                if label_a >= label_b:
                    return relation
        # following condition is not allowed in any of the 3 relations
        elif time_a.start > time_b.start:
            return relation

        if label_a != label_b:
            # epsilon amount of overlap is allowed?
            if time_a.end - epsilon <= time_b.start:
                relation = Relation.Follows
            elif (time_a.start <= time_b.start) and (time_a.end + epsilon >= time_b.end):
                relation = Relation.Contains
            elif (time_a.start < time_b.start) and (time_a.end - epsilon < time_b.end) and (time_a.end - time_b.start >= min_overlap - epsilon):
                relation = Relation.Overlaps
        else:
            # same relation but not the same span? Is it even possible?
            if time_a.end - epsilon <= time_b.start:
                relation = Relation.Follows
        return relation

    def checkDistance(IPI, value,  mindist, maxdist):
        if len(IPI) == 0:
            return True
        # [-1][1]: last row, second column
        distance = value - IPI[-1][1]
        if mindist <= distance <= maxdist:
            return True
        return False

    # I'm a bit confused about this part, what is SR and ESR actually?
    # we check 3 variables here except minSeason
    def calculateSR_ESR(sID_list, maxper, minPS, mindist, maxdist):
        # instead of detecting the SUP set first and separating NearSUP sets from it like the algorithm, we are doing it in every step
        # so we are going into the trouble of calculating NearSUP, distInterval and things since calculating SUP is already taking a lot of work

        # sort H granules
        sID_list = sorted(sID_list)
        Esr = 0
        IPI = []
        # what's the point of adding a -1 at the last
        sID_list.append(-1)
        # ultimately candidateESR and candidate are the same thing, candidate= prevCandidateESR
        candidate = None
        # SUP set
        candidateESR = None
        
        # value is an iterator, [value] is the value pointed by the iterator
        for value in sID_list:
            # only for the first element in candidateESR
            if candidateESR is None:
                # initializing
                candidateESR = [value]
            else:  
                perESR = (value) - (candidateESR[-1])
                if  0<perESR <= maxper:
                    # last entry satisfies maxPeriod
                    candidateESR.append(value)
                else:
                    # calculating maxSeason, all the candidateESR add up to the SUP set
                    Esr += math.floor(len(candidateESR)/minPS)
                    # instead of appending, start a new list
                    candidateESR = [value]
            # following condition will be true for the first iteration
            # also when doesn't satisfy distInterval below
            if candidate is None:
                if Utils.checkDistance(IPI, value, mindist, maxdist):
                    # initializing
                    candidate = [value]
                continue            
            per = (value) - (candidate[-1])
            
            if  0<per <= maxper:
                candidate.append(value)
            else:
                ps = len(candidate)
                # if greater than minDensity
                if ps >= minPS:
                    # add a new row to IPI
                    IPI.append((candidate[0], candidate[-1]))
                # difference between last entry of the last season and current iterator
                if Utils.checkDistance(IPI, value, mindist, maxdist):
                    # initializing
                    candidate = [value]
                else:
                    candidate = None 
        # number of seasons that satisfy distInterval
        Sr = len(IPI) 
        return Sr, Esr
    # get seasons, only seasons satisfying maxPeriod and minDensity
    def getIPI_PS(sID_list, maxper, minPS): 
        sID_list = sorted(sID_list) 
        result = []
        sID_list.append(-1)
        # initializing with first H value
        candidate = [sID_list[0]] 
        # starting from 2nd H value
        for value in sID_list[1:]:
            per = (value) - (candidate[-1])
            if 0<per <= maxper:
                candidate.append(value)
            else:
                ps = len(candidate)
                # exceeds minDensity
                if ps >= minPS:
                    # can't be string
                    # result[str((candidate[0], candidate[-1]))] = ps
                    result.append(candidate.copy())
                candidate = [value]   
        return result

    def checkArgs(epsilon, minoverlap, maxpatternsize):
        if 0 > epsilon:
            raise ValueError("epsilon must be larger or equal 0")
        if 0 >= minoverlap:
            raise ValueError("minoverlap must be larger 0")
        if 2 > maxpatternsize:
            if maxpatternsize != -1:
                raise ValueError("maxpatternsize must be larger or equal 2")

    