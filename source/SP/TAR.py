from SP.FTPMining import FTPMining
from Model import *
import csv
import os
import json

def TAR(config):
    datapath = config.Dataset
    sid = set()
    database = []
    with open(datapath, "r") as csvseq:
        reader = csv.reader(csvseq)
        csvseq.seek(0, 0)
        for row in reader:
            database.append(row)
            # how many H granules are there
            sid.add(int(row[-1]))


    name = config.Name
    resultDir = config.ResultsDir
    maxPatternSize = config.MaxPatternSize
    savePatterns = config.SavePatterns
    epsilon = config.Epsilon
    minOverlap = config.MinOverlap
    maxper = config.MaxPer
    minSR = config.MinSR 
    minPS = config.MinPS
    details = config.InDetails
    mindist = config.Mindist
    maxdist = config.Maxdist
    # ====================================
    minper=config.Minper
    minconf=config.Minconf
    minoccur=config.Minoccur
    granularity_of_G=config.G_granularity
    fineness=config.Fineness

    FTP = FTPMining(database, maxper, minSR, minPS, epsilon, minOverlap, mindist, maxdist, maxPatternSize,minper,minconf,minoccur,granularity_of_G,fineness)

    FTP.Find1Event()
    FTP.Find2FrequentPatterns()
    FTP.FindKFrequentPatterns()
    # ====================================================================== Our code starts here
    FTP.Find1PeriodicPatterns()
    # FTP.FindKPeriodicPatterns()

    if not os.path.exists(resultDir): 
        os.makedirs(resultDir) 
     
    if savePatterns:  
        for level in FTP.Nodes:
            list_json_obj = [] 
            for id in FTP.Nodes[level]: 
                node = FTP.Nodes[level][id]
                if details:
                    json_object = node.to_dict(FTP.EventTable, maxper, minPS, FTP.EventInstanceTable)
                else: 
                    json_object = node.to_dict(FTP.EventTable, maxper, minPS)
                if json_object:
                    list_json_obj.append(json_object)

            with open(os.path.join(resultDir, '{}.json'.format(level)), 'w') as outfile:
                json.dump(list_json_obj, outfile)
                 