import os
import argparse
import pandas as pd
from sklearn import metrics

class MutualInformation:
    def __init__(self, pathInput, pathOutput):
        self.pathInput = pathInput
        self.pathOutput = pathOutput
        self.listMI = []

    def CalculateNormalizedMI(self):
        mi_df = pd.read_csv(self.pathInput)
        columns_name = mi_df.columns
        for att1 in columns_name:
            listAnMI = [att1]
            column1 = mi_df[att1].to_numpy()
            mi1 = metrics.mutual_info_score(column1, column1)
            if mi1 == 0.0:
                 mi1=1
            for att2 in columns_name:
                if att1 == att2:
                    nmi2 = 1.0
                else:
                    column2 = mi_df[att2].to_numpy()
                    mi2 = metrics.mutual_info_score(column1, column2)
                    nmi2 = mi2 / mi1
                listAnotherAtt = [att2, nmi2]
                listAnMI.append(listAnotherAtt)
            self.listMI.append(listAnMI)

    def SaveListMI(self):
        dirname = os.path.dirname(self.pathOutput)
        if dirname:
            if not os.path.exists(dirname):
                os.makedirs(dirname)
        wfile = open(self.pathOutput, 'a')
        for f in self.listMI:
            wfile.writelines(str(f) + '\n')
        wfile.close()

parser = argparse.ArgumentParser(description="MI")
parser.add_argument('-i', '--input', help='input', required=True)
parser.add_argument('-o', '--output', help='output', required=True)
args = vars(parser.parse_args())
pathInput = args['input']
pathOutput = args['output']

obMI = MutualInformation(pathInput,pathOutput)
obMI.CalculateNormalizedMI()
obMI.SaveListMI()