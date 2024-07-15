from SP import *
from Model import Utils
import argparse
import os

parser = argparse.ArgumentParser(description='STPM')
parser.add_argument('-name', '--name', help='Name of experiment', default='Experiment')
parser.add_argument('-i', '--input', help='input path', required=True)
parser.add_argument('-o', '--output', help='output path', required=True)
parser.add_argument('-maxper', '--maxper', help='maximum period', default=100)
parser.add_argument('-minsr', '--minsr', help='recurrence threshold', default=0.5)
parser.add_argument('-minden', '--minden', help='periodic support threshold', default=0.5)
parser.add_argument('-e', '--epsilon', help='epsilon', default=0)
parser.add_argument('-mo', '--minoverlap', help='min overlap', default=1)
parser.add_argument('-mps', '--maxpatternsize', help='max pattern size', default=4)
parser.add_argument('-mindist', '--mindistance', help='min distance', default=1)
parser.add_argument('-maxdist', '--maxdistance', help='max distance', default=100)
# ==========================================
parser.add_argument('-minper', '--minper', help='minimum period', default=2)
parser.add_argument('-minconf', '--minconf', help='minimum confidence', default=0.6)
parser.add_argument('-minoccur', '--minoccur', help='minimum occurence', default=3)
parser.add_argument('-granularity_of_G', '--granularity_of_G', help='granularity_of_G', default=3600)
parser.add_argument('-fineness', '--fineness', help='fineness', default=24)


parser.add_argument('-nosave', '--nosavepattern', help='save pattern to file', default=False, action='store_true')
parser.add_argument('-details', '--details', help='write to file in details mode', default=True, action='store_true')
args = vars(parser.parse_args())


dataset = args['input']
maxper = float(args['maxper'])
minSR = float(args['minsr'])
minPS = float(args['minden'])
epsilon = float(args['epsilon'])
minoverlap = float(args['minoverlap'])
name = args['name']
max_pattern_size = int(args['maxpatternsize'])
mindist = int(args['mindistance'])
maxdist = int(args['maxdistance'])
# ===================================
minper=int(args['minper'])
minconf=float(args['minconf'])
minoccur=int(args['minoccur'])
granularity_of_G=int(args['granularity_of_G'])
fineness=int(args['fineness'])

save_patterns = not bool(args['nosavepattern'])
in_details = bool(args['details'])
output = os.path.join(args['output'], '{}_maxper{}_minden{}_minSR{}'.format(name, maxper, minPS, minSR))


Utils.checkArgs(epsilon, minoverlap, max_pattern_size)

config = TARConf(name, maxper, minSR, minPS, dataset, max_pattern_size, epsilon, minoverlap, mindist, maxdist, output, save_patterns, in_details,minper,minconf,minoccur,granularity_of_G,fineness)

TAR(config)
