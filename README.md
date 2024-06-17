# Mining Seasonal Temporal Patterns in Time Series

## Paper abstract
Very large time series are increasingly available from an ever wider range of IoT-enabled sensors, from which significant insights can be obtained through mining temporal patterns from them. A useful type of patterns found in many real-world applications exhibits periodic occurrences, and is thus called seasonal temporal pattern (STP). Compared to regular patterns, mining seasonal temporal patterns is more challenging since traditional measures such as support and confidence do not capture the seasonality characteristics. Further, the anti-monotonicity property does not hold for STPs, and thus, resulting in an exponential search space. This paper presents our Frequent Seasonal Temporal Pattern Mining from Time Series (FreqSTPfTS) solution providing: (1) The first solution for seasonal temporal pattern mining (STPM) from time series that can mine STP at different data granularities. (2) The STPM algorithm that uses efficient data structures and two pruning techniques to reduce the search space and speed up the mining process. (3) An approximate version of STPM that uses mutual information, a measure of data correlation, to prune unpromising time series from the search space. (4) An extensive experimental evaluation showing that STPM outperforms the baseline in runtime and memory consumption, and can scale to big datasets. The approximate STPM is up to an order of magnitude faster and less memory consuming than the baseline, while maintaining high accuracy. 

## Prerequisites
- Python 3.7.3 (or later)
- numpy 1.18.2 (or later)
- pandas 1.0.3 (or later)
- sklearn 0.22.2 (or later)

## To run E-STPM 
```
python3 ESTPM.py -i path_input -o path_output -maxper maxPeriod_threshold -minden minDensity_threshold -minsr minSeason_threshold -e epsilon -mindist minDistance -maxdist maxDistance -mps max_patternLength
```

path_input: path of sequence database  
path_output: path of folder contains output  
maxPeriod_threshold: the maximum period threshold  
minDensity_threshold: the minimum density threshold  
minSeason_threshold: the minimum season threshold  
epsilon: the epsilon value  
minDistance: the minimum distance  
maxDistance: the maximum distance  
max_patternLength: the maximum length of a pattern  

Example: 
```
python3 ESTPM.py -i data/SequenceDBRE.csv -o result/ -maxper 8 -minden 7 -minsr 8 -e 4 -mindist 90 -maxdist 270 -mps 3
```

## To run A-STPM
### 1. Compute Mutual Information  
```
python3 MI.py -i path_input -o path_output
```

path_input: path of symbolic database  
path_output: path of output file contains MI values  

Example: 
``` 
python3 MI.py -i data/SymbolicDBRE.csv -o MIResult/MI
```

### 2. Compute mu values for each pair of minSeason and minDensity
```
python3 MU.py -i path_input -o path_output -s season_threshold -d density_threshold -n sequenceDB_size
```

path_input: path of symbolic database  
path_output: path of output file contains MU values corresponding with the minSeason and minDensity thresholds  
season_threshold: the minimum season threshold  
density_threshold: the minimum density threshold  
sequenceDB_size: the size of the sequence database  

Example: 
```
python3 MU.py -i data/SymbolicDBRE.csv -o MUResult/MU -s 8 -d 7 -n 1460
```

### 3. Run approximate STPM
```
python3 ASTPM.py -i path_input -o path_output -imi path_MI -imu path_MU -maxper maxPeriod_threshold -minden minDensity_threshold -minsr minSeason_threshold -e epsilon -mindist minDistance -maxdist maxDistance -mps max_patternLength
```

path_input: path of sequence database  
path_output: path of folder contains output  
path_MI: path of file MI  
path_MU: path of file MU  
maxPeriod_threshold: the maximum period threshold  
minDensity_threshold: the minimum density threshold  
minSeason_threshold: the minimum season threshold  
epsilon: the epsilon value  
minDistance: the minimum distance  
maxDistance: the maximum distance  
max_patternLength: the maximum length of a pattern  

Example: 
```
python3 ASTPM.py -i data/SequenceDBRE.csv -o result/ -imi MIResult/MI -imu MUResult/MU -maxper 8 -minden 7 -minsr 8 -e 4 -mindist 90 -maxdist 270 -mps 3
```
