class MITARConf:
    def __init__(self, name, maxper, minSR, minPS, dataset, MIs, MUs, max_pattern_size, epsilon, min_overlap, mindist, maxdist, output_path, save_patterns, in_details):
        self.Name = name
        self.Dataset = dataset
        self.MIs = MIs
        self.MUs = MUs
        self.ResultsDir = output_path
        self.MaxPatternSize = max_pattern_size
        self.SavePatterns = save_patterns
        self.Epsilon = epsilon
        self.MinOverlap = min_overlap
        self.Mindist = mindist
        self.Maxdist = maxdist
        self.MaxPer = maxper
        self.MinSR = minSR
        self.MinPS = minPS
        self.InDetails = in_details
