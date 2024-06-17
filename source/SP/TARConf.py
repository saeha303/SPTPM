class TARConf:
    def __init__(self, name, maxper, minSR, minPS, dataset, max_pattern_size, epsilon, min_overlap, mindist, maxdist, output_path, save_patterns, in_details,minper,minconf,minoccur,granularity_of_G,fineness):
        self.Name = name
        self.Dataset = dataset
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

        # =================================
        self.Minper=minper
        self.Minoccur=minoccur
        self.Minconf=minconf
        self.G_granularity=granularity_of_G
        self.Fineness=fineness