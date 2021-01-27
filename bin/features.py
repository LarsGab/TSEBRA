#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# features.py: Classes for creating a list of features for a transcript
# ==============================================================

class Edge_features:
    def __init__(self, tx1, tx2, evi):
        # tx1, tx2, evi are lists of lists of gtf lines of all introns,
        # start- and stop-codons for two overlapping transcripts and relavent evidence

        # sets of ids for all introns and start-/stop-codons of tx1, tx2
        # id: 'chr_strand_hinttype_start_end'
        self.tx1 = self.__init_set__(tx1)
        self.tx2 = self.__init_set__(tx2)
        # list of ids of hints with evidence
        # can include duplicates, as an intron can be
        # supported by mulitple sources
        self.evi = self.__init_list__(evi)
        self.evi_set = set(self.evi)

        # feature vector specifies the support of introns, start/stop codons
        # in tx1 and tx2 by the other transcript or the evidence
        # In the followinf is n := symbol for intersection of sets,
        #   - := symbol for the relative complement (difference) of sets
        #   tx1, tx2 sets of ids
        #   evi: list of ids of hints with evidence (with duplicates)
        #   evi_set: set of ids of evidence hints (no duplicates)
        #   (evi_set = set(evi))
        # feature_vector:= [
        #       |tx1 n evi_set| / |tx1| ,
        #       1 - |evi - tx1| / |evi| ,
        #       |tx2 n evi_set| / |tx2| ,
        #       1 - |evi - tx2| / |evi|
        # ]
        self.feature_vector = self.__create_feature_vector__()

    def __init_set__(self, list):
        id_set = set()
        for l in list:
            id_set.add('{}_{}_{}_{}_{}'.format(l[0], l[6], l[2], l[3], l[4]))
        return id_set

    def __init_list__(self, list):
        id_list = []
        for l in list:
            id_list.append('{}_{}_{}_{}_{}'.format(l[0], l[6], l[2], l[3], l[4]))
        return id_list

    def __create_feature_vector__(self):
        vec = []
        size_tx1 = len(self.tx1)
        size_tx2 = len(self.tx2)
        size_evi = len(self.evi)
        if size_tx1 == 0:
            vec.append(1)
        else:
            vec.append(len(self.tx1.intersection(self.evi_set)) / size_tx1)
        if size_evi == 0:
            vec.append(0)
        else:
            vec.append(1.0 - (len([l for l in self.evi if l not in self.tx1]) / size_evi))
        if size_tx2 == 0:
            vec.append(1)
        else:
            vec.append(len(self.tx2.intersection(self.evi_set)) / size_tx2)
        if size_evi == 0:
            vec.append(0)
        else:
            vec.append(1.0 - (len([l for l in self.evi if l not in self.tx2]) / size_evi))

        return vec
