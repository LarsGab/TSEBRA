#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# features.py: Handles the features for a transcript
# ==============================================================

class Node_features:
    def __init__(self, tx, evi, hint_source_weight={'P' : 0.1, 'E' : 10, 'C' : 5,  'M' : 1}):

        self.sw = hint_source_weight
        self.evi_list = {'intron' : [], 'start_codon' : [], 'stop_codon': []}
        self.numb_introns = 0
        self.__init_hints__(tx, evi)
        # feature vector specifies the support of
        # introns, start/stop codons for a transcript
        # self.feature_vector[0] : (supported introns by evidence of tx) / (number of introns in tx)
        # self.feature_vector[1] : (supported start/stop codons by evidence of tx) / (number of \
        # start/stop codons in tx)
        # self.feature_vector[2] : sum of multiplicities of intron evidence for tx
        # self.feature_vector[3] : sum of multiplicities of start/stop codon evidence for tx
        # self.feature_vector[4] : 1 if tx is from anno_pref, 0 otherwise
        self.feature_vector = self.create_feature_vec()

    def __init_hints__(self, tx, evi):
        cds_len = 0
        for type in ['intron', 'start_codon', 'stop_codon']:
            for line in tx.transcript_lines[type]:
                hint = evi.get_hint(line[0], line[3], line[4], line[2], \
                    line[6])
                if hint:
                    self.evi_list[type].append(hint)
        if tx.transcript_lines['intron']:
            self.numb_introns = len(tx.transcript_lines['intron'])

    def create_feature_vec(self):
        return [self.relative_support(['intron'], self.numb_introns), \
                self.relative_support(['start_codon', 'stop_codon'], 2.0), \
                self.absolute_support(['intron']), \
                self.absolute_support(['start_codon', 'stop_codon'])]

    def relative_support(self, gene_feature_types, abs_numb):
        # fraction of gene_feature_types that are supported by hints
        if abs_numb > 0:
            hint_numb = 0
            for type in gene_feature_types:
                hint_numb += len(self.evi_list[type])
            return hint_numb / abs_numb
        return 1

    def absolute_support(self, gene_feature_types):
        score = 0.0
        for type in gene_feature_types:
            for hint in self.evi_list[type]:
                for src in hint.keys():
                    score += self.sw[src] * hint[src]
        return score

    def get_features(self):
        return self.feature_vector
