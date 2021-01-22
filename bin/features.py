#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# features.py: Classes for creating a list of features for a transcript
# ==============================================================
class Node_features:
    def __init__(self, tx, evi, anno_pref):

        # feature vector specifies the support of
        # introns, start/stop codons for a transcript
        # self.feature_vector[0] : (supported introns by evidence of tx) / (number of introns in tx)
        # self.feature_vector[1] : (supported start/stop codons by evidence of tx) / (number of \
        # start/stop codons in tx)
        # self.feature_vector[2] : sum of multiplicities of intron evidence for tx
        # self.feature_vector[3] : sum of multiplicities of start/stop codon evidence for tx
        # self.feature_vector[4] : 1 if tx is from anno_pref, 0 otherwise        
        self.feature_vector = self.__init_features__(tx, evi, anno_pref)

    def __init_features__(self, tx, evi, pref):
        f = [None]*5

        evi_list = {'intron' : [], 'start_codon' : [], 'stop_codon': []}
        cds_len = 0
        for type in ['intron', 'start_codon', 'stop_codon']:
            for line in tx.transcript_lines[type]:
                mult = evi.get_hint(line[0], line[3], line[4], line[2], \
                    line[6])
                if mult > 0:
                    evi_list[type].append(mult)


        if tx.transcript_lines['intron']:
            f[0] = len(evi_list['intron']) / len(tx.transcript_lines['intron'])
        else:
            f[0] = 1.0
        f[2] = sum(evi_list['intron'])
        f[1] = (len(evi_list['start_codon']) + len(evi_list['stop_codon'])) / 2
        f[3] = sum(evi_list['start_codon']) + sum(evi_list['stop_codon'])
        if tx.source_anno == pref:
            f[4] = 1.0
        else:
            f[4] = 0.0
        return f

    def get_features(self):
        return self.feature_vector
