#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# features.py: Classes for creating a list of features for a transcript
# ==============================================================
class Node_features:
    def __init__(self, tx, evi, anno_pref, hint_source_weight={'P' : 0.1, 'E' : 5, 'C' : 0.5,  'M' : 1}):

        self.sw = hint_source_weight
        # evi_list[gene structure type] = [{hint.src : hint.multi}]
        self.evi_list = {'intron' : [], 'start_codon' : [], 'stop_codon': []}
        self.numb_introns = 0
        self.pref = None
        self.__init_hints__(tx, evi, anno_pref)
        # feature vector specifies the support of
        # introns, start/stop codons for a transcript
        # self.feature_vector[0] : (supported introns by evidence of tx) / (number of introns in tx)
        # self.feature_vector[1] : (supported start/stop codons by evidence of tx) / (number of \
        # start/stop codons in tx)
        # self.feature_vector[3] : sum of multiplicities of intron evidence for tx
        # self.feature_vector[4] : sum of multiplicities of start/stop codon evidence for tx
        # self.feature_vector[6] : 1 if tx is from anno_pref, 0 otherwise
        self.feature_vector = self.create_feature_vec()

    def __init_hints__(self, tx, evi, anno_pref):
        cds_len = 0
        for type in ['intron', 'start_codon', 'stop_codon']:#, 'CDS']:
            for line in tx.transcript_lines[type]:
                '''
                if line[2] == 'CDS':
                    cds_len += line[4] - line[3] + 1
                    cds_parts = evi.get_cds_parts(line[0], line[3], line[4], line[7])
                    if cds_parts:
                        cds_parts = sorted(cds_parts, key=lambda c:c[0])
                        start = cds_parts[0][0]
                        end = cds_parts[0][1]
                        for c in cds_parts:
                            if c[0] < line[3]:
                                c[0] = line[3]
                            if c[1] > line[4]:
                                c[1] = line[4]
                            f[5] += c[1] - c[0] + 1
                            if end >= c[0]:
                                end = c[1]
                            else:
                                f[2] += end - start + 1
                                start = c[0]
                                end = c[1]
                            if not start == cds_parts[-1][0]:
                                f[2] += end - start + 1
                else:
                '''
                hint = evi.get_hint(line[0], line[3], line[4], line[2], \
                    line[6])
                if hint:
                    self.evi_list[type].append(hint)
        if tx.transcript_lines['intron']:
            self.numb_introns = len(tx.transcript_lines['intron'])
        self.pref = tx.source_anno == anno_pref

    def create_feature_vec(self):
        return [self.relative_support(['intron'], self.numb_introns), \
                self.relative_support(['start_codon', 'stop_codon'], 2.0), \
                self.absolute_support(['intron']), \
                self.absolute_support(['start_codon', 'stop_codon']), \
                self.preferred_anno() \
                ]


    def relative_support(self, gene_feature_types, abs_numb):
        # fraction of gene_feature_types that are supported by hints
        if abs_numb > 0:
            hint_numb = 0
            for type in gene_feature_types:
                hint_numb += len(self.evi_list[type])
            return hint_numb / abs_numb
        return 1.0

    def absolute_support(self, gene_feature_types):
        score = 0.0
        for type in gene_feature_types:
            for hint in self.evi_list[type]:
                for src in hint.keys():
                    score += self.sw[src] * hint[src]
        return score

    def preferred_anno(self):
        if self.pref:
            return 1.0
        return 0.0

    def get_features(self):
        return self.feature_vector
