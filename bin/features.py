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
        # self.feature_vector[3] : sum of multiplicities of intron evidence for tx
        # self.feature_vector[4] : sum of multiplicities of start/stop codon evidence for tx
        # self.feature_vector[6] : 1 if tx is from anno_pref, 0 otherwise
        self.feature_vector = self.__init_features__(tx, evi, anno_pref)

    def __init_features__(self, tx, evi, pref):
        f = [None]*7
        f[2] = 0.0
        f[5] = 0.0

        evi_list = {'intron' : [], 'start_codon' : [], 'stop_codon': []}
        cds_len = 0
        for type in ['intron', 'start_codon', 'stop_codon', 'CDS']:
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
                mult = evi.get_hint(line[0], line[3], line[4], line[2], \
                    line[6])
                if mult > 0:
                    evi_list[type].append(mult)


        if tx.transcript_lines['intron']:
            f[0] = len(evi_list['intron']) / len(tx.transcript_lines['intron'])
        else:
            f[0] = 1.0
        f[3] = sum(evi_list['intron'])
        f[1] = (len(evi_list['start_codon']) + len(evi_list['stop_codon'])) / 2
        f[4] = sum(evi_list['start_codon']) + sum(evi_list['stop_codon'])
        if tx.source_anno == pref:
            f[6] = 1.0
        else:
            f[6] = 0.0
        #f[2] /= cds_len
        #f[5] /= cds_len
        return f

    def get_features(self):
        return self.feature_vector


'''
class Edge_features:
    def __init__(self, tx1, tx2, evi):
        # tx1, tx2, evi are lists of lists of gtf lines of all introns,
        # start- and stop-codons for two overlapping transcripts and
        # relavent evidence

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
        #       |tx1 n tx2 n evi_set| / |tx1| ,
        #       |(tx1 n tx2) - evi_set| / |tx1| ,
        #       |(tx1 n evi_set) - tx2| / |tx1| ,
        #       |(tx1 - tx2) - evi_set| / |tx1| ,
        #       |evi - tx1| / |evi| ,
        #       |tx1 n tx2 n evi_set| / |tx2| ,
        #       |(tx1 n tx2) - evi_set| / |tx2| ,
        #       |(tx2 n evi_set) - tx1| / |tx2| ,
        #       |(tx2 - tx1) - evi_set| / |tx2| ,
        #       |evi - tx2| / |evi|
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
        tx1_inter_tx2 = self.tx1.intersection(self.tx2)
        tx1_inter_tx2_inter_evi = tx1_inter_tx2.intersection(self.evi_set)
        #print(self.tx1)
        #print(self.evi)
        #print(size_evi)
        if size_tx1 == 0:
            vec += [0,0,1,1]
        else:
            vec.append(len(tx1_inter_tx2_inter_evi) / size_tx1)
            new_subset = tx1_inter_tx2.difference(self.evi_set)
            vec.append(len(new_subset) / size_tx1)
            new_subset = self.tx1.intersection(self.evi_set).difference(self.tx2)
            vec.append(len(new_subset) / size_tx1)
            new_subset = self.tx1.difference(self.evi_set).difference(self.tx2)
            vec.append(len(new_subset) / size_tx1)
        if size_evi == 0:
            vec.append(0)
        else:
            new_subset = [l for l in self.evi if l not in self.tx1]
            vec.append(len(new_subset) / size_evi)
        if size_tx2 == 0:
            vec += [0,0,1,1]
        else:
            vec.append(len(tx1_inter_tx2_inter_evi) / size_tx2)
            new_subset = tx1_inter_tx2.difference(self.evi_set)
            vec.append(len(new_subset) / size_tx2)
            new_subset = self.tx2.intersection(self.evi_set).difference(self.tx1)
            vec.append(len(new_subset) / size_tx2)
            new_subset = self.tx2.difference(self.evi_set).difference(self.tx1)
            vec.append(len(new_subset) / size_tx2)
        if size_evi == 0:
            vec.append(0)
        else:
            new_subset = [l for l in self.evi if l not in self.tx2]
            vec.append(len(new_subset) / size_evi)

        return vec


class Features:
    def __init__(self):
        self.name2class = {'numb_introns' : Float_feature, \
            'transcript_length' : Float_feature, 'intron_length' : Float_feature, \
            'fraction_intron_leng' : Float_feature, 'source' : String_feature}
        self.features = {}

    def add(self, f_name, value, pref = None):
        arg = [value]
        if not pref == None:
            arg.append(pref)
        if f_name in self.features.keys():
            self.features[f_name] = self.name2class[f_name](*arg)
        elif f_name in self.feature_names():
            self.features.update({f_name : self.name2class[f_name](*arg)})
        else:
            eprint("Can't add {}. Feature not available.".format(f_name))

    def get_value(self, f_name):
        if f_name in self.features.keys():
            return self.features[f_name].value
        elif f_name in self.feature_names():
            eprint("Can't get {}. Feature has not been added yet.".format(f_name))
        else:
            eprint("Can't get {}. Feature not available.".format(f_name))

    def get_value_list(self, f_names = []):
        if not f_names:
            f_names = self.features.keys()
        result = []
        for k in f_names:
            result.append(self.features[k].value)
        return result

    def get_dist_list(self, f_names = []):
        if not f_names:
            f_names = self.features.keys()
        result = []
        for k in f_names:
            result.append(self.features[k].dist())
        return result

    def feature_names(self):
        return self.name2class.keys()


class Float_feature:
    def __init__(self, value, pref=0):
        self.value = value
        self.pref = pref

    def dist(self):
        return abs(self.value - self.pref)

class String_feature:
    def __init__(self, value, pref=''):
        self.value = source
        self.pref = pref

    def dist(self):
        if self.value == self.pref:
            return 0
        return 1
'''
