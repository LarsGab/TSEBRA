#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# evdence.py: classes to handle the hints from multiple hintfiles
# ==============================================================
from bisect import bisect_left

class NotGtfFormat(Exception):
    pass

class AttributeMissing(Exception):
    pass

class Hint:
    # data structure for one hint
    def __init__(self, line):
        # hint types currently used for the descision rule
        allowed_types = ['intron', 'start', 'stop', 'CDSpart']
        line = line.split('\t')
        if not len(line) == 9:
            raise NotGtfFormat('File not in gtf Format. Error at line: {}'.format(line))
        self.chr, self.source_program, self.type, self.start, self.end, \
            self.score, self.strand, self.phase, attribute = line
        self.start = int(self.start)
        self.end = int(self.end)

        try:
            self.src = attribute.split('src=')[1].split(';')[0]
        except IndexError:
            raise AttributeMissing('Source of Hint is missing in line {}.'.format(line))

        self.mult = ''
        if 'mult=' in attribute:
            self.mult = attribute.split('mult=')[1].split(';')[0]
        else:
            self.mult= '1'

        self.pri = ''
        if 'pri=' in attribute:
            self.pri = attribute.split('pri=')[1].split(';')[0]

        if self.type == 'stop_codon':
            self.type = 'stop'
        elif self.type == 'start_codon':
            self.type = 'start'

    def hint2list(self):
        attribute = ['src=' + self.src]
        if int(self.mult) > 1:
            attribute.append('mult={}'.format(self.mult))
        if self.pri:
            attribute.append('pri={}'.format(self.pri))
        return [self.chr, self.source_program, self.type, self.start, self.end, \
            self.score, self.strand, self.phase, ';'.join(attribute)]

class Hintfile:
    # data strucure for a gff file with hints
    def __init__(self, path):
        self.path = path

        # dictonary containing evidence
        # self.hints[chromosom_id] = [Hints()]
        self.hints = {}
        self.read_file()

    def read_file(self):
        with open(self.path, 'r') as file:
            prothint = file.read().split('\n')
        for line in prothint:
            if not line:
                continue
            if line[0] == '#':
                continue
            new_hint = Hint(line)
            if not new_hint.chr in self.hints.keys():
                self.hints.update({new_hint.chr : []})
            self.hints[new_hint.chr].append(new_hint)

class Evidence:
    # data structure for hints
    def __init__(self):
        # hint_keys[chr][start_end_type_strand][src] = multiplicity
        self.hint_keys = {}
        self.cds_parts = {}
        self.cds_start = {}
        self.cds_end = {}
        self.cds_keys_start = {}
        self.cds_keys_end = {}

    def add_hintfile(self, path_to_hintfile):
        # read hintfile
        hintfile = Hintfile(path_to_hintfile)
        for chr in hintfile.hints.keys():
            if chr not in self.hint_keys.keys():
                self.hint_keys.update({chr : {}})
                self.cds_parts.update({chr : []})
                self.cds_start.update({chr : []})
                self.cds_keys_start.update({chr : []})
                self.cds_end.update({chr : []})
                self.cds_keys_end.update({chr : []})
            for hint in hintfile.hints[chr]:
                if hint.type == 'CDSpart':
                    self.cds_parts[chr].append(hint)
                else:
                    new_key = '{}_{}_{}_{}'.format(hint.start, hint.end, \
                        hint.type, hint.strand)
                    if not new_key in self.hint_keys[chr].keys():
                        self.hint_keys[chr].update({new_key : {}})
                    if not hint.src in self.hint_keys[chr][new_key].keys():
                        self.hint_keys[chr][new_key].update({hint.src : 0})
                    self.hint_keys[chr][new_key][hint.src] += int(hint.mult)
            i = 0
            for h in self.cds_parts[chr]:
                self.cds_start[chr].append([h.start, i])
                self.cds_end[chr].append([h.end, i])
                i+=1
            self.cds_start[chr] = sorted(self.cds_start[chr], key=lambda h:h[0])
            self.cds_end[chr] = sorted(self.cds_end[chr], key=lambda h:h[0])
            self.cds_keys_start[chr] = [h[0] for h in self.cds_start[chr]]
            self.cds_keys_end[chr] = [h[0] for h in self.cds_end[chr]]

    def get_cds_parts(self, chr, start, end, phase):
        result = []
        # indices of hints already added to result
        visited = []
        if chr in self.cds_parts.keys():
            start = int(start)
            end = int(end)
            index1 = bisect_left(self.cds_keys_start[chr], start)
            if index1 < len(self.cds_keys_start[chr]):
                index2 = self.cds_start[chr][index1][1]
                while self.cds_parts[chr][index2].start <= end:
                    if self.cds_parts[chr][index2].phase == phase:
                        result.append([self.cds_parts[chr][index2].start, self.cds_parts[chr][index2].end])
                        visited.append(index2)
                    index1 += 1
                    if index1 == len(self.cds_keys_start[chr]):
                        break
                    index2 = self.cds_start[chr][index1][1]

            index1 = bisect_left(self.cds_keys_end[chr], start)
            if index1 < len(self.cds_keys_end[chr]):
                index2 = self.cds_end[chr][index1][1]
                while self.cds_parts[chr][index2].start <= end:
                    if self.cds_parts[chr][index2].phase == phase and index2 not in visited:
                        result.append([self.cds_parts[chr][index2].start, self.cds_parts[chr][index2].end])
                    index1 += 1
                    if index1 == len(self.cds_keys_end[chr]):
                        break
                    index2 = self.cds_end[chr][index1][1]
        return result

    def get_hint(self, chr, start, end, type, strand):
        if type == 'start_codon':
            type = 'start'
        elif type == 'stop_codon':
            type = 'stop'
        key = '{}_{}_{}_{}'.format(start, end, type, strand)
        if chr in self.hint_keys.keys():
            if key in self.hint_keys[chr].keys():
                return self.hint_keys[chr][key]
        return {}
