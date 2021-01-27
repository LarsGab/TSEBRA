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

        if self.type == 'stop':
            self.type = 'stop_codon'
        elif self.type == 'start':
            self.type = 'start_codon'

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
        self.start_key = {}
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
        for k in self.hints.keys():
            self.hints[k] = sorted(self.hints[k], key=lambda h:h.start)
            self.start_key.update({k : [h.start for h in self.hints[k]]})

    def hints_in_range(self, start, end, chr):
        # find all hints between start and end
        result = []
        if chr in self.hints.keys():
            start = int(start)
            end = int(end)
            index = bisect_left(self.start_key[chr], start)
            if index < len(self.hints[chr]):
                while self.hints[chr][index].start < end:
                    if self.hints[chr][index].end <= end:
                        result.append(self.hints[chr][index].hint2list())
                    index += 1
                    if index >= len(self.hints[chr]):
                        break
        return result
