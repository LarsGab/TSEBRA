#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# evdence.py: Handles the extrinsic evidence from the hintfiles
# ==============================================================
import csv

class NotGtfFormat(Exception):
    pass

class AttributeMissing(Exception):
    pass

class Hint:
    """
        Class handling the data structures and methods for a hint
    """
    def __init__(self, line):
        """
            Create a hint from a gff line. The line has to include 'src=' as
            an attribute in the last column. Only introns, start/stop codons
            are used.

            Args:
                line (list(str)): GFF line for one hint from extrinsic evidence.
        """
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
        """
            Returns:
                line (list(str)): GFF line for the hint.
        """
        attribute = ['src=' + self.src]
        if int(self.mult) > 1:
            attribute.append('mult={}'.format(self.mult))
        if self.pri:
            attribute.append('pri={}'.format(self.pri))
        return [self.chr, self.source_program, self.type, self.start, self.end, \
            self.score, self.strand, self.phase, ';'.join(attribute)]

class Hintfile:
    """
        Class handling the data structures and methods for a hintfile
    """
    def __init__(self, path):
        """
            Args:
                path (str): Path to the hintfile.
        """
        # dictonary containing evidence
        # self.hints[chromosom_id] = [Hints()]
        self.hints = {}
        self.src = set()
        self.read_file(path)

    def read_file(self, path):
        """
            Read a gff file with intron or start/stop codon hints
            and create a dict of Hints.
        """
        #
        with open(path, 'r') as file:
            hints_csv = csv.reader(file, delimiter='\t')
            for line in hints_csv:
                if line[0][0] == '#':
                    continue
                new_hint = Hint(line)
                if not new_hint.chr in self.hints.keys():
                    self.hints.update({new_hint.chr : []})
                self.hints[new_hint.chr].append(new_hint)
                self.src.add(new_hint.src)

class Evidence:
    """
        Class handling the data structures and methods for extrinsic evidence
        from one or more hintfiles.
    """
    def __init__(self):
        # hint_keys[chr][start_end_type_strand][src] = multiplicity
        self.hint_keys = {}
        self.src = set()

    def add_hintfile(self, path_to_hintfile):
        """
            Read hintfile
        """
        # read hintfile
        hintfile = Hintfile(path_to_hintfile)
        self.src = self.src.union(hintfile.src)
        for chr in hintfile.hints.keys():
            if chr not in self.hint_keys.keys():
                self.hint_keys.update({chr : {}})
            for hint in hintfile.hints[chr]:
                new_key = '{}_{}_{}_{}'.format(hint.start, hint.end, \
                    hint.type, hint.strand)
                if not new_key in self.hint_keys[chr].keys():
                    self.hint_keys[chr].update({new_key : {}})
                if not hint.src in self.hint_keys[chr][new_key].keys():
                    self.hint_keys[chr][new_key].update({hint.src : 0})
                self.hint_keys[chr][new_key][hint.src] += int(hint.mult)

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
