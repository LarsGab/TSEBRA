from bisect import bisect_left

class NotGtfFormat(Exception):
    pass
class AttributeMissing(Exception):
    pass

hint_source_weight = {'P' : 5, 'E' : 0.1, 'C' : 2,  'M' : 1}

class Hint:
    def __init__(self, line):
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
        if self.src not in hint_source_weight.keys():
            print('!!!!!' +  self.src)


    def hint2list(self):
        attribute = ['src=' + self.src]
        if self.mult:
            attribute.append('mult={}'.format(self.mult))
        if self.pri:
            attribute.append('pri={}'.format(self.pri))
        return [self.chr, self.source_program, self.type, self.start, self.end, \
            self.strand, self.score, self.phase, ';'.join(attribute)]

class Hintfile:
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
                        print(len(self.hints[chr]))
                        print(index)
                        break
        return result

class Evidence:
    def __init__(self):
        # hint_keys[chr][start_end_type_strand] = multiplicity
        self.hint_keys = {}

    def add_hintfile(self, path_to_hintfile):
        hintfile = Hintfile(path_to_hintfile)
        for chr in hintfile.hints.keys():
            if chr not in self.hint_keys.keys():
                self.hint_keys.update({chr : {}})
            for hint in hintfile.hints[chr]:
                new_key = '{}_{}_{}_{}'.format(hint.start, hint.end, \
                    hint.type, hint.strand)
                if new_key in self.hint_keys[chr].keys():
                    self.hint_keys[chr][new_key] += hint_source_weight[hint.src] * int(hint.mult)
                else:
                    self.hint_keys[chr].update({new_key : hint_source_weight[hint.src] * int(hint.mult)})

    def get_hint(self, chr, start, end, type, strand):
        if type == 'start_codon':
            type = 'start'
        elif type == 'stop_codon':
            type = 'stop'
        key = '{}_{}_{}_{}'.format(start, end, type, strand)
        if chr in self.hint_keys.keys():
            if key in self.hint_keys[chr].keys():
                return self.hint_keys[chr][key]
        return 0
