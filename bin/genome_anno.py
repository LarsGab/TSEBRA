#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# Class structure for a genome annotation file
# ==============================================================
import os
import multiprocessing as mp
import features as ft

callback = []

class transcript:
    #data structures and methods for a transcript

    def __init__(self, id, gene_id, chr, source_anno):
        self.id = id
        self.chr = chr
        self.gene_id = gene_id
        self.transcript_lines = {}
        self.gtf = []
        self.prfl_match = []
        self.source_anno = source_anno
        self.start = -1
        self.end = -1

        #FEATURES
        #self.feature_names = ['numb_introns', 'transcript_length', 'intron_length', \
        #    'fraction_intron_leng']
        self.features = ft.Features()
        #self.numb_introns = 0
        #self.transcript_length = 0
        #self.intron_length = 0
        #self.fraction_intron_leng = -1.0

    def add_line(self, line):
        if line[2] not in self.transcript_lines.keys():
            self.transcript_lines.update({line[2] : []})
        self.transcript_lines[line[2]].append(line)

    def find_introns(self):
        if not 'intron' in self.transcript_lines.keys():
            self.transcript_lines.update({'intron' : []})
            key = ''
            if 'exon' in self.transcript_lines.keys():
                key = 'exon'
            elif 'CDS' in self.transcript_lines.keys():
                key = 'CDS'
            if key:
                exon_lst = []
                for line in self.transcript_lines[key]:
                    exon_lst.append(line)
                exon_lst = sorted(exon_lst, key=lambda e:e[0])
                if self.id =='anno1_7789_t':
                    print(exon_lst)
                for i in range(1, len(exon_lst)):
                    intron = []
                    intron += exon_lst[i][0:2]
                    intron.append('intron')
                    intron.append(exon_lst[i-1][4] + 1)
                    intron.append(exon_lst[i][3] - 1)
                    intron += exon_lst[i][5:]
                    if self.id =='anno1_7789_t':
                        print(exon_lst[i])
                        print(intron)
                    self.transcript_lines['intron'].append(intron)

    def get_gtf(self):
        #returns the annotation of the transcript in gtf
        gtf = []
        for k in self.transcript_lines.keys():

        #self.transcript_line[8] = self.id
        #gtf = ['\t'.join(list(map(str, self.transcript_line)))]
            for g in self.transcript_lines[k]:
                g[8] = "transcript_id \"{}\"; gene_id \"{}\"".format(self.id, self.gene_id)
                gtf.append('\t'.join(list(map(str, g))))
        return('\n'.join(gtf))

    def get_transcript_length(self):
        if 'transcript' in self.transcript_lines.keys():
            self.start = self.transcript_lines['transcript'][0][3]
            self.end = self.transcript_lines['transcript'][0][4]
        else:
            for k in self.transcript_lines.keys():
                for line in self.transcript_lines[k]:
                    if line[3] < self.start or self.start < 0:
                        self.start = line[3]
                    if line[4] > self.end:
                        self.end = line[4]
        return self.end - self.start + 1

    def get_intron_length(self):
        length = 0
        if not 'intron' in self.transcript_lines.keys():
            self.find_introns()
        for line in self.transcript_lines['intron']:
            length += line[4] - line[3] + 1
        return length

    def add_features(self):
        #for k in self.feature_names:
            #self.features.update({k : 0.0})
        if not 'intron' in self.transcript_lines.keys():
            self.find_introns()

        #self.features.update({'numb_introns' : \
        #    ft.Intron_numb(len(self.transcript_lines['intron']))})
        self.features.add('numb_introns', len(self.transcript_lines['intron']), 2)
        self.features.add('transcript_length', self.get_transcript_length(), 0)
        self.features.add('intron_length', self.get_intron_length(), 0)
        self.features.add('fraction_intron_leng', self.features.get_value('intron_length') \
                / self.features.get_value('transcript_length'), 1)

    def print_features(self):
        print(self.id)
        print(self.features.get_value_list())


class gene:
    #data structures and methods for a gene

    def __init__(self, id, chr, source_anno):
        self.id = id
        self.chr = chr
        self.gtf_line = []
        self.transcript = {}
        self.score = 0
        self.source_anno = source_anno

    def transcript_update(self, id):
        if not id in self.transcript.keys():
            self.transcript.update({ id : transcript(id, self.source_anno)})

    def gtf(self):
        #returns the annotation of the transcript in gtf
        if self.gtf_line:
            self.gtf_line[8] = self.id
        gtf = ['\t'.join(list(map(str, self.gtf_line)))]
        for k in self.transcript.keys():
            gtf.append(self.transcript[k].get_gtf(self.id))
        return('\n'.join(gtf))

    def add_features(self):
        for k in self.transcript.keys():
            self.transcript[k].add_features()

    def print_features(self):
        for k in self.transcript.keys():
            self.transcript[k].print_features()

class Anno:
    #data structures and methods for one genome annotation file
    def __init__(self, path, id):
        self.id = id
        self.genes = {}
        self.gene_gtf = {}
        self.transcripts = {}
        self.path = path

    def addGtf(self):
        with open (self.path, 'r') as file:
            file_lines = file.readlines()
        for line in file_lines:
            if not (line[0] == '#' or line == '\n'):
                line = line.strip('\n').split('\t')
                print(line)
                line[3] = int(line[3])
                line[4] = int(line[4])
                if line[2] == 'gene':
                    gene_id = line[8]
                    self.genes_update(gene_id)
                    if not gene_id in self.gene_gtf.keys():
                        self.gene_gtf.update({gene_id : line})
                    else:
                        print('ERROR, gene_id not unique: {}'.format(gene_id))
                elif line[2] == 'transcript':
                    gene_id = line[8].split('.')[0]
                    transcript_id = line[8]
                    self.genes_update(gene_id, transcript_id)
                    self.transcript_update(transcript_id, gene_id, line[0])
                    self.transcripts[transcript_id].add_line(line)
                else:
                    gene_id = line[8].split('gene_id "')[1].split('";')[0]
                    transcript_id = line[8].split('transcript_id "')[1].split('";')[0]
                    self.genes_update(gene_id, transcript_id)
                    self.transcript_update(transcript_id, gene_id, line[0])
                    self.transcripts[transcript_id].add_line(line)

    def genes_update(self, gene_id, transcript_id=''):
        if not gene_id in self.genes.keys():
            self.genes.update({ gene_id : []})
        if transcript_id and transcript_id not in self.genes[gene_id]:
            self.genes[gene_id].append(transcript_id)

    def transcript_update(self, t_id, g_id, chr):
        if not t_id in self.transcripts.keys():
            self.transcripts.update({ t_id : transcript(t_id, g_id, chr, self.id)})

    def print_gtf(self):
        c = 0
        for k in self.genes.keys():
            if k in self.gene_gtf.keys():
                print(self.gene_gtf[k])
            for t_id in self.genes[k]:
                print(self.transcripts[t_id].get_gtf())
            c+=1
            if c == 200:
                break
    def print_features(self):
        c = 0
        for k in self.transcripts.keys():
            self.transcripts[k].print_features()
            c+=1
            if c == 200:
                break
    def add_features(self, threads):
        global callback
        job_results = []
        pool = mp.Pool(threads)
        for k in self.transcripts.keys():
            #job(p)
            self.feature_job(k)
            #pool.apply_async(self.feature_job, (k,), callback=callback_mp)
            #r = pool.apply_async(self.feature_job, (k,), callback=callback_mp)
            #job_results.append(r)
        #for r in job_results:
            #r.wait()
        pool.close()
        pool.join()

        for c in callback:
            self.transcripts[c[1]] = c[0]
        callback = []

    def feature_job(self, k):
        self.transcripts[k].add_features()
        return [self.transcripts[k], k]

    def change_id(self, new_id):
        self.id = new_id
        for k in self.transcripts.keys():
            self.transcripts.source_anno = self.id

    def get_transcript_list(self):
        return list(self.transcripts.values())
#gathering results for multiprocessing
def callback_mp(result):
    global callback
    callback.append(res)
