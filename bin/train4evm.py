#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# parameter training for PrEvCo
# ==============================================================
import argparse
import sys
import os
import multiprocessing as mp
import subprocess as sp
import itertools
import numpy as np
import csv

from genome_anno import Anno
from overlap_graph import Graph
from evidence import Evidence

class EvalError(Exception):
    pass
class GeneSetError(Exception):
    pass

class ModeError(Exception):
    pass

class Score:
    def __init__(self, tp, fn, fp):
        self.tp = int(tp)
        self.fn = int(fn)
        self.fp = int(fp)

    def sens(self):
        if self.fn + self.tp > 0:
            return self.tp / (self.fn + self.tp)
        return 0

    def spec(self):
        if self.fp + self.tp > 0:
            return self.tp / (self.fp + self.tp)
        return 0

    def f1(self):
        if self.fn + self.tp + self.fp > 0:
            return self.tp / (1/2 * (self.fn + self.fp) + self.tp)
        return 0

bin = os.path.dirname(os.path.realpath(__file__))
anno = ''
pseudo = ''
threads = mp.cpu_count()
workdir = ''
eval_level = ['trans', 'cds']
quiet = False
top_para = []
pred_result = []
scores = []
parameter = {'P' : 0, 'E' : 0, 'C' : 0,  'M' : 0, \
    'intron_support' : 0, 'stasto_support' : 0, \
    'e_1' : 0, 'e_2' : 0, 'e_3' : 0, 'e_4' : 0, 'e_5' : 0}
para_label = ['P', 'E', 'C', 'M', 'intron_support', 'stasto_support', \
    'e_1', 'e_2', 'e_3', 'e_4', 'e_5']
partition_graphs = []
train = False
def main():
    global partition_graphs, workdir, train
    args = parseCmd()
    workdir = os.path.abspath(args.workdir)
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    partition_list = []
    with open(args.partition_list, 'r') as file:
        part = csv.reader(file, delimiter='\t')
        for p in part:
            if p:
                partition_list.append(p)

    for part in partition_list:
        partition_graphs.append(build_graph(part[3]))
    if args.mode == 'train':
        train = True
        para_train()
    elif args.mode == 'test':
        parameter = [1,10,0,0,0.75,1,0,0,10,0,0]
        l = ['Braker 1:', 'Braker 2:', 'Protein:', 'Transcript:']
        with open('{}/parameter.out'.format(workdir), 'w+') as file:
            for i in range(0, len(l)):
                file.write('{}\t{}\n'.format(l[i], parameter[i]))
        print(job(parameter))
    else:
        raise ModeError('Mode must be either "train" or "test"')
    #grouped_para_train(graph, evi, anno)

def para_train():
    para_level_sw = [[0,1,5,10,20]] * 2
    para_level_supp = [[0.0, 0.25, 0.5, 0.75, 1.0], [0.0, 0.5, 1.0]]
    para_level_epsi = [[0.0, 0.25, 0.5, 0.75, 1.0], [0.0, 0.5, 1.0], [0, 5, 10, 20], [0,5,10, 20], [0,1]]
    para_list = list(itertools.product(para_level_sw[0], para_level_sw[1], [0], [0], \
                                       para_level_supp[0], para_level_supp[1], para_level_epsi[0], \
                                       para_level_epsi[1], para_level_epsi[2], para_level_epsi[3], para_level_epsi[4]))
    print(len(para_list))
    multi_jobs(para_list)
    print(top_para[-1][0])

def random_para_train(graph, evi, anno, size):
    para_list = np.array(    [np.random.rand(size) * 50,\
                                np.random.rand(size) * 50, \
                                [0] * size, \
                                [0] * size, \
                                np.random.rand(size), \
                                np.random.randint(3, size=(size)) / 4,
                                np.random.rand(size), \
                                np.random.rand(size), \
                                np.random.rand(size) * 20, \
                                np.random.rand(size) * 20, \
                                np.random.randint(1, size=(size))])
    multi_jobs(para_list.T, graph, evi, anno,'it1')

def grouped_para_train(graph, evi, anno):
    para_level_sw = [0,1,5,10]
    para_level_supp = [[0, 0.5, 1], [0, 0.5, 1]]
    para_level_epsi = [[0, 0.1, 0.5], \
                        [0, 1, 2, 3], [0, 1]]
    leng = 0
    default = [[1,1,0,0], [0,0], [0,0,0,0,0]]
    for i in range(1, 6):
        new_para = list(itertools.product(para_level_supp[0], para_level_supp[1]))
        para_list = [default[0] + list(p)  + default[2] for p in new_para]
        print(len(para_list))
        leng += len(para_list)
        multi_jobs(para_list, graph, evi, anno,'it{}_1'.format(i))
        default[1] = list(map(float, top_para[-1][0].split('_')[4:6]))

        new_para = list(itertools.product(para_level_sw, repeat=2))
        para_list = [list(p) + [0,0] + default[1] + default[2] for p in new_para]
        #para_list = list(itertools.product(para_level_sw, para_level_sw, para_level_sw, para_level_sw, \
            #para_level_supp[0], para_level_supp[1], para_level_epsi[0], para_level_epsi[0], \
            #para_level_epsi[1], para_level_epsi[1], para_level_epsi[2]))
        print(len(para_list))
        leng += len(para_list)
        multi_jobs(para_list, graph, evi, anno,'it{}_2'.format(i))
        default[0] = list(map(float, top_para[-1][0].split('_')[:4]))


        new_para = list(itertools.product(para_level_epsi[0], para_level_epsi[0], \
            para_level_epsi[1], para_level_epsi[1], para_level_epsi[2]))
        para_list = [default[0] + default[1] + list(p) for p in new_para]
        print(len(para_list))
        leng += len(para_list)
        multi_jobs(para_list, graph, evi, anno,'it{}_3'.format(i))
        default[2] = list(map(float, top_para[-1][0].split('_')[6:]))

    print(leng)

def build_graph(exec_dir):
    os.chdir(exec_dir)
    #print(exec_dir)
    gene_set = []
    # read gene prediciton files
    c = 1
    for g in ['braker1.gtf', 'braker2.gtf']:
        #print(g)
        gene_set.append(Anno(g, 'braker{}'.format(c)))
        gene_set[-1].addGtf()
        gene_set[-1].norm_tx_format()
        c += 1

    evi = Evidence()
    for h in ['braker_pasa.gff', 'braker_protein.gff']:
        evi.add_hintfile(h)

    # detect overlapping transcripts
    # two transcript overlap, if there is overlap in the cds
    graph = Graph(gene_set, anno_pref=2, verbose=0)
    graph.build()
    return [graph, evi, gene_set, exec_dir]

def multi_jobs(p_list):
    global top_para, scores
    job_results = []
    scores = []
    '''
    for p in p_list[:5]:
        #print(job(p))
        collector(job(p))
    '''
    pool = mp.Pool(threads)
    for p in p_list:
        r = pool.apply_async(job, (p,), callback=collector)
        job_results.append(r)
    for r in job_results:
        r.wait()
    pool.close()
    pool.join()

    with open ('{}/full_scores.out'.format(workdir), 'w+') as file:
        file.write('\n'.join(['\t'.join(map(str,s)) for s in scores]))
    with open ('{}/top.out'.format(workdir), 'w+') as file:
        file.write('\n'.join(['\t'.join(map(str,s)) for s in top_para]))


def collector(r):
    global scores, top_para
    scores.append(r)
    if len(top_para) < 10:
        top_para.append(r)
        top_para =sorted(top_para, key=lambda t:t[1])
    elif r[1] > top_para[0][1]:
        top_para[0] = r
        top_para = sorted(top_para, key=lambda t:t[1])

def collector_test(r):
    global pred_result
    pred_result.append(r)

def job(p):
    job_dir = '{}/{}'.format(workdir, '_'.join(map(str, p)))
    os.mkdir(job_dir)
    results = []
    if train:
        for part in partition_graphs:
            results.append(prediction(p, part, job_dir))
    else:
        #for part in partition_graphs:
            #collector_test(prediction(p, part, job_dir))

        job_results = []
        pool = mp.Pool(mp.cpu_count())

        for part in partition_graphs:
            r = pool.apply_async(prediction, (p, part, job_dir), callback=collector_test)
            job_results.append(r)
        for r in job_results:
            r.wait()
        pool.close()
        pool.join()

        results = pred_result
    f1_score = 0
    sens = 0
    spec = 0
    for m in eval_level:
        list = []
        for r in results:
            list.append(r[m])
        score = sum_score_lst(list)
        f1_score += 1/2 * score.f1()
        sens += 1/2 * score.sens()
        spec += 1/2 * score.spec()
    if train:
        cmd = 'rm -r {}'.format(job_dir)
        sp.call(cmd, shell=True)
    return ['_'.join(map(str, p)), f1_score, sens ,spec]

def prediction(p, data, job_dir):
    part_name = data[3].split('/')[-1]
    for i in range(0, len(para_label)):
        parameter[para_label[i]] = p[i]
    data[0].para = parameter
    data[0].add_node_features(data[1])
    combined_prediction = data[0].get_decided_graph()
    combined_gtf = ''
    for gene_set in data[2]:
        combined_gtf += gene_set.get_subset_gtf(combined_prediction[gene_set.id])
        combined_gtf += '\n'
    combined_gtf = combined_gtf.strip('\n')
    combined_gtf_path = '{}/{}.gtf'.format(job_dir, part_name)
    with open(combined_gtf_path, 'w+') as file:
        file.write(combined_gtf)

    score = {}
    if not combined_gtf:
        count = count_trans_cds('{}/annot.gtf'.format(data[3]))
        for m in eval_level:
            score.update({m : Score(0,count[m],0)})
        return score

    for m in eval_level:
        if os.stat('{}/annot.gtf'.format(data[3])).st_size == 0:
            fp = 0
            if m == 'trans':
                for c in combined_prediction:
                    fp += len(c)
            elif m == 'cds':
                for line in combined_gtf.split('\n'):
                    if line.split('\t')[2].lower() == 'cds':
                        fp += 1
            score.update({m : Score(0,0,fp)})
        else:
            cmd = 'sort -k1,1 -k4,4n -k5,5n -o {} {}'.format(combined_gtf_path, combined_gtf_path)
            sp.call(cmd, shell=True)
            cmd = '{}/dev/compare_intervals_exact.pl --f1 {}/annot.gtf --f2 {} --pseudo {}/pseudo.gff3 --{}'.format(bin, \
                data[3], combined_gtf_path, data[3], m)
            #print(cmd)
            q = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
            stdout, stderr = q.communicate()
            if stderr.decode():
                sys.stderr.write('Error in predicition {} with: {}'.format(job_dir, stderr.decode()))
                score.update({m : Score(0,0,0)})
            else:
                out = [o.split('\t') for o in stdout.decode().split('\n') if o]
                #print(out)
                score.update({m : Score(out[0][1], out[0][2], out[1][2])})
    return score

def sum_score_lst(list):
    score = Score(0,0,0)
    for l in list:
        score.tp += l.tp
        score.fn += l.fn
        score.fp += l.fp
    return score

def count_trans_cds(file_path):
    cds = 0
    tx =[]
    with open(file_path, 'r') as file:
        gtf = csv.reader(file, delimiter='\t')
        for line in gtf:
            id = get_attribute(line[8], 'transcript_id')
            if id not in tx:
                tx.append(id)
            if line[2].lower() == 'cds':
                cds += 1
    return {'trans' : len(tx), 'cds' : cds}

def get_attribute(attributes, a_name):
    expression = a_name + '\s"([^";]+)'
    return re.search(expression, attributes).groups()[0]

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--partition_list', type=str,
        help='')
    parser.add_argument('--workdir', type=str,
        help='')
    parser.add_argument('--mode', type=str,
        help='test or train')

    return parser.parse_args()

if __name__ == '__main__':
    main()
