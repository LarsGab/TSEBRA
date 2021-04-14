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
import copy
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
quiet = False
top_para = []
pred_result = []
scores = []
eval_level = []
parameter = {'P' : 0, 'E' : 0, 'C' : 0,  'M' : 0, \
    'intron_support' : 0, 'stasto_support' : 0, \
    'e_1' : 0, 'e_2' : 0, 'e_3' : 0, 'e_4' : 0}
para_label = ['P', 'E', 'C', 'M', 'intron_support', 'stasto_support', \
    'e_1', 'e_2', 'e_3', 'e_4']
partition_graphs = []
train = False
def main():
    global partition_graphs, workdir, train, eval_level
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
        eval_level = ['trans', 'cds']
        train = True

        para_train()

    elif args.mode == 'test':
        eval_level = ['gene', 'trans', 'cds']

        if args.parameter:
            split = [p for p in args.parameter.split('_') if p]
            print(split)
            parameter = list(map(float, split))
        else:
            parameter = [1.0,1.0,0.0,0.0,0.675,0.5,0.0,0.0,10.0,0.0,1.0]

        para_level_tab = []
        for i in range(0, len(para_label)):
            para_level_tab.append([para_label[i], parameter[i]])
        write_tab(para_level_tab, '{}/para_level.out'.format(workdir))

        eval = job(parameter)

        write_tab(eval, '{}/{}.eval'.format(workdir, args.parameter))
    else:
        raise ModeError('Mode must be either "train" or "test"')
    #grouped_para_train(graph, evi, anno)

def para_train():
    global scores, top_para
    para_level_sw = [[0.0, 1.0, 10.0]] * 2
    para_level_supp = [[0.0, 0.5, 1.0], [0.0, 0.5, 1.0]]
    para_level_epsi = [[0.0, 0.5, 1.0], [0.0, 0.5, 1.0], [0.0, 10.0, 20.0], [0.0, 10.0, 20.0]]
    para_level = para_level_sw + [[0.0]]*2 + para_level_supp + para_level_epsi
    maxiter = 2
    para_level = [[p[0]] for p in para_level]
    for i in range(0, maxiter):
        prefix = 'it{}'.format(i+1)
        para_level_tab = []
        for i in range(0, len(para_label)):
            para_level_tab.append([para_label[i], para_level[i]])
        write_tab(para_level_tab, '{}/{}_para_level.out'.format(workdir, prefix))

        para_list = list(itertools.product(*para_level))
        print(len(para_list))
        multi_jobs(para_list, prefix)
        print(top_para)
        best_para = list(map(float, top_para[-1][0].split('_')))
        print(best_para)

        for j in [4, 5, 6, 7]:
            para_level[j] = get_new_para_level_0_1(para_level[j], para_level[j].index(best_para[j]))
        for j in [0, 1, 8, 9]:
            para_level[j] = get_new_para_level(para_level[j], para_level[j].index(best_para[j]))
        print(para_level)
        scores = []
        top_para = []

def write_tab(tab, out_path):
    with open(out_path, 'w+') as file:
        table = csv.writer(file, delimiter='\t')
        for line in tab:
            table.writerow(line)

def get_new_para_level_0_1(old_level, best):
    if old_level[best] in [0,1]:
        new = set([old_level[best]])
        sign = int((-1) ** old_level[best])
        diff = 1/4 * (old_level[best + sign] - old_level[best] )
        for i in [1,3]:
            new.add(i * diff + old_level[best])
        return sorted(list(new))
    else:
        return get_new_para_level(old_level, best)
def get_new_para_level(old_level, best):
    factor = 2
    new = set([old_level[best]])
    i = 0
    if old_level[0] == 0 and best == 0:
        i = 1
    old_level = [old_level[i]/factor] + old_level \
        + [old_level[-1]*factor]
    old_level = sorted(old_level)
    best += 1
    for i in [-1, 0]:
        j = best + i
        new.add((old_level[j] + old_level[j+1])/2)
    return sorted(list(new))

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

    graph = Graph(gene_set, verbose=0)
    graph.build()
    return [graph, evi, gene_set, exec_dir]

def multi_jobs(p_list, prefix):
    global top_para, scores
    job_results = []
    scores = []
    '''
    for p in p_list[:55]:
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

    write_tab(scores, '{}/{}_full_scores.out'.format(workdir, prefix))
    write_tab(top_para, '{}/{}_top.out'.format(workdir, prefix))

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
    test_result = []
    for m in eval_level:
        list = []
        for r in results:
            list.append(r[m])
        score = sum_score_lst(list)
        test_result.append(['{}_F1'.format(m), score.f1()])
        test_result.append(['{}_Sn'.format(m), score.sens()])
        test_result.append(['{}_Sp'.format(m), score.spec()])
        f1_score += 1/2 * score.f1()
        sens += 1/2 * score.sens()
        spec += 1/2 * score.spec()
    if train:
        cmd = 'rm -r {}'.format(job_dir)
        sp.call(cmd, shell=True)
        return ['_'.join(map(str, p)), f1_score, sens ,spec]
    return test_result

def prediction(p, data, job_dir):
    graph = copy.deepcopy(data[0])
    part_name = data[3].split('/')[-1]

    for i in range(0, len(para_label)):
        parameter[para_label[i]] = p[i]

    graph.para = parameter
    graph.add_node_features(data[1])
    combined_prediction = graph.get_decided_graph()
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
            score.update({m : Score(0, count[m], 0)})
        return score

    for m in eval_level:
        if os.stat('{}/annot.gtf'.format(data[3])).st_size == 0:
            fp = 0
            if m in ['trans', 'gene']:
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
    tx = []
    gene = []
    with open(file_path, 'r') as file:
        gtf = csv.reader(file, delimiter='\t')
        for line in gtf:
            id = get_attribute(line[8], 'transcript_id')
            if id not in tx:
                tx.append(id)
            id = get_attribute(line[8], 'gene_id')
            if id not in gene:
                gene.append(id)
            if line[2].lower() == 'cds':
                cds += 1
    return {'gene' : len(gene), 'trans' : len(tx), 'cds' : cds}

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
    parser.add_argument('--parameter', help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
