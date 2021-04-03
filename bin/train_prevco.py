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
class EvalError(Exception):
    pass

bin = os.path.dirname(os.path.realpath(__file__))
gtf = []
ref_anno = ''
pseudo = ''
threads = mp.cpu_count()
hintfiles = []
workdir = ''
pref = ''
eval_level = ['trans', 'cds']
quiet = False
top_para = []
scores = []
parameter = {'P' : 0, 'E' : 0, 'C' : 0,  'M' : 0, \
    'intron_support' : 0, 'stasto_support' : 0, \
    'e_1' : 0, 'e_2' : 0, 'e_3' : 0, 'e_4' : 0, 'e_5' : 0}
para_label = ['P', 'E', 'C', 'M', 'intron_support', 'stasto_support', \
    'e_1', 'e_2', 'e_3', 'e_4', 'e_5']
def main():
    from genome_anno import Anno
    from overlap_graph import Graph
    from evidence import Evidence
    anno = []
    args = parseCmd()
    init(args)

    # read gene prediciton files
    c = 1
    for g in gtf:
        anno.append(Anno(g, 'anno{}'.format(c)))
        anno[-1].addGtf()
        anno[-1].norm_tx_format()
        c += 1

    evi = Evidence()
    for h in hintfiles:
        evi.add_hintfile(h)

    # detect overlapping transcripts
    # two transcript overlap, if there is overlap in the cds
    graph = Graph(anno, anno_pref=pref, verbose=0)
    graph.build()
    no_epsi_train(graph, evi, anno)
    #grouped_para_train(graph, evi, anno)



def no_epsi_train(graph, evi, anno):
    para_level_sw = [[0, 1], list(range(0, 20))]
    para_level_supp = [[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], [0.0, 0.5, 1.0]]
    para_list = list(itertools.product(para_level_sw[0], para_level_sw[1], [0], [0], \
        para_level_supp[0], para_level_supp[1], [0], [0], [0], [0], [0]))
    print(len(para_list))
    multi_jobs(para_list, graph, evi, anno,'it1')
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



def multi_jobs(p_list, graph, evi, anno, name):
    global top_para, scores
    job_results = []
    scores = []
    '''
    for p in p_list[:2]:
        print(p)
        collector(job(p,graph, evi, anno))
    '''
    pool = mp.Pool(threads)
    for p in p_list:
        r = pool.apply_async(job, (p,graph, evi, anno), callback=collector)
        job_results.append(r)
    for r in job_results:
        r.wait()
    pool.close()
    pool.join()

    with open ('{}/{}_full_scores.out'.format(workdir, name), 'w+') as file:
        file.write('\n'.join(['\t'.join(map(str,s)) for s in scores]))
    with open ('{}/{}_top.out'.format(workdir, name), 'w+') as file:
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

def job(p, graph, evi, anno):
    for i in range(0, len(para_label)):
        parameter[para_label[i]] = p[i]
    graph.para = parameter
    graph.add_node_features(evi)
    combined_prediction = graph.get_decided_graph()
    combined_gtf = ''
    for a in anno:
        combined_gtf += a.get_subset_gtf(combined_prediction[a.id])
        combined_gtf += '\n'
    combined_gtf = combined_gtf.strip('\n')
    combined_gtf_path = '{}/{}.gtf'.format(workdir, '_'.join(map(str, p)))
    with open(combined_gtf_path, 'w+') as file:
        file.write(combined_gtf)
    f1_score = 0
    for level in eval_level:
        cmd = 'python3 {}/dev/f1_score.py --pred {} --anno {} --pseudo {} --mode={}'.format(bin, \
            combined_gtf_path, ref_anno, pseudo, level)
        q = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = q.communicate()
        if stderr.decode():
            raise EvalError(stderr.decode())
        f1_score += 1/2 * float(stdout.decode().strip('\n'))
    del graph
    del evi
    del anno
    cmd = 'rm {}'.format(combined_gtf_path)
    sp.call(cmd, shell=True)
    return ['_'.join(map(str, p)), f1_score]

def init(args):
    global gtf, hintfiles, threads, workdir, ref_anno, pseudo, pref, quiet
    if args.gtf:
        gtf = args.gtf.split(',')
    if args.hintfiles:
        hintfiles = args.hintfiles.split(',')
    if args.anno:
        ref_anno = args.anno
    if args.pseudo:
        pseudo = args.pseudo
    if args.threads:
        threads = args.threads
    if args.workdir:
        workdir = args.workdir
        if not os.path.exists(workdir):
            os.mkdir(workdir)
    if args.pref:
        pref = 'anno{}'.format(args.pref)
    if args.quiet:
        quiet = True

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='PrEvCo: gene Predcition and extrinsic Evidence Combiner')
    parser.add_argument('-p', '--pref', type=int, required=True,
        help='Index (>=1) of the preferred gene prediction source file in the gene prediciton list.')
    parser.add_argument('-g', '--gtf', type=str, required=True,
        help='List (separated by commas) of gene prediciton files in gtf .\n(gene_pred1.gtf,gene_pred2.gtf,gene_pred3.gtf)')
    parser.add_argument('-e', '--hintfiles', type=str, required=True,
        help='List (separated by commas) of files containing extrinsic evidence in gff.\n(hintsfile1.gff,hintsfile2.gtf,3.gtf)')
    parser.add_argument('-w', '--workdir', type=str, required=True,
        help='Working directory.')
    parser.add_argument('-a', '--anno', type=str, required=True,
        help='Reference annotation with the correct gene structure.')
    parser.add_argument('-t', '--threads', type=int, required=True,
        help='Number of threads used for multiprocessing.')
    parser.add_argument('--pseudo', type=str, required=True,
        help='Pseudo genes in the annotation.')
    parser.add_argument('-q', '--quiet', action='store_true',
        help='No standard output.')
    return parser.parse_args()

if __name__ == '__main__':
    main()
