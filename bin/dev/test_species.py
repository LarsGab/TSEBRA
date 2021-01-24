#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# test_species.py: automated testing for multiple species and testsettings for BRAKER1 and BRAKER2
# ==============================================================
import argparse
import os
import subprocess as sp
import numpy as np
import multiprocessing as mp

class EvaluationError(Exception):
    pass

sw = {}
combiner_bin = os.path.dirname(os.path.realpath(__file__))

header = ['gene Sn', 'gene Sp', 'trans Sn', 'trans Sp', 'cds Sn', 'cds Sp', 'tx per gene']
test_order = ['Arabidopsis_thaliana_species_excluded', 'Arabidopsis_thaliana_family_excluded',\
            'Arabidopsis_thaliana_order_excluded', 'Caenorhabditis_elegans_species_excluded', \
            'Caenorhabditis_elegans_family_excluded', 'Caenorhabditis_elegans_order_excluded', \
            'Danio_rerio_order_excluded', 'Drosophila_melanogaster_species_excluded', \
            'Drosophila_melanogaster_family_excluded', 'Drosophila_melanogaster_order_excluded', \
            'Medicago_truncatula_order_excluded', 'Solanum_lycopersicum_order_excluded']
full_eval = {}
summary_eval = {}

def main():
    global combiner_bin, sw
    args = parseCmd()
    if args.combiner:
        combiner_bin = args.combiner
    if args.sw:
        sw = args.sw
    braker2_level = ['species_excluded', 'family_excluded', 'order_excluded']

    if args.species:
        species_list = args.species.split(',')
    else:
        with open(args.data + '/species.tab', 'r') as file:
            species_list = file.read().split('\n')
    species_list = [s for s in species_list if s]

    # [...,[braker_list, hintfile_list, out, test id, species path],..]
    param = []
    for species in species_list:
        species_path = "{}/{}".format(args.data, species)
        for level in braker2_level:
            braker2 = "{}/braker2/{}/braker_fixed.gtf".format(species_path, level)
            if os.path.exists(braker2):
                test_id = "{}_{}".format(species, level)
                braker1 = "{}/braker1/braker_{}.gtf".format(species_path, level)
                evidence = "{}/braker1/hintsfile_{}.gff".format(species_path, level) \
                        + ",{}/braker2/{}/hintsfile.gff".format(species_path, level)
                out = "{}/{}_{}".format(args.out, species, level)
                param.append([braker1 + ',' + braker2, evidence, out, test_id, species_path])
    print(param)

    #for p in param:
        #job(p)
    '''
    job_results = []
    pool = mp.Pool(mp.cpu_count())
    for p in param:
        pool.apply_async(job, (p,))
    pool.close()
    pool.join()
    '''
    for p in param:
        evaluation(p)
        
    write_full_eval()
    write_summary_eval()

def job(para):
    combine(para[0], para[1], para[2] + ".gtf")

def evaluation(para):
    global full_eval, summary_eval
    gtf2ucsc(para[2] + ".gtf", para[2] + "_ucsc.gtf", para[3])
    accuracies = eval("{}/anno/".format(para[4]), para[2] + ".gtf")
    tx_gene = tx_per_gene(para[2] + ".gtf")
    txt = [a[0] for a in accuracies]
    txt.append(tx_gene[2][0])
    print(header)
    print(txt)
    if not header == txt:
        raise EvaluationError('Accuracy assessment output for {}'.format(para[3]))
    txt = [a[1] for a in accuracies]
    txt.append(str(round(float(tx_gene[2][1]), 2)))
    full = [para[3]] + txt
    txt = list(map(float, txt))
    summary = [para[3], sum(txt[2:-1])/4, txt[-1]]
    if para[3] in full_eval.keys() or para[3] in summary_eval.keys():
        raise EvaluationError('{} already in evaluation dictionarys.'.format(para[3]))
    full_eval.update({para[3] : full})
    summary_eval.update({para[3] : summary})

def combine(braker, evidence, out):
    # run the combiner
    cmd = "{}/../prevco.py --gtf {} --hintfiles {} --out {} --sw {} -q -p 2".format(combiner_bin, braker, \
        evidence, out, sw)
    print(cmd)
    sp.call(cmd, shell=True)

def gtf2ucsc(gtf, out, name):
    color = list(map(str, list(np.random.choice(range(256), size=3))))
    cmd = "{}/gtf2ucsc.py --gtf {} --out {} --name {} --mode augustus --color {}".\
        format(combiner_bin, gtf, out, name , ','.join(color))
    sp.call(cmd, shell=True)

def eval(anno, prediction):
    cmd = "{}/compute_accuracies.sh {}/annot.gtf {}/pseudo.gff3 {} gene trans cds".\
        format(combiner_bin, anno, anno, prediction)
    print(cmd)
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if stderr.decode():
        raise EvaluationError(stderr.decode())
    stdout = stdout.decode()
    return [s.split('\t') for s in stdout.split('\n') if s]

def tx_per_gene(gtf):
    cmd = "{}/tx_per_gene.py --gtf {}".format(combiner_bin, gtf)
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if stderr.decode():
        raise EvaluationError(stderr.decode())
    stdout = stdout.decode()
    stdout = [s for s in stdout.split('\n') if s]
    return [s.split('\t') for s in stdout[-3:]]

def write_full_eval():
    full_eval_out = '# Mode\t{}\n'.format(header)
    for id in test_order:
        full_eval_out += '# {}\n'.format('\t'.join(full_eval[id]))
    full_eval_out += '\\begin{table}[h]\n\\centering\n\\begin{tabular}{p{4cm}||c|c||c|c||c|c||c}\n'
    full_eval_out += 'Mode&{}\\\\\n\\hline\\hline\n'.format('&'.join(header.split('\t')))
    for id in test_order:
        full_eval_out += '{}\\\\\n\\hline\n'.format('&'.join(line))
    full_eval_out += '\\end{tabular}\n\\end{table}'

    full_eval_out = full_eval_out.replace('_', ' ')
    with open(args.out + 'full_evaluation.txt', 'w+') as file:
        file.write(full_eval_out)

def write_summary_eval():
    summary_list = []
    for id in test_order:
        summary_list.append(summary_eval[id])
    summary_out = '# {}\n'.format('\t'.join(map(str,[s[0] for s in summary_list])))
    summary_out += '# {}\n'.format('\t'.join(map(str,[s[1] for s in summary_list])))
    summary_out += '&'.join(map(str,[round(s[1], 2) for s in summary_list]))
    summary_out = summary_out.replace('_', ' ')
    with open(args.out + 'summary_eval.txt', 'w+') as file:
        file.write(summary_out)
def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--data', type=str,
        help='')
    parser.add_argument('--out', type=str,
        help='')
    parser.add_argument('--sw', type=str,
        help='P,E,C,M')
    parser.add_argument('--species', type=str,
        help='')
    parser.add_argument('--combiner', type=str,
        help='')

    return parser.parse_args()

if __name__ == '__main__':
    main()
