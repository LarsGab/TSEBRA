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

combiner_bin = os.path.dirname(os.path.realpath(__file__))

header = ['gene_Sn', 'gene_Sp', 'trans_Sn', 'trans_Sp', 'cds_Sn', 'cds_Sp', 'tx_per_gene']
test_order = ['Arabidopsis_thaliana_species_excluded', \
        'Arabidopsis_thaliana_family_excluded', \
        'Arabidopsis_thaliana_order_excluded', \
        'Bombus_terrestris_species_excluded', \
        'Bombus_terrestris_order_excluded', \
        'Caenorhabditis_elegans_species_excluded', \
        'Caenorhabditis_elegans_family_excluded', \
        'Caenorhabditis_elegans_order_excluded', \
        'Danio_rerio_order_excluded', \
        'Drosophila_melanogaster_species_excluded', \
        'Drosophila_melanogaster_family_excluded', \
        'Drosophila_melanogaster_order_excluded', \
        'Medicago_truncatula_order_excluded', \
        'Parasteatoda_tepidariorum_order_excluded', \
        'Populus_trichocarpa_order_excluded', \
        'Rhodnius_prolixus_order_excluded', \
        'Tetraodon_nigroviridis_order_excluded', \
        'Xenopus_tropicalis_order_excluded']

full_eval = {}
summary_eval = {}
cmd_lst_path = ''
parameter = ''
def main():
    global combiner_bin, cmd_lst_path, parameter
    args = parseCmd()
    if args.config:
        parameter = args.config
    if args.combiner:
        combiner_bin = args.combiner
    braker2_level = ['species_excluded', 'family_excluded', 'order_excluded']

    if args.species:
        species_list = args.species.split(',')
    else:
        species_list = []
        with open('{}/species.tab'.format(args.data), 'r') as file:
            for s in file.read().split('\n'):
                if s:
                    if not s[0] == '#':
                        species_list.append(s)

    cmd_lst_path = args.out + '/cmd.lst'
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    full_eval = []
    summary_eval = []
    header = ''

    # [...,[braker_list, hintfile_list, out, test id, species path],..]
    param = []
    for species in species_list:
        species_path = "{}/{}".format(args.data, species)
        for level in braker2_level:
            braker2 = "{}/braker2/{}/braker_fixed.gtf".format(species_path, level)
            if os.path.exists(braker2):
                test_id = "{}_{}".format(species, level)
                braker1 = "{}/braker1/braker_fixed.gtf".format(species_path)
                #evidence = "{}/braker1/hintsfile.gff".format(species_path) \
                        #+ ",{}/braker2/{}/hintsfile.gff".format(species_path, level)
                evidence = '{}/varus/braker_evm_hints_varus.gff'.format(species_path) \
                    + ',{}/protHint/{}/braker_evm_hints_spaln.gff'.format(species_path, level)
                out = "{}/{}_{}".format(args.out, species, level)
                param.append([braker1 + ',' + braker2, evidence, out, test_id, species_path])

    pool = mp.Pool(mp.cpu_count())
    for p in param:
        pool.apply_async(job, (p,))
    pool.close()
    pool.join()

    job_results = []
    pool = mp.Pool(mp.cpu_count())
    for p in param:
        r = pool.apply_async(evaluation, (p,), callback=collector)
        job_results.append(r)
    for r in job_results:
        r.wait()
    pool.close()
    pool.join()

    write_full_eval(args.out)
    write_summary_eval(args.out)

def job(para):
    combine(para[0], para[1], para[2] + ".gtf")

def evaluation(para):
    gtf2ucsc(para[2] + ".gtf", para[2] + "_ucsc.gtf", para[3])
    accuracies = eval("{}/anno/".format(para[4]), para[2] + ".gtf")
    tx_gene = tx_per_gene(para[2] + ".gtf")
    txt = [a[0] for a in accuracies]
    txt.append(tx_gene[2][0].strip(':'))
    if not header == txt:
        raise EvaluationError('Accuracy assessment output for {}'.format(para[3]))
    txt = [a[1] for a in accuracies]
    txt.append(str(round(float(tx_gene[2][1]), 2)))
    full = [para[3]] + txt
    txt = list(map(float, txt))
    summary = [para[3], sum(txt[2:-1])/4, txt[-1]]
    if para[3] in full_eval.keys() or para[3] in summary_eval.keys():
        raise EvaluationError('{} already in evaluation dictionarys.'.format(para[3]))
    return [para[3], full, summary]

def collector(result):
    global full_eval, summary_eval
    full_eval.update({result[0] : result[1]})
    summary_eval.update({result[0] : result[2]})

def combine(braker, evidence, out):
    # run combiner
    cmd = "{}/../prevco.py --gtf {} --hintfiles {} --out {} -c {} -q -p 2".format(combiner_bin, braker, \
        evidence, out, parameter)
    with open(cmd_lst_path, 'a+') as file:
        file.write(cmd + '\n')
    sp.call(cmd, shell=True)

def gtf2ucsc(gtf, out, name):
    color = list(map(str, list(np.random.choice(range(256), size=3))))
    cmd = "{}/gtf2ucsc.py --gtf {} --out {} --name {} --mode augustus --color {}".\
        format(combiner_bin, gtf, out, name , ','.join(color))
    with open(cmd_lst_path, 'a+') as file:
        file.write(cmd + '\n')
    sp.call(cmd, shell=True)

def eval(anno, prediction):
    cmd = "{}/compute_accuracies.sh {}/annot.gtf {}/pseudo.gff3 {} gene trans cds".\
        format(combiner_bin, anno, anno, prediction)
    with open(cmd_lst_path, 'a+') as file:
        file.write(cmd + '\n')
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if stderr.decode():
        raise EvaluationError(stderr.decode())
    stdout = stdout.decode()
    return [s.split('\t') for s in stdout.split('\n') if s]

def tx_per_gene(gtf):
    cmd = "{}/tx_per_gene.py --gtf {}".format(combiner_bin, gtf)
    with open(cmd_lst_path, 'a+') as file:
        file.write(cmd + '\n')
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if stderr.decode():
        raise EvaluationError(stderr.decode())
    stdout = stdout.decode()
    stdout = [s for s in stdout.split('\n') if s]
    return [s.split('\t') for s in stdout[-3:]]

def write_full_eval(out):
    full_eval_out = '# Mode\t{}\n'.format('\t'.join(header))
    for id in test_order:
        full_eval_out += '# {}\n'.format('\t'.join(full_eval[id]))
    full_eval_out += '\\begin{table}[h]\n\\centering\n\\begin{tabular}{p{4cm}||c|c||c|c||c|c||c}\n'
    full_eval_out += 'Mode&{}\\\\\n\\hline\\hline\n'.format('&'.join(header))
    for id in test_order:
        full_eval_out += '{}\\\\\n\\hline\n'.format('&'.join(full_eval[id]))
    full_eval_out += '\\end{tabular}\n\\end{table}'

    full_eval_out = full_eval_out.replace('_', ' ')
    with open(out + '/full_evaluation.txt', 'w+') as file:
        file.write(full_eval_out)

def write_summary_eval(out):
    summary_list = []
    for id in test_order:
        summary_list.append(summary_eval[id])
    summary_out = '# {}\n'.format('\t'.join(map(str,[s[0] for s in summary_list])))
    summary_out += '# {}\n'.format('\t'.join(map(str,[s[1] for s in summary_list])))
    summary_out += '&'.join(map(str,[round(s[1], 2) for s in summary_list]))
    summary_out = summary_out.replace('_', ' ')
    with open(out + '/summary_eval.txt', 'w+') as file:
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
    parser.add_argument('--config', type=str,
        help='config file')
    parser.add_argument('--species', type=str,
        help='')
    parser.add_argument('--combiner', type=str,
        help='')

    return parser.parse_args()

if __name__ == '__main__':
    main()
