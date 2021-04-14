import argparse
import itertools
import subprocess as sp
import os
import multiprocessing as mp

species = 'Arabidopsis_thaliana'
t_setting = 'family_excluded'
pref = '2'
default_sw = ['0.1', '10', '5', '1']
default_support = ['0.75', '1']
default_epsilon = ['0', '0.5', '25', '0']
label = ['P', 'E', 'C', 'M', 'intron_support', 'stasto_support', 'e_1', 'e_2', 'e_3', 'e_4']
prevco = ''
braker = ''
hints = ''
anno = ''
def main():
    global prevco, braker, hints, anno
    args = parseCmd()
    file_list = []
    prevco = args.prevco
    braker = '{}/{}/braker1/braker_fixed.gtf'.format(args.data, species)
    braker += ',{}/{}/braker2/{}/braker_fixed.gtf'.format(args.data, species, t_setting)
    hints = '{}/{}/braker1/hintsfile.gff'.format(args.data, species)
    hints += ',{}/{}/braker2/{}/hintsfile.gff'.format(args.data, species, t_setting)
    anno = '{}/{}/anno'.format(args.data, species)
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    if args.mode == 1:
        para_level_sw = ['0', '0.1', '0.5', '1', '5', '10', '20', '50']
        para_list = list(itertools.product(para_level_sw, repeat=4))
        for para in para_list:
            id = '_'.join(para)
            path = args.out + '/' + id
            os.mkdir(path)
            full_para = list(para) + list(default_support) + list(default_epsilon)
            with open(path + '/para.cfg', 'w+') as file:
                for i in range(0, len(label)):
                    file.write(label[i] + ' ' + full_para[i] + '\n')
            file_list.append(path)
    elif args.mode == 2:
        para_level_supp = [['0', '0.1', '0.25', '0.5', '0.6', '0.75', '0.8', '0.9', '1'], ['0', '0.5', '1']]
        para_list = list(itertools.product(para_level_supp[0], para_level_supp[1]))
        for para in para_list:
            id = '_'.join(para)
            path = args.out + '/' + id
            os.mkdir(path)
            full_para = list(default_sw) + list(para) + list(default_epsilon)
            with open(path + '/para.cfg', 'w+') as file:
                for i in range(0, len(label)):
                    file.write(label[i] + ' ' + full_para[i] + '\n')
            file_list.append(path)
    elif args.mode == 3:
        para_level_epsi = [['0', '0.01', '0.1', '0.2', '0.3', '0.5'], \
                            ['0', '5', '10', '25', '50'], ['0']]
        para_list = list(itertools.product(para_level_epsi[0], para_level_epsi[0], \
            para_level_epsi[1], para_level_epsi[1]))
        for para in para_list:
            id = '_'.join(para)
            path = args.out + '/' + id
            os.mkdir(path)
            full_para = default_sw + default_support + list(para)
            with open(path + '/para.cfg', 'w+') as file:
                for i in range(0, len(label)):
                    file.write(label[i] + ' ' + full_para[i] + '\n')
            file_list.append(path)

    pool = mp.Pool(mp.cpu_count())
    for f in file_list:
        pool.apply_async(job, (f,))
    pool.close()
    pool.join()




def job(p):
    cmd = '{}/bin/prevco.py -c {}/para.cfg -g {} -e {} -q -o {}/combined.gtf'.format(prevco, p, braker, hints, p)
    print(cmd)
    sp.call(cmd, shell=True)
    cmd = "{}/bin/dev/compute_accuracies.sh {}/annot.gtf {}/pseudo.gff3 {}/combined.gtf trans cds".\
        format(prevco, anno, anno, p)
    print(cmd)
    q = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = q.communicate()
    if stderr.decode():
        raise EvaluationError(stderr.decode())
    stdout = stdout.decode()
    with open('{}/full_eval.out'.format(p), 'w+') as file:
        file.write(stdout)
    accuracies = [s.split('\t') for s in stdout.split('\n') if s]
    txt = [a[1] for a in accuracies]
    txt = list(map(float, txt))
    summary = str((txt[2]*txt[3])/(txt[2]+txt[3]) + (txt[0]*txt[1])/(txt[0]+txt[1]))
    with open('{}/summary_eval.out'.format(p), 'w+') as file:
        file.write(summary)

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--data', type=str,
        help='')
    parser.add_argument('--mode', type=int,
        help='1: sw, 2: support, 3: epsilon')
    parser.add_argument('--out', type=str,
        help='')
    parser.add_argument('--prevco', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
