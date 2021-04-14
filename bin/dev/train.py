import argparse
import itertools
import subprocess as sp
import os
import multiprocessing as mp

class EvalError(Exception):
    pass

dev = os.path.dirname(os.path.realpath(__file__))
braker = ''
hints = ''
anno = ''
pseudo = ''
label = ['P', 'E', 'C', 'M', 'intron_support', 'stasto_support', 'e_1', 'e_2', 'e_3', 'e_4', 'e_5']
pref = 2
scores = []
top10 = []
workdir = ''
modes = ['trans', 'cds']
def main():
    global braker, hints, anno, workdir, pseudo
    args = parseCmd()
    braker = args.pred
    hints = args.hints
    anno = args.anno
    pseudo = args.pseudo
    workdir = args.workdir
    para_level_sw = ['0', '0.01', '0.1', '0.5' '1', '5', '10', '25', '50', '100']
    para_level_supp = [['0','0.25', '0.5', '0.75', '1'], ['0', '0.5', '1']]
    para_level_epsi = [['0', '0.01', '0.1', '0.25', '0.5'], \
                        ['0', '5', '50'], ['0', '1']]
    leng = 0
    default = [['1','1','0','0'], ['0','0'], ['0','0','0','0','0']]
    for i in range(1, 6):
        new_para = list(itertools.product(para_level_supp[0], para_level_supp[1]))
        para_list = [default[0] + list(p)  + default[2] for p in new_para]
        print(len(para_list))
        leng += len(para_list)
        multi_jobs(para_list, 'it{}_1'.format(i))
        default[1] = top10[-1][0].split('_')[4:6]

        new_para = list(itertools.product(para_level_sw, repeat=2))
        para_list = [list(p) + ['0','0'] + default[1] + default[2] for p in new_para]
        #para_list = list(itertools.product(para_level_sw, para_level_sw, para_level_sw, para_level_sw, \
            #para_level_supp[0], para_level_supp[1], para_level_epsi[0], para_level_epsi[0], \
            #para_level_epsi[1], para_level_epsi[1], para_level_epsi[2]))
        print(len(para_list))
        leng += len(para_list)
        multi_jobs(para_list, 'it{}_2'.format(i))
        default[0] = top10[-1][0].split('_')[:4]


        new_para = list(itertools.product(para_level_epsi[0], para_level_epsi[0], \
            para_level_epsi[1], para_level_epsi[1], para_level_epsi[2]))
        para_list = [default[0] + default[1] + list(p) for p in new_para]
        print(len(para_list))
        leng += len(para_list)
        multi_jobs(para_list, 'it{}_3'.format(i))
        default[2] = top10[-1][0].split('_')[6:]

    print(leng)
def multi_jobs(p_list, name):
    global top10, scores
    print(p_list)
    job_results = []
    scores = []

    pool = mp.Pool(mp.cpu_count())
    for p in p_list:
        r = pool.apply_async(job, (p,), callback=collector)
        job_results.append(r)
    for r in job_results:
        r.wait()
    pool.close()
    pool.join()

    with open ('{}/{}_full_scores.out'.format(workdir, name), 'w+') as file:
        file.write('\n'.join(['\t'.join(map(str,s)) for s in scores]))
    with open ('{}/{}_top10.out'.format(workdir, name), 'w+') as file:
        file.write('\n'.join(['\t'.join(map(str,s)) for s in top10]))


def collector(r):
    global scores, top10
    scores.append(r)
    if len(top10) < 10:
        top10.append(r)
        top10 =sorted(top10, key=lambda t:t[1])
    elif r[1] > top10[0][1]:
        top10[0] = r
        top10 =sorted(top10, key=lambda t:t[1])

def job(para):
    cfg_file = '{}/{}.cfg'.format(workdir, '_'.join(para))
    combined_gtf = '{}/{}.gtf'.format(workdir, '_'.join(para))
    with open(cfg_file, 'w+') as file:
        for i in range(0, len(label)):
            file.write(label[i] + ' ' + para[i] + '\n')
    cmd = '{}/../prevco.py -c {} -g {} -e {} -q -p {} -o {}'.format(dev, \
        cfg_file, braker, hints, pref, combined_gtf)
    sp.call(cmd, shell=True)
    f1_score = 0
    for m in modes:
        cmd = 'python3 {}/f1_score.py --pred {} --anno {} --pseudo {} --mode={}'.format(dev, \
            combined_gtf, anno, pseudo, m)
        q = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = q.communicate()
        if stderr.decode():
            raise EvalError(stderr.decode())
        f1_score += 1/2 * float(stdout.decode().strip('\n'))
    cmd = 'rm {}'.format(cfg_file)
    sp.call(cmd, shell=True)
    cmd = 'rm {}'.format(combined_gtf)
    sp.call(cmd, shell=True)
    return ['_'.join(para), f1_score]



def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--pred', type=str,
        help='list separated by ","')
    parser.add_argument('--hints', type=str,
        help='list separated by ","')
    parser.add_argument('--anno', type=str,
        help='')
    parser.add_argument('--pseudo', type=str,
        help='')
    parser.add_argument('--workdir', type=str,
        help='')

    return parser.parse_args()

if __name__ == '__main__':
    main()
