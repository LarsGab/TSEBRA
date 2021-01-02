import argparse
import os
import subprocess as sp
import random

def main():
    args = parseCmd()
    flank_reg = 500
    temp_dir = os.path.abspath(os.path.dirname(__file__)) + '/temp/'
    compare_script = os.path.abspath(os.path.dirname(__file__)) \
        + '/../../program/braker/scripts/compare_intervals_exact.pl'
    pseudo = os.path.abspath(os.path.dirname(__file__)) \
        + '/../example/anno/pseudo.gff3'
    os.mkdir(temp_dir)
    sort_braker1 = temp_dir + 'sort_braker1.gtf'
    sort_braker2 = temp_dir + 'sort_braker2.gtf'
    sort_anno = temp_dir + 'sort_anno.gtf'
    sp.call('sort -k1,1 -k4,4n -k5,5n {} > {}'.format(args.braker1, sort_braker1), shell=True)
    sp.call('sort -k1,1 -k4,4n -k5,5n {} > {}'.format(args.braker2, sort_braker2), shell=True)
    sp.call('sort -k1,1 -k4,4n -k5,5n {} > {}'.format(args.anno, sort_anno), shell=True)
    out_braker1 = '{}/keys_braker1.lst'.format(args.out_dir)
    cmd = '{} --pseudo {} --f1 {} --f2 {} -f3 {} --trans --shared13 --out {}'.format(compare_script, \
            pseudo, sort_braker1, sort_braker2, sort_anno, out_braker1)
    sp.call(cmd, shell=True)
    out_braker2 = '{}/keys_braker2.lst'.format(args.out_dir)
    cmd = '{} --pseudo {} --f1 {} --f2 {} -f3 {} --trans --shared23 --out {}'.format(compare_script, \
            pseudo, sort_braker1, sort_braker2, sort_anno, out_braker2)
    print(cmd)
    sp.call(cmd, shell=True)

    os.remove(sort_braker1)
    os.remove(sort_braker2)
    os.remove(sort_anno)
    os.rmdir(temp_dir)

    key_list = []
    with open(out_braker1, 'r') as file:
        key_list += file.read().split('\n')
    with open(out_braker2, 'r') as file:
        key_list += file.read().split('\n')
    sample_keys = random.sample(key_list, 100)
    #create link
    #https://genome-euro.ucsc.edu/s/LarsGab/combiner?db=dm6&position=chr2L%3A6001-11000
    links = []
    link_template = "https://genome-euro.ucsc.edu/s/LarsGab/combiner?db=dm6&position=chr{}%3A{}-{}"
    for key in sample_keys:
        if key:
            if not key[0] == '#':
                key = [k for k in key.split(' ') if k]
                chr = key[0].split('_')[0]
                start = int(key[0].split('_')[1])
                end = int(key[-1].split('_')[2])
                start -= int((end - start)/6)
                end += int((end - start)/6)
                links.append(link_template.format(chr, start, end))
    with open(args.out_dir + '/ucsc_links.txt', 'w+') as file:
        file.write('\n'.join(links))



def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--braker2', type=str,
        help='')
    parser.add_argument('--braker1', type=str,
        help='')
    parser.add_argument('--anno', type=str,
        help='')
    parser.add_argument('--out_dir', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
