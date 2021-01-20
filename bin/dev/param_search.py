import argparse
import itertools
import subprocess as sp
import os

species = 'Medicago_truncatula'
def main():
    args = parseCmd()
    para_level = ['0', '1']
    para_list = list(itertools.product(para_level, repeat=4))
    i = 1
    for para in para_list:
        print(para)
        out = '{}/test{}/'.format(args.out, i)
        i += 1
        os.mkdir(out)
        with open(out + 'para.txt', 'w+') as file:
            file.write('# P,E,C,M\n{}'.format(para))
        para = ','.join(para)
        cmd = 'python3 /home/lars/work/combiner/bin/test_species.py ' \
            + '--data {} --out {} --sw {} --species {}'.format(args.data, out, para, species)
        sp.call(cmd, shell=True)

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
    return parser.parse_args()

if __name__ == '__main__':
    main()
