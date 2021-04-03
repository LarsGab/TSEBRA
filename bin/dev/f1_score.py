import argparse
import subprocess as sp
import os
class EvalError(Exception):
    pass

modes = ['gene', 'trans', 'cds']
dev = os.path.dirname(os.path.realpath(__file__))
def main():
    args = parseCmd()
    cmd = 'sort -k1,1 -k4,4n -k5,5n -o {}_sort {}'.format(args.pred, args.pred)
    sp.call(cmd, shell=True)
    cmd = '{}/compare_intervals_exact.pl --f1 {} --f2 {}_sort --pseudo {} --{}'.format(dev, args.anno, \
        args.pred, args.pseudo, args.mode)
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if stderr.decode():
        raise EvalError(stderr.decode())
    stdout = stdout.decode().split('\n')
    fn_fp = 0
    for line in stdout:
        if line:
            line = line.split('\t')
            tp = int(line[1])
            fn_fp += int(line[2])
    f1 = tp /(tp + 1/2 * fn_fp)
    cmd = 'rm {}_sort'.format(args.pred)
    sp.call(cmd, shell=True)
    print(f1)

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--anno', type=str,
        help='')
    parser.add_argument('--pseudo', type=str,
        help='')
    parser.add_argument('--pred', type=str,
        help='')
    parser.add_argument('--mode', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
