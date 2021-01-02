import argparse
import subprocess as sp
import os
augustus_script = os.path.abspath(os.path.dirname(__file__)) + "/../../intron-ppx/Augustus/scripts/augustus2browser.pl"

def main():
    args = parseCmd()
    if args.mode == 'rnaseq':
        cfg = {'intron' : {1        0.34       M       1  1e+100      RM       1       1       E       1     1e6       W       1       1       P       1 1'

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--gtf', type=str,
        help='')
    parser.add_argument('--out', type=str,
        help='')
    parser.add_argument('--name', type=str,
        help='')
    parser.add_argument('--color', type=str,
        help='')
    parser.add_argument('--mode', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
