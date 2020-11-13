#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# main script for genome combiner
# ==============================================================
import argparse
import genome_anno as ga


gtf = ''
def main():
    args = parseCmd()
    init(args)

    #read gtf files
    anno = ga.genome(gtf, 'anno1')
    anno.addGtf()
    anno.print_gtf()

    #compute features

    #detect overaps

    #combine gtfs

    #write result

def init(args):
    global gtf
    if args.gtf:
        gtf = args.gtf

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--gtf', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
