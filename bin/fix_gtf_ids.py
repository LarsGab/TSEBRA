#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# Fixes an transcript and gene id error, where transcripts/genes have
# the same ID on different chromosomes or strands.
# ==============================================================
import sys
import os
import argparse

class FormatError(Exception):
    pass

class Chr:
    def __init__(self):
        self.genes = {}
        self.txs = {}

def start2int(line):
    line[3] = int(line[3])
    return line

def main():
    # replace gene/tx oldID with chr_strand_oldID
    args = parseCmd()
    result = ''
    with open(args.gtf, 'r') as file:
        for line in file.readlines():
            line = line.split('\t')
            if len(line) == 9:
                if line[2] in ['gene', 'transcript']:
                    continue
                id_prefix = line[0] + line[6]
                transcript_id = line[8].split('transcript_id "')[1].split('";')[0]
                temp = line[8].split('transcript_id "')
                line[8] = '{}transcript_id "{}_{}";{}'.format(temp[0], id_prefix, transcript_id, '";'.join(temp[1].split('";')[1:]))
                gene_id = line[8].split('gene_id "')[1].split('";')[0]
                temp = line[8].split('gene_id "')
                line[8] = '{}gene_id "{}_{}";{}'.format(temp[0], id_prefix, gene_id, '";'.join(temp[1].split('";')[1:]))
                result += '\t'.join(line)
    with open(args.out, 'w+') as file:
        file.write(result)

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
    return parser.parse_args()

if __name__ == '__main__':
    main()
