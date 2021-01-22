#!/usr/bin/env python3
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
    args = parseCmd()
    with open(args.gtf, 'r') as file:
        gtf = file.readlines()
    chr = {}
    print(len(gtf))

    for line in gtf:
        line = line.split('\t')
        chr_id = line[0] + line[6]
        if chr_id not in chr.keys():
            chr.update({chr_id : Chr()})
        if line[2] == 'gene':
            gene_id = line[8]
            gene_id = gene_id.strip('\n')
            if not gene_id in chr[chr_id].genes.keys():
                chr[chr_id].genes.update({gene_id : [line]})
            else:
                raise FormatError('GeneID {} not unique in chromosome {}.'.format(gene_id, chr_id))
        else:
            if line[2] == 'transcript':
                transcript_id = line[8]
            else:
                transcript_id = line[8].split('transcript_id "')[1].split('";')[0]
            transcript_id = transcript_id.strip('\n')
            if not transcript_id in chr[chr_id].txs.keys():
                chr[chr_id].txs.update({transcript_id : [line]})
            else:
                chr[chr_id].txs[transcript_id].append(line)

    result = {}
    for k in chr.keys():
        for g_key in chr[k].genes.keys():
            if not g_key in result:
                result.update({g_key : chr[k].genes[g_key]})
            else:
                chr[k].genes[g_key][8] = k + '_' + g_key
                result.update({k + '_' + g_key : chr[k].genes[g_key]})
        for t_key in chr[k].txs.keys():
            if len(chr[k].txs[t_key]) == 1 and chr[k].txs[t_key][0][2] == 'transcript':
                print(t_key)
                continue
            if not t_key in result:
                result.update({t_key : chr[k].txs[t_key]})
            else:
                new_key = k + '_' + t_key
                for line in chr[k].txs[t_key]:
                    if line[2] == 'transcript':
                        line[8] = new_key
                    else:
                        temp = line[8].split('transcript_id "')
                        line[8] = '{}transcript_id "{}";{}'.format(temp[0], new_key, '";'.join(temp[1].split('";')[1:]))
                result.update({k + '_' + t_key : chr[k].txs[t_key]})
    out = []
    for value in result.values():
        out += value
    out = list(map(start2int, out))
    print(len(out))
    out = sorted(out, key=lambda o : int(o[3]))
    out = ['\t'.join(map(str,o)) for o in out]
    out = ''.join(out)
    with open(args.out, 'w+') as file:
        file.write(out)
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
