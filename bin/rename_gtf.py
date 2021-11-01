#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# Rename the transcripts and genes of a GTF file.
# ==============================================================
import argparse
import os
import csv
class FileNotFound(Exception):
    pass

def main():
    args = parseCmd()
    from genome_anno import Anno

    args = parseCmd()

    if not os.path.exists(args.gtf):
        raise FileNotFound('File not found: {}'.format(args.gtf))
    prefix = ''
    if args.prefix:
        prefix = args.prefix


    anno = Anno(args.gtf, id='')
    anno.addGtf()
    anno.norm_tx_format()
    anno.find_genes()
    tx_tab = anno.rename_tx_ids(prefix)
    anno.write_anno(args.out)
    if args.translation_tab:
        with open(args.translation_tab, 'w+') as file:
            out_writer = csv.writer(file, delimiter='\t', quotechar = "|", lineterminator = '\n')
            for line in tx_tab:
                out_writer.writerow(line)

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='Renames the transcripts and genes of a GTF file.')
    parser.add_argument('--gtf', type=str, required=True,
        help='Path to a gene prediciton file in GTF format, for example the output of TSEBRA.')
    parser.add_argument('--prefix', type=str,
        help='The string is added as a prefix to all transcript and gene IDs.')
    parser.add_argument('--translation_tab', type=str,
        help='Writes the translation table for old transcript IDs to new transcript IDs to the given file path.')
    parser.add_argument('--out', type=str, required=True,
        help='Path to the output file.')
    return parser.parse_args()

if __name__ == '__main__':
    main()
