#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# PrEvCo: gene Predcition and extrinsic Evidence Combiner
# ==============================================================
import argparse
import sys
import os

gtf = []
anno = []
hintfiles = []
hints = []
graph = None
out = ''
pref = 'braker2'
v = 0
quiet = False
hint_source_weight = {'P' : 5, 'E' : 0.1, 'C' : 2,  'M' : 1}

def main():
    '''
    Overview:
    1. Read gene predicitions from .gtf files
    2. Read Evidence from .gff files
    3. Detect overlapping transcripts
    4. Create feature vector (for a feature list see features.py) for all transcripts
    5. Compare the feature vectors of all pairs of overlapping transcripts
    6. Exclude transcripts based on the 'decision rule' and 5.
    7. Create combined gene predicitions (all transcripts that weren't excluded)
    '''
    from genome_anno import Anno
    from overlap_graph import Graph
    from evidence import Evidence
    global anno, graph

    args = parseCmd()
    init(args)

    # read gene prediciton files
    c = 1
    if v > 0:
        print(gtf)

    for g in gtf:
        if not quiet:
            print('### READING GTF')
        #!!! change braker to anno
        anno.append(Anno(g, 'braker{}'.format(c)))
        anno[-1].addGtf()
        anno[-1].norm_tx_format()
        c += 1

    evi = Evidence(hint_source_weight)
    for h in hintfiles:
        if not quiet:
            print('### READING EVIDENCE')
        evi.add_hintfile(h)

    # detect overlapping transcripts
    # two transcript overlap, if there is overlap in the cds
    graph = Graph(anno)
    if not quiet:
        print('### BUILD OVERLAP GRAPH')
    graph.build()

    # add features
    if not quiet:
        print('### ADD FEATURES')
    graph.add_node_features(evi)

    # apply decision rule to exclude a set of transcripts
    if not quiet:
        print('### APPLY DECISION RULE')
    combined_prediction = graph.get_decided_graph()
    combined_gtf = ''

    if v > 0:
        print(combined_prediction.keys())
        for a in anno:
            print('Numb_tx in {}: {}'.format(a.id, len(combined_prediction[a.id])))

    # write result to output file
    for a in anno:
        combined_gtf += a.get_subset_gtf(combined_prediction[a.id])
        combined_gtf += '\n'
    combined_gtf = combined_gtf.strip('\n')
    with open(out, 'w+') as file:
        file.write(combined_gtf)

def init(args):
    global gtf, hintfiles, threads, hint_source_weight, out, v, pref, quiet
    if args.gtf:
        gtf = args.gtf.split(',')
    if args.hintfiles:
        hintfiles = args.hintfiles.split(',')
    if args.sw:
        sw = args.sw.split(',')
        i = 0
        for key in ['P', 'E', 'C', 'M']:
            hint_source_weight[key] = float(sw[i])
            i += 1
    if args.out:
        out = args.out
    if args.verbose:
        v = args.verbose
    if args.pref:
        pref = 'braker{}'.format(args.pref)
    if args.quiet:
        quiet = True

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='PrEvCo: gene Predcition and extrinsic Evidence Combiner')
    parser.add_argument('--sw', type=str,
        help='P,E,C,M')
    parser.add_argument('-v', '--verbose', type=int,
        help='')
    parser.add_argument('-p', '--pref', type=int, required=True,
        help='Index (>=1) of the preferred gene prediction source file in the gene prediciton list.')
    parser.add_argument('-g', '--gtf', type=str, required=True,
        help='List (separated by commas) of gene prediciton files in gtf .\n(gene_pred1.gtf,gene_pred2.gtf,gene_pred3.gtf)')
    parser.add_argument('-e', '--hintfiles', type=str, required=True,
        help='List (separated by commas) of files containing extrinsic evidence in gff.\n(hintsfile1.gff,hintsfile2.gtf,3.gtf)')
    parser.add_argument('-o', '--out', type=str, required=True,
        help='Outputfile for the combined gene prediciton in gtf.')
    parser.add_argument('-q', '--quiet', action='store_true',
        help='No standard output.')
    return parser.parse_args()

if __name__ == '__main__':
    main()
