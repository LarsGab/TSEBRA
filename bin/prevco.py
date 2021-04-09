#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# PrEvCo: gene Predcition and extrinsic Evidence Combiner
# ==============================================================
import argparse
import sys
import os

class ConfigFileError(Exception):
    pass

gtf = []
anno = []
hintfiles = []
graph = None
out = ''
#pref = ''
v = 0
quiet = False
parameter = {'P' : 0, 'E' : 0, 'C' : 0,  'M' : 0, \
    'intron_support' : 0, 'stasto_support' : 0, \
    'e_1' : 0, 'e_2' : 0, 'e_3' : 0, 'e_4' : 0}#, 'e_5' : 0}

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

    if v > 0:
        print(gtf)

    # read gene prediciton files
    c = 1
    for g in gtf:
        if not quiet:
            print('### READING GTF')
        anno.append(Anno(g, 'anno{}'.format(c)))
        anno[-1].addGtf()
        anno[-1].norm_tx_format()
        c += 1

    evi = Evidence()
    for h in hintfiles:
        if not quiet:
            print('### READING EVIDENCE')
        evi.add_hintfile(h)

    # detect overlapping transcripts
    # two transcript overlap, if there is overlap in the cds
    #graph = Graph(anno, anno_pref=pref, para=parameter, verbose=v)
    graph = Graph(anno, para=parameter, verbose=v)
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

def set_parameter(cfg_file):
    global parameter
    with open(cfg_file, 'r') as file:
        for line in file.readlines():
            line = line.strip('\n')
            if not line:
                continue
            if not line[0] == '#':
                line = line.split(' ')
                if line[0] not in parameter.keys():
                    raise ConfigFileError('{} is not a valid parameter.'.format(line[0]))
                parameter[line[0]] = float(line[1])

def init(args):
    global gtf, hintfiles, threads, hint_source_weight, out, v, quiet#, pref
    if args.gtf:
        gtf = args.gtf.split(',')
    if args.hintfiles:
        hintfiles = args.hintfiles.split(',')
    if args.cfg:
        cfg_file = args.cfg
    else:
        cfg_file = os.path.dirname(os.path.realpath(__file__)) + '/../config/default.cfg'
    set_parameter(cfg_file)
    if args.out:
        out = args.out
    if args.verbose:
        v = args.verbose
    #if args.pref:
        #pref = 'anno{}'.format(args.pref)
    if args.quiet:
        quiet = True

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='PrEvCo: gene Predcition and extrinsic Evidence Combiner')
    parser.add_argument('-c', '--cfg', type=str,
        help='List of parameter settings, if not set default parameters are used.')
    parser.add_argument('-v', '--verbose', type=int,
        help='')
    #parser.add_argument('-p', '--pref', type=int, required=True,
        #help='Index (>=1) of the preferred gene prediction source file in the gene prediciton list.')
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
