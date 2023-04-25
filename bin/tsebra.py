#!/usr/bin/env python3
# ==============================================================
# author: Lars Gabriel
#
# TSEBRA: Transcript Selector for BRAKER
# ==============================================================
import argparse
import sys
import os
import csv

class ConfigFileError(Exception):
    pass

class GeneSetMissing(Exception):
    pass

gtf = []
enforce_tx = []
anno = []
hintfiles = []
graph = None
out = ''
v = 0
quiet = False
filter_sing_exon = False
ignore_tx_phase = False
scores_tab = ''
parameter = {'intron_support' : 0, 'stasto_support' : 0, \
    'e_1' : 0, 'e_2' : 0, 'e_3' : 0, 'e_4' : 0}

def main():
    """
        Overview:

        1. Read gene predicitions from .gtf files.
        2. Read Evidence from .gff files.
        3. Detect overlapping transcripts.
        4. Create feature vector (for a list of all features see features.py)
           for all transcripts.
        5. Compare the feature vectors of all pairs of overlapping transcripts.
        6. Exclude transcripts based on the 'transcript comparison rule' and 5.
        7. Remove Transcripts with low evidence support.
        8. Create combined gene predicitions (all transcripts that weren't excluded).
    """

    from genome_anno import Anno
    from overlap_graph import Graph
    from evidence import Evidence

    global anno, graph, parameter

    args = parseCmd()
    init(args)

    if v > 0:
        print(gtf)

    # read gene prediciton files
    c = 1
    keep = []
    
    for g in gtf:
        if not quiet:
            sys.stderr.write(f'### READING GENE PREDICTION: [{g}]\n')
        anno.append(Anno(g, f'anno{c}'))
        anno[-1].addGtf()
        anno[-1].norm_tx_format()
        c += 1
    for g in enforce_tx:
        if not quiet:
            sys.stderr.write(f'### READING GENE PREDICTION: [{g}]\n')
        anno.append(Anno(g, f'anno{c}'))
        anno[-1].addGtf()
        anno[-1].norm_tx_format()
        keep.append(f'anno{c}')
        c += 1
    
    # read hintfiles
    evi = Evidence()
    for h in hintfiles:
        if not quiet:
            sys.stderr.write(f'### READING EXTRINSIC EVIDENCE: [{h}]\n')
        evi.add_hintfile(h)
    for src in evi.src:
        if src not in parameter.keys():
            sys.stderr.write(f'ConfigError: No weight for src={src}, it is set to 1\n')
            parameter.update({src : 1})

    # create graph with an edge for each unique transcript
    # and an edge if two transcripts overlap
    # two transcripts overlap if they share at least 3 adjacent protein coding nucleotides
    
    graph = Graph(anno, para=parameter, keep_tx=keep, filter_single=filter_sing_exon, ignore_phase=ignore_tx_phase, verbose=v)
    if not quiet:
        sys.stderr.write('### BUILD OVERLAP GRAPH\n')
    graph.build()

    # add features
    if not quiet:
        sys.stderr.write('### ADD FEATURES TO TRANSCRIPTS\n')
    graph.add_node_features(evi)

    # apply decision rule to exclude a set of transcripts
    if not quiet:
        sys.stderr.write('### SELECT TRANSCRIPTS\n')
    combined_prediction = graph.get_decided_graph()

    if v > 0:
        sys.stderr.write(str(combined_prediction.keys()) + '\n')
        for a in anno:
            sys.stderr.write('Numb_tx in {}: {}\n'.format(a.id, len(combined_prediction[a.id])))

    # write result to output file
    if not quiet:
        sys.stderr.write('### WRITE COMBINED GENE PREDICTION\n')
    combined_anno = Anno('', 'combined_annotation')
    for a in anno:
        txs = a.get_subset([t[0] for t in combined_prediction[a.id]])
        for id, new_gene_id in combined_prediction[a.id]:
            txs[id].set_gene_id(new_gene_id)
        combined_anno.add_transcripts(txs, a.id + '.')
    combined_anno.find_genes()
    combined_anno.write_anno(out)

    if scores_tab:
        if not quiet:
            sys.stderr.write('### WRITE TRANSCRIPT SCORES\n')
        tab_out = [['### TX_ID','intron_support', 'stasto_support', 's1', 's2', 's3', 's4']]
        for node in graph.nodes.values():
            tab_out += [[node.id] + list(node.feature_vector)]
        write_csv(scores_tab, tab_out)

    if not quiet:
        sys.stderr.write('### FINISHED\n\n')
        sys.stderr.write('### The combined gene prediciton is located at {}.\n'.format(\
            out))

def set_parameter(cfg_file):
    """
        Read parameters from the cfg file and store them in parameter.

        Args:
            cfg_file (str): Path to configuration file.
    """
    global parameter
    with open(cfg_file, 'r') as file:
        cfg = csv.reader(file, delimiter=' ')
        for line in cfg:
            if not line[0][0] == '#':
                if line[0] not in parameter.keys():
                    parameter.update({line[0] : None})
                parameter[line[0]] = float(line[1])

def write_csv(out_path, tab):
    """
        Write table to out_path.
        Args:
            (str) : path to the output file
            (list) : table  
    """
    with open(out_path, 'w+') as file:
        out_writer = csv.writer(file, delimiter='\t', quotechar = "|", lineterminator = '\n')
        for line in tab:
            out_writer.writerow(line)

def init(args):
    global gtf, hintfiles, threads, hint_source_weight, out, enforce_tx, v, scores_tab, filter_sing_exon, ignore_tx_phase, quiet
    if args.gtf:
        gtf = args.gtf.split(',')
    if args.keep_gtf:
        enforce_tx = args.keep_gtf.split(',')
    if not args.keep_gtf and not args.gtf:
        raise GeneSetMissing('At least one gene set has to be provided '\
            + 'either with --gtf or --kepp_all!')
    if args.hintfiles:
        hintfiles = args.hintfiles.split(',')
    if args.cfg:
        cfg_file = args.cfg
    else:
        cfg_file = os.path.dirname(os.path.realpath(__file__)) + '/../config/default.cfg'
    set_parameter(cfg_file)
    if args.score_tab:
        scores_tab = args.score_tab
    if args.filter_single_exon_genes:
        filter_sing_exon = args.filter_single_exon_genes
    if args.ignore_tx_phase:
        ignore_tx_phase = args.ignore_tx_phase
    if args.out:
        out = args.out
    if args.verbose:
        v = args.verbose
    if args.quiet:
        quiet = True

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='TSEBRA: Transcript Selector for BRAKER\n\n' \
        + 'TSEBRA combines gene predictions by selecing ' \
        + 'transcripts based on their extrisic evidence support.')
    parser.add_argument('-g', '--gtf', type=str,
        help='List (separated by commas) of gene prediciton files in gtf.\n' \
            + '(e.g. gene_pred1.gtf,gene_pred2.gtf,gene_pred3.gtf)')
    parser.add_argument('-k', '--keep_gtf', type=str,
        help='List (separated by commas) of gene prediciton files in gtf.\n' \
            + 'These gene sets are used the same way as other inputs, but TSEBRA '\
            + 'ensures that all transcripts from these gene sets are included in the output.')
    parser.add_argument('-e', '--hintfiles', type=str,
        help='List (separated by commas) of files containing extrinsic evidence in gff.\n' \
            + '(e.g. hintsfile1.gff,hintsfile2.gtf,3.gtf)')
    parser.add_argument('-c', '--cfg', type=str,
        help='Configuration file that sets the parameter for TSEBRA. ' \
            + 'You can find the recommended parameter at config/default.cfg.')
    parser.add_argument('--filter_single_exon_genes', action='store_true',
        help='Filter out all single-exon genes out that are not' \
            + ' supported by at least one start- or stop-codon hint.')
    parser.add_argument('--ignore_tx_phase', action='store_true',
        help='Ignore the phase of transcripts while detecting clusters ' \
            + 'of overlapping transcripts.')
    parser.add_argument('-s', '--score_tab', type=str,
        help='Prints the transcript scores as a table to the specified file.')
    parser.add_argument('-o', '--out', type=str, required=True,
        help='Outputfile for the combined gene prediciton in gtf.')
    parser.add_argument('-q', '--quiet', action='store_true',
        help='Quiet mode.')
    parser.add_argument('-v', '--verbose', type=int,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()