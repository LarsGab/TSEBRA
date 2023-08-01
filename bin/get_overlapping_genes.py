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
parameter = {'intron_support' : 0, 'stasto_support' : 0, \
    'e_1' : 0, 'e_2' : 0, 'e_3' : 0, 'e_4' : 0}
cfg_file = os.path.dirname(os.path.realpath(__file__)) + '/../config/braker3.cfg'
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
#     init(args)
    set_parameter(cfg_file)
    if v > 0:
        print(gtf)
    tx_keys = []
    # read gene prediciton files
    c = 1
    keep = []
    for g in [args.geneset1, args.geneset2]:
        tx_keys.append([])
        if not quiet:
            sys.stderr.write(f'### READING GENE PREDICTION: [{g}]\n')
        anno.append(Anno(g, f'anno{c}'))
        anno[-1].addGtf()
        anno[-1].norm_tx_format()
        keep.append(f'anno{c}')
        for tx in anno[-1].transcripts.values():
            cds = tx.get_type_coords('CDS', False)
            key = ['_'.join(list(map(str, c_1))) for c_1 in cds]
            tx_keys[-1].append(key)
        c+=1
        
    

    # read hintfiles
    evi = Evidence()

    # create graph with an edge for each unique transcript
    # and an edge if two transcripts overlap
    # two transcripts overlap if they share at least 3 adjacent protein coding nucleotides
    graph = Graph(anno, para=parameter, keep_tx=keep, verbose=v)
    if not quiet:
        sys.stderr.write('### BUILD OVERLAP GRAPH\n')
    graph.build()

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
    
    out_only_g1 = []
    out_only_g2 = []
    out_overlap_g1 = []
    out_overlap_g2 = []
    
    gene_gtf = sorted(combined_anno.gene_gtf.values(), key=lambda g: (g[0],g[3],g[4]))
    for gene in gene_gtf:
        gtf_gene = [[],[]]
        current_anno_sources = set([])
#         gtf_gene.append(gene)
        for tx_id in combined_anno.genes[gene[8]]:
            n_id = f'{combined_anno.transcripts[tx_id].source_anno};{".".join(tx_id.split(".")[1:])}'
#             gtf_gene += combined_anno.transcripts[tx_id].get_gtf()
#             current_anno_sources = current_anno_sources.union(graph.nodes[n_id].gene_sets)
            cds = combined_anno.transcripts[tx_id].get_type_coords('CDS', False)
            key = ['_'.join(list(map(str, c_1))) for c_1 in cds]
        
            for i, k in enumerate(tx_keys):
                if key in k:
                    gtf_gene[i].append(gene)
                    gtf_gene[i] += combined_anno.transcripts[tx_id].get_gtf()                    
            
#         print(current_anno_sources)
#         print(gtf_gene)
        if gtf_gene[0] and gtf_gene[1]:
            print(current_anno_sources, 'A')
            out_overlap_g1 += gtf_gene[0]
            out_overlap_g2 += gtf_gene[1]
        elif gtf_gene[0]:
            out_only_g1 += gtf_gene[0]
        elif gtf_gene[1]:
            out_only_g2 += gtf_gene[1]
        else:
            print(current_anno_sources)
    
    
    for i,j in zip([out_only_g1,out_only_g2,out_overlap_g1,out_overlap_g2], 
                   [f'{args.out}_only_g1', f'{args.out}_only_g2', f'{args.out}_overlap_g1',f'{args.out}_overlap_g2']):
        with open(j, 'w+') as file:
            out_writer = csv.writer(file, delimiter='\t', quotechar = "|", lineterminator = '\n')
            for line in i:
                out_writer.writerow(line)
    

def set_parameter(cfg_file):
    """
        read parameters from the cfg file and store them in parameter.

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

def init(args):
    global gtf, hintfiles, threads, hint_source_weight, out, enforce_tx, v, quiet
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
        cfg_file = os.path.dirname(os.path.realpath(__file__)) + '/../config/braker3.cfg'
    set_parameter(cfg_file)
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
    parser = argparse.ArgumentParser(description='Input: Two gtf files; Output: 3 GTF files with overlapping/not overlapping genes.')
    parser.add_argument('-g1', '--geneset1', type=str,
        help='') 
    parser.add_argument('-g2', '--geneset2', type=str,
        help='') 
    parser.add_argument('-o', '--out', type=str, required=True,
        help='')
    parser.add_argument('-q', '--quiet', action='store_true',
        help='Quiet mode.')
    parser.add_argument('-v', '--verbose', type=int,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
