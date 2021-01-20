#!/usr/bin/env python3
import argparse
import sys
import os

def main():
    sys.path.remove('/home/lars/work/prevco/bin/dev')
    sys.path.append('/home/lars/work/prevco/bin')
    from genome_anno import Anno
    from overlap_graph import Graph
    args = parseCmd()
    gtf = [Anno(args.gtf, 'anno1')]
    gtf[0].addGtf()
    gtf[0].norm_tx_format()
    graph = Graph(gtf)
    graph.build()
    print(len(graph.anno.keys()))
    components =  graph.connected_components()
    print(len(components))
    tx_per_gene = 0
    i = 0
    gtf_out = []
    for c in components:
        i += 1
        tx_per_gene += len(c)
        for t in c:
            id = t.split(';')[1]
            gtf[0].transcripts[id].gene_id = 'g{}'.format(i)
            gtf_out.append(gtf[0].transcripts[id].get_gtf())
    print('Numb_genes:\t{}\nNumb_txs:\t{}\ntx_per_gene:\t{}'.format(\
        len(components), tx_per_gene, tx_per_gene/len(components)))
    #with open(args.out, 'w+') as file:
        #file.write('\n'.join(gtf_out))

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--gtf', type=str,
        help='')
    #parser.add_argument('--out', type=str,
        #help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
