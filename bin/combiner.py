#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# main script for genome combiner
# ==============================================================
import argparse
import genome_anno as ga
import overlap_graph as og

gtf = ''
threads = 1
def main():
    args = parseCmd()
    init(args)

    #read gtf files
    anno = []
    c = 1
    #print(gtf)
    for g in gtf:
        anno.append(ga.Anno(g, 'anno_{}'.format(c)))
        anno[-1].addGtf()
        c += 1

    #compute features
        anno[-1].add_features(threads)
        #anno[-1].print_gtf()
        #anno[-1].print_features()

    #detect overlaps
    graph = og.Graph(anno)
    graph.build()
    graph.connected_components()
    
    #combine gtfs

    #write result

def init(args):
    global gtf, threads
    if args.gtf:
        gtf = args.gtf.split(',')
    if args.threads:
        threads = args.threads

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--gtf', type=str,
        help='')
    parser.add_argument('--threads', type=int,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
