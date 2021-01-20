#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# main script for genome combiner
# ==============================================================
import argparse
import sys
import os
import numpy as np

gtf = []
threads = 1
anno = []
hintfiles = []
hints = []
graph = None
out = ''
hint_source_weight = {'P' : 5, 'E' : 0.1, 'C' : 2,  'M' : 1}
#feature_names = ['numb_introns' ,'transcript_length', 'intron_length', \
#    'fraction_intron_leng']

def main():
    combiner_bin = os.path.dirname(os.path.realpath(__file__))
    sys.path.remove(combiner_bin)
    sys.path.append(combiner_bin + '/../')
    from combiner.bin.genome_anno import Anno
    from combiner.bin.overlap_graph import Graph
    from combiner.bin.evidence import Evidence
    global anno, graph

    args = parseCmd()
    init(args)

    #read gtf files
    c = 1
    print(gtf)
    for g in gtf:
        print('READING GTF')
        anno.append(Anno(g, 'braker{}'.format(c)))
        anno[-1].addGtf()
        anno[-1].norm_tx_format()
        c += 1

    evi = Evidence(hint_source_weight)
    for h in hintfiles:
        print('READING HINTS')
        evi.add_hintfile(h)

    #compute features
        #anno[-1].add_transcript_features(threads)
        #anno[-1].print_gtf()
        #anno[-1].print_features()

    #detect overlaps
    graph = Graph(anno)
    graph.build()
    graph.add_node_features(evi)
    combined_prediction = graph.get_decided_graph()
    combined_gtf = ''
    print(combined_prediction.keys())
    for a in anno:
        print('Numb_tx in {}: {}'.format(a.id, len(combined_prediction[a.id])))
        combined_gtf += a.get_subset_gtf(combined_prediction[a.id])
        combined_gtf += '\n'
    combined_gtf = combined_gtf.strip('\n')
    with open(out, 'w+') as file:
        file.write(combined_gtf)



    '''
    components = graph.connected_components()

    #combine gtfs
    combined_components = []
    for c in components:
        combined_components.append(decide_component(c, ))
    print(combined_components)
    #write result
    '''

def decide_component(component):
    subtree_ranks = {}
    subtree_features = []
    print(component)
    print(graph.component_to_no_edge_subtree(component))
    for c in graph.component_to_no_edge_subtree(component):
        feature_dist = np.zeros(len(feature_names))
        for node in c.split(';'):
            if node:
                tx = graph.get_transcript_from_node(node)
                feature_dist += np.array(tx.feature_dist_list(feature_names))
        feature_dist = feature_dist / len(c.split(';'))
        subtree_features.append([c, feature_dist])
        subtree_ranks.update({c :[]})
    for i in range(0, len(feature_names)):
        subtree_features = sorted(subtree_features, key=lambda f:f[1][i])
        for j in range(0, len(subtree_features)):
            if not (j > 0 and subtree_features[j][1][i] == subtree_features[j-1][1][i]):
                rank = j
            subtree_ranks[subtree_features[j][0]].append(rank)
    result = []
    for k in subtree_ranks.keys():
        result.append([k, np.array(subtree_ranks[k]).mean()])
    result = sorted(result, key=lambda r:r[1])
    print(subtree_ranks)
    print(subtree_features)
    print(result)
    print('\n')
    return result[0][0]

def init(args):
    global gtf, hintfiles, threads, hint_source_weight, out
    if args.gtf:
        gtf = args.gtf.split(',')
    if args.hintfiles:
        hintfiles = args.hintfiles.split(',')
    if args.threads:
        threads = args.threads
    if args.sw:
        sw = args.sw.split(',')
        i = 0
        for key in ['P', 'E', 'C', 'M']:
            hint_source_weight[key] = float(sw[i])
            i += 1
        print(hint_source_weight)
    if args.out:
        out = args.out
def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    #os.path.abspath(os.path.dirname(__file__)))


    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--gtf', type=str,
        help='')
    parser.add_argument('--hintfiles', type=str,
        help='')
    parser.add_argument('--threads', type=int,
        help='')
    parser.add_argument('--sw', type=str,
        help='P,E,C,M')
    parser.add_argument('--out', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
