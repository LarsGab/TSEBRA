#!/usr/bin/env python3
import os
import sys
import pytest

sys.path.append('/home/lars/work/prevco/bin')

from genome_anno import Anno
from overlap_graph import Graph, Node
from evidence import Hintfile

testDir = os.path.abspath(os.path.dirname(__file__))
example_files = testDir + '/graph/'

def compare_lists(list1, list2):
    assert len(list1) == len(list2)
    list1 = [set(l) for l in list1]
    list2 = [set(l) for l in list2]
    for element in list1:
        assert element in list2

def test_example_1():
    result = [['anno1;t1', 'anno2;t1', 'anno1;t2', 'anno2;t2'], ['anno1;t3'], ['anno2;t3']]
    anno1 = Anno(example_files + '/ex1_anno1.gtf', 'anno1')
    anno1.addGtf()
    anno1.norm_tx_format()
    anno2 = Anno(example_files + '/ex1_anno2.gtf', 'anno2')
    anno2.addGtf()
    graph = Graph([anno1, anno2])
    graph.build()
    component_list = graph.connected_components()
    compare_lists(result, component_list)

def test_example_2():
    result = [['anno2;t1', 'anno1;t1'], ['anno2;t2']]
    anno1 = Anno(example_files + '/ex2_anno1.gtf', 'anno1')
    anno1.addGtf()
    anno1.norm_tx_format()
    anno2 = Anno(example_files + '/ex2_anno2.gtf', 'anno2')
    anno2.addGtf()
    anno2.norm_tx_format()
    graph = Graph([anno1, anno2])
    graph.build()
    print(graph.print_nodes())
    component_list = graph.connected_components()
    compare_lists(result, component_list)

def test_example_3():
    result = [['anno1;t1', 'anno2;t1'], ['anno2;t2']]
    anno1 = Anno(example_files + '/ex3_anno1.gtf', 'anno1')
    anno1.addGtf()
    anno1.norm_tx_format()
    anno2 = Anno(example_files + '/ex3_anno2.gtf', 'anno2')
    anno2.addGtf()
    graph = Graph([anno1, anno2])
    graph.build()
    component_list = graph.connected_components()
    compare_lists(result, component_list)
'''
def test_edge_features():
    result = [[0.0, 0.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 1.0, 1.0], \
        [0.0, 0.5, 0.0, 0.5, 0, 0.0, 0.5, 0.0, 0.5, 0], \
        [1.0, 0.0, 0.0, 0.0, 14/31, 0.5, 0.0, 0.4, 0.1, 0.0], \
        [0.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.5, 0.5, 0.0, 0.0], \
        [0.4, 0.0, 0.5, 0.1, 0.0, 1.0, 0.0, 0.0, 0.0, 17/31], \
        [0.5, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0]]

    anno1 = Anno(example_files + '/ex_feature_anno1.gtf', 'anno1')
    anno1.addGtf()
    anno1.norm_tx_format()
    anno2 = Anno(example_files + '/ex_feature_anno2.gtf', 'anno2')
    anno2.addGtf()
    anno2.norm_tx_format()
    graph = Graph([anno1, anno2])
    graph.build()
    hints = [Hintfile(example_files + '/ex_feature_hint1.gff'), \
        Hintfile(example_files + '/ex_feature_hint2.gff')]
    graph.add_edge_features(hints)
    for e in graph.edges.values():
        for r in e.features.values():
            assert r in result
'''
