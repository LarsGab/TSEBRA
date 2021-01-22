#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# Graph for transcripts of multiple genome annotations
# It can detect overlapping transcripts
# Add a feature vector to each node
# Compare nodes with the 'decision rule'
# ==============================================================
from features import Node_features


class Edge:
    def __init__(self, tx1_id, tx2_id):
        self.node1 = tx1_id
        self.node2 = tx2_id
        self.node_to_remove = None

class Node:
    def __init__(self, a_id, t_id):
        self.id = '{};{}'.format(a_id, t_id)
        self.transcript_id = t_id
        self.anno_id = a_id
        self.component_id = None

        # dict of edge_ids of edges that are incident
        # self.edge_to[id of incident Node] = edge_id
        self.edge_to = {}
        self.feature_vector = [None] * 7

class Graph:
    def __init__(self, genome_anno_lst, anno_pref='braker2', verbose=0):
        # self.nodes['anno;txid'] = Node(anno, txid)
        self.nodes = {}

        # self.edges['ei'] = Edge()
        self.edges = {}

        # self.anno[annoid] = Anno()
        self.anno = {}

        # preferred annotation source
        self.anno_pref = anno_pref

        # list of connected graph components
        self.component_list = []

        # subset of all transcripts that weren't excluded by the decision rule
        self.decided_graph = []

        # dict of duplicate genome annotation ids to new ids
        self.duplicates = {}
        # init annotations, check fpr duplicate ids
        self.init_anno(genome_anno_lst)

        # variables for verbose mode
        self.v = verbose
        self.f = [[],[],[],[],[],[],[]]
        self.ties = 0

    def init_anno(self, genome_anno_lst):
        # make sure that the genome_anno ids are unique
        counter = 0
        for ga in genome_anno_lst:
            if ga.id in self.anno.keys():
                counter += 1
                new_id = "duplicate.anno.{}".format(counter)
                self.duplicates.update({new_id : ga.id})
                ga.change_id(new_id)
            self.anno.update({ga.id : ga})

    def __tx_from_key__(self, key):
        # get transcript of a node
        anno_id, tx_id = key.split(';')
        return self.anno[anno_id].transcripts[tx_id]

    def build(self):
        # build graph
        # put all transcripts of a chromosome in one list and sort it by start coordinates
        # create vertex in graph for each transcript
        # tx_start_end[chr] = [tx_id, coord, id for start or end]
        # for every tx one element for start and one for end
        tx_start_end = {}
        for k in self.anno.keys():
            for tx in self.anno[k].get_transcript_list():
                if tx.chr not in tx_start_end.keys():
                    tx_start_end.update({tx.chr : []})
                key = '{};{}'.format(tx.source_anno, \
                    tx.id)
                self.nodes.update({key : Node(tx.source_anno, \
                    tx.id)})
                tx_start_end[tx.chr].append([key, tx.start, 0])
                tx_start_end[tx.chr].append([key, tx.end, 1])

        # detect overlapping nodes
        edge_count = 0
        for chr in tx_start_end.keys():
            tx_start_end[chr] = sorted(tx_start_end[chr], key=lambda t:(t[1], t[2]))
            open_intervals = []
            for interval in tx_start_end[chr]:
                if interval[2] == 0:
                    open_intervals.append(interval[0])
                else:
                    open_intervals.remove(interval[0])
                    for match in open_intervals:
                        tx1 = self.__tx_from_key__(interval[0])
                        tx2 = self.__tx_from_key__(match)
                        if self.compare_tx_cds(tx1, tx2):
                            new_edge_key = 'e{}'.format(edge_count)
                            edge_count += 1
                            self.edges.update({new_edge_key : Edge(interval[0], match)})
                            self.nodes[interval[0]].edge_to.update({match : new_edge_key})
                            self.nodes[match].edge_to.update({interval[0] : new_edge_key})


    def compare_tx_cds(self, tx1, tx2):
        # check for overlapping cds in two txs
        coords = []
        coords += tx1.get_cds_coords()
        coords += tx2.get_cds_coords()
        coords = sorted(coords, key=lambda c:c[0])
        for i in range(1, len(coords)):
            if coords[i][0] <= coords[i-1][1]:
                return True
        return False

    def print_nodes(self):
        for k in self.nodes.keys():
            print(self.nodes[k].id)
            print(self.nodes[k].transcript_id)
            print(self.nodes[k].anno_id)
            print(self.nodes[k].edge_to.keys())
            print('\n')

    def connected_components(self):
        # returns list of connected components of the graph and adds component id to vertices
        visited = []
        self.component_list = []
        component_index = 0
        for key in list(self.nodes.keys()):
            component = [key]
            if key in visited:
                continue
            visited.append(key)
            not_visited = list(self.nodes[key].edge_to.keys())
            component += not_visited
            while not_visited:
                next_node = not_visited.pop()
                visited.append(next_node)
                new_nodes = [n for n in self.nodes[next_node].edge_to.keys() if n not in component]
                not_visited += new_nodes
                component += new_nodes
            self.component_list.append(component)
            component_index += 1
            for node in component:
                self.nodes[node].component_id = 'g_{}'.format(component_index)
        return self.component_list

    def add_node_features(self, evi):
        for key in self.nodes.keys():
            tx = self.__tx_from_key__(key)
            new_node_feature = Node_features(tx, evi, self.anno_pref)
            self.nodes[key].feature_vector = new_node_feature.get_features()

    def decide_node(self, edge):
        # decision rule
        # compare the feature vectors f of two txs
        # compare elements of f from f[0] to f[7],
        # a tx is excluded from the 'decided graph',if it has a smaller value during a comparison
        # in this case the decision makes no other comparison
        n1 = self.nodes[edge.node1]
        n2 = self.nodes[edge.node2]
        for i in range(0,5):
            if n1.feature_vector[i] > n2.feature_vector[i]:
                self.f[i].append(n2.id)
                return n2.id
            elif n1.feature_vector[i] < n2.feature_vector[i]:
                self.f[i].append(n1.id)
                return n1.id
        return None

    def decide_component(self, component):
        # return all ids of vertices of a graph component, that weren't excluded by the decision rule
        result = component.copy()
        for node_id in component:
            for e_id in self.nodes[node_id].edge_to.values():
                node_to_remove = self.edges[e_id].node_to_remove
                if node_to_remove:
                    if node_to_remove in result:
                        result.remove(node_to_remove)
        return result

    def decide_graph(self):
        # applies the decision rule to all pairs of connected vertices and
        # excludes all transcript with no evidencesupport
        # returns node.id for all nodes in final prediction
        for key in self.edges.keys():
            self.edges[key].node_to_remove = self.decide_node(self.edges[key])
            #edge.decision = self.decide_edge(edge)
        self.decided_graph = []
        if not self.component_list:
            self.connected_components()
        for component in self.component_list:
            if len(component) > 1:
                self.decided_graph += self.decide_component(component)
            else:
                self.decided_graph += component

    def get_decided_graph(self):
        # returns the result of decide graph as a dict with
        # result[anno_id] = [[tx_ids, new_gene_id]]
        # a connected component in decided graph is a gene in the output
        if not self.decided_graph:
            self.decide_graph()
        result = {}
        for key in self.anno.keys():
            result.update({key : []})
        for node in self.decided_graph:
            anno_id, tx_id = node.split(';')
            result[anno_id].append([tx_id, self.nodes[node].component_id])
        if self.v > 0:
            print('NODES: {}'.format(len(self.nodes.keys())))
            f = list(map(set, self.f))
            print('f1: {}'.format(len(f[0])))
            u = f[0]
            print('f2: {}'.format(len(f[1])))
            print('f2/f1: {}'.format(len(f[1].difference(u))))
            u = u.union(f[1])
            print('f3: {}'.format(len(f[2])))
            print('f3/f2/f1: {}'.format(len(f[2].difference(u))))
            u = u.union(f[2])
            print('f4: {}'.format(len(f[3])))
            print('f4/f3/f2/f1: {}'.format(len(f[3].difference(u))))
            u = u.union(f[3])
            print('f5: {}'.format(len(f[4])))
            print('f5/f4/f3/f2/f1: {}'.format(len(f[4].difference(u))))
            u = u.union(f[4])
            print('f6: {}'.format(len(f[5])))
            print('f6/...: {}'.format(len(f[5].difference(u))))
            u = u.union(f[5])
            print('f7: {}'.format(len(f[6])))
            print('f7/...: {}'.format(len(f[6].difference(u))))
        return result
