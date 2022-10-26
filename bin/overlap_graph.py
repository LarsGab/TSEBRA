#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# Graph for transcripts of multiple genome annotations.
# It can detect overlapping transcripts.
# Add a feature vector to each node.
# Compare nodes with the 'decision rule'.
# ==============================================================
from features import Node_features
import numpy as np


class Edge:
    """
        Class handling an edge in the overlap graph.
    """
    def __init__(self, n1_id, n2_id):
        """
            Args:
                n1_id (str): Node ID from overlap graph
                n2_id (str): Node ID from overlap graph
        """
        self.node1 = n1_id
        self.node2 = n2_id
        self.node_to_remove = None

class Node:
    """
        Class handling a node that represents a transcript in the overlap graph.
    """
    def __init__(self, a_id, t_id):
        """
            Args:
                a_id (str): Annotation ID of the transcript from Anno object
                t_id (str): Transcript ID from Transcrpt object
        """
        self.id = '{};{}'.format(a_id, t_id)
        self.transcript_id = t_id
        self.is_in_ref_anno = 0.0
        # ID of original annotation/gene prediction
        self.anno_id = a_id
        # unique ID for a cluster of overlapping transcripts
        self.component_id = None

        # dict of edge_ids of edges that are incident
        # self.edge_to[id of incident Node] = edge_id
        self.edge_to = {}
        self.feature_vector = [None] * 4
        self.evi_support = False
        self.enforce = False
        self.gene_sets = set()

class Graph:
    """
        Overlap graph that can detect and filter overlapping transcripts.
    """
    def __init__(self, genome_anno_lst, para, keep_tx=[], verbose=0):
        """
            Args:
                genome_anno_lst (list(Anno)): List of Anno class objects
                                              containing genome annotations.
                para (dict(float)): Dictionary for parameter used for filtering of transcripts.
                verbose (int): Verbose mode if verbose >0 .
        """
        # self.nodes['anno;txid'] = Node(anno, txid)
        self.nodes = {}

        # self.edges['ei'] = Edge()
        self.edges = {}

        # self.anno[annoid] = Anno()
        self.anno = {}
        
        self.max_features = np.zeros(6)

        # list of connected graph components
        self.component_index = 0
        self.component_list = []

        # subset of all transcripts that weren't removed by the transcript comparison rule
        self.decided_graph = []

        # dict of duplicate genome annotation ids to new ids
        self.duplicates = {}

        # variables for verbose mode
        self.v = verbose
        self.f = [[],[],[],[]]
        self.ties = 0

        # parameters for decision rule
        self.para = para

        # list of transcript set names that are enforced
        self.keep_tx = keep_tx

        # init annotations, check for duplicate ids
        self.init_anno(genome_anno_lst)

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
        """
            Gets a transcript of a node.

            Args:
                key (str): ID of a node as 'anno_id;tx_id'

            Returns:
                (Transcript): Transcript class object with id = tx_id
                              from Anno() with id = anno_id
        """
        anno_id, tx_id = key.split(';')
        return self.anno[anno_id].transcripts[tx_id]

    def build(self):
        """
            Builds the overlap graph for >=1 Anno() objects.
            Each node of the graph represents a unique transcript from any annotation.
            Two nodes have an edge if their transcripts overlap.
            Two transcripts overlap if they share at least 3 adjacent protein coding nucleotides.
        """

        # tx_start_end[chr] = [tx_id, coord, id for start or end]
        # for every tx one element for start and one for end
        # this dict is used to check for overlapping transcripts
        tx_start_end = {}
        # check for duplicate txs, list of ['start_end_strand']
        unique_tx_keys = {}

        for k in self.anno.keys():
            for tx in self.anno[k].get_transcript_list():
                key = f'{tx.source_anno};{tx.id}'
                if tx.chr not in tx_start_end.keys():
                    tx_start_end.update({tx.chr : []})
                    unique_tx_keys.update({tx.chr : {}})
                unique_key = '{}_{}_{}'.format(tx.start, tx.end, tx.strand)
                if unique_key in unique_tx_keys[tx.chr].keys():
                    check = False
                    coords = tx.get_type_coords('CDS')
                    for t in unique_tx_keys[tx.chr][unique_key]:
                        if coords == t.get_type_coords('CDS'):
                            check = True
                            break
                    if check:
                        if tx.source_anno in self.keep_tx:
                            self.nodes[f'{t.source_anno};{t.id}'].enforce = True
                        self.nodes[f'{t.source_anno};{t.id}'].gene_sets.add(t.source_anno)
                        continue
                else:
                    unique_tx_keys[tx.chr].update({unique_key : []})
                unique_tx_keys[tx.chr][unique_key].append(tx)
                self.nodes.update({key : Node(tx.source_anno, \
                    tx.id)})
                self.nodes[f'{tx.source_anno};{tx.id}'].gene_sets.add(tx.source_anno)
                if tx.source_anno in self.keep_tx:
                    self.nodes[key].enforce = True
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
        """
            Check if two transcripts share at least 3 adjacent protein
            coding nucleotides on the same strand and reading frame.

            Args:
                tx1 (Transcript): Transcript class object of first transcript
                tx2 (Transcript): Transcript class object of second transcript

            Returns:
                (boolean): TRUE if they overlap and FALSE otherwise
        """
        if not tx1.strand == tx2.strand:
            return False
        tx1_coords = tx1.get_type_coords('CDS')
        tx2_coords = tx2.get_type_coords('CDS')
        for phase in ['0', '1', '2']:
            coords = []
            coords += tx1_coords[phase]
            coords += tx2_coords[phase]
            coords = sorted(coords, key=lambda c:c[0])
            for i in range(1, len(coords)):
                if coords[i-1][1] - coords[i][0] > 1:
                    return True
        return False

    def add_reference_anno_label(self, ref_anno):
        """
            Sets the value of is_in_ref_anno for each node to 1
            if the coding sequence of the corresponding transcript matches the
            coding sequence of a transcript in the reference anno

            Args:
                ref_anno (Anno): Anno() obeject of reference annotation
        """
        def get_cds_keys(tx):
            keys = [tx.chr, tx.strand] + [str(c[0]) + '_' + str(c[1]) \
                for c in tx.get_type_coords('CDS', frame=False)]
            return keys
        ref_anno_keys = []
        ref_anno_cds = []
        for tx in ref_anno.transcripts.values():
            cds_keys = get_cds_keys(tx)
            ref_anno_cds += cds_keys
            ref_anno_keys.append('_'.join(cds_keys))
        ref_anno_cds = set(ref_anno_cds)
        ref_anno_keys = set(ref_anno_keys)
        
        false_cds_keys = set([])
        correct_cds_keys = set([])
        numb_correct_tx = 0
        for n in self.nodes:
            self.nodes[n].is_in_ref_anno = 0.0
            c_keys = get_cds_keys(self.__tx_from_key__(n))
            if '_'.join(c_keys) in ref_anno_keys:
                self.nodes[n].is_in_ref_anno = 1.0
            


    def print_nodes(self):
        # prints all nodes of the graph (only used for development)
        for k in self.nodes.keys():
            print(self.nodes[k].id)
            print(self.nodes[k].transcript_id)
            print(self.nodes[k].anno_id)
            print(self.nodes[k].edge_to.keys())
            print('\n')

    def connected_components(self):
        """
            Compute all clusters of connected transcripts.
            A cluster is connected component of the graph.
            Adds component IDs to nodes.

            Returns:
                (list(list(str))): Lists of list of all node IDs of a component.
        """
        visited = []
        self.component_list = []
        self.component_index = 0
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
            self.component_index += 1
            for node in component:
                self.nodes[node].component_id = 'g_{}'.format(self.component_index)
        return self.component_list

    def add_node_features(self, evi):
        """
            Compute for all nodes the feature vector based on the evidence support by evi.

            Args:
                evi (Evidence): Evidence class object with all hints from any source.
        """
        self.max_features = np.zeros(6, float)
        all_features = []
        for key in self.nodes.keys():
            tx = self.__tx_from_key__(key)
            new_node_feature = Node_features(tx, evi, self.para)
            self.nodes[key].feature_vector = np.array(new_node_feature.get_features())
            self.max_features = np.maximum(self.nodes[key].feature_vector, 
                                          self.max_features)
            all_features.append(self.nodes[key].feature_vector)
            
                
        std = np.std(np.array(all_features)[:,2:], axis=0)
        mean = np.mean(np.array(all_features)[:,2:], axis=0)
        for key in self.nodes.keys():
            tx = self.__tx_from_key__(key)
#             print(self.nodes[key].feature_vector, mean, std)
            self.nodes[key].feature_vector[2:] -= mean
            self.nodes[key].feature_vector[2:] /= std
#             print(self.nodes[key].feature_vector, '\n\n')
            if self.nodes[key].feature_vector[0] >= self.para['intron_support'] \
                or self.nodes[key].feature_vector[1] >= self.para['stasto_support']:
                self.nodes[key].evi_support = True
            if self.nodes[key].feature_vector[0] >= 0.9999 and \
                self.nodes[key].feature_vector[1] == 1:
                self.nodes[key].evi_support = True

    def decide_edge(self, edge, iter_range = range(0,6)):
        """
            Apply transcript comparison rule to two overlapping transcripts

            Args:
                edge (Edge): edge between two transcripts

            Returns:
                (str): node ID of the transcript that is marked for removal
        """
        n1 = self.nodes[edge.node1]
        n2 = self.nodes[edge.node2]
        tx1 = self.__tx_from_key__(n1.id)
        tx2 = self.__tx_from_key__(n2.id)
#         if (not n1.evi_support) or (not n2.evi_support):
#             return None
        
        if len(tx1.transcript_lines['intron']) == 0 or \
            len(tx2.transcript_lines['intron']) == 0:
            iter_range = [1,3]
        for i in iter_range:
            diff = n1.feature_vector[i] - n2.feature_vector[i]
            #print(diff)
            if diff > self.para[f'e_{i+1}']:
#                 self.f[i].append(n2.id)
                return n2.id
            elif diff < (-1 * self.para[f'e_{i+1}']):
#                 self.f[i].append(n1.id)
                return n1.id
            
        if len(tx1.transcript_lines['intron']) == 0 and \
            len(tx2.transcript_lines['intron']) == 0:
            if tx1.start >= tx2.start and tx1.end <= tx2.end:
                return n1.id
            elif tx1.start <= tx2.start and tx1.end >= tx2.end:
                return n2.id
        
        return None

    def decide_component(self, component):
        """
            Applies transcript comparison rule to all transcripts of one component
            and returns the node IDs of all transcripts that are not removed by
            a comparison.

            Args:
                component (list(str)): List of node IDs

            Returns:
                (list(str)): Filtered subset of component list.
        """
        # return all ids of vertices of a graph component, that weren't excluded by the decision rule
        result = component.copy()
        for node_id in component:
            for e_id in self.nodes[node_id].edge_to.values():
                node_to_remove = self.edges[e_id].node_to_remove
                if node_to_remove:
                    if node_to_remove in result and \
                        not self.nodes[node_to_remove].enforce:
                        result.remove(node_to_remove)
            if node_id in result and not self.nodes[node_id].evi_support and \
                not self.nodes[node_id].enforce:
                result.remove(node_id)
        new_components = [[]]
        visited = []
        for k, n_id in enumerate(result):
            if n_id not in visited:
                if k > 0:
                    self.component_index += 1
                not_visited = [n_id]
                while not_visited:
                    n2_id = not_visited.pop()
                    visited.append(n2_id)
                    new_components[-1].append(n2_id)
                    not_visited += [n for n in self.nodes[n2_id].edge_to \
                        if n in result and n not in not_visited + visited]
                    if k > 0:
                        self.nodes[n2_id].component_id = f'g_{self.component_index}'
        return result

    def decide_graph(self):
        """
            Create list of connected components of the graph and apply the
            transcript comparison rule to all components.
        """
        for key in self.edges.keys():
            self.edges[key].node_to_remove = self.decide_edge(self.edges[key])
        self.decided_graph = []
        if not self.component_list:
            self.connected_components()
        for component in self.component_list:
            if len(component) > 1:
                self.decided_graph += self.decide_component(component)
            elif self.nodes[component[0]].evi_support \
                or self.nodes[component[0]].enforce:
                self.decided_graph += component

    def get_decided_graph(self):
        """
            Filter graph with the transcript comparison rule.
            Then, remove all transcripts with low evidence support and
            compute the subset of transcripts that are included in the
            combined gene prediciton.

            Returns:
                (dict(list(list(str))): Dictionary with transcript IDs and new
                gene IDs of all transcripts included in the combined gene prediciton
                for all input annotations

        """
        if not self.decided_graph:
            self.decide_graph()
        # result[anno_id] = [[tx_ids, new_gene_id]]
        result = {}
        for key in self.anno.keys():
            result.update({key : []})
        for node in self.decided_graph:
            # if self.nodes[node].evi_support or self.nodes[node].enforce:
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

        return result
