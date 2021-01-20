#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# Data structure for a graph that detects overlapping transcripts of multiple genome annotations
# ==============================================================
from combiner.bin.features import Edge_features, Node_features
#from features import Edge_features
class Edge:
    def __init__(self, tx1_id, tx2_id):
        self.node1 = tx1_id
        self.node2 = tx2_id
        self.features = {}
        self.decision = None
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
        self.evi_support = False

class Graph:
    def __init__(self, genome_anno_lst, anno_pref='braker2'):
        # self.nodes['anno;txid'] = Node(anno, txid)
        self.nodes = {}
        # self.edges['ei'] = Edge()
        self.edges = {}
        # self.anno[annoid] = Anno()
        self.anno = {}
        self.duplicates = {}

        self.component_list = []
        self.decided_graph = []
        self.init_anno(genome_anno_lst)
        self.anno_pref = anno_pref
        self.f = [[],[],[],[],[],[],[]]
        self.ties = 0

    def init_anno(self, genome_anno_lst):
        #make sure that the genome_anno ids are unique
        counter = 0
        for ga in genome_anno_lst:
            if ga.id in self.anno.keys():
                counter += 1
                new_id = "duplicate.anno.{}".format(counter)
                self.duplicates.update({new_id : ga.id})
                ga.change_id(new_id)
            self.anno.update({ga.id : ga})

    def get_transcript_from_node(self, node_id):
        return self.anno[self.nodes[node_id].anno_id].transcripts[self.nodes[node_id].transcript_id]

    def __tx_from_key__(self, key):
        #print(key)
        anno_id, tx_id = key.split(';')
        return self.anno[anno_id].transcripts[tx_id]

    def build(self):
        # build graph
        # put all transcripts of a chromosome in one list and sort the
        # lists by start coordinates
        print('SORTING')
        # transcripts[chromosom_id] = [list of Transcripts()]
        #transcripts = {}

        # transcripts_end[chromosom_id] = [[index of transcripts[chr], transcripts[chr][index]]]
        #transcripts_end = {}


        # tx_start_end[chr] = [tx_id, coord, id for start or end]
        # for every tx one entry for start and one for end
        tx_start_end = {}
        for k in self.anno.keys():
            for tx in self.anno[k].get_transcript_list():
                if tx.chr not in tx_start_end.keys():
                    tx_start_end.update({tx.chr : []})
                    #transcripts.update({tx.chr : []})
                #transcripts[tx.chr].append(tx)
                key = '{};{}'.format(tx.source_anno, \
                    tx.id)
                self.nodes.update({key : Node(tx.source_anno, \
                    tx.id)})
                tx_start_end[tx.chr].append([key, tx.start, 1])
                tx_start_end[tx.chr].append([key, tx.end, 0])

        #detect overlaps
        edge_count = 0
        for chr in tx_start_end.keys():
            tx_start_end[chr] = sorted(tx_start_end[chr], key=lambda t:(t[1], t[2]))
            open_intervals = []
            for interval in tx_start_end[chr]:
                if interval[2] == 1:
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
        print('FINISHED BUILDING')

    def compare_tx_cds(self, tx1, tx2):
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
        #returns list of connected components
        visited = []
        self.component_list
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

    def component_to_no_edge_subtree(self, component):
        return self.no_edge_subtree('', component)

    def no_edge_subtree(self, prefix, node_list):
        result = [prefix + n for n in node_list]
        for i in range(0, len(node_list)-1):
            new_node_list = [n for n in node_list[i+1:] if n not in \
                self.nodes[node_list[i]].edge_to.keys()]
            result += self.no_edge_subtree(result[i] + ';', new_node_list)
        return result

    def add_node_features(self, evi):
        print('# ADD NODE FEATURES')
        for key in self.nodes.keys():
            tx = self.__tx_from_key__(key)
            new_node_feature = Node_features(tx, evi, self.anno_pref)
            self.nodes[key].feature_vector = new_node_feature.get_features()
            if self.nodes[key].feature_vector[0] > 0  or self.nodes[key].feature_vector[1] > 0:
                self.nodes[key].evi_support = True

    def add_edge_features(self, hintfiles):
        # hintfiles list of Hintfiles()
        print('ADD EDGE FEATURES')
        for key in self.edges.keys():
            anno1, tx1_id = self.edges[key].node1.split(';')
            tx1 = self.anno[anno1].transcripts[tx1_id]
            #print(self.edges[key].node2.split(';'))
            anno2, tx2_id = self.edges[key].node2.split(';')
            tx2 = self.anno[anno2].transcripts[tx2_id]

            if tx1.start < tx2.start:
                start = tx1.start
            else:
                start = tx2.start
            if tx1.end > tx2.end:
                end = tx1.end
            else:
                end = tx2.end
            evi = []
            for h_file in hintfiles:
                evi += h_file.hints_in_range(start, end, tx1.chr)
            type_dict = {'intron' : ['intron'], 'start_stop' : \
                ['start_codon', 'stop_codon']}
            for type_key in type_dict.keys():
                tx1_lst = []
                tx2_lst = []
                evi_lst = []
                for type in type_dict[type_key]:
                    if type in tx1.transcript_lines.keys():
                        tx1_lst += tx1.transcript_lines[type]
                    if type in tx2.transcript_lines.keys():
                        tx2_lst += tx2.transcript_lines[type]
                    evi_lst += [e for e in evi if e[2] == type]
                self.edges[key].features.update({type_key : Edge_features(tx1_lst, \
                    tx2_lst, evi_lst).feature_vector})
        #anno1;3053_t
        #print(self.edges[e_id].features)
    def decide_node(self, edge):
        n1 = self.nodes[edge.node1]
        n2 = self.nodes[edge.node2]
        # print(n1.feature_vector)
        # print(n2.feature_vector)

        for i in range(0,7):
            if n1.feature_vector[i] > n2.feature_vector[i]:
                self.f[i].append(n2.id)
                return n2.id
            elif n1.feature_vector[i] < n2.feature_vector[i]:
                self.f[i].append(n1.id)
                return n1.id
        return None


    def decide_edge(self, edge):
        if 'intron' in edge.features.keys():
            vec = edge.features['intron']
            if (vec[0] + vec[2]) > (vec[5] + vec[7]):
                self.f[0].append(edge.node2)
                return [edge.node1]
            if (vec[0] + vec[2]) < (vec[5] + vec[7]):
                self.f[0].append(edge.node1)
                return [edge.node2]

        vec = edge.features['start_stop']
        if (vec[0] + vec[2]) > (vec[5] + vec[7]):
            self.f[1].append(edge.node2)
            return [edge.node1]
        if (vec[0] + vec[2]) < (vec[5] + vec[7]):
            self.f[1].append(edge.node1)
            return [edge.node2]
        if 'intron' in edge.features.keys():
            vec = edge.features['intron']
            if vec[4] < vec[9]:
                self.f[2].append(edge.node2)
                return [edge.node1]
            if vec[4] > vec[9]:
                self.f[2].append(edge.node1)
                return [edge.node2]

        vec = edge.features['start_stop']
        if vec[4] < vec[9]:
            self.f[3].append(edge.node2)
            return [edge.node1]
        if vec[4] > vec[9]:
            self.f[3].append(edge.node1)
            return [edge.node2]
        #anno1 = node1.split(';')[0]
        #anno2 = node2.split(';')[0]
        if edge.node1.split(';')[0] == 'braker2':
            self.f[4].append(edge.node2)
            return [edge.node1]
        elif edge.node2.split(';')[0] == 'braker2':
            self.f[4].append(edge.node1)
            return [edge.node2]
        return [edge.node1, edge.node2]

    def decide_component(self, component):
        result = component.copy()
        #visited_edges = []
        for node_id in component:
            for e_id in self.nodes[node_id].edge_to.values():
                #if e_id in visited_edges:
                    #continue
                #visited_edges.append(e_id)
                node_to_remove = self.edges[e_id].node_to_remove
                if node_to_remove:
                    if node_to_remove in result:
                        result.remove(node_to_remove)
                '''
                node_to_remove = [node for node in [edge.node1, edge.node2] \
                    if node not in edge.decision]
                for n in node_to_remove:
                    if n in result:
                        result.remove(n)
                '''
        return result

    def decide_graph(self):
        print('# DECIDE EDGES')
        for key in self.edges.keys():
            self.edges[key].node_to_remove = self.decide_node(self.edges[key])
            #edge.decision = self.decide_edge(edge)
            #ar = ['anno2;g482.t1', 'anno1;g479.t1']
            #if edge.node1 in ar and edge.node2 in ar:
                #print(edge.decision)
                #print(edge.features)
        print('# CREATE COMBINATION')
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
        if not self.decided_graph:
            self.decide_graph()
        result = {}
        for key in self.anno.keys():
            result.update({key : []})
        for node in self.decided_graph:
            if self.nodes[node].evi_support:
                anno_id, tx_id = node.split(';')
                result[anno_id].append([tx_id, self.nodes[node].component_id])
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
