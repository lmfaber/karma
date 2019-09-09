import itertools
import os

import matplotlib.pyplot as plt
import networkx as nx
from contig import Contig
from logs import logger
import numpy as np


class ReadGraph(nx.Graph):

    def __init__(self, incoming_graph_data=None, **attr):
        super().__init__(incoming_graph_data, **attr)

        self.original_contigs = []
        self.mcl_cluster = []

    @classmethod
    def from_contigs(cls, contigs):
        """Creates initital graph with a list of contig obejcts.
        
        Arguments:
            contigs {[type]} -- [description]
        """
        G = nx.Graph()
        logger.debug('Initial graph calculation.')
        for (first_contig, second_contig) in itertools.combinations(contigs, 2):
            # Normalize the weight:
            # (overlap/len(C1))(overlap/len(C2)) / 2
            overlap = len(first_contig.readset.intersection(second_contig.readset))
            mapped_reads_contig_A = len(first_contig.readset)
            mapped_reads_contig_B = len(second_contig.readset)

            try:
                weight = ( (overlap/mapped_reads_contig_A) + (overlap/mapped_reads_contig_B) ) / 2
            except ZeroDivisionError:
                weight = 0
            
            if weight > 0:
                G.add_edge(first_contig.name, second_contig.name, weight = weight)
            else:
                G.add_nodes_from([first_contig.name, second_contig.name])
        return cls(incoming_graph_data=G)

    def set_original_contigs(self, original_contigs):
        self.original_contigs = original_contigs

    @classmethod
    def from_equivalence_classes(cls, equivalence_class_file, sequences_from_fasta):
        """Builds the ReadGraph from a Salmon equivalence class file.
        
        Arguments:
            equivalence_class_file {str} -- Salmon equivalence file
        
        Returns:
            [ReadGraph] -- [description]
        """
        # Create a hash to identify the contig names from the equivalence class numbers
        with open(equivalence_class_file, 'r') as reader:
            no_of_contigs = int(reader.readline())
            _ = reader.readline()
            contig_hash = {}
            for i in range(no_of_contigs):
                contig = reader.readline().rstrip('\n')
                contig_hash[str(i)] = contig
            eq_classes = [line.rstrip('\n') for line in reader.readlines()]
        assert no_of_contigs == len(contig_hash)

        # Count the total number of mapped reads for each contig
        max_number_of_mapped_reads = {contig: 0 for contig in contig_hash.values()}
        for line in eq_classes:
            eq_size, *contig_ids, count = line.split('\t')
            count = int(count)
            for contig_id in contig_ids:
                contig_id = contig_hash[contig_id]
                max_number_of_mapped_reads[contig_id] += count
        assert no_of_contigs == len(max_number_of_mapped_reads.keys())

        # Count how many reads are shared by the contigs.
        G = nx.Graph()
        G.add_nodes_from(contig_hash.values())

        for line in eq_classes:
            eq_size, *contig_ids, count = line.split('\t')
            count = int(count)
            if eq_size == '1':
                continue
            else:
                for (A, B) in itertools.combinations(contig_ids, 2):
                    A = contig_hash[A]
                    B = contig_hash[B]
                    if G.has_edge(A, B):
                        existing_edge = G.get_edge_data(A, B)
                        old_weight = existing_edge['weight']
                        new_weight = old_weight + count
                        G.add_edge(A, B, weight = new_weight)
                    else:
                        G.add_edge(A, B, weight = count)
        assert len(G.nodes()) == no_of_contigs

        # Normalize and update weights
        weighted_graph = nx.Graph()
        weighted_graph.add_nodes_from(contig_hash.values())
        for A, B, data in G.edges(data = True):
            shared_reads = data['weight']
            logger.debug(f'Shared Reads: {A} / {B} -- {shared_reads}')
            if shared_reads == 0:
                weighted_graph.add_nodes_from([A, B])
            else:
                total_reads_A = max_number_of_mapped_reads[A]
                total_reads_B = max_number_of_mapped_reads[B]
                normalized_weight = ( (shared_reads / total_reads_A) + (shared_reads / total_reads_B) ) / 2
                weighted_graph.add_edge(A, B, weight = normalized_weight)

        assert len(weighted_graph.nodes()) == no_of_contigs

        # Add sequences that are not present in the equivalence class file
        original_sequence_names = set([name.lstrip('>') for name in sequences_from_fasta.keys()])
        difference = original_sequence_names.difference(set(weighted_graph.nodes())) # Missing contigs
        for missing_node in difference: # Add missing contigs
            weighted_graph.add_node(missing_node)
        assert len(weighted_graph.nodes()) == len(sequences_from_fasta), "The read graph has not enough nodes. Maybe Salmon could couldn't add all contigs to a equivalence class"

        return cls(incoming_graph_data=weighted_graph)

    def get_unconnected_nodes(self):
        """[summary]
        
        Returns:
            [type] -- [description]
        """
        unconnected_nodes = []
        for node in self.nodes():
            if len(list(nx.all_neighbors(self, node))) == 0:
                unconnected_nodes.append(node)
        return unconnected_nodes
    
    def get_connected_nodes(self):
        """ Returns all connected nodes in a graph as list """
        connected_nodes = []
        for node in self.nodes():
            if len(list(nx.all_neighbors(self, node))) != 0:
                connected_nodes.append(node)
        return connected_nodes

    def __calculate_node_weights(self):
        """
        For all nodes:
        Adds the weights for each connection for a node in the graph.
        Returns a dictionary.
        """
        node_weights = {}
        for node in self.nodes():
            edges = self.edges(node, data=True)

            node_weight = 0
            for _, _, w in edges:
                node_weight += w['weight']
            node_weights[node] = node_weight
        logger.debug(f'Node weights: {node_weights}')
        return node_weights

    def update_graph(self, contigs):
        """Updates the graph with new given contig objs. And trims the graph again.
        
        Arguments:
            contigs {[type]} -- [description]
        """
        logger.debug('Updating graph.')
        logger.debug(f'{self.original_contigs}')
        for (existing_contig, new_contig) in itertools.product(self.original_contigs, contigs):
            overlap = len(existing_contig.readset.intersection(new_contig.readset))
            # Normalize the weight:
            # (overlap/len(C1))(overlap/len(C2)) / 2
            mapped_reads_contig_A = len(existing_contig.readset)
            mapped_reads_contig_B = len(new_contig.readset)

            try:
                weight = ( (overlap/mapped_reads_contig_A) + (overlap/mapped_reads_contig_B) ) / 2 
            except ZeroDivisionError:
                weight = 0

            if weight > 0:
                self.add_edge(existing_contig.name, new_contig.name, weight = weight)
            else:
                self.add_node(new_contig.name)

    def save_graph(self, filename):
        """
        Saves the graph as abc-file. E.g.
        A B 20
        A C 10
        """
        # Add empty edges with weight of zero for mcl clustering to work
        G = nx.Graph(self)
        for (A, B) in itertools.combinations(self.nodes, 2):
            if not G.has_edge(A, B):
                G.add_edge(A, B, weight=0)
        
        if os.path.isfile(filename):
            os.remove(filename)
        nx.write_weighted_edgelist(G, filename)
    
    def get_contigs_not_in_mcl_cluster(self, mcl_cluster_file, original_contigs=None):
        """
        Remove all contigs that don't have a connection to the original contigs.
        Returns all contigs that either:
            1. are not in a group with the original contigs.
            2. Dont have a connection to another contig
        """

        # The last group only contains unlabeled contigs to every  sequence is a original sequence.
        if original_contigs == None:
            original_cluster = set([contig.name for contig in self.original_contigs])
        else:
            original_cluster = set(original_contigs)

        not_connected_to_original_contigs = []

        with open(mcl_cluster_file, 'r') as cluster_reader:
            for line in cluster_reader:
                line = line.rstrip('\n')
                line = set(line.split('\t'))

                if len(original_cluster.intersection(line)) == 0:
                    not_connected_to_original_contigs += list(line)
                else:
                    self.mcl_cluster.append(list(line))
        if len(self.mcl_cluster) == 0:
            logger.debug('MCL CLUSTER IS EMPTY')
            self.mcl_cluster = [[seq] for seq in original_cluster]

        return not_connected_to_original_contigs

    def get_contigs_not_in_mcl_cluster_stdout(self, stdout, original_contigs=None):
        """
        Remove all contigs that don't have a connection to the original contigs.
        Returns all contigs that either:
            1. are not in a group with the original contigs.
            2. Dont have a connection to another contig
        """

        # The last group only contains unlabeled contigs to every  sequence is a original sequence.
        if original_contigs == None:
            original_cluster = set([contig.name for contig in self.original_contigs])
        else:
            original_cluster = set(original_contigs)

        not_connected_to_original_contigs = []

        # Remove last row and split at new lines
        stdout_iterator = (line for line in stdout.split('\n')[:-1])
        for line in stdout_iterator:
            line = line.rstrip('\n')
            line = set(line.split('\t'))

            if len(original_cluster.intersection(line)) == 0:
                not_connected_to_original_contigs += list(line)
            else:
                self.mcl_cluster.append(list(line))
        # if len(self.mcl_cluster) == 0:
        #     logger.debug('MCL CLUSTER IS EMPTY')
        #     self.mcl_cluster = [[seq] for seq in original_cluster]

        return not_connected_to_original_contigs

    def calculate_representative_sequences(self, lowest=False):
        """
        Calculates the sum weight of each node to select a representative sequence per group. A sequence that has a high score of summed weights should mean that it has the highest similarity to all other sequences, thus being chosen as 'consensus sequence'.
        To further increase the diversity and information of the assembly, the sequence with the lowest score is also chosen. The thought here is, that this sequence is the farthest away from all other sequences, so it will raise the information level of the assembly.
        Just a thought.
        """

        node_weights = self.__calculate_node_weights()

        # For each cluster take the sequence with the highest and lowest weight for now.
        representative_sequences = []
        for i, cluster in enumerate(self.mcl_cluster):

            sub_node_weight = dict((k, node_weights[k]) for k in cluster)
            max_sequence_name = max(sub_node_weight, key=sub_node_weight.get)
            representative_sequences.append( f'>{max_sequence_name}' )
            if lowest:
                min_sequence_name = min(sub_node_weight, key=sub_node_weight.get)
                representative_sequences.append( f'>{min_sequence_name}' )
        
        return representative_sequences

    def nodes_list(self):
        """
        Returns the nodes in the graph.
        """
        return(list(self.nodes()))

    def edge_list(self):
        """ Return the edge list as a binary string for stdin for mcl. """
        edges = (f"{A} {B} {data['weight']}" for A, B, data in self.edges(data=True))
        return '\n'.join(edges).encode('utf-8')

    def calc_distance_between_subgraphs(self, nodes_A, nodes_B):
        import itertools
        weight = 0
        for A, B in itertools.product(nodes_A, nodes_B):
            if self.has_edge(A, B):
                weight += self[A][B]['weight']
        return weight


    def draw_graph(self, output_file, draw):
        """Draws a graph without edges of weight zero.
        
        Arguments:
            G {graph} -- Networkx graph.
            output_file {str} -- Output file.
        """
        if not draw:
            return
        else:
            if os.path.exists(output_file):
                os.remove(output_file)
            # figsize is intentionally set small to condense the graph
            fig, ax = plt.subplots(figsize=(10,10))
            margin=0.33
            fig.subplots_adjust(margin, margin, 1.-margin, 1.-margin)
            ax.axis('equal')

            pos = nx.circular_layout(self)

            node_list = self.nodes()
            n = len(node_list)
            angle_dict = {}
            for i, node in zip(range(n), node_list):
                theta = 2.0 * np.pi * i/n
                angle_dict[node] = theta

            # Nodes
            elarge = [(u, v) for (u, v, d) in self.edges(data=True) if d['weight'] > 0]
            esmall = [(u, v) for (u, v, d) in self.edges(data=True) if d['weight'] == 0.0]
            nx.draw_networkx_nodes(self, pos=pos, node_size=500, ax = ax)

            # Edges
            nx.draw_networkx_edges(self, pos, edgelist=elarge,width=2,ax = ax)
            nx.draw_networkx_edges(self, pos, edgelist=esmall,width=0,ax = ax)

            # Edge labels of non zero weights
            edge_labels = {}
            for u,v in self.edges():
                if self[u][v]['weight'] != 0:
                    edge_labels[(u,v)] = round(self[u][v]['weight'], 3)
            nx.draw_networkx_edge_labels(self, pos, label_pos = 0.4, edge_labels = edge_labels,ax = ax)

            # Draw labels
            labels = {key:key for key in node_list}
            description = nx.draw_networkx_labels(self, pos=pos, labels=labels,ax = ax)

            # Correct the position of the Node names
            r = fig.canvas.get_renderer()
            trans = plt.gca().transData.inverted()
            for node, t in description.items():
                bb = t.get_window_extent(renderer=r)
                bbdata = bb.transformed(trans)
                orientation = ''

                if 5/6 * np.pi >= angle_dict[node] >= 1/6 * np.pi:
                    orientation = 'TOP'
                elif 7/6 * np.pi >= angle_dict[node] >= 5/6 * np.pi:
                    orientation = 'LEFT'
                elif 11/6 * np.pi >= angle_dict[node] >= 7/6 * np.pi:
                    orientation = 'DOWN'
                else:
                    orientation = 'RIGHT'

                if orientation == 'TOP' or orientation == 'DOWN':
                    radius = 1.2+bbdata.height/1.1
                elif orientation == 'RIGHT' or orientation == 'LEFT':
                    radius = 1.15+bbdata.width/2.

                offset = 0
                if (orientation == 'TOP' or orientation == 'DOWN') and (1/2 * np.pi >= angle_dict[node] or angle_dict[node] >= 3/2 * np.pi): # right
                    offset = 2.2*(len(node) / 100)
                elif (orientation == 'TOP' or orientation == 'DOWN') and (1/2 * np.pi <= angle_dict[node] or angle_dict[node] <= 3/2 * np.pi): # left
                    offset = -2.2 * (len(node) / 100)

                position = (radius*np.cos(angle_dict[node])+offset,radius* np.sin(angle_dict[node]))   

                t.set_position(position)
                t.set_clip_on(False)
            plt.axis('off')
            plt.savefig(output_file, format='SVG')
            plt.close()
