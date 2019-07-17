from Bio import SeqIO
import os
import networkx as nx
import itertools
import subprocess
import matplotlib.pyplot as plt
import glob
import shutil
from cmd import Cmd
from copy import deepcopy

from logs import logger
from contig import Contig


def makedir(path):
    if not os.path.exists(path):
        os.mkdir(path)

def basename_woe(path):
    """ Basename without extension """
    return os.path.splitext(os.path.basename(path))[0]


class ReadGraph():

    def __init__(self, labeled, unlabeled, mapping_dir, no, graph_dir, mcl_dir, draw = False, threads = 4):
        # Initialized variables
        self.labeled = labeled
        self.unlabeled = unlabeled
        self.list_of_contigs = self.labeled + self.unlabeled

        self.threads = threads
        self.mapping_dir = mapping_dir
        self.cluster_number = no
        self.graph_dir = graph_dir
        self.mcl_dir = mcl_dir
          

        # Created in the Class.
        self.contig_objs = []
        self.G = nx.Graph()

        self.cluster_file = None
        self.mcl_cluster = []

        # Information for main program:
        self.representative_sequences = []

        if draw:
            self.drawing_routine()
        else:
            self.load_contig_info()
            self.create_graph()
            self.mcl_clustering()
            self.extract_cluster()
            self.calculate_representative_sequences()


    def load_contig_info(self):
        """
        Reads the sam file for each contig and saves the needed information in a Contig object.
        """
        logger.debug(f'Creating new contig object...')
        for contig_name in self.list_of_contigs:
            sam_file = f'{self.mapping_dir}/{contig_name}.sam'
            with open(sam_file, 'r') as sam_reader:
                new_contig = Contig(contig_name)
                for line in sam_reader:
                    infos = line.split('\t')
                    name = infos[2]
                    read = infos[0]
                    position = infos[3]
                    # TODO: Check if this can be left out. Should be possible, because mapping took place on single contigs.
                    #if name == new_contig.name:
                    #    new_contig.add_read(read, position)
                    new_contig.add_read(read, position)
                self.contig_objs.append(new_contig)

    def create_graph(self):
        logger.info('Calculate graph.')
        for (first_contig, second_contig) in itertools.combinations(self.contig_objs, 2):
            overlap = len(first_contig.readset.intersection(second_contig.readset))
            # Normalize the weight:
            # (overlap/len(C1))(overlap/len(C2)) / 2
            mapped_reads_contig_A = len(first_contig.readset)
            mapped_reads_contig_B = len(first_contig.readset)

            try:
                weight = ( (overlap/mapped_reads_contig_A) + (overlap/mapped_reads_contig_B) ) / 2
                self.G.add_edge(first_contig.name, second_contig.name, weight = weight)
            except ZeroDivisionError:
                weight = 0


    def mcl_clustering(self):
        """
        MCL: The input is then a file or stream in which each line encodes an edge in terms of two labels 
         (the 'A' and the 'B') and a numerical value (the 'C'), all separated by white space.
        A B 20
        A C 10
        The output is then a file where each line is a cluster of tab-separated labels.
        """
        logger.debug('Save cluster as file')

        # TODO: CARE!!! clusters are not in order, so i need a file check that is good. Currently deleting the file and creating a new one, so the problem doesnt occur
        abc_file = f'{self.graph_dir}/graph_{self.cluster_number}.abc'
        if os.path.isfile(abc_file):
            os.remove(abc_file)
        nx.write_weighted_edgelist(self.G, abc_file)

        logger.debug('Clustering...')
        mcl_output_file = f'{self.mcl_dir}/cluster_{self.cluster_number}.mcl'
        self.cluster_file = mcl_output_file
        if os.path.isfile(mcl_output_file):
            os.remove(mcl_output_file)
        mcl = Cmd(f'mcl {abc_file} --abc -o {mcl_output_file} -te {self.threads} -resource 4 -V all')
        mcl.run()

    def extract_cluster(self):
        """
        Extracts only the clusters from the mcl output, that were in the originally kmer based cluster. All other clusters are then readded to the unlabeled cluster.
        """
        original_cluster = set(self.labeled)

        contigs_with_connections = []

        new_unlabeled_list = []
        with open(self.cluster_file, 'r') as cluster_reader:
            for line in cluster_reader:
                line = line.rstrip('\n')
                line = set(line.split('\t'))
                if len(original_cluster.intersection(line)) != 0:
                    contigs_with_connections += list(line)
                    self.mcl_cluster.append(list(line))
                else:
                    new_unlabeled_list += list(line)
        #logger.debug(f'Original clusters: {self.labeled}')
        logger.debug(f'Clusters extracted: {contigs_with_connections}')
        #logger.debug(f'Labels not clustered: {new_unlabeled_list}')
        ## Remove unnecessary contigs from graph
        for contig in new_unlabeled_list:
            self.G.remove_node(contig)
        
        self.unlabeled = new_unlabeled_list

    def calculate_representative_sequences(self):
        """
        Calculates the sum weight of each node to select a representative sequence per group. A sequence that has a high score of summed weights should mean that it has the highest similarity to all other sequences, thus being chosen as 'consensus sequence'.
        To further increase the diversity and information of the assembly, the sequence with the lowest score is also chosen. The thought here is, that this sequence is the farthest away from all other sequences, so it will raise the information level of the assembly.
        Just a thought.
        """

        node_weights = {}
        for node in self.G.nodes():
            edges = self.G.edges(node, data=True)

            node_weight = 0
            for _, _, w in edges:
                node_weight += w['weight']

            node_weights[node] = node_weight

        # For each cluster take the sequence with the highest weight for now. TODO lowest weight?
        for i, cluster in enumerate(self.mcl_cluster):
            logger.debug(f'{i}: {cluster}')
            sub_node_weight = dict((k, node_weights[k]) for k in cluster)
            max_sequence_name = max(sub_node_weight, key=sub_node_weight.get)
            self.representative_sequences.append( f'>{max_sequence_name}' )
        
    def drawing_routine(self):
        
        self.load_contig_info()
        self.create_graph()

        # Draw before clustering
        graph_visual_dir = f'{self.graph_dir}/visual'
        makedir(graph_visual_dir)

        plt.figure(figsize=(20,10))
        plt.subplot(121)
        self.draw_graph(self.G, self.unlabeled, output_file)
        
        self.mcl_clustering()
        self.extract_cluster()

        # Draw after clustering
        plt.subplot(122)
        self.draw_graph(self.G, self.unlabeled, output_file)

        output_file = f'{graph_visual_dir}/cluster_{self.cluster_number}'
        plt.savefig(output_file, format='SVG')
        plt.close()

        self.calculate_representative_sequences()

    def draw_graph(self, G, nodes_to_remove):
        """Draws a graph without edges of weight zero.
        
        Arguments:
            G {graph} -- Networkx graph.
            output_file {str} -- Output file.
        """
        # Copy the graph
        G = nx.Graph()
        G = deepcopy(self.G)

        for contig in nodes_to_remove:
            if contig in G.nodes():
                G.remove_node(contig)

        elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > 0]
        esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] == 0.0]


        pos = nx.circular_layout(G)
        # nodes
        nx.draw_networkx_nodes(G, pos, node_size=500)
        # edges
        nx.draw_networkx_edges(G, pos, edgelist=elarge,width=2)
        nx.draw_networkx_edges(G, pos, edgelist=esmall,width=0)
        # labels
        nx.draw_networkx_labels(G, pos)
        # Edge labels of non zero weights:
        edge_labels = {}
        for u,v in G.edges():
            if G[u][v]['weight'] != 0:
                edge_labels[(u,v)] = round(G[u][v]['weight'], 3)

        nx.draw_networkx_edge_labels(G, pos, label_pos = 0.3, edge_labels = edge_labels)
        plt.axis('off')
        plt.tight_layout()

    def calculate_alignments(self):
        # TODO
        alignment_dir = f'{RESULT_DIR}/alignments'
        makedir(alignment_dir)

        with open(mcl_output_file, 'r') as reader, open(fastaFile, 'r') as fastaReader:
            original_fasta_sequences = SeqIO.to_dict(SeqIO.parse(fastaReader, 'fasta'))
            for j, line in enumerate(reader, 1):
                with open(f'{alignment_dir}/cluster_{cluster_no}_{j}.fa', 'w') as writer:
                    line = line.split()
                    for contig in line:
                        SeqIO.write(original_fasta_sequences[contig], writer, 'fasta')
                os.system(f'mafft --quiet --auto {alignment_dir}/cluster_{cluster_no}_{j}.fa > {alignment_dir}/cluster_{cluster_no}_{j}.aln')
        
