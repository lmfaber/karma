#!/usr/bin/env python3
from numba.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning, NumbaWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaWarning)

import os
from cmd import Cmd
from collections import OrderedDict
import sys
import itertools

from annotate import Dammit
from cmd_parser import args
from contig import Contig
from hisat2 import Hisat2
from kmer import KmerClustering
from logs import logger
from mcl import Mcl
from read_graph import ReadGraph
from writer import Writer
from salmon import Salmon
import networkx as nx
from multiprocessing import Pool


def makedir(path):
    """ Create path if it doesnt exist"""
    if not os.path.exists(path):
        os.makedirs(path)

def basename_woe(path):
    """ Basename without extension """
    return os.path.splitext(os.path.basename(path))[0]

def read_fasta_file(fasta_file):
    """Reads a fasta file and saves the entries as OrderedDict from collections
    
    Arguments:
        fasta_file {str} -- Fasta file to read.
    """
    logger.info('Reading fasta file.')
    sequences = OrderedDict()
    sequence_names = []
    with open(fasta_file, 'r') as reader:
        seq_name = reader.readline().rstrip('\n').split(' ')[0]
        sequence = ''
        for line in reader:
            if line.startswith('>'):
                sequences[seq_name] = sequence
                seq_name = line.rstrip('\n').split(' ')[0]
                sequence = ''
            else:
                sequence += line.rstrip('\n')
        sequences[seq_name] = sequence
    logger.debug(f'Read {len(sequences)} sequences in total.')
    return(sequences)
    
def flatten(lst):
	return sum( ([x] if not isinstance(x, list) else flatten(x) for x in lst), [] )

def save_options(file_path):
    """ Saves all parameters set in a file. """
    if os.path.exists(file_path):
        logger.warning('Parameter file will be overwritten')
    with open(file_path, 'w') as writer:
        writer.write(f"Command: {' '.join(sys.argv)}\n\n")
        for entry in vars(args):
            writer.write(f"{entry}: {getattr(args, entry)}\n")

def create_lookup_dict(clusters_with_subcluster, sequences):
    """ 
    Creates a lookup dictionary for the kmer + mcl clustered sequences. Saves the original cluster number to which each subgroup belongs and the actual cluster.
    """
    mcl_subclusters = {}
    index = 0
    for cl_no, cluster in enumerate(clusters_with_subcluster, 1):
        for mcl_cluster in cluster:
            mcl_subclusters[index] = {'previous_cluster': cl_no, 'mcl_subcluster': mcl_cluster}
            index += 1

    # Check if dictionary is okay.
    length_of_dict = sum( [len(subcluster['mcl_subcluster']) for subcluster in mcl_subclusters.values()] )
    logger.debug(f"length of dict: {length_of_dict}, {len(sequences)}")
    assert length_of_dict == len(sequences), 'The creation of the lookup dictionary went wrong.'
    return mcl_subclusters

def calc_connections_between_mcl_subclusters(mcl_subclusters, weight_cutoff=0):
    """ 
    Calculates if there is a connections between subclusters and the according weight. If there is a connection, those two groups are meant to be connected later. Returns a list of tuples containing each two subclusters, that share reads.
    """
    mcl_groups_to_combine = []
    for index_A, index_B in itertools.combinations(mcl_subclusters, 2):
        nodes_A = mcl_subclusters[index_A]['mcl_subcluster']
        nodes_B = mcl_subclusters[index_B]['mcl_subcluster']
    
        weight = 0
        for A, B in itertools.product(nodes_A, nodes_B):
            if full_graph.has_edge(A, B):
                weight += full_graph[A][B]['weight']
                if weight > weight_cutoff:
                    mcl_groups_to_combine.append([index_A, index_B])
    return mcl_groups_to_combine

def remove_already_added_clusters(from_dict, remove):
    """
    Removes the clusters that are connected and were already added to the final output from the lookup dictionary.
    """
    for mcl in remove:
        from_dict.pop(mcl)
    return from_dict

def combine_connected_subclusters(mcl_subclusters, mcl_groups_to_combine):
    """ Creates a graph to identify all connected subclusters with each other and returns the nested array of clusters with their respective subclusters. """
    mcl_cluster_connection_graph = nx.Graph()
    mcl_cluster_connection_graph.add_edges_from(mcl_groups_to_combine)

    connected_subclusters = []
    for new_cluster in nx.k_edge_subgraphs(mcl_cluster_connection_graph, k=1):
        connected_subclusters.append( [mcl_subclusters[a]['mcl_subcluster'] for a in new_cluster] )
    assert len(list(nx.k_edge_subgraphs(mcl_cluster_connection_graph, k=1))) == len(connected_subclusters), 'Combining went wrong'
    return connected_subclusters

def add_remaining_kmer_based_clusters(mcl_subclusters):
    """ Restores the remaining kmer based clusters as they were before. """
    remaining_cluster = []
    first_cl_no = -1
    first = True
    for key, value in mcl_subclusters.items():
        current_orig_cluster = value['previous_cluster']
        if first_cl_no != current_orig_cluster:
            first_cl_no = value['previous_cluster']
            if first:
                first = False
            else:
                remaining_cluster.append(new_subgroup)
            new_subgroup = []
        if first_cl_no == current_orig_cluster:
            new_subgroup.append(value['mcl_subcluster'])
    remaining_cluster.append(new_subgroup)
    return remaining_cluster


if __name__ == "__main__":

    ####             ####
    ## output folders ###
    ####             ####
    BASENAME = basename_woe(args.FASTA_FILE)

    kmer_dir = f'{args.OUTPUT_DIR}/kmer'               # Kmer based clustering
    mapping_dir = f'{args.OUTPUT_DIR}/salmon/mapping'  # Sam files from mapping
    index_dir = f'{args.OUTPUT_DIR}/salmon/index'      # Index files for mapping
    graph_dir = f'{args.OUTPUT_DIR}/graphs'            # abc-graph
    graph_visual_dir = f'{graph_dir}/visual'           # plotted graphs
    dammit_dir = f'{args.OUTPUT_DIR}/dammit'           # de novo annotation results

    # Create all relevant output directories
    directories = [kmer_dir, mapping_dir, index_dir]
    if args.DRAW:
        directories.append(graph_visual_dir)
        directories.append(graph_dir)
    if args.ANNOTATE:
        directories.append(dammit_dir)
    [makedir(dir) for dir in directories]

    save_options(f"{args.OUTPUT_DIR}/parameters.txt")

    sequences = read_fasta_file(args.FASTA_FILE)

    ####             ####
    ## kmer clustering ##
    ####             ####
    logger.info('Starting kmer based clustering.')

    k = KmerClustering(sequences = sequences, 
                        output_dir = kmer_dir, 
                        kmer_size=args.KMER_SIZE, 
                        threads = args.THREADS)

    k.run(neighbors=args.N_NEIGHBORS, 
            components=args.N_COMPONENTS, 
            dist=args.MIN_DIST, 
            r_state=args.RANDOM_STATE,
            min_cluster_size=args.MIN_CLUSTER_SIZE)

    labeled_contigs = k.clusters
    unlabeled_contigs = k.unlabeled_cluster[0]

    logger.debug(f'LABELED: {labeled_contigs}')
    logger.debug(f'NOT LABELED: {unlabeled_contigs}')

    ####     ####
    ## Mapping ##
    ####     ####
    if args.R:
        READS = [args.R]
    else:
        READS = [args.R1, args.R2]

    logger.info('Salmon mapping...')
    salmon = Salmon(input_file=args.FASTA_FILE, 
                    output_dir=mapping_dir, 
                    index_name=index_dir, 
                    threads=args.THREADS)
    salmon.build_index()
    salmon.run(reads=READS)

    eqc_file = f'{mapping_dir}/aux_info/eq_classes.txt'

    # Create full graph and take subgraph for clustering
    logger.info('Build full graph...')
    full_graph = ReadGraph.from_equivalence_classes(eqc_file, sequences)
    # full_graph.save_graph(f"{output_dir}/full_graph_edges.txt")

    ## Extract all subclusters, remove single contigs and save remaining contigs in a list.
    logger.info('Trim kmer based clusters...')
    trimmed_clusters = []
    logger.debug(f'LABELED CONTIGS:  {labeled_contigs}')
    logger.debug(f'UNLABELED CONTIGS : {unlabeled_contigs}')
    logger.debug(flatten(labeled_contigs))
    assert len(flatten(labeled_contigs) + flatten(unlabeled_contigs)) == len(sequences), "LOST A SEQUENCE :("
    logger.debug(f'SEQUENCES : {len(sequences)}')
    i = 1
    logger.debug(f'full graph nodes: {full_graph.nodes}')
    for  cluster in labeled_contigs:
        logger.debug(f'= = = = = Cluster {i} = = = = =')
        unlabeled_before = len(unlabeled_contigs)

        cluster_graph = ReadGraph(full_graph.subgraph(cluster))
        logger.debug(cluster)
        logger.debug(cluster_graph.nodes())
        assert len(cluster_graph.nodes()) == len(cluster), "Cluster graph not equal the cluster group"

        # Draw
        # logger.info(f"NODES: {cluster_graph.nodes}")
        cluster_graph.draw_graph(f'{graph_visual_dir}/cluster_{i}_A.svg', draw = args.DRAW)

        contigs_to_remove = cluster_graph.get_unconnected_nodes()
        cluster_graph.remove_nodes_from(contigs_to_remove)

        remaining_nodes = cluster_graph.nodes_list()
        
        unlabeled_contigs += contigs_to_remove

        # assert len(flatten(unlabeled_contigs)) == unlabeled_before + len(flatten(contigs_to_remove))
        if len(remaining_nodes) > 0:
            trimmed_clusters.append(remaining_nodes)
            i += 1

    logger.debug(f'LABELED: {trimmed_clusters}')
    logger.debug(f'LABELED: {unlabeled_contigs}')
    logger.debug(len(flatten(trimmed_clusters)))
    logger.debug(len(flatten(unlabeled_contigs)))
    logger.debug(len(sequences))
    assert len(flatten(trimmed_clusters) + flatten(unlabeled_contigs)) == len(sequences), "LOST A SEQUENCE :("


    ## Update the graphs with the unlabeled group.
    ## Remove nodes/contigs that:
    ## 1. Are not connected to any other node.
    ## 2. Are not but remove all contigs that are in a mcl-cluste that doesn't contain any original sequences.

    logger.info('Redefine kmer based clusters...')
    cluster_representative_sequences = []
    clusters_with_subcluster = []
    for i, cluster in enumerate(trimmed_clusters, 1):
    # TODO: TEST WITHOUT TRIMMING FROM THE BEGINNING
    # for i, cluster in enumerate(labeled_contigs, 1):
        logger.debug(f"= = = = = CLUSTER {i} = = = = =")
        logger.debug(f"Original Cluster: {cluster}")
        cluster_with_unlabeled = cluster + unlabeled_contigs

        cluster_graph = ReadGraph(full_graph.subgraph(cluster_with_unlabeled))
        assert len(cluster_graph.nodes) == len(cluster_with_unlabeled)

        # Draw
        # logger.info(f"NODES: {cluster_graph.nodes}")
        cluster_graph.draw_graph(f'{graph_visual_dir}/cluster_{i}_B.svg', draw = args.DRAW)

        mcl = Mcl(threads=args.THREADS, inflation=args.INFLATION)
        mcl_result = mcl.run_pipe(graph_file = cluster_graph.edge_list())
        
        logger.debug(f"MCL RESULTS: {mcl_result}")

        # Extract the groups from mcl that contain original sequences.
        unconnected_contigs = cluster_graph.get_unconnected_nodes()

        non_mcl_cluster_contigs = cluster_graph.get_contigs_not_in_mcl_cluster_stdout(mcl_result, original_contigs=cluster)

        assert len(flatten(cluster_graph.mcl_cluster)) != 0, f"Cluster {i} has no mcl clusters:\n {cluster}"
        cluster_graph.remove_nodes_from(non_mcl_cluster_contigs)
        assert len(flatten(non_mcl_cluster_contigs)) + len(cluster_graph.nodes) == len(cluster_with_unlabeled)
        cluster_graph.remove_nodes_from(unconnected_contigs)

        unlabeled_contigs = unconnected_contigs + non_mcl_cluster_contigs

        cluster_representative_sequences += cluster_graph.calculate_representative_sequences(lowest=args.LOWEST)

        clusters_with_subcluster.append(cluster_graph.mcl_cluster)
        logger.debug(f"MCL CLUSTER: {cluster_graph.mcl_cluster}")
        
        # Draw
        # logger.info(f"NODES: {cluster_graph.nodes}")
        cluster_graph.draw_graph(f'{graph_visual_dir}/cluster_{i}_C.svg', draw = args.DRAW)

        logger.debug(f'SUBCLUSTER: {(clusters_with_subcluster)}')
        logger.debug(f'UNLABELED: {(unlabeled_contigs)}')
        logger.debug('')
    

    logger.debug(len(flatten(clusters_with_subcluster)))
    logger.debug(len(flatten(unlabeled_contigs)))
    logger.debug(len(sequences))

    assert len(flatten(clusters_with_subcluster)) + len(flatten(unlabeled_contigs)) == len(sequences)

    ## Create graph and perform mcl for the unlabeled group.
    cluster_graph = ReadGraph(full_graph.subgraph(unlabeled_contigs))
    # DRAW
    cluster_graph.draw_graph(f'{graph_visual_dir}/cluster_unlabeled.svg', draw = args.DRAW)

    mcl = Mcl(threads=args.THREADS, inflation=args.INFLATION) 
    mcl_result = mcl.run_pipe(graph_file = cluster_graph.edge_list())

    cluster_graph.get_contigs_not_in_mcl_cluster_stdout(mcl_result, original_contigs=unlabeled_contigs)
    logger.debug(f"LAST GRAPH, mcl: {cluster_graph.mcl_cluster}")
    logger.debug(f"LAST GRAPH, unconncected: {cluster_graph.get_unconnected_nodes()}")
    logger.debug(f"mcl result: {mcl_result}")

    # Add those remaining contigs that from a cluster
    for mcl in cluster_graph.mcl_cluster:
        clusters_with_subcluster.append([mcl])
    
    # Add those remaining contigs that do NOT form a cluster.
    leftover = [[[singles]] for singles in cluster_graph.get_unconnected_nodes()]
    if len(cluster_graph.get_unconnected_nodes()) != 0:
        clusters_with_subcluster += leftover

    logger.debug(f"len last cluster mcl: {len(cluster_graph.mcl_cluster)} ")
    # logger.debug(f"Clusters with subcluster: {clusters_with_subcluster}")
    logger.debug(f"len leftover: {len(leftover)}")
    logger.debug(len(flatten(clusters_with_subcluster)))
    logger.debug(len(sequences))
    assert len(flatten(clusters_with_subcluster)) == len(sequences)

    clustered_sequence_names = []
    unlabeled_representative_sequences = cluster_graph.calculate_representative_sequences(lowest=args.LOWEST)
    leftover_names = [f'>{name[0][0]}' for name in leftover]
    clustered_sequence_names = unlabeled_representative_sequences + cluster_representative_sequences 
    clustered_sequence_names += leftover_names

    ####                ####
    ## Rearrange clusters ##
    ####                ####
    if args.REARRANGE:
        logger.info('Rearrange clusters.')
        # logger.info(clusters_with_subcluster)

        mcl_subclusters = create_lookup_dict(clusters_with_subcluster, sequences)
        mcl_groups_to_combine = calc_connections_between_mcl_subclusters(mcl_subclusters, weight_cutoff=args.THRESHOLD)

        # New object that stores the new clusters with subclusters in a nested array list.
        new_cluster_subcluster = []

        # Combine those clusters that belong together.
        new_cluster_subcluster += combine_connected_subclusters(mcl_subclusters, mcl_groups_to_combine)

        # Remove already added clusters from lookup dict
        mcl_subclusters = remove_already_added_clusters(from_dict=mcl_subclusters, remove=set(flatten(mcl_groups_to_combine)) )

        # Add remaining clusters according to the kmer based clustering.
        new_cluster_subcluster += add_remaining_kmer_based_clusters(mcl_subclusters)

        assert(len(flatten(new_cluster_subcluster)) == len(sequences)), 'rearraning groups went wrong.'
    else:
        new_cluster_subcluster = clusters_with_subcluster


    ####    ####
    ## Output ##
    ####    ####

    # Write fasta file
    fasta_output_file = f'{args.OUTPUT_DIR}/{BASENAME}.fa'
    sequences_to_write = {name: sequences[name] for name in clustered_sequence_names}
    writer = Writer(fasta_output_file)
    writer.write_fasta(sequences_to_write)

    # Write clstr file
    clstr_output_file = f'{args.OUTPUT_DIR}/{BASENAME}.clstr'
    writer = Writer(clstr_output_file)
    writer.write_clstr(new_cluster_subcluster, clustered_sequence_names)

    reduction_percent = ( len(clustered_sequence_names) / len(sequences) ) * 100
    logger.info(f'Reduced to ~{ round(reduction_percent, 2) } % of original sequences.')

    ####    ####
    ## Dammit ##
    ####    ####
    if args.ANNOTATE:
        logger.info('Starting annotation.')

        labels_after_kmer_clustering = k.unlabeled_cluster + k.clusters
        labels_after_read_graph = new_cluster_subcluster    

        # logger.info(f'LABELS AFTER KMER CLUSTERING: {labels_after_kmer_clustering}')
        # logger.info(f'LABELS AFTER READ GRAPH     : {labels_after_read_graph}')

        dammit = Dammit(first_file = args.FASTA_FILE, second_file = fasta_output_file, output_dir = dammit_dir, kmer_clusters = labels_after_kmer_clustering, karma_clusters = labels_after_read_graph, representative_sequences = clustered_sequence_names, threads = args.THREADS)
        dammit.update_database(database_dir=args.DATABASE_DIR, busco_group=args.BUSCO_GROUP)
        dammit.run()
        dammit.postprocessing()

        eval_output_file = f'{args.OUTPUT_DIR}/{BASENAME}.csv'
        dammit.save(eval_output_file)
