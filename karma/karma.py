#!/usr/bin/env python3
import glob
import os
import pickle
import subprocess
from cmd import Cmd
from collections import OrderedDict

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

def save_contigs_in_separate_files(sequences, output_folder):
    makedir(contig_dir)
    for name, sequence in sequences.items():
        file_name = f"{contig_dir}/{name.lstrip('>')}.fa"
        with open(file_name, 'w') as fasta_writer:
            fasta_writer.write(f'{name}\n')
            fasta_writer.write(sequence)

def read_fasta_file(fasta_file):
    """Reads a fasta file and saves the entries as OrderedDict from collections
    
    Arguments:
        fasta_file {str} -- Fasta file to read.
    """
    logger.info('Reading fasta file.')
    sequences = OrderedDict()
    sequence_names = []
    with open(fasta_file, 'r') as reader:
        seq_name = reader.readline().rstrip('\n')
        sequence = ''
        for line in reader:
            if line.startswith('>'):
                sequences[seq_name] = sequence
                seq_name = line.rstrip('\n')
                sequence = ''
            else:
                sequence += line.rstrip('\n')
        sequences[seq_name] = sequence
    logger.debug(f'Read {len(sequences)} sequences in total.')
    return(sequences)

def runner(contig_name, mapper):
    logger.debug(f'Mapping and saving {contig_name}')
    mapper.build_index()
    sam_info = mapper.run(READS)
    contig_obj = Contig(contig_name)
    contig_obj.load_from_iterator(sam_info)

    contig_pickle_file = f'{contig_dir}/{contig_name}.pickle'
    with open(contig_pickle_file, 'wb') as writer:
        pickle.dump(contig_obj, writer)
    return (contig_name, contig_obj)

def build_mapping_objects(n, total, contig_name, mapping_threads):
    logger.debug(f'{n}/{total}')
    # If mapping was already done load the pickle file.
    contig_pickle_file = f'{contig_dir}/{contig_name}.pickle'

    if os.path.exists(contig_pickle_file):
        # TODO: Better method for checking the file.
        with open(contig_pickle_file, 'rb') as reader:
            contig_obj = pickle.load(reader)
            mapping_obj = None
    else:
        contig_file = f'{contig_dir}/{contig_name}.fa'
        index_name = f'{index_dir}/{contig_name}'
        mapping_obj = Hisat2(input_file=contig_file, index_name=index_name, threads=mapping_threads)
        contig_obj = None
    return (contig_name, contig_obj, mapping_obj)

def perform_mappings(contigs, contig_dir, index_dir, threads):
    """Create all contig objects either trough mapping or loading from existing files.
    
    Arguments:
        contig_paths {[type]} -- [description]
        index_dir {[type]} -- [description]
        mapping_dir {[type]} -- [description]
    """
    contig_names = [name.lstrip('>') for name in contigs.keys()]
    contig_objs = {}
    mappings = []

    if threads <= 5:
        pool_threads = 1
        mapping_threads = 1
        logger.warning("This may take forever. Please use more threads...")
    else:
        # pool_threads = 5
        # mapping_threads = (threads // pool_threads)
        pool_threads = 10
        mapping_threads = 2

    total = len(contig_names)
    with Pool(threads) as pool:
        arguments = ((i, total, contig_name, mapping_threads) for (i, contig_name) in enumerate(contig_names))
        for (contig_name, contig_obj, mapping_obj) in pool.starmap(build_mapping_objects, arguments):
            if contig_obj:
                contig_objs[contig_name] = contig_obj
            if mapping_obj:
                mappings.append(mapping_obj)

    logger.info(f"Perform a whopping amount of {len(mappings)} mappings. Here we go.")
    with Pool(pool_threads) as pool:
        arguments = ((con, mapp) for (con, mapp) in zip(contig_names, mappings))
        for (contig_name, contig_obj) in pool.starmap(runner, arguments):
            contig_objs[contig_name] = contig_obj

    return contig_objs
    


if __name__ == "__main__":
    
    ####             ####
    ## output folders ###
    ####             ####
    BASENAME = basename_woe(args.FASTA_FILE)
    contig_dir = f'{args.OUTPUT_DIR}{BASENAME}/contigs'          # Saved single contigs
    kmer_dir = f'{args.OUTPUT_DIR}{BASENAME}/kmer'               # Kmer based clustering
    mapping_dir = f'{args.OUTPUT_DIR}{BASENAME}/hisat2/mapping'  # Sam files from mapping
    index_dir = f'{args.OUTPUT_DIR}{BASENAME}/hisat2/index'      # Index files for mapping
    # index_dir = f'/dev/shm/lasse/index'      # Index files for mapping
    graph_dir = f'{args.OUTPUT_DIR}{BASENAME}/graphs'            # abc-graph
    graph_visual_dir = f'{graph_dir}/visual'                     # plotted graphs
    mcl_dir = f'{args.OUTPUT_DIR}{BASENAME}/mcl'                 # mcl clustering results
    output_dir = f'{args.OUTPUT_DIR}{BASENAME}/karma'                             # actual results
    dammit_dir = f'{args.OUTPUT_DIR}{BASENAME}/dammit'           # de novo annotation results



    sequences = read_fasta_file(args.FASTA_FILE)

    # Split all contigs in separate file.
    save_contigs_in_separate_files(sequences, contig_dir)

    # Create all output directories
    directories = [kmer_dir, mapping_dir, index_dir, graph_dir, graph_visual_dir, mcl_dir, output_dir, dammit_dir]
    [makedir(dir) for dir in directories]

    ####             ####
    ## kmer clustering ##
    ####             ####
    logger.info('Starting kmer based clustering.')

    cluster_output = f'{args.OUTPUT_DIR}{BASENAME}/kmer/{BASENAME}.txt'

    k = KmerClustering(sequences = sequences, output_file = cluster_output, threads = args.THREADS, kmer_size=args.KMER_SIZE)
    labeled_contigs = k.clusters
    unlabeled_contigs = k.unlabeled_cluster[0]

    logger.debug(f'LABELED: {labeled_contigs}')
    logger.debug(f'NOT LABELED: {unlabeled_contigs}')


    ####     ####
    ## Mapping ##
    ####     ####
    logger.info('Starting mapping...')
    if args.R:
        READS = [args.R]
    else:
        READS = [args.R1, args.R2]

    if args.MAPPING_METHOD == 'salmon':
        # Perfrom salmon mapping.
        logger.info('Use Salmon...')
        salmon = Salmon(args.FASTA_FILE, mapping_dir, index_dir, args.THREADS)
        salmon.build_index()
        salmon.run(READS)

        eqc_file = f'{mapping_dir}/aux_info/eq_classes.txt'

        # Create full graph and take subgraph for clustering
        logger.info('Build full graph...')
        full_graph = ReadGraph.from_equivalence_classes(eqc_file)

        ## Extract all subclusters, remove single contigs and save remaining contigs in a list.
        logger.info('Trim kmer based clusters...')
        trimmed_clusters = []
        for i, cluster in enumerate(labeled_contigs, 1):       
            sub_graph = full_graph.get_subgraph(cluster)
            cluster_graph = ReadGraph(sub_graph)

            # Draw
            if args.DRAW:
                draw_filename = f'{graph_visual_dir}/cluster_{i}_A.svg'
                cluster_graph.draw_graph(draw_filename)

            contigs_to_remove = cluster_graph.get_unconnected_nodes()
            trimmed_cluster = cluster_graph.get_connected_nodes()
            
            unlabeled_contigs += contigs_to_remove
            trimmed_clusters.append(trimmed_cluster)

        ## Update the graphs with the unlabeled group.
        ## Remove nodes/contigs that:
        ## 1. Are not connected to any other node.
        ## 2. Are not but remove all contigs that are in a mcl-cluste that doesn't contain any original sequences.

        logger.info('Redefine kmer based clusters...')
        final_cluster = []
        cluster_representative_sequences = []
        clusters_with_subcluster = []
        for i, cluster in enumerate(trimmed_clusters, 1):

            clusters_with_unlabeled = cluster + unlabeled_contigs

            sub_graph = full_graph.get_subgraph(clusters_with_unlabeled)
            cluster_graph = ReadGraph(sub_graph)

            # Draw
            if args.DRAW:
                draw_filename = f'{graph_visual_dir}/cluster_{i}_B.svg'
                cluster_graph.draw_graph(draw_filename)

            abc_file = f'{graph_dir}/graph_{i}'
            cluster_graph.save_graph(abc_file)
            mcl_output_file = f'{mcl_dir}/graph_{i}'
            mcl = Mcl()
            mcl.run(graph_file = abc_file, output_file = mcl_output_file, threads=args.THREADS)

            # Extract the groups from mcl that contain original sequences.
            unconnected_contigs = cluster_graph.get_unconnected_nodes()
            non_mcl_cluster_contigs = cluster_graph.get_contigs_not_in_mcl_cluster(mcl_output_file, original_contigs=cluster)
            cluster_graph.remove_nodes(non_mcl_cluster_contigs)
            cluster_graph.remove_nodes(unconnected_contigs)

            unlabeled_contigs = unconnected_contigs + non_mcl_cluster_contigs
            # unlabeled_contigs = non_mcl_cluster_contigs

            logger.debug(f'2.Unlabeled contigs: {unlabeled_contigs}')

            representative_sequences = cluster_graph.calculate_representative_sequences()
            cluster_representative_sequences += representative_sequences
            logger.debug(f'MCL CLUSTER: {cluster_graph.mcl_cluster}')
            clusters_with_subcluster.append(cluster_graph.mcl_cluster)
            # Draw
            if args.DRAW:
                draw_filename = f'{graph_visual_dir}/cluster_{i}_C.svg'
                cluster_graph.draw_graph(draw_filename)

            final_cluster.append(cluster_graph.get_nodes())

        logger.info(f'Clusters are : {final_cluster}')


        ## Create graph and perform mcl for the unlabeled group.
        unlabeled_sub_graph = full_graph.get_subgraph(unlabeled_contigs)
        cluster_graph = ReadGraph(unlabeled_sub_graph)

        draw_filename = f'{graph_visual_dir}/cluster_unlabeled.svg'
        cluster_graph.draw_graph(draw_filename)

        mcl_output_file = f'{mcl_dir}/graph_unlabeled'
        mcl = Mcl()
        mcl.run(graph_file = abc_file, output_file = mcl_output_file, threads=args.THREADS)

        cluster_graph.get_contigs_not_in_mcl_cluster(mcl_output_file, original_contigs = unlabeled_contigs)
        unlabeled_representative_sequences = cluster_graph.calculate_representative_sequences()
        clusters_with_subcluster.append(cluster_graph.mcl_cluster)
        logger.debug(f'UNLABELED_CLUSTER: {cluster_graph.mcl_cluster}')

        clustered_sequence_names = []
        clustered_sequence_names = unlabeled_representative_sequences + cluster_representative_sequences

        logger.debug(f'clustered_sequence_names: {clustered_sequence_names}')
        logger.debug(f'clusters_with_subcluster: {clusters_with_subcluster}')

        logger.info('Salmon nice')


    else:

        # Perfom read mapping on all Contigs.
        contig_objs = perform_mappings(sequences, contig_dir, index_dir, args.THREADS)

        ####      ####
        ## Graphing ##
        ####      ####
        logger.info('Building graphs and redefine kmer clusters.')
        logger.debug(f'First group of labeled contigs {labeled_contigs[0]}')
        logger.debug(f'Unlabeled contigs: {unlabeled_contigs}')

        ## Create all cluster-graphs with the contig objects and remove the ones that don't have read overlap
        ## It is important that all unused contigs get removed before the graphs are updated.
        trimmed_clusters = []
        for i, cluster in enumerate(labeled_contigs, 1):
            relevant_contigs = [contig_objs[contig] for contig in cluster]
            cluster_graph = ReadGraph.from_contigs(relevant_contigs)

            
            # Draw
            draw_filename = f'{graph_visual_dir}/cluster_{i}_A.svg'
            cluster_graph.draw_graph(draw_filename)

            contigs_to_remove = cluster_graph.get_unconnected_nodes()
            cluster_graph.remove_nodes(contigs_to_remove)

            remaining_contigs = [contig_objs[contig] for contig in list(set(cluster).difference(set(contigs_to_remove)))]
            cluster_graph.set_original_contigs(remaining_contigs)
            
            unlabeled_contigs += contigs_to_remove
            trimmed_clusters.append(cluster_graph)


        ## Update the graphs with the unlabeled group.
        ## Remove nodes/contigs that:
        ## 1. Are not connected to any other node.
        ## 2. Are not but remove all contigs that are in a mcl-cluste that doesn't contain any original sequences.

        final_cluster = []
        cluster_representative_sequences = []
        clusters_with_subcluster = []


        for i, cluster_graph in enumerate(trimmed_clusters, 1):
            # get unlabeled contig objects
            relevant_contigs = [contig_objs[contig] for contig in unlabeled_contigs]
            logger.debug(f"Unlabeled contigs: {unlabeled_contigs}")
            cluster_graph.update_graph(relevant_contigs)

            abc_file = f'{graph_dir}/graph_{i}'
            cluster_graph.save_graph(abc_file)
            mcl_output_file = f'{mcl_dir}/graph_{i}'
            mcl = Mcl()
            mcl.run(graph_file = abc_file, output_file = mcl_output_file, threads=args.THREADS)

            # Extract the groups from mcl that contain original sequences.
            unconnected_contigs = cluster_graph.get_unconnected_nodes()
            non_mcl_cluster_contigs = cluster_graph.get_contigs_not_in_mcl_cluster(mcl_output_file)
            cluster_graph.remove_nodes(non_mcl_cluster_contigs)
            cluster_graph.remove_nodes(unconnected_contigs)

            unlabeled_contigs = unconnected_contigs + non_mcl_cluster_contigs


            logger.debug(f'2.Unlabeled contigs: {unlabeled_contigs}')

            representative_sequences = cluster_graph.calculate_representative_sequences()
            cluster_representative_sequences += representative_sequences
            logger.debug(f'MCL CLUSTER: {cluster_graph.mcl_cluster}')
            clusters_with_subcluster.append(cluster_graph.mcl_cluster)
            # Draw
            draw_filename = f'{graph_visual_dir}/cluster_{i}_B.svg'
            cluster_graph.draw_graph(draw_filename)

            final_cluster.append(cluster_graph.get_nodes())

        logger.info(f'Clusters are : {final_cluster}')
        
        ## Create graph and perform mcl for the unlabeled group.

        logger.debug(f'1: unlabeled contigs {unlabeled_contigs}')

        relevant_contigs = [contig_objs[contig] for contig in unlabeled_contigs]
        cluster_graph = ReadGraph.from_contigs(relevant_contigs)

        draw_filename = f'{graph_visual_dir}/cluster_unlabeled.svg'
        cluster_graph.draw_graph(draw_filename)
                
        abc_file = f'{graph_dir}/graph_unlabeled'
        cluster_graph.save_graph(abc_file)
        mcl_output_file = f'{mcl_dir}/graph_unlabeled'
        mcl = Mcl()
        mcl.run(graph_file = abc_file, output_file = mcl_output_file, threads=args.THREADS)

        cluster_graph.get_contigs_not_in_mcl_cluster(mcl_output_file, original_contigs = unlabeled_contigs)
        unlabeled_representative_sequences = cluster_graph.calculate_representative_sequences()

        clusters_with_subcluster.append(cluster_graph.mcl_cluster)
        logger.debug(f'UNLABELED_CLUSTER: {cluster_graph.mcl_cluster}')

        clustered_sequence_names = []
        clustered_sequence_names = unlabeled_representative_sequences + cluster_representative_sequences

        logger.debug(f'clustered_sequence_names: {clustered_sequence_names}')
        logger.debug(f'clusters_with_subcluster: {clusters_with_subcluster}')

    ####    ####
    ## Output ##
    ####    ####

    # Write fasta file
    fasta_output_file = f'{output_dir}/{BASENAME}.fa'
    sequences_to_write = {name: sequences[name] for name in clustered_sequence_names}
    writer = Writer(fasta_output_file)
    writer.write_fasta(sequences_to_write)

    # Write clstr file
    clstr_output_file = f'{output_dir}/{BASENAME}.clstr'
    writer = Writer(clstr_output_file)
    writer.write_clstr(clusters_with_subcluster, clustered_sequence_names)

    reduction_percent = ( len(clustered_sequence_names) / len(sequences) ) * 100
    logger.info(f'Reduced to ~{ round(reduction_percent, 2) } % of original sequences.')

    ####    ####
    ## Dammit ##
    ####    ####
    if args.ANNOTATE:
        logger.info('Starting annotation.')

        labels_after_kmer_clustering = k.unlabeled_cluster + k.clusters
        labels_after_read_graph = clusters_with_subcluster    

        logger.info(f'LABELS AFTER KMER CLUSTERING: {labels_after_kmer_clustering}')
        logger.info(f'LABELS AFTER READ GRAPH     : {labels_after_read_graph}')

        dammit = Dammit(first_file = args.FASTA_FILE, second_file = fasta_output_file, output_dir = dammit_dir, kmer_clusters = labels_after_kmer_clustering, karma_clusters = labels_after_read_graph, representative_sequences = clustered_sequence_names, threads = args.THREADS)
        dammit.update_database(database_dir=args.DATABASE_DIR, busco_group=args.BUSCO_GROUP)
        dammit.run()
        dammit.postprocessing()

        eval_output_file = f'{output_dir}/{BASENAME}.csv'
        dammit.save(eval_output_file)
