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

if __name__ == "__main__":
        
    BASENAME = basename_woe(args.FASTA_FILE)
    sequences = read_fasta_file(args.FASTA_FILE)

    ####             ####
    ## output folders ###
    ####             ####
    contig_dir = f'{args.OUTPUT_DIR}{BASENAME}/contigs'          # Saved single contigs
    kmer_dir = f'{args.OUTPUT_DIR}{BASENAME}/kmer'               # Kmer based clustering
    mapping_dir = f'{args.OUTPUT_DIR}{BASENAME}/hisat2/mapping'  # Sam files from mapping
    index_dir = f'{args.OUTPUT_DIR}{BASENAME}/hisat2/index'      # Index files for mapping
    graph_dir = f'{args.OUTPUT_DIR}{BASENAME}/graphs'            # abc-graph
    graph_visual_dir = f'{graph_dir}/visual'                     # plotted graphs
    mcl_dir = f'{args.OUTPUT_DIR}{BASENAME}/mcl'                 # mcl clustering results
    output_dir = f'{BASENAME}/karma'                             # actual results
    dammit_dir = f'{args.OUTPUT_DIR}{BASENAME}/dammit'           # de novo annotation results

    directories = [kmer_dir, mapping_dir, index_dir, graph_dir, graph_visual_dir, mcl_dir, output_dir, dammit_dir]
    for directory in directories:
        makedir(directory)

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

    # Split all contigs in separate file. Perfom read mapping on all Contigs.

    save_contigs_in_separate_files(sequences, contig_dir)

    ####     ####
    ## Mapping ##
    ####     ####
    if args.R:
        READS = [args.R]
    else:
        READS = [args.R1, args.R2]

    logger.info('Mapping all contigs')
    # Create all contig objects:
    all_contigs = labeled_contigs + [unlabeled_contigs]

    contig_objs = {}
    for cluster in all_contigs:
        # Map the the reads to the contigs in each cluster
        logger.debug(f'CLUSTER: {cluster}')

        logger.debug('Mapping to contigs')
        for contig_name in cluster:
            # Create contig object
            contig_pickle_file = f'{contig_dir}/{contig_name}.pickle'
            if os.path.isfile(contig_pickle_file):
                logger.debug('Reading pickle file.')
                with open(contig_pickle_file, 'rb') as reader:
                    contig = pickle.load(reader)
            else:

                contig_file = f'{contig_dir}/{contig_name}.fa'
                index_name = f'{index_dir}/{contig_name}'
                output_sam = f'{mapping_dir}/{contig_name}.sam'
                hisat2 = Hisat2(input_file=contig_file, output_file=output_sam, index_name=index_name, threads=args.THREADS)
                hisat2.build_index()
                hisat2.run(READS)

                logger.debug('Creating new contig from sam files')
                contig = Contig(contig_name)
                contig.load_contig_info_from_sam(output_sam)
                with open(contig_pickle_file, 'wb') as writer:
                    pickle.dump(contig, writer)

                # Remove index and sam files. Don't remove index files because it will just 
                if args.KEEP_SAM == False:
                    os.remove(output_sam)
                    # for f in glob.glob(f'{index_name}.*.ht2'):
                    #     os.remove(f)

            contig_objs[contig_name] = contig


    ####      ####
    ## Graphing ##
    ####      ####
    logger.info('Building graphs and redefine kmer clusters.')
    logger.debug(f'First group of labeled contigs {labeled_contigs[0]}')
    logger.debug(f'Unlabeled contigs: {unlabeled_contigs}')

    ## Create all cluster-graphs with the contig objects and remove the ones that don't have read overlap
    ## It is important that all unused contigs get removed before the graphs are updated.
    trimmed_labeled_contigs = []
    for i, cluster in enumerate(labeled_contigs, 1):
        relevant_contigs = []
        for contig in cluster:
            relevant_contigs.append(contig_objs[contig])
        
        cluster_graph = ReadGraph()
        removed_contigs = cluster_graph.create_graph(relevant_contigs)
        
        # Draw
        draw_filename = f'{graph_visual_dir}/cluster_{i}_A.svg'
        cluster_graph.draw_graph(draw_filename)

        cluster_graph.trim()

        unlabeled_contigs += removed_contigs
        trimmed_labeled_contigs.append(cluster_graph)

    ## Update the graphs with the unlabeled group.
    ## Remove nodes/contigs that:
    ## 1. Are not connected to any other node.
    ## 2. Are not but remove all contigs that are in a mcl-cluste that doesn't contain any original sequences.

    final_cluster = []
    cluster_representative_sequences = []
    clusters_with_subcluster = []
    for i, cluster in enumerate(trimmed_labeled_contigs, 1):

        # get unlabeled contig objects
        relevant_contigs = []
        logger.debug(f"Unlabeled contigs: {unlabeled_contigs}")
        for contig in unlabeled_contigs:
            relevant_contigs.append(contig_objs[contig])

        cluster.update_graph(relevant_contigs)

        abc_file = f'{graph_dir}/graph_{i}'
        cluster.save_graph(abc_file)

        mcl_output_file = f'{mcl_dir}/graph_{i}'
        mcl = Mcl()
        mcl.run(graph_file = abc_file, output_file = mcl_output_file, threads=args.THREADS)

        # Extract the groups from mcl that contain original sequences.
        not_connected = cluster.extract_mcl_clusters(mcl_output_file)
        cluster.remove_contigs(not_connected)
        trimmed = cluster.trim()

        unlabeled_contigs = not_connected + trimmed


        logger.debug(f'2.Unlabeled contigs: {unlabeled_contigs}')

        representative_sequences = cluster.calculate_representative_sequences()
        cluster_representative_sequences += representative_sequences
        logger.debug(f'MCL CLUSTER: {cluster.mcl_cluster}')
        clusters_with_subcluster.append(cluster.mcl_cluster)
        # Draw
        draw_filename = f'{graph_visual_dir}/cluster_{i}_B.svg'
        cluster.draw_graph(draw_filename)

        final_cluster.append(cluster.get_nodes())

    logger.info(f'Clusters are : {final_cluster}')
    

    ## Create graph and perform mcl for the unlabeled group.
    relevant_contigs = []
    logger.debug(f'1: unlabeled contigs {unlabeled_contigs}')
    for contig in unlabeled_contigs:
        relevant_contigs.append(contig_objs[contig])

        cluster_graph = ReadGraph()
        cluster_graph.create_graph(relevant_contigs)

        draw_filename = f'{graph_visual_dir}/cluster_unlabeled.svg'
        cluster_graph.draw_graph(draw_filename)
                
        abc_file = f'{graph_dir}/graph_unlabeled'
        cluster_graph.save_graph(abc_file)

        mcl_output_file = f'{mcl_dir}/graph_unlabeled'
        mcl = Mcl()
        mcl.run(graph_file = abc_file, output_file = mcl_output_file, threads=args.THREADS)

        cluster_graph.extract_mcl_clusters(mcl_output_file, original_contigs = unlabeled_contigs)
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
    logger.info('Starting annotation.')

    labels_after_kmer_clustering = k.unlabeled_cluster + k.clusters
    labels_after_read_graph = clusters_with_subcluster    

    logger.info(f'LABELS AFTER KMER CLUSTERING: {labels_after_kmer_clustering}')
    logger.info(f'LABELS AFTER READ GRAPH     : {labels_after_read_graph}')

    dammit = Dammit(first_file = args.FASTA_FILE, second_file = fasta_output_file, output_dir = dammit_dir, database_dir = args.DATABASE_DIR, busco_group=args.BUSCO_GROUP, kmer_clusters = labels_after_kmer_clustering, karma_clusters = labels_after_read_graph, representative_sequences = clustered_sequence_names, threads = args.THREADS)
    dammit.postprocessing()

    eval_output_file = f'{output_dir}/{BASENAME}.csv'
    dammit.save(eval_output_file)


