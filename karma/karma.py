import os
import glob 

from logs import logger
from cmd_parser import args
from kmer import KmerClustering
from cmd import Cmd
from hisat2 import Hisat2
from read_graph import ReadGraph

import subprocess


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

if __name__ == "__main__":
    
    logger.info('Starting kmer based clustering.')
    k = KmerClustering(args.FASTA_FILE, threads = args.THREADS, kmer_size=args.KMER_SIZE)
    k.run()

    sequences = k.sequences
    labeled_contigs = k.clusters
    unlabeled_contigs = k.unlabeled_cluster[0]

    # Split all contigs in separate file
    # Perfom read mapping on all Contigs
    if args.R:
        READS = [args.R]
    else:
        READS = [args.R1, args.R2]
    
    BASENAME = basename_woe(args.FASTA_FILE)

    contig_dir = f'{BASENAME}/contigs'
    save_contigs_in_separate_files(sequences, contig_dir)

    ## Mapping
    logger.info('Lean back and relax. Mapping could take a while.')
    mapping_dir = f'{BASENAME}/hisat2/mapping'
    index_dir = f'{BASENAME}/hisat2/index'

    makedir(mapping_dir)
    makedir(index_dir)

    for contig_file in glob.glob(f'{contig_dir}/*.fa'):
        contig_basename = basename_woe(contig_file)
        index_name = f'{index_dir}/{contig_basename}'
        output_sam = f'{mapping_dir}/{contig_basename}.sam'
        hisat2 = Hisat2(input_file=contig_file, output_file=output_sam, index_name=index_name, threads=args.THREADS)
        hisat2.build_index()
        hisat2.run(READS)

    ## Graphing
    logger.info('Build graphs and redefine groups.')


    logger.debug(labeled_contigs[0])
    logger.debug(unlabeled_contigs)


    graph_dir = f'{BASENAME}/graphs'
    makedir(graph_dir)
    mcl_dir = f'{BASENAME}/mcl'
    makedir(mcl_dir)

    clustered_sequence_names = []
    # Create graph for all Clusters from KmerClustering including the unlabeled cluster
    for i, cluster in enumerate(labeled_contigs):

        cluster_graph = ReadGraph(labeled=cluster, unlabeled=unlabeled_contigs, mapping_dir = mapping_dir, no=i, graph_dir=graph_dir, mcl_dir=mcl_dir)

        #cluster_graph.draw_graph()
        cluster_graph.mcl_clustering()
        unlabeled_contigs = cluster_graph.extract_cluster()
        INFOS = cluster_graph.calculate_representative_sequences()
        clustered_sequence_names.append(INFOS)
        logger.debug(INFOS)
    
    reduction_percent = ( len(clustered_sequence_names) / len(sequences) ) * 100
    logger.info(f'Reduced to { reduction_percent } % of original sequences.')

