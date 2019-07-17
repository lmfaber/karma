import glob
import os
import subprocess
from cmd import Cmd

from cmd_parser import args
from annotate import Dammit
from hisat2 import Hisat2
from kmer import KmerClustering
from logs import logger
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

if __name__ == "__main__":
    
    BASENAME = basename_woe(args.FASTA_FILE)

    logger.info('Starting kmer based clustering.')
    k = KmerClustering(args.FASTA_FILE, threads = args.THREADS, kmer_size=args.KMER_SIZE)
    k.run()

    sequences = k.sequences
    labeled_contigs = k.clusters
    unlabeled_contigs = k.unlabeled_cluster[0]


    # Split all contigs in separate file
    # Perfom read mapping on all Contigs
    contig_dir = f'{BASENAME}/contigs'
    save_contigs_in_separate_files(sequences, contig_dir)

    ####     ####
    ## Mapping ##
    ####     ####
    if args.R:
        READS = [args.R]
    else:
        READS = [args.R1, args.R2]

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

    ####      ####
    ## Graphing ##
    ####      ####
    logger.info('Build graphs and redefine groups.')
    logger.debug(f'First group of labeled contigs {labeled_contigs[0]}')
    logger.debug(f'Unlabeled contigs: {unlabeled_contigs}')


    graph_dir = f'{BASENAME}/graphs'
    makedir(graph_dir)
    mcl_dir = f'{BASENAME}/mcl'
    makedir(mcl_dir)

    clustered_sequence_names = []
    cluster_with_subcluster = []
    # Create graph for all Clusters from KmerClustering including the unlabeled cluster
    for i, cluster in enumerate(labeled_contigs, 1):

        cluster_graph = ReadGraph(labeled=cluster, unlabeled=unlabeled_contigs, mapping_dir = mapping_dir, no=i, graph_dir=graph_dir, mcl_dir=mcl_dir, draw=True)

        unlabeled_contigs = cluster_graph.unlabeled
        representative_sequences = cluster_graph.representative_sequences

        cluster_with_subcluster.append(cluster_graph.mcl_cluster)
        clustered_sequence_names += representative_sequences
        logger.debug(representative_sequences)
    # Join unlabeled sequences into the clustered sequences:
    clustered_sequence_names += unlabeled_contigs
    for seq in unlabeled_contigs:
        cluster_with_subcluster.append(seq)
    
    ####    ####
    ## Output ##
    ####    ####
    output_dir = f'{BASENAME}/karma'
    makedir(output_dir)

    # Write fasta file
    fasta_output_file = f'{output_dir}/{BASENAME}.fa'
    sequences_to_write = {name: sequences[name] for name in clustered_sequence_names}
    writer = Writer(fasta_output_file)
    writer.write_fasta(sequences_to_write)

    # Write clstr file
    clstr_output_file = f'{output_dir}/{BASENAME}.clstr'
    writer = Writer(clstr_output_file)
    writer.write_clstr(cluster_with_subcluster, clustered_sequence_names)

    reduction_percent = ( len(clustered_sequence_names) / len(sequences) ) * 100
    logger.info(f'Reduced to { reduction_percent } % of original sequences.')


    ####    ####
    ## Dammit ##
    ####    ####
    logger.info('Perform de novo annotation.')
    dammit_dir = f'{BASENAME}/dammit'
    makedir(dammit_dir)
    dammit = Dammit(first_file = args.FASTA_FILE, second_file = fasta_output_file, output_dir = dammit_dir, database_dir = args.DATABASE_DIR, busco_group=args.BUSCO_GROUP)
