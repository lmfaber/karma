import os
import sys
from collections import Counter, OrderedDict
from multiprocessing import Pool

import hdbscan
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import umap
from logs import logger


def cdhit_like_output(sequence_ids, sequences, cluster_ids, probabilities, filename='cluster'):
    """    Creates a cd-hit-est like output for the clusters from umap and hdbscan

    Arguments:
        sequence_ids {[type]} -- [description]
        cluster {[type]} -- [description]

    Keyword Arguments:
        filename {str} -- Output filename (default: {'cluster'})
    """

    with open(f'{filename}.clstr', 'w') as writer, open(f'{filename}_no_classfication.clstr', 'w') as no_cluster_writer:
        unique_cluster = sorted(set(cluster_ids))
        logger.debug(f'Write {len(unique_cluster) - 1} cluster.')
        for cluster_id, cl in enumerate(unique_cluster):
            internal_cluster_counter = 0
            if cl == -1:
                no_cluster_writer.write(f'>Cluster {cluster_id}\n')
            else:
                writer.write(f'>Cluster {cluster_id}\n')
            clustered_sequences = []
            for sequence_id, sequence, cluster_no, probability in zip(sequence_ids, sequences, cluster_ids, probabilities):
                if cluster_no == cl:
                    clustered_sequences.append( (len(sequence), sequence_id, round(probability*100.0, 2)) )
            sorted_by_length = sorted(clustered_sequences, key = lambda x: x[0])

            for sequence_length, sequence_id, probability in sorted_by_length:
                # Either write to no cluster file or the cluster file
                if cl == -1:
                    no_cluster_writer.write(f'{internal_cluster_counter}\t{sequence_length}nt, {sequence_id}... at {probability}%\n')
                else:
                    writer.write(f'{internal_cluster_counter}\t{sequence_length}nt, {sequence_id}... at {probability}%\n')
                internal_cluster_counter += 1

def kmers_of_seq(sequence, kmer_size, palindromic=False):
    """
    Returns all kmers of length n from a given sequence.
    
    Arguments:
        sequence {str} -- Input sequence.
    
    Keyword Arguments:
        kmer_size {int} -- Kmer size
    """
    for i in range(len(sequence)-kmer_size+1):
        if not palindromic:
            yield(sequence[i:i+kmer_size])
        else:
            tmp_seq = sequence[i:i+kmer_size]
            if is_palindrome(tmp_seq):
                yield(tmp_seq)

def is_palindrome(sequence):
    """
    Checks if sequence in palindrome.
    
    Arguments:
        sequence {str} -- Input string.
    """
    return sequence == sequence[::-1]

def complement_sequence(sequence):
    """Complements a DNA sequence
    
    Arguments:
        sequence {str} -- DNA string
    """
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return( ''.join([complement_dict[N] for N in sequence]) )

def reverse_complement_sequence(sequence):
    """Reverse complements a DNA sequence
    
    Arguments:
        sequence {str} -- DNA string
    """
    compl = complement_sequence(sequence)
    return( compl[::-1] )

def count_kmer_occurence(sequences, kmer_size):
    """Counts actually appearing kmers in list of sequences and returns a list of Counter objects.
    
    Arguments:
        sequences {dictionary} -- Dictionary where the sequences are the values.
    """
    logger.debug('Count ALL kmers per sequence...')

    # Save them in a ordered way, otherwise the clustering wont work.
    counts = []
    for sequence in sequences.values():
        c = Counter()

        if kmer_size == '5p6':
            # 5-mers
            k = 5
            for i in range(len(sequence)-k+1):
                c[sequence[i:i+k]] += 1

            # 6p-mers
            k = 6
            for i in range(len(sequence)-k+1):
                tmp_seq = sequence[i:i+k]
                if is_palindrome(tmp_seq):
                    c[sequence[i:i+k]] += 1
            counts.append(c)

        else:
            for i in range(len(sequence)-kmer_size+1):
                c[sequence[i:i+kmer_size]] += 1
            counts.append(c)

    assert len(counts) == len(sequences), 'Something went wrong while counting kmers. Exiting...'
    
    return counts 

def mask_list(list_to_mask, mask):

    labeled = []
    unlabeled = []
    for i in set(mask):
        positions = (j for j, k in enumerate(mask) if k == i)
        if i == -1:
            unlabeled.append( [list(list_to_mask.keys())[c].lstrip('>') for c in positions] )
        else:
            labeled.append( [list(list_to_mask.keys())[c].lstrip('>') for c in positions] )
    return(labeled, unlabeled)

class KmerClustering():

    def __init__(self, sequences, output_file, kmer_size, threads = 2,):
        self.sequences = sequences
        self.output_file = output_file
        self.threads = threads
        self.kmer_size = kmer_size

        self.clusters = []
        self.unlabeled_cluster = []
        self.kmers = None
        self.sorted_kmer_set = set()

        if os.path.isfile(self.output_file):
            logger.info(f'Read from previous calculation: {self.output_file}')
            self.read_clusters()
        else:
            self.run()
        self.__fix_fasta_headers()
        
    def __fix_fasta_headers(self):
        logger.debug(f'unlabeled: {self.unlabeled_cluster}')
        self.unlabeled_cluster = [[name.split(' ')[0] for name in self.unlabeled_cluster[0]]]
        renamed_cluster = []
        for cluster in self.clusters:
            renamed_cluster.append([name.split(' ')[0] for name in cluster])
        assert len(renamed_cluster) == len(self.clusters), 'Something went wrong while renaming fasta header. See if removing spaces from fasta headers solves the problem.'
        self.clusters = renamed_cluster

    def fill_array_for_contig(self, *args):
        """
        Calculates all positions in the array that have to be changed. Returns a tuple of three. (row, col, no. of kmer occurrence)
        """
        row = args[0]
        contig = args[1] # Counter object
        length = args[2]

        posCounts = []
        for kmer, count in contig.items():
            col = self.kmers[kmer]
            # Normalize kmer appearance to length of sequence
            # posCounts.append((row, col, count/length))
            posCounts.append((row, col, count))
        return(posCounts)

    def save_groups_to_file(self):
        """
        Saves the names of contigs from labeled clusters and unlabeled clusters into one file. Names are tab separated.
        The first line contains the unlabeled contigs.
        """

        with open(self.output_file, 'w') as cluster_writer:
            cluster_writer.write('\t'.join(self.unlabeled_cluster[0]) + '\n')
            for cluster in self.clusters:
                cluster_writer.write('\t'.join(cluster) + '\n')

    def read_clusters(self):
        """
        Reads existing cluster from previous runs.
        """
        with open(self.output_file, 'r') as cluster_reader:
            self.unlabeled_cluster = [cluster_reader.readline().rstrip('\n').split('\t')]
            for line in cluster_reader:
                self.clusters.append(line.rstrip('\n').split('\t'))

    def extract_kmers(self, sequences, kmer_size):
        """Extracts all kmers from a list of sequences.
        
        Arguments:
            sequences {dictionary} -- Receives a dictionary with values as sequences.
        """
        logger.info('Extracting kmers from contigs.')
        kmer_set = set()

        if kmer_size == '5p6':
            # Add 5-mers
            for sequence in sequences.values():
                for k in kmers_of_seq(sequence, 5, palindromic=False):
                    kmer_set.add(k)

            # Add palindromic 6-mers
            for sequence in sequences.values():
                for k in kmers_of_seq(sequence, 6, palindromic=True):
                    kmer_set.add(k)

        else:
            logger.info(f'Accounting only {kmer_size}-mers.')
            for sequence in sequences.values():
                for k in kmers_of_seq(sequence, kmer_size, palindromic=False):
                    kmer_set.add(k)


        # logger.debug(f'No of {KMER_SIZE}-mers: {len(kmer_set)}')
        self.sorted_kmer_set = sorted(kmer_set)
        kmer_set.clear()

        kmer_dict = {}
        for i, kmer in enumerate(self.sorted_kmer_set):
            kmer_dict[kmer] = i
        logger.debug(f'Dict has {len(kmer_dict)} entries.')
        return( kmer_dict )

    def run(self):
        ## 2. Extract kmers from sequences. Using sets. Dictionary! Used in fill_array for  function
        self.kmers = self.extract_kmers(self.sequences, kmer_size=self.kmer_size)

        # Initialize empty numpy array:
        logger.debug('Initialize empty numpy array...')
        kmer_profile = np.zeros(shape=(len(self.sequences), len(self.kmers)), dtype=np.float)

        ## 3. Count all accourences of kmers in each contig
        counts = count_kmer_occurence(self.sequences, kmer_size=self.kmer_size)

        lengths_of_sequences = (len(seq) for seq in self.sequences)

        logger.debug(f'Calculate position and counts with {self.threads} cores.')
        pcount = []

        with Pool(self.threads) as pool:
            for i in pool.starmap(self.fill_array_for_contig, ((row, countObj, length) for row, (countObj, length) in enumerate(zip(counts, lengths_of_sequences)) ) ):
                pcount += i

        logger.debug(f'Size of pcount: {sys.getsizeof(pcount)}')

        logger.info('Change entries in numpy array...')
        for i in pcount:
            row, col, count =  i[0], i[1], i[2]
            kmer_profile[row][col] = count

        # Sanity check: All colums should at least have one non zero value per column
        for n in range(kmer_profile.shape[1]):
            try:
                assert np.any(kmer_profile[:,n]), 'At least one column has all zero values...'
            except AssertionError:
                logger.error(f'Values of column {n} are all zero, which should not be the case.')
                logger.debug(f'Problematic kmer: {self.sorted_kmer_set[n]}\n{kmer_profile[:,n]}')
                exit(1)
        self.sorted_kmer_set.clear()

        logger.debug(f'KMER-PROFILE - Size: {sys.getsizeof(kmer_profile)}, Shape: {kmer_profile.shape}')


        ## 4. dimension reduction with umap
        logger.info('Dimension reduction with UMAP.')
        neighbors = 2
        components = 10
        dist = 0
        reducer = umap.UMAP(n_neighbors=neighbors, n_components=components, min_dist=dist, random_state=42).fit_transform(kmer_profile)

        ## 5. HDBSCAN
        min_cluster_size = 2
        logger.info(f'Clustering with HDBSCAN.\n\tMin cluster size: {min_cluster_size}')
        clusterer = hdbscan.HDBSCAN(min_cluster_size = min_cluster_size).fit(reducer)

        self.clusters, self.unlabeled_cluster = mask_list(self.sequences, clusterer.labels_)
        logger.debug(f'Cluster labels: {clusterer.labels_}')

        ## Save to file 
        self.save_groups_to_file()


        # bname = os.path.splitext(os.path.basename(inFile))[0]
        # f = f'{bname}_allow_single_cluster_k{KMER_SIZE}_min{min_cluster_size}'
        # cdhit_like_output(sequence_names, sequences, clusterer.labels_, clusterer.probabilities_, filename=f)
