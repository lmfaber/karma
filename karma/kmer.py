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


class KmerClustering:
    def __init__(self, sequences, output_dir, kmer_size, threads):
        self.sequences = sequences
        self.output_dir = output_dir
        self.output_eval = f"{self.output_dir}/eval.txt"
        self.output_file = f"{self.output_dir}/cluster.txt"

        self.threads = threads
        self.kmer_size = kmer_size

        self.clusters = []
        self.unlabeled_cluster = []
        self.kmers = None
        self.sorted_kmer_set = set()

    @staticmethod
    def __mask_list(list_to_mask, mask):
        """ Return only specifid items from a list according to a mask. """
        labeled = []
        unlabeled = []
        for i in set(mask):
            positions = (j for j, k in enumerate(mask) if k == i)
            if i == -1:
                unlabeled.append(
                    [list(list_to_mask.keys())[c].lstrip(">") for c in positions]
                )
            else:
                labeled.append(
                    [list(list_to_mask.keys())[c].lstrip(">") for c in positions]
                )
        return (labeled, unlabeled)

    @staticmethod
    def __is_palindrome(sequence):
        """
        Checks if sequence in palindrome.

        Arguments:
            sequence {str} -- Input string.
        """
        return sequence == sequence[::-1]

    def __count_kmer_occurence(self, sequences, kmer_size):
        """Counts actually appearing kmers in list of sequences and returns a list of Counter objects.

        Arguments:
            sequences {dictionary} -- Dictionary where the sequences are the values.
        """
        logger.debug("Count ALL kmers per sequence...")

        # Save them in a ordered way, otherwise the clustering wont work.
        counts = []
        for sequence in sequences.values():
            c = Counter()

            if kmer_size == "5p6":
                # 5-mers
                k = 5
                for i in range(len(sequence) - k + 1):
                    c[sequence[i : i + k]] += 1

                # 6p-mers
                k = 6
                for i in range(len(sequence) - k + 1):
                    tmp_seq = sequence[i : i + k]
                    if self.__is_palindrome(tmp_seq):
                        c[sequence[i : i + k]] += 1
                counts.append(c)

            else:
                for i in range(len(sequence) - kmer_size + 1):
                    c[sequence[i : i + kmer_size]] += 1
                counts.append(c)

        assert len(counts) == len(
            sequences
        ), "Something went wrong while counting kmers. Exiting..."

        return counts

    def __fix_fasta_headers(self):
        """ Renames fasta header. Splits the sequence name at the first space. """
        logger.debug(f"unlabeled: {self.unlabeled_cluster}")
        self.unlabeled_cluster = [
            [name.split(" ")[0] for name in self.unlabeled_cluster[0]]
        ]
        renamed_cluster = []
        for cluster in self.clusters:
            renamed_cluster.append([name.split(" ")[0] for name in cluster])
        assert len(renamed_cluster) == len(
            self.clusters
        ), "Something went wrong while renaming fasta header. See if removing spaces from fasta headers solves the problem."
        self.clusters = renamed_cluster

    def fill_array_for_contig(self, *args):
        """
        Calculates all positions in the array that have to be changed. Returns a tuple of three. (row, col, no. of kmer occurrence)
        """
        row = args[0]
        contig = args[1]  # Counter object
        length = args[2]

        posCounts = []
        for kmer, count in contig.items():
            col = self.kmers[kmer]
            # Normalize kmer appearance to length of sequence
            posCounts.append((row, col, count / length))
            # posCounts.append((row, col, count))
        return posCounts

    def __save_groups_to_file(self):
        """
        Saves the names of contigs from labeled clusters and unlabeled clusters into one file. Names are tab separated.
        The first line contains the unlabeled contigs.
        """

        with open(self.output_file, "w") as cluster_writer:
            cluster_writer.write("\t".join(self.unlabeled_cluster[0]) + "\n")
            for cluster in self.clusters:
                cluster_writer.write("\t".join(cluster) + "\n")

    def __read_clusters(self):
        """
        Reads existing cluster from previous runs.
        """
        with open(self.output_file, "r") as cluster_reader:
            self.unlabeled_cluster = [
                cluster_reader.readline().rstrip("\n").split("\t")
            ]
            for line in cluster_reader:
                self.clusters.append(line.rstrip("\n").split("\t"))

    def __extract_kmers(self, sequences, kmer_size):
        """Extracts all kmers from a list of sequences.

        Arguments:
            sequences {dictionary} -- Receives a dictionary with values as sequences.
        """
        logger.info("Extracting kmers from contigs.")
        kmer_set = set()

        if kmer_size == "5p6":
            # Add 5-mers
            for sequence in sequences.values():
                for k in self.__kmers_of_seq(sequence, 5, palindromic=False):
                    kmer_set.add(k)

            # Add palindromic 6-mers
            for sequence in sequences.values():
                for k in self.__kmers_of_seq(sequence, 6, palindromic=True):
                    kmer_set.add(k)

        else:
            logger.info(f"Accounting only {kmer_size}-mers.")
            for sequence in sequences.values():
                for k in self.__kmers_of_seq(sequence, kmer_size, palindromic=False):
                    kmer_set.add(k)

        self.sorted_kmer_set = sorted(kmer_set)
        kmer_set.clear()

        kmer_dict = {}
        for i, kmer in enumerate(self.sorted_kmer_set):
            kmer_dict[kmer] = i
        logger.debug(f"Dict has {len(kmer_dict)} entries.")
        return kmer_dict

    def __kmers_of_seq(self, sequence, kmer_size, palindromic=False):
        """
        Returns all kmers of length n from a given sequence.

        Arguments:
            sequence {str} -- Input sequence.

        Keyword Arguments:
            kmer_size {int} -- Kmer size
        """
        for i in range(len(sequence) - kmer_size + 1):
            if palindromic:
                tmp_seq = sequence[i : i + kmer_size]
                if self.__is_palindrome(tmp_seq):
                    yield (tmp_seq)
            else:
                yield (sequence[i : i + kmer_size])

    def __calc_kmer_profile(self):
        """ Calculate kmer profile. """
        ## 2. Extract kmers from sequences. Using sets. Dictionary! Used in fill_array for  function
        self.kmers = self.__extract_kmers(self.sequences, kmer_size=self.kmer_size)

        # Initialize empty numpy array:
        logger.debug("Initialize empty numpy array...")
        kmer_profile = np.zeros(
            shape=(len(self.sequences), len(self.kmers)), dtype=np.float
        )

        ## 3. Count all accourences of kmers in each contig
        counts = self.__count_kmer_occurence(self.sequences, kmer_size=self.kmer_size)

        lengths_of_sequences = (len(seq) for seq in self.sequences)

        logger.debug(f"Calculate position and counts with {self.threads} cores.")
        pcount = []

        with Pool(self.threads) as pool:
            for i in pool.starmap(
                self.fill_array_for_contig,
                (
                    (row, countObj, length)
                    for row, (countObj, length) in enumerate(
                        zip(counts, lengths_of_sequences)
                    )
                ),
            ):
                pcount += i

        logger.debug("Change entries in numpy array...")
        for i in pcount:
            row, col, count = i[0], i[1], i[2]
            kmer_profile[row][col] = count

        # Sanity check: All colums should at least have one non zero value per column
        for n in range(kmer_profile.shape[1]):
            try:
                assert np.any(
                    kmer_profile[:, n]
                ), "At least one column has all zero values..."
            except AssertionError:
                logger.error(
                    f"Values of column {n} are all zero, which should not be the case."
                )
                logger.debug(
                    f"Problematic kmer: {self.sorted_kmer_set[n]}\n{kmer_profile[:,n]}"
                )
                exit(1)

        # Sanity check: All rows should not be all zeros
        for n in range(kmer_profile.shape[0]):
            try:
                assert np.any(kmer_profile[n, :])
            except AssertionError:
                logger.error(
                    f"Values of row {n} are all zero, which should not be the case."
                )
                exit(1)
        self.sorted_kmer_set.clear()

        logger.debug(
            f"KMER-PROFILE - Size: {sys.getsizeof(kmer_profile)}, Shape: {kmer_profile.shape}"
        )
        return kmer_profile

    def __write_eval_information(self, **kwargs):
        """ Writes the parameters and the mean probability for the cluster results in the evaluation file. """
        with open(f"{self.output_eval}", "w") as writer:
            header = "\t".join(list(kwargs.keys()))
            writer.write(header + "\n")
            values = "\t".join([str(a) for a in list(kwargs.values())])
            writer.write(values + "\n")

    def run(self, neighbors, components, dist, r_state, min_cluster_size):
        if os.path.isfile(self.output_file):
            logger.info(f"Read from previous calculation: {self.output_file}")
            self.__read_clusters()
        else:
            ## 1. Calculate kmer profile for each sequence
            logger.info("Calculate kmer profiles.")
            kmer_profile = self.__calc_kmer_profile()

            ## 2. Reduce the dimension with UMAP
            logger.info("Dimension reduction with UMAP.")
            reducer = umap.UMAP(
                n_neighbors=neighbors,
                n_components=components,
                min_dist=dist,
                random_state=r_state,
            ).fit_transform(kmer_profile)

            ## 3. Perform clustering using HDBSCAN
            logger.info(
                f"Perform clustering with HDBSCAN. (min_cluster_size: {min_cluster_size})"
            )
            if min_cluster_size == 1:
                clusterer = hdbscan.HDBSCAN(allow_single_cluster=True).fit(reducer)
            else:
                clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size).fit(
                    reducer
                )

            self.clusters, self.unlabeled_cluster = self.__mask_list(
                self.sequences, clusterer.labels_
            )

            ## Save to file
            self.__save_groups_to_file()

            no_unlabeled = list(clusterer.labels_).count(-1)
            no_groups = max(clusterer.labels_) + 1
            mean_probability = np.mean(clusterer.probabilities_)
            self.__write_eval_information(
                kmer_size=self.kmer_size,
                n_neighbors=neighbors,
                n_components=components,
                min_dist=dist,
                random_state=r_state,
                min_cluster_size=min_cluster_size,
                unlabeled=no_unlabeled,
                no_groups=no_groups,
                mean_probability=mean_probability,
            )

        self.__fix_fasta_headers()
