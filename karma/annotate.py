import os
from cmd import Cmd

import numpy as np
import pandas as pd
from dammit.fileio.gff3 import GFF3Parser
from logs import logger


class Dammit:
    """ Perform de novo annotation with dammit. """

    def __init__(self, first_file, second_file, output_dir, kmer_clusters, karma_clusters, representative_sequences, threads = 6):
        self.transcriptomes = {'before': first_file, 'after': second_file}
        self.output_dir = output_dir
        self.kmer_clusters = kmer_clusters
        self.karma_clusters = karma_clusters
        self.threads = threads
        self.representative_sequences = set([a.lstrip('>') for a in representative_sequences])

        self.database_dir = ''
        self.busco_group = ''

        self.gff_files = {}
        self.namemaps = {}

    def update_database(self, database_dir, busco_group):
        """
        Updates the dammit database.
        """
        logger.info('Update dammit database.')
        self.database_dir = database_dir
        self.busco_group = busco_group
        database = Cmd(f'dammit databases --install --n_threads {self.threads} --database-dir {self.database_dir} --busco-group {self.busco_group}')
        database.run()

    def run(self):
        """
        Executes the dammit annotation for the original and reduced fasta file.
        """
        logger.info('Run dammit annotation.')
        for name, transcriptome in self.transcriptomes.items():

            output_dir = f'{self.output_dir}/{name}'
            annotation_file = f'{output_dir}/{os.path.basename(transcriptome)}.dammit.gff3'
            self.gff_files[name] = annotation_file
            namemap_file = f'{output_dir}/{os.path.basename(transcriptome)}.dammit.namemap.csv'
            self.namemaps[name] = namemap_file
            if not os.path.exists(annotation_file) and not os.path.exists(namemap_file):
                dammit = Cmd(f'dammit annotate {transcriptome} -o {output_dir} --database-dir {self.database_dir} --busco-group {self.busco_group} --n_threads {self.threads}')
                dammit.run()

    def postprocessing(self):
        """
        Calculates some basic information from the two annotations.
        """
        logger.debug('Postprocessing dammit output.')
        
        before = self.gff_files['before'] 
        namemap_before = self.namemaps['before'] 

        after = self.gff_files['after'] 
        namemap_after = self.namemaps['after']

        BEFORE = self.gene_dict(before, namemap_before, suffix='before')
        AFTER = self.gene_dict(after, namemap_after, suffix = 'after')

        logger.debug(f'GENE DICT BEFORE: {BEFORE}')
        logger.debug(f'GENE DICT AFTER : {AFTER}')


        self.lost_genes = self.calculate_lost_genes(BEFORE, AFTER)
        self.removed_sequences_per_gene = self.removed_sequences_per_gene(BEFORE, AFTER)
        clusters_per_gene = self.calculate_clusters_per_gene(BEFORE, AFTER)
        sequences_per_gene = self.calculate_sequences_per_gene(BEFORE, AFTER)
        self.clusters_and_sequences_per_gene = self.combine_dicts(clusters_per_gene, sequences_per_gene)

    def calculate_sequences_per_gene(self, gene_dict_before, gene_dict_after):
        """Calculates how many sequences are per gene in the assembly.
        
        Arguments:
            gene_dict_before {[type]} -- [description]
            gene_dict_after {[type]} -- [description]
        
        Returns:
            [type] -- [description]
        """
        A = {}
        for gene, transcripts in gene_dict_before.items():
            A[gene] = len(transcripts)

        B = {}
        for gene, transcripts in gene_dict_after.items():

            B[gene] = len(transcripts)
        
        # Combine the dicts
        combined_seqs_per_cluster = {}
        for gene, no_of_sequences in A.items():

            clusters = []
            clusters.append(no_of_sequences)

            if gene in B:
                clusters.append(B[gene])
            else:
                clusters.append(0)
            combined_seqs_per_cluster[gene] = clusters
        return combined_seqs_per_cluster

    def calculate_clusters_per_gene(self, gene_dict_before, gene_dict_after):
        """
        Calculates how many clusters were generated per gene.
        
        Arguments:
            gene_dict_before {dict} -- dictionary containing gene/cluster as key/value
        
        Returns:
            dict -- key/value: gene/no_of_clusters_before/number-of-cluster-after/no-of-seq-per-gene
        """

        clusters_per_gene_A = {}
        for gene, gene_sequences in gene_dict_before.items():
            # A
            cluster_counter_A = 0
            already_checked_clusters_A = []

            for cluster_gene_A in gene_sequences:
                for i, cluster in enumerate(self.kmer_clusters):
                    if i in already_checked_clusters_A:
                        break
                    if cluster_gene_A in cluster:
                        cluster_counter_A += 1
                        already_checked_clusters_A.append(i)
                clusters_per_gene_A[gene] = cluster_counter_A

        # Flatten the clusters from read graph based clustering
        flattened_clusters = []
        for i in self.karma_clusters:
            for j in i:
                flattened_clusters.append(j)
        logger.info(f'FLATTENED CLUSTER: {flattened_clusters}')

        clusters_per_gene_B = {}
        for gene, gene_sequences in gene_dict_after.items():
            cluster_counter_B = 0
            already_checked_clusters_B = []

            for cluster_gene_B in gene_sequences:
                #B(after)
                for i, cluster in enumerate(flattened_clusters):
                    if i in already_checked_clusters_B:
                        break
                    if cluster_gene_B in cluster:
                        cluster_counter_B += 1
                        already_checked_clusters_B.append(i)
                clusters_per_gene_B[gene] = cluster_counter_B

        logger.debug('CLUSTERS PER GENE A _ B')
        logger.debug(clusters_per_gene_A)
        logger.debug(clusters_per_gene_B)

        clusters_per_gene = {}
        for gene, no_of_cluster in clusters_per_gene_A.items():
            clusters = []
            clusters.append(no_of_cluster)
            if gene in clusters_per_gene_B:
                clusters.append(clusters_per_gene_B[gene])
            else:
                clusters.append(0)
            
            clusters_per_gene[gene] = clusters

        return clusters_per_gene

    def calculate_lost_genes(self, gene_dict_before, gene_dict_after):
        """Calculates the lost genes. Genes that were detected before kmer based clustering, but not after graph clustering.
        
        Arguments:
            gene_dict_before {dict} -- [description]
            gene_dict_after {dict} -- [description]
        
        Returns:
            tuple -- (number_of_lost_genes, lost_genes)
        """
        genes_detected_before_clustering = set(gene_dict_before.keys())
        genes_detected_after_clustering = set(gene_dict_after.keys())
        lost_genes = genes_detected_before_clustering.difference(genes_detected_after_clustering)
        number_of_lost_genes = len(lost_genes)
        if number_of_lost_genes == 0:
            lost_genes = '-'
        return (number_of_lost_genes, lost_genes)

    def removed_sequences_per_gene(self, gene_dict_before, gene_dict_after):
        """
        Calculates the min/max/mean of removed sequences per gene cluster. Compares before graph clustering and after.
        """
        removed_sequences = {'min': None,
                            'max': None,
                            'mean': None}

        removed_per_cluster = []
        for (key_b, value_b), (key_a, value_a) in zip(gene_dict_before.items(), gene_dict_after.items()):
            b = set(value_b)
            a = set(value_a)
            no_removed_sequences = len(b.difference(a))
            removed_per_cluster.append( no_removed_sequences ) 

        removed_sequences['min'] = min(removed_per_cluster)
        removed_sequences['max'] = max(removed_per_cluster)
        removed_sequences['mean'] = np.sum(removed_per_cluster)/len(gene_dict_before.keys())
        return removed_sequences

    def gene_dict(self, gff, namemap, suffix):
        """Creates a dictionary of genes/transcripts.
        
        Arguments:
            gff {str} -- gff file from dammit
            namemap {str} -- namemap file from dammit
        """

        annotations = GFF3Parser(filename=gff).read()
        names = annotations.sort_values(by=['seqid', 'score'], ascending=True).query('score < 1e-05').drop_duplicates(subset='seqid')[['seqid', 'Name']]
        new_file = names.dropna(axis=0,how='all')

        name_conversion = {}
        with open(namemap, 'r') as namemap_reader:
            namemap_reader.readline()
            for line in namemap_reader:
                line = line.rstrip('\n').split(',')
                original = line[0]
                new = line[1]
                name_conversion[new] = original

        old_transcript_names = [name_conversion[a] for a in new_file['seqid']]
        new_file['transcript_id'] = old_transcript_names

        new_file = new_file.sort_values(by=['Name', 'transcript_id'])
        # dataframe_to_save = new_file[['transcript_id', 'Name']]
        new_file.to_csv(f'{self.output_dir}/genes_{suffix}.csv', sep='\t', columns=['transcript_id', 'Name'] )

        unique_genes = set(new_file['Name'])

        contigs_per_gene = {}
        for gene in unique_genes:
            a =  new_file.loc[new_file['Name'] == gene]
            transcripts = []
            for c in a['transcript_id']:
                transcripts.append(c)
            contigs_per_gene[gene] = transcripts

        return(contigs_per_gene)

    def combine_dicts(self, a, b):
        """Combines two dictionaries and adds both values to the new dictionary.
        
        Arguments:
            a {[type]} -- [description]
            b {[type]} -- [description]
        
        Returns:
            [type] -- [description]
        """
        keys_a = set(a.keys())
        keys_b = set(b.keys())
        keys = keys_a.union(keys_b)

        new_dict = {}
        for key in keys:
            if key in a and key in b:
                new_dict[key] = a[key] + b[key]
            elif key in a:
                new_dict[key] = a[key]
            elif key in b:
                new_dict[key] = b[key]

        return new_dict

    def save(self, output):
        """Saves all calculates metrics to a file.
        
        Arguments:
            output {[type]} -- [description]
        """

        logger.info('= = = = = = = = =')
        with open(output, 'w') as writer:
            lost_gene_line = f'Lost genes:\t{self.lost_genes[0]}\t{self.lost_genes[1]}'
            writer.write(lost_gene_line + '\n')
            logger.info(lost_gene_line)
            
            header = 'Removed sequences per gene'
            writer.write(header + '\n')
            logger.info(header)
            for (key, value) in self.removed_sequences_per_gene.items():
                line = f'{key}\t{value}'
                writer.write(line + '\n')
                logger.info(line)

            header = 'Gene\tNo. of clusters(kmer)\tNo. of clusters(read graph)\tSequences per gene(kmer)\tSequences_per_gene(read graph)'
            writer.write(header + '\n')
            logger.info(header)
            for (key, value) in self.clusters_and_sequences_per_gene.items():
                line = key + '\t' + '\t'.join((str(i) for i in value))
                writer.write(line + '\n')
                logger.info(line)
