class Writer:

    def __init__(self, file_name):
        self.file_name = file_name

    def write_fasta(self, sequences):
        """
        Sequences to write in dictionary.
        
        Arguments:
            sequences {dict} -- key/value: header/sequence
        """
        n = 80
        with open(self.file_name, 'w') as fasta_writer:
            for header, sequence in sequences.items():
                # Write header
                fasta_writer.write(header + '\n')
                # Write sequence
                for part in [sequence[i:i+n] for i in range(0, len(sequence), n)]:
                    fasta_writer.write(part + '\n')

    def write_clstr(self, clusters, representative_sequences):
        """
        Writes a clstr-like ouput just as cdhit-est. Mark the representative sequence with a star.
        """
        with open(f'{self.file_name}', 'w') as clstr_writer:
            for id, cluster in enumerate(clusters, 1):
                clstr_writer.write(f'>Cluster {id}\n')
                internal_cluster_counter = 0

                for sub_id, subcluster in enumerate(cluster, 1):
                    for seq in subcluster:
                        output_line = f'{internal_cluster_counter}_{sub_id}\t{seq}'
                        if f'>{seq}' in representative_sequences:
                            output_line += '\t*'
                        output_line += '\n'
                        
                        clstr_writer.write(output_line)
                        internal_cluster_counter += 1
