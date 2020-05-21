class Writer:
    """A helper class to write fastas and clstr files."""

    def __init__(self, file_name):
        self.file_name = file_name

    def write_fasta(self, sequences: dict) -> None:
        """ Writes fasta sequences from a dictionary.

        Args:
            sequences: Sequences to write to a file. Key/Value: Header/Sequence

        """
        line_length = 80
        with open(self.file_name, "w") as fasta_writer:
            for header, sequence in sequences.items():
                # Write header
                fasta_writer.write(header + "\n")
                # Write sequence
                split_starts = range(0, len(sequence), line_length)
                for part in [sequence[i : i + line_length] for i in split_starts]:
                    fasta_writer.write(part + "\n")

    def write_clstr(self, clusters: list, representative_sequences: list) -> None:
        """ Writes a clstr-like ouput just as cdhit-est.
        Mark the representative sequence with a star.

        Args:
            clusters: All clusters.
            representative_sequences: Representative sequences to the cluster list.

        """
        with open(f"{self.file_name}", "w") as clstr_writer:
            for cluster_id, cluster in enumerate(clusters, 1):
                clstr_writer.write(f">Cluster {cluster_id}\n")
                internal_cluster_counter = 0

                for sub_id, subcluster in enumerate(cluster, 1):
                    for seq in subcluster:
                        output_line = f"{internal_cluster_counter}_{sub_id}\t{seq}"
                        if f">{seq}" in representative_sequences:
                            output_line += "\t*"
                        output_line += "\n"

                        clstr_writer.write(output_line)
                        internal_cluster_counter += 1
