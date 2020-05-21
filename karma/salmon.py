from cmd import Cmd
import os
from logs import logger


class Salmon:
    """ Salmon mapping. """

    def __init__(self, input_file, output_dir, index_name, threads):
        which = Cmd("which salmon")
        which.run()

        self.input_file = input_file
        self.output_dir = output_dir
        self.index_name = index_name
        self.threads = threads

    def build_index(self) -> None:
        """ Builds the salmon index.

        """
        logger.debug("Build salmon index.")
        # TODO: Implement check to avoid duplicate runs
        indexing = Cmd(
            f"salmon index -p {self.threads} -t {self.input_file} -i {self.index_name} --keepDuplicates"
        )
        indexing.run()

    def run(self, reads: list) -> None:
        """ Run the salmon mapping with the given reads.

        Args:
            reads: List of reads. Either paired end or single end.

        """
        logger.debug("Perform salmon mapping.")
        if not os.path.exists(f"{self.output_dir}/aux_info/eq_classes.txt"):
            if len(reads) == 1:  # single end reads
                salmon = Cmd(
                    f"salmon quant --libType A --validateMappings --dumpEq -p {self.threads} -i {self.index_name} --unmatedReads {reads[0]} -o {self.output_dir}"
                )
            elif len(reads) == 2:  # paired end reads
                salmon = Cmd(
                    f"salmon quant --libType A --validateMappings --dumpEq -p {self.threads} -i {self.index_name} -1 {reads[0]} -2 {reads[1]} -o {self.output_dir}"
                )
            salmon.run()
        else:
            logger.info("Skipping mapping.")
