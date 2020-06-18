from logs import logger
from cmd import Cmd
import glob


class Hisat2:
    """ Hisat2 mapping. """

    def __init__(self, input_file: str, index_name: str, threads: int):
        which = Cmd('which hisat2')
        which.run()

        # Check if indexing already run
        self.index_build_has_run = True if len(glob.glob(f'{index_name}.*.ht2')) == 8 else False

        self.input_file = input_file
        self.index_name = index_name
        self.threads = threads

    def build_index(self) -> None:
        """Build the Hisat2 index."""
        if not self.index_build_has_run:
            logger.debug('Build Hisat2 index.')
            indexing = Cmd(f'hisat2-build -q -p {self.threads} {self.input_file} {self.index_name}')
            indexing.run()
            self.index_build_has_run = True
        else:
            logger.debug('Skipping index building.')

    def run(self, reads: str):
        """Run the Hisat2 mapping with the given reads.

        Arguments:
            reads: Reads to perform mapping with.
        """
        logger.debug('Perform Hisat2 mapping.')
        if len(reads) == 1: # single end reads
            hisat = Cmd(f'hisat2 -q --threads {self.threads} -k 1 -x {self.index_name} -U {reads[0]} --no-unal | \
                        samtools view --threads {self.threads} -hS -F 4 -q 1 -O SAM')
        elif len(reads) == 2: # paired end reads
            hisat = Cmd(f'hisat2 -q --threads {self.threads} -k 1 -x {self.index_name} -1 {reads[0]} -2 {reads[1]} --no-unal | \
                            samtools view --threads {self.threads} -hS -F 4 -q 1 -O SAM')
        hisat.run()
        self.mapping_has_run = True
        return (entry for entry in hisat.stdout.split('\n')[:-1] if not entry.startswith('@'))
