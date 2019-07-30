from logs import logger
from cmd import Cmd
import glob
import os

class Hisat2:

    def __init__(self, input_file, output_file, index_name, threads):

        hisat2 = Cmd('which hisat2')
        hisat2.run()
        if hisat2.status != 0:
            logger.error('Hisat2 not found in path.')
            exit(1)

        # Check if indexing already run
        self.index_build_has_run = True if len(glob.glob(f'{index_name}.*.ht2')) == 8 else False

        # Check if mapping already run
        self.mapping_has_run = True if os.path.isfile(output_file) else False

        self.input_file = input_file
        self.output_file = output_file
        self.index_name = index_name
        self.threads = threads

    def build_index(self):
        if not self.index_build_has_run:
            logger.debug('Build index.')
            indexing = Cmd(f'hisat2-build -q -p {self.threads} {self.input_file} {self.index_name}')
            indexing.run()
            self.index_build_has_run = True
        else:
            logger.debug('Skipping index building.')

    def run(self, reads):
        if not self.mapping_has_run:
            logger.debug('Starting mapping...')
            if len(reads) == 1: # single end reads
                hisat = Cmd(f'hisat2 -q -x {self.index_name} -U {reads[0]} | \
                                samtools view -hS -F 4 -q 1 | \
                                samtools sort | \
                                samtools view -o {self.output_file} ')
            elif len(reads) == 2: # paired end reads
                hisat = Cmd(f'hisat2 -q --threads {self.threads} -k 1 -x {self.index_name} -1 {reads[0]} -2 {reads[1]} | \
                                samtools view -hS -F 4 -q 1 | \
                                samtools sort | \
                                samtools view -o {self.output_file}')
            hisat.run()
            self.mapping_has_run = True
        else:
            logger.debug('Skipping mapping.')
