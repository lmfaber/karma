from logs import logger
from cmd import Cmd
import glob
import os

class Salmon:
    """ Salmon mapping. """

    def __init__(self, input_file, output_dir, index_name, threads):
        which = Cmd('which salmon')
        which.run()

        self.input_file = input_file
        self.output_dir = output_dir
        self.index_name = index_name
        self.threads = threads
    
    def build_index(self):
        """ Build the salmon index. """
        logger.debug('Build salmon index.')
        indexing = Cmd(f'salmon index -p {self.threads} -t {self.input_file} -i {self.index_name} --keepDuplicates')
        indexing.run()

    def run(self, reads):
        """ Run the salmon mapping with the given reads. """
        logger.debug('Perform salmon mapping.')
        if len(reads) == 1: # single end reads
            salmon = Cmd(f'salmon quant --libType A --validateMappings --dumpEq -p {self.threads} -i {self.index_name} --unmatedReads {reads[0]} -o {self.output_dir}')
        elif len(reads) == 2: # paired end reads
            salmon = Cmd(f'salmon quant --libType A --validateMappings --dumpEq -p {self.threads} -i {self.index_name} -1 {reads[0]} -2 {reads[1]} -o {self.output_dir}')
        salmon.run()
