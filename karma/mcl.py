import glob
import os
from cmd import Cmd

from logs import logger


class Mcl:
    """ Run MCL clustering on a given graph file. """
    def __init__(self):
        mcl = Cmd('which mcl')
        mcl.run()

    def run(self, graph_file, output_file, threads):
        """
        MCL: The input is then a file or stream in which each line encodes an edge in terms of two labels 
         (the 'A' and the 'B') and a numerical value (the 'C'), all separated by white space.
        A B 20
        A C 10
        The output is then a file where each line is a cluster of tab-separated labels.
        """
        logger.debug('MCL clustering...')
        if os.path.exists(output_file):
            os.remove(output_file)
        mcl = Cmd(f'mcl {graph_file} --abc -o {output_file} -te {threads} -resource 4 -V all')
        mcl.run()
