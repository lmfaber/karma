import argparse
import logging
from logs import logger

parser = argparse.ArgumentParser(description='Description.')
parser.add_argument('-i', '--input', dest='FASTA_FILE', metavar='STR', help='')
parser.add_argument('-k', '--kmer', dest='KMER_SIZE', metavar='STR', default = None, help='')
parser.add_argument('-1', '--left', dest='R1', metavar='STR', help='')
parser.add_argument('-2', '--right', dest='R2', metavar='STR', help='')
parser.add_argument('-s', '--single', dest='R', metavar='STR', default=None, help='')
parser.add_argument('-o', '--output', dest='OUTPUT_DIR', metavar='STR', help='')
parser.add_argument('-t', '--threads', dest='THREADS', metavar='INT', default=4, help='Threads to use. (default: 6)')
parser.add_argument("-l", "--log", dest="LOGLEVEL", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default = logging.INFO, help="Set the logging level")
args = parser.parse_args()

# Set loggging level
if args.LOGLEVEL != logging.INFO:
    logger.setLevel(args.LOGLEVEL)
    logger.info(f'Loglevel set to {args.LOGLEVEL}.')
