import argparse
from logs import logger
import os

parser = argparse.ArgumentParser(description='Description.')
parser.add_argument('-i', '--input', dest='FASTA_FILE', required=True, metavar='STR', help='Input fasta file to cluster.')
parser.add_argument('-k', '--kmer', dest='KMER_SIZE', metavar='STR', default = None, help='Kmer size to perform kmer based clustering. default: 5p6-mer')
parser.add_argument('-1', '--left', dest='R1', metavar='STR', default=None, help='Left reads.')
parser.add_argument('-2', '--right', dest='R2', metavar='STR', default=None, help='Right reads.')
parser.add_argument('-s', '--single', dest='R', metavar='STR', default=None, help='Single end reads.')
parser.add_argument('-o', '--output', dest='OUTPUT_DIR', metavar='STR', default=None, help='Output directory.')
parser.add_argument('-t', '--threads', dest='THREADS', type=int, metavar='INT', default=6, help='Threads to use. (default: 6)')

parser.add_argument('-m', '--mapping', dest='MAPPING_METHOD', default='salmon', choices=['salmon', 'hisat'], help='Choose mapping method. Salmon is a lot faster, but may be inaccurate. Hisat is more precise but takes A LOT of time. (default: salmon)')
parser.add_argument('-d', '--draw', dest='DRAW', action='store_true', default=False, help='Draw the graphs before and after read based clustering. May take a lot of time. (default: False)')
parser.add_argument('--keep-sam', dest='KEEP_SAM', action='store_true', default=False, help='Keep sam files. Takes a lot of space. default: False')
parser.add_argument("-l", "--log", dest="LOGLEVEL", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default = 'INFO', help="Set the logging level.")

# Annotation options
parser.add_argument('--annotate', dest='ANNOTATE', action='store_true', default=False, help='Do you want to annotate? Because thats how you annotate.')
parser.add_argument('--busco-group', dest='BUSCO_GROUP', metavar='STR', choices=['actinobacteria', 'actinopterygii',
'alveolata_stramenophiles', 'arthropoda', 'ascomycota', 'aves', 'bacillales', 'bacteria',
'bacteroidetes', 'basidiomycota', 'betaproteobacteria', 'clostridia', 'cyanobacteria',
'deltaepsilonsub', 'dikarya', 'diptera', 'embryophyta', 'endopterygota',
'enterobacteriales', 'euarchontoglires', 'eukaryota', 'eurotiomycetes', 'firmicutes',
'fungi', 'gammaproteobacteria', 'hymenoptera', 'insecta', 'lactobacillales',
'laurasiatheria', 'mammalia', 'metazoa', 'microsporidia', 'nematoda', 'pezizomycotina',
'proteobacteria', 'protists', 'rhizobiales', 'saccharomyceta', 'saccharomycetales',
'sordariomyceta', 'spirochaetes', 'tenericutes', 'tetrapoda', 'vertebrata'], default='metazoa', help='(default: metazoa)')
parser.add_argument('--database-dir', dest='DATABASE_DIR', metavar='STR', help='Directory with already built dammit database.')

args = parser.parse_args()

# Set loggging level
if args.LOGLEVEL != 'INFO':
    logger.setLevel(args.LOGLEVEL)
    logger.info(f'Loglevel set to {args.LOGLEVEL}.')

# Set output directory
if args.OUTPUT_DIR:
    args.OUTPUT_DIR = os.path.abspath(args.OUTPUT_DIR) + '/'
else:
    args.OUTPUT_DIR = ''

# Check if input files are valid
all_input_files = [args.FASTA_FILE]
if args.R:
    all_input_files.append(args.R)
if args.R1:
    all_input_files.append(args.R1)
if args.R2:
    all_input_files.append(args.R2)
for infile in all_input_files:
    if not os.path.exists(infile):
        logger.error(f'{infile} is not existing. Please check input.')
        exit(1)


logger.info(args.OUTPUT_DIR)
