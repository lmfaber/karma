import argparse
from logs import logger
import os


parser = argparse.ArgumentParser(
    description="karma. Using Kmers And Read MAppings to optimize (de novo) transcriptome assemblies.",
    formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(
        prog, max_help_position=40
    ),
)
# KARMA
group_karma = parser.add_argument_group("KARMA", "Necessary arguments for karma.")
group_karma.add_argument(
    "-i",
    "--input",
    dest="FASTA_FILE",
    required=True,
    metavar="STR",
    help="Input fasta file to cluster.",
)
group_karma.add_argument(
    "-1", "--left", dest="R1", metavar="STR", default=None, help="Left reads."
)
group_karma.add_argument(
    "-2", "--right", dest="R2", metavar="STR", default=None, help="Right reads."
)
group_karma.add_argument(
    "-s", "--single", dest="R", metavar="STR", default=None, help="Single end reads."
)
group_karma.add_argument(
    "-o",
    "--output",
    dest="OUTPUT_DIR",
    metavar="STR",
    default="karma_output",
    help="Output directory.",
)
group_karma.add_argument(
    "-t",
    "--threads",
    dest="THREADS",
    type=int,
    metavar="INT",
    default=6,
    help="Threads to use.",
)
group_karma.add_argument(
    "-l",
    "--log",
    dest="LOGLEVEL",
    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    default="INFO",
    help="Set the logging level.",
)

group_karma.add_argument(
    "--lowest",
    dest="LOWEST",
    action="store_true",
    default=False,
    help="Representative sequence selection. Per default karma will select the sequence with the highest shared read information in one cluster. With this option the sequence with the lowest shared read information will be selected as well. Increasing redundancy but may increase overall information.",
)
group_karma.add_argument(
    "--rearrange",
    dest="REARRANGE",
    action="store_true",
    default=False,
    help="After kmer and read graph based clustering the clusters can be rearranged based on shared reads.",
)
group_karma.add_argument(
    "--threshold",
    dest="THRESHOLD",
    type=int,
    default=0,
    help="Threshold for rearranging groups.",
)
group_karma.add_argument(
    "-d",
    "--draw",
    dest="DRAW",
    action="store_true",
    default=False,
    help="Draw the graphs before and after read based clustering. May take a lot of time.",
)

# KMER
group_kmer = parser.add_argument_group(
    "KMER",
    "Kmer based arguments required for UMAP (Dimension reduction) and HDBSCAN (clustering algorithm). See \nhttps://umap-learn.readthedocs.io/en/latest/parameters.html#\nhttps://hdbscan.readthedocs.io/en/latest/parameter_selection.html",
)
group_kmer.add_argument(
    "-k",
    "--kmer",
    dest="KMER_SIZE",
    metavar="STR",
    default="5p6",
    help="Kmer size to perform kmer based clustering.",
)
# UMAP
group_kmer.add_argument(
    "--n_neighbors",
    dest="N_NEIGHBORS",
    type=int,
    default=2,
    help="UMAP: This parameter controls how UMAP balances local versus global structure in the data.",
)
group_kmer.add_argument(
    "--n_components",
    dest="N_COMPONENTS",
    type=int,
    default=10,
    help="UMAP: allows the user to determine the dimensionality of the reduced dimension space we will be embedding the data into",
)
group_kmer.add_argument(
    "--min_dist",
    dest="MIN_DIST",
    type=int,
    default=0,
    help="UMAP: Minimum distance between two data points.",
)
group_kmer.add_argument(
    "--random_state",
    dest="RANDOM_STATE",
    type=int,
    default=42,
    help="UMAP: Initialization value for reproducible results.",
)
# HDBSCAN
group_kmer.add_argument(
    "--min_cluster_size",
    dest="MIN_CLUSTER_SIZE",
    type=int,
    default=2,
    help="HDBSCAN: Minimum cluster size.",
)

# MCL clustering options
group_mcl = parser.add_argument_group(
    "MCL", "MCL related arguments. See https://micans.org/mcl/"
)
group_mcl.add_argument(
    "--inflation",
    dest="INFLATION",
    type=float,
    metavar="INT",
    default=2,
    help="Sets the main inflation value to INT. This value is the main handle for affecting cluster granularity. Between 1.2 and 6.",
)

# Dammit annotation options
group_dammit = parser.add_argument_group(
    "Dammit",
    "Parameters for de novo gene annotation. See https://github.com/dib-lab/dammit",
)
group_dammit.add_argument(
    "--annotate",
    dest="ANNOTATE",
    action="store_true",
    default=False,
    help="Switch for annotating.",
)
group_dammit.add_argument(
    "--busco-group",
    dest="BUSCO_GROUP",
    metavar="STR",
    choices=[
        "actinobacteria",
        "actinopterygii",
        "alveolata_stramenophiles",
        "arthropoda",
        "ascomycota",
        "aves",
        "bacillales",
        "bacteria",
        "bacteroidetes",
        "basidiomycota",
        "betaproteobacteria",
        "clostridia",
        "cyanobacteria",
        "deltaepsilonsub",
        "dikarya",
        "diptera",
        "embryophyta",
        "endopterygota",
        "enterobacteriales",
        "euarchontoglires",
        "eukaryota",
        "eurotiomycetes",
        "firmicutes",
        "fungi",
        "gammaproteobacteria",
        "hymenoptera",
        "insecta",
        "lactobacillales",
        "laurasiatheria",
        "mammalia",
        "metazoa",
        "microsporidia",
        "nematoda",
        "pezizomycotina",
        "proteobacteria",
        "protists",
        "rhizobiales",
        "saccharomyceta",
        "saccharomycetales",
        "sordariomyceta",
        "spirochaetes",
        "tenericutes",
        "tetrapoda",
        "vertebrata",
    ],
    default="metazoa",
    help="Which BUSCO group to use. Should be chosen based on the organism being annotated.",
)
group_dammit.add_argument(
    "--database-dir",
    dest="DATABASE_DIR",
    metavar="STR",
    help="Directory to store databases. Existing databases will not be overwritten.",
)

args = parser.parse_args()

# Set loggging level
if args.LOGLEVEL != "INFO":
    logger.setLevel(args.LOGLEVEL)
    logger.info(f"Loglevel set to {args.LOGLEVEL}.")

# Set output directory
args.OUTPUT_DIR = os.path.abspath(args.OUTPUT_DIR)

# Set correct kmer dataype
if args.KMER_SIZE != "5p6":
    args.KMER_SIZE = int(args.KMER_SIZE)

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
        logger.error(f"{infile} is not existing. Please check input.")
        exit(1)
