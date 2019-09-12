# karma
Using Kmers And Read MAppings to optimize (de novo) transcriptome assemblies.

## Installation
[![Install with conda](https://anaconda.org/lmfaber/karma/badges/installer/conda.svg)](https://anaconda.org/lmfaber/karma)

```
conda install -c bioconda -c conda-forge -c lmfaber karma
```

## Parameter
```
usage: karma.py [-h] -i STR [-1 STR] [-2 STR] [-s STR] [-o STR] [-t INT]
                [-l {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--lowest]
                [--rearrange] [--threshold THRESHOLD] [-d] [-k STR]
                [--n_neighbors N_NEIGHBORS] [--n_components N_COMPONENTS]
                [--min_dist MIN_DIST] [--random_state RANDOM_STATE]
                [--min_cluster_size MIN_CLUSTER_SIZE] [--inflation INT]
                [--annotate] [--busco-group STR] [--database-dir STR]

karma. Using Kmers And Read MAppings to optimize (de novo) transcriptome
assemblies.

optional arguments:
  -h, --help                            show this help message and exit

KARMA:
  Necessary arguments for karma.

  -i STR, --input STR                   Input fasta file to cluster. (default:
                                        None)
  -1 STR, --left STR                    Left reads. (default: None)
  -2 STR, --right STR                   Right reads. (default: None)
  -s STR, --single STR                  Single end reads. (default: None)
  -o STR, --output STR                  Output directory. (default:
                                        karma_output)
  -t INT, --threads INT                 Threads to use. (default: 6)
  -l {DEBUG,INFO,WARNING,ERROR,CRITICAL}, --log {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                                        Set the logging level. (default: INFO)
  --lowest                              Representative sequence selection. Per
                                        default karma will select the sequence
                                        with the highest shared read
                                        information in one cluster. With this
                                        option the sequence with the lowest
                                        shared read information will be
                                        selected as well. Increasing
                                        redundancy but may increase overall
                                        information. (default: False)
  --rearrange                           After kmer and read graph based
                                        clustering the clusters can be
                                        rearranged based on shared reads.
                                        (default: False)
  --threshold THRESHOLD                 Threshold for rearranging groups.
                                        (default: 0)
  -d, --draw                            Draw the graphs before and after read
                                        based clustering. May take a lot of
                                        time. (default: False)

KMER:
  Kmer based arguments required for UMAP (Dimension reduction) and HDBSCAN
  (clustering algorithm). See https://umap-
  learn.readthedocs.io/en/latest/parameters.html#
  https://hdbscan.readthedocs.io/en/latest/parameter_selection.html

  -k STR, --kmer STR                    Kmer size to perform kmer based
                                        clustering. (default: 5p6)
  --n_neighbors N_NEIGHBORS             UMAP: This parameter controls how UMAP
                                        balances local versus global structure
                                        in the data. (default: 2)
  --n_components N_COMPONENTS           UMAP: allows the user to determine the
                                        dimensionality of the reduced
                                        dimension space we will be embedding
                                        the data into (default: 10)
  --min_dist MIN_DIST                   UMAP: Minimum distance between two
                                        data points. (default: 0)
  --random_state RANDOM_STATE           UMAP: Initialization value for
                                        reproducible results. (default: 42)
  --min_cluster_size MIN_CLUSTER_SIZE   HDBSCAN: Minimum cluster size.
                                        (default: 2)

MCL:
  MCL related arguments. See https://micans.org/mcl/

  --inflation INT                       Sets the main inflation value to INT.
                                        This value is the main handle for
                                        affecting cluster granularity. Between
                                        1.2 and 6. (default: 2)

Dammit:
  Parameters for de novo gene annotation. See https://github.com/dib-
  lab/dammit

  --annotate                            Switch for annotating. (default:
                                        False)
  --busco-group STR                     Which BUSCO group to use. Should be
                                        chosen based on the organism being
                                        annotated. (default: metazoa)
  --database-dir STR                    Directory to store databases. Existing
                                        databases will not be overwritten.
                                        (default: None)


```

## Example Command
`python karma.py -i file.fasta -1 left.fastq -2 right.fastq -o output_dir -t 30 --annotate --busco-group nematoda --database-dir database_dir`
If you do not wish to annotate leave out `--annotate`, `--busco-group` and `--database_dir`.
