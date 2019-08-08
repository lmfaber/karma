# karma
Kmers And Read MAppings are used to optimize (de novo) transcriptome assemblies.


```
usage: karma.py [-h] -i STR [-k STR] [-1 STR] [-2 STR] [-s STR] [-o STR]
                [-t INT] [-m {salmon,hisat}] [-d] [--keep-sam]
                [-l {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--annotate]
                [--busco-group STR] [--database-dir STR]

Description.

optional arguments:
  -h, --help            show this help message and exit
  -i STR, --input STR   Input fasta file to cluster.
  -k STR, --kmer STR    Kmer size to perform kmer based clustering. default:
                        5p6-mer
  -1 STR, --left STR    Left reads.
  -2 STR, --right STR   Right reads.
  -s STR, --single STR  Single end reads.
  -o STR, --output STR  Output directory.
  -t INT, --threads INT
                        Threads to use. (default: 6)
  -m {salmon,hisat}, --mapping {salmon,hisat}
                        Choose mapping method. Salmon is a lot faster, but may
                        be inaccurate. Hisat is more precise but takes A LOT
                        of time. (default: salmon)
  -d, --draw            Draw the graphs before and after read based
                        clustering. May take a lot of time. (default: False)
  --keep-sam            Keep sam files. Takes a lot of space. default: False
  -l {DEBUG,INFO,WARNING,ERROR,CRITICAL}, --log {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the logging level.
  --annotate            Do you want to annotate? Because thats how you
                        annotate.
  --busco-group STR     (default: metazoa)
  --database-dir STR    Directory with already built dammit database.

```

## Example Command
`python karma.py -i file.fasta -1 left.fastq -2 right.fastq -o output_dir -t 30 --annotate --busco-group nematoda --database-dir database_dir`
If you do not wish to annotate leave out `--annotate`, `--busco-group` and `--database_dir`.
