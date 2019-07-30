conda activate karma
python3 karma.py -i /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/cel_flux.PE.fasta -1 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/cel_flux/fastp/cel_flux_1.fastq -2 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/cel_flux/fastp/cel_flux_2.fastq --busco-group nematoda --database-dir /mnt/prostlocal/lasse/kmer_bases_cluster/dammit-db -l DEBUG --keep-sam 


python3 karma2.py -i /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/cel_flux.PE.fasta -1 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/cel_flux/fastp/cel_flux_1.fastq -2 /mnt/prostlocal/lasse/snakemake_pipelines/multiAssembly/output/cel_flux/fastp/cel_flux_2.fastq --busco-group nematoda --database-dir /mnt/prostlocal/lasse/kmer_bases_cluster/dammit-db -l DEBUG --keep-sam 
