from cmd import Cmd
from dammit.fileio.gff3 import GFF3Parser
import pandas as pd
from logs import logger
from collections import Counter

def basename_woe(path):
    """ Basename without extension """
    return os.path.splitext(os.path.basename(path))[0]

class Dammit:

    def __init__(self, first_file, second_file, output_dir, database_dir, busco_group, threads = 6):
        self.transcriptomes = {'before': first_file, 'after': second_file}
        self.output_dir = output_dir
        self.database_dir = database_dir

        if not busco_group:
            self.busco_group = 'metazoa'
        else:
            self.busco_group = busco_group

        self.threads = threads

        self.namemaps = ['/home/lasse/Schreibtisch/karma/karma/cel_flux.PE/dammit/before/cel_flux.PE.fasta.dammit.namemap.csv', '/home/lasse/Schreibtisch/karma/karma/cel_flux.PE/dammit/after/cel_flux.PE.fa.dammit.namemap.csv']
        self.gff_files = ['/home/lasse/Schreibtisch/karma/karma/cel_flux.PE/dammit/before/cel_flux.PE.fasta.dammit.gff3', '/home/lasse/Schreibtisch/karma/karma/cel_flux.PE/dammit/after/cel_flux.PE.fa.dammit.gff3']
        # Run
        #self.update_database()
        #self.run()
        self.postprocessing()

    def run(self):
        for name, transcriptome in self.transcriptomes.items():
            output_dir = f'{self.output_dir}/{name}'
            gff_outfile = f'basename_woe(transcriptome).fasta.dammit.gff3'
            if not os.path.exists(f'{output_dir}/{gff_outfile}'):
                dammit = Cmd(f'dammit annotate {transcriptome} -o {output_dir} --database-dir {self.database_dir} --busco-group {self.busco_group} --n_threads {self.threads}')
                dammit.run()

    def update_database(self):
        database = Cmd(f'dammit databases --install --database-dir {self.database_dir} --busco-group {self.busco_group}')
        database.run()

    def postprocessing(self):
        for i, (namemap, gff_file) in enumerate(zip(self.namemaps, self.gff_files)):

            name_conversion = {}
            with open(namemap, 'r') as namemap_reader:
                namemap_reader.readline()
                for line in namemap_reader:
                    line = line.rstrip('\n').split(',')
                    original = line[0]
                    new = line[1]
                    name_conversion[new] = original

            annotations = GFF3Parser(filename=gff_file).read()
            names = annotations.sort_values(by=['seqid', 'score'], ascending=True).query('score < 1e-05').drop_duplicates(subset='seqid')[['seqid', 'Name']]
            new_file = names.dropna(axis=0,how='all')

            old_transcript_names = [name_conversion[a] for a in new_file['seqid']]
            new_file['transcript_id'] = old_transcript_names

            new_file = new_file.sort_values(by=['Name', 'transcript_id'])
            new_file.to_csv(f'{self.output_dir}/{i}_GENE_IDS.csv')

        before = pd.read_csv(f'{self.output_dir}/0_GENE_IDS.csv', header = 0)
        after = pd.read_csv(f'{self.output_dir}/1_GENE_IDS.csv', header = 0)

        names_before = before['Name']
        names_after = after['Name']

        # lost transcripts:
        transcripts_before = set(names_before)
        transcripts_after = set(names_after)
        lost_transcripts = transcripts_before.difference(transcripts_after)
        logger.info(f'Lost transcripts: {lost_transcripts}')

        # reduction of multiple transcripts
        counter_before = Counter(names_before)
        counter_after = Counter(names_after)

        total_counts_before = sum(counter_before.values())
        total_counts_after = sum(counter_after.values())
        logger.info(f'Reduced from {total_counts_before} by {total_counts_before-total_counts_after}')

        differences = 0
        for a in counter_before.items():
            current = a[0]
            c1 = counter_before[current]
            c2 = counter_after[current]
            if c2 != 0:
                differences += c1-c2
        differences = differences / len(counter_before.items())
        logger.info(f'Reduced by ~ {differences} sequences per cluster.')

        # annotations.to_csv('ANNOTATIONS.csv')
