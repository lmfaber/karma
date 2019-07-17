from cmd import Cmd

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

        # Run
        self.update_database()
        self.run()

    def run(self):
        for name, transcriptome in self.transcriptomes.items():
            output_dir = f'{self.output_dir}/{name}'
            dammit = Cmd(f'dammit annotate {transcriptome} -o {output_dir} --database-dir {self.database_dir} --busco-group {self.busco_group} --n_threads {self.threads}')
            dammit.run()

    def update_database(self):
        database = Cmd(f'dammit databases --install --database-dir {self.database_dir} --busco-group {self.busco_group}')
        database.run()

    def postprocessing(self):
        pass

