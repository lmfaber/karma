from logs import logger

class Contig():
    
    def __init__(self, name):
        self.name = name
        # self.reads = {}
        self.readset = set()

    def add_read(self, name, position):
        # self.reads[name] = position
        self.readset.add(name)

    def has_read(self, name):
        return(name in self.reads)

    def load_contig_info_from_sam(self, sam_file):
        """
        Reads the sam file for each contig and saves the needed information in a Contig object.
        TODO: Save contig objects so the program doesnt have to read all sam files again
        """
        logger.debug(f'Creating new contig object...')
        with open(sam_file, 'r') as sam_reader:
            for line in sam_reader:
                read, _, name, position, *_ = line.split('\t')
                # TODO: Check if this can be left out. Should be possible, because mapping took place on single contigs.
                #if name == new_contig.name:
                #    new_contig.add_read(read, position)
                self.add_read(read, position)
    
    def load_from_iterator(self, sam_infos):
        """Read the infos from iterator in RAM.
        
        Arguments:
            sam_infos {[type]} -- [description]
        """
        for line in sam_infos:
            read, _, name, position, *_ = line.split('\t')
            self.add_read(read, position)
            