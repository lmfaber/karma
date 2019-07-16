class Contig():
    
    def __init__(self, name):
        self.name = name
        self.reads = {}
        self.readset = set()

    def add_read(self, name, position):
        self.reads[name] = position
        self.readset.add(name)

    def has_read(self, name):
        return(name in self.reads)