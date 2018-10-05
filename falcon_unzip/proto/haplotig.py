class Haplotig:
    def __init__(self, name = '', phase = ('', -1, 0), seq = None, edges = None, path = None, labels = None, cstart = None, cend = None):
        self.name = name
        self.phase = phase
        self.seq = [] if seq == None else seq
        self.edges = [] if edges == None else edges
        self.path = [] if path == None else path
        # Position of the haplotig on the collapsed primary contig (2-asm-falcon/p_ctg.fa).
        self.cstart = -1 if cstart == None else cstart
        self.cend = -1 if cend == None else cend
        self.labels = {} if labels == None else labels    # Additional user info, such as coordinates, annotations, and various attributes.
