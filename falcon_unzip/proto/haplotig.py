class Haplotig:
    def __init__(self, name, phase, seq, edges, path, labels, cstart, cend):
        self.name = name            # ''
        self.phase = phase          # (ctg_id, -1, 0)
        self.seq = seq              # ''
        self.edges = edges          # list()
        self.path = path            # list()
        # Position of the haplotig on the collapsed primary contig (2-asm-falcon/p_ctg.fa).
        self.cstart = cstart        # -1
        self.cend = cend            # -1
        # Additional user info, such as coordinates, annotations, and various attributes.
        self.labels = labels        # dict()
