
class PhasingBlock:
    def __init__(self, line):
        self.block_id, self.start, self.end_inclusive, self.block_len, self.num_snps, self.avg_snp_dist = [None, None, None, None, None, None]
        if line:
            self.parse(line)

    def parse(self, line):
        line = line.strip()
        if len(line) == 0: return
        if line[0] != 'P': return
        sl = line.split()
        self.block_id, self.start, self.end_inclusive, self.block_len, self.num_snps, self.avg_snp_dist = [int(val) for val in sl[1:6]] + [float(sl[6])]

def load_phasing_blocks(in_path):
    blocks = {}
    with open(in_path) as fp:
        for line in fp:
            line = line.strip()
            if len(line) == 0: continue
            if line[0] != 'P': continue
            new_block = PhasingBlock(line)
            blocks[new_block.block_id] = new_block

def load_rid_to_phase(rid_phase_map):
    arid2phase = {}
    with open(rid_phase_map) as f:
        for row in f:
            row = row.strip().split()
            arid2phase[row[0]] = (row[1], row[2], row[3])  # ctg_id, phase_blk_id, phase_id
    assert arid2phase, 'Empty rid_phase_map: {!r}'.format(rid_phase_map)
    return arid2phase
