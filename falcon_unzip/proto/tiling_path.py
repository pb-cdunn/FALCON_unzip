from falcon_kit import io

class TilingPathEdge:
    def __init__(self, split_line = None):
        self.ctg_id, self.v, self.w, self.b, self.e, self.score, self.identity = None, None, None, None, None, None, None
        self.parsed = False
        self.split_line = None
        if split_line:
            self.set_from(split_line)
            self.split_line = split_line
    def set_from(self, split_line):
        self.parsed = False
        self.ctg_id = split_line[0]
        self.v = split_line[1]
        self.w = split_line[2]
        self.b = int(split_line[4])
        self.e = int(split_line[5])
        self.score = float(split_line[6])
        self.identity = float(split_line[7])
        self.parsed = True
    # def __str__(self):
    #     seq_id = self.w.split(':')[1]
    #     return '%s %s %s %s %d %d %d %.2f' % (self.ctg_id, self.v, self.w, seq_id, self.b, self.e, int(self.score), self.identity)


class TilingPath:
    def __init__(self, edge_list, contig_sequence_len = None):
        self.edges = edge_list  # These are TilingPathEdge objects.

        offset = 0

        # If the total contig sequence len is known, use that to
        # calculate the length of the first read (in case proper
        # contigs are specified). This is needed to offset the coordinates
        # which can be calculated from the tiling path.
        if contig_sequence_len != None:
            _, tiling_len = calc_node_coords(edge_list)
            assert(contig_sequence_len >= tiling_len)
            offset = contig_sequence_len - tiling_len   # This is the length of the first node.

        # The self.coords is a dict: self.coords[v] = coordinate_on_contig
        self.coords, self.contig_len = calc_node_coords(edge_list, offset)

        if contig_sequence_len != None:
            assert(self.contig_len == contig_sequence_len)

        self.v_to_edge = {}
        self.w_to_edge = {}
        for i in xrange(len(self.edges)):
            e = self.edges[i]
            self.v_to_edge[e.v] = i
            self.w_to_edge[e.w] = i

    def dump_as_split_lines(self):
        ret = []
        for e in self.edges:
            ret.append(e.split_line)
        return ret

    def get_subpath(self, start_coord, end_coord):
        """
        If end_coord is <= 0, then the entire suffix is taken.
        """
        if len(self.edges) == 0:
            return []

        if end_coord <= 0:
            end_coord = self.contig_len
        assert(start_coord <= end_coord)
        # end_coord -= 1  # Make the end inclusive.
        # start_node = None
        # end_node = None
        start_edge = None
        end_edge = None
        if start_coord < self.coords[self.edges[0].v]:
            start_edge = 0
        if end_coord <= self.coords[self.edges[0].v]:
            end_edge = 1
        for i in xrange(len(self.edges)):
            e = self.edges[i]
            if start_coord >= self.coords[e.v] and start_coord < self.coords[e.w]:
                start_edge = i
            if end_coord > self.coords[e.v] and end_coord <= self.coords[e.w]:
                end_edge = i + 1
        if end_coord >= self.coords[self.edges[-1].w]:
            end_edge = len(self.edges)
        assert(start_edge != None and end_edge != None)

        # Since the start_coord and end_coord can end within an edge,
        # we return the position in the final contigas.

        new_start_coord = start_coord - self.coords[self.edges[start_edge].v]
        new_end_coord = end_coord - self.coords[self.edges[start_edge].v]
        new_path = self.edges[start_edge:end_edge]

        new_path = [val.split_line for val in new_path]

        return new_path, new_start_coord, new_end_coord

    # def load_from_split_lines(self, split_lines):

def calc_node_coords(tiling_path, offset=0):
    """
    For a single tiling path (tiling_path is a list
    of edges for a particular contig) calculates the
    genomic coordinate of every node in the path.
    In case there are cycles in the tiling path,
    the existing node's coordinate will be overwritten.
    Offset refers to the length of the first node. If
    not specified, the contig length should not
    consider the length of the first node.
    """
    if not tiling_path:
        return {}, 0
    coord_map = {}
    contig_len = 0
    edge0 = tiling_path[0]
    coord_map[edge0.v] = offset
    for edge in tiling_path:
        if edge.v not in coord_map:
            raise Exception(
                'Tiling path is not in sorted order. Node "{v!r}" does not yet have an assigned coordinate.'.format(v=edge.v))
        coord = coord_map[edge.v]
        coord += abs(int(edge.b) - int(edge.e))
        coord_map[edge.w] = coord
        contig_len = max(contig_len, coord)
    return coord_map, contig_len

def load_tiling_paths(tp_file, whitelist_seqs, contig_lens, contig_prefix = None):
    """
    Parameters:
        whitelist_seqs - a dict or a set object containing contig IDs to load. If None, no filter will be applied.
        contig_lens -   if a dict with contig sequence lengths is specified, the difference between the
                        contig len and the length obtained from the tiling path will be used to offset
                        the tiling path coordinates.
    """
    tiling_path_edges = {}
    counter = io.FilePercenter(tp_file)
    with open(tp_file) as fp:
        for line in fp:     # Example row: "0 000000007:B 000000005:B 000000005 9 0 1980 99.95"
            counter(len(line))
            line = line.strip()
            if len(line) == 0: return
            sl = line.split()
            new_edge = TilingPathEdge(sl)
            if contig_prefix != None and new_edge.ctg_id.startswith(contig_prefix) == False:
                continue
            ctg_id = new_edge.ctg_id
            # The whitelist filter.
            if whitelist_seqs != None and (ctg_id in whitelist_seqs) == False:
                continue
            tiling_path_edges.setdefault(ctg_id, [])
            tiling_path_edges[ctg_id].append(new_edge)
    del counter
    # Convert the flat lists to objects for easier comprehention.
    tiling_paths = {}
    for ctg_id, edges in tiling_path_edges.iteritems():
        ctg_len = None
        if contig_lens != None:
            ctg_len = contig_lens[ctg_id]
        tiling_paths[ctg_id] = TilingPath(edges, ctg_len)
    return tiling_paths

def convert_split_lines_to_tiling_path(path_split_lines, contig_sequence_len = None):
    edges = []
    for sl in path_split_lines:
        edges.append(TilingPathEdge(sl))
    return TilingPath(edges, contig_sequence_len)

def find_a_ctg_placement(p_paths, a_paths):
    placement = {}
    for a_ctg_id, a_tp in a_paths.iteritems():
        first_node = a_tp.edges[0].v
        last_node = a_tp.edges[-1].w
        p_ctg_id = a_ctg_id.split('-')[0].split('_')[0]
        p_tp = p_paths[p_ctg_id]
        start, end = p_tp.coords[first_node], p_tp.coords[last_node]
        placement.setdefault(p_ctg_id, {})
        placement[p_ctg_id][a_ctg_id] = (start, end, p_ctg_id, a_ctg_id, first_node, last_node)
    return placement
