from falcon_kit.FastaReader import FastaReader
import os
import re
import networkx as nx
import falcon_unzip.proto.sam2m4 as sam2m4

RCMAP = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))

def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def reverse_end(node_id):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end

def revcmp_seq(seq):
    return "".join([RCMAP[c] for c in seq[::-1]])

def load_sg_seq(all_read_ids, fasta_fn):
    seqs = {}
    # load all p-read name into memory
    f = FastaReader(fasta_fn)
    for r in f:
        if r.name not in all_read_ids:
            continue
        seqs[r.name] = r.sequence.upper()
    assert seqs, 'No sg_seqs. Maybe empty fasta input? {!r}'.format(fasta_fn)
    return seqs

def load_all_seq(fasta_fn):
    seqs = {}
    f = FastaReader(fasta_fn)
    for r in f:
        seqs[r.name.split()[0]] = r.sequence.upper()
    return seqs

def path_to_seq(preads, path, with_first_read):
    ret = ''

    if len(path) == 0:
        return ret

    if with_first_read:
        ctg_id, v, w, wrid, sp, tp = path[0][0:6]
        vrid, vorient = v.split(':')
        ret += preads[vrid] if vorient == 'E' else "".join([RCMAP[c] for c in preads[vrid][::-1]])

    for edge in path:
        ctg_id, v, w, wrid, sp, tp = edge[0:6]
        sp, tp = int(sp), int(tp)
        ret += preads[wrid][sp:tp] if sp < tp else "".join([RCMAP[c] for c in preads[wrid][sp:tp:-1]])

    ret = ''.join([str(val) for val in ret])

    return ret

def write_haplotig_path_and_edges(haplotigs, haplotigs_path_path, haplotigs_edges_path):
    with open(haplotigs_path_path, "w") as fp_h_tig_path, \
            open(haplotigs_edges_path, "w") as fp_edges:
        for h_name, haplotig in haplotigs.iteritems():
            for edge in haplotig.edges:
                fp_edges.write(' '.join([str(val) for val in edge]) + '\n')
            for edge in haplotig.path:
                fp_h_tig_path.write(' '.join([str(val) for val in edge]) + '\n')

def write_haplotigs(haplotigs, haplotig_seqs_out_path, fp_proto_log, hack_qnames_for_blasr=False):
    num_htigs = 0
    with open(haplotig_seqs_out_path, "w") as fp_h_tig_fa:
        # sorted_hnames = sorted(haplotigs.keys(), key = lambda x: (x.split('-')[-1]))
        # for h_name in sorted_hnames:
        #     haplotig = haplotigs[h_name]
        for h_name, haplotig in haplotigs.iteritems():
            num_htigs += 1
            fp_proto_log('[IS] Writing {}, h_name = {}'.format(num_htigs, h_name))
            qname = haplotig.name
            qseq = ''.join(haplotig.seq)
            if hack_qnames_for_blasr == True:
                qname += '/%d/0_%d' % (num_htigs, len(qseq))
            fp_h_tig_fa.write('>%s\n%s\n' % (qname, qseq))

def write_diploid_groups(diploid_groups, out_path, fp_proto_log):
    haplotigs = {}
    for pb, phases in diploid_groups.iteritems():
        # Sanity check.
        if len(phases.keys()) != 2:
            continue
        # `phases` is a dict, where key is the phase ID, and value is [haplotig, aln]
        hap0, aln0 = phases[0]
        hap1, aln1 = phases[1]
        haplotigs[hap0.name] = hap0
        haplotigs[hap1.name] = hap1
    write_haplotigs(haplotigs, out_path, fp_proto_log)

def extract_weakly_unphased_haplotig_paths(graph, weight_param='weight'):
    best_paths = []
    for sub_hg in nx.weakly_connected_component_subgraphs(graph):
        best_path = extract_unphased_haplotig_paths(sub_hg)
        best_paths.append(best_path)
    return best_paths

def extract_unphased_haplotig_paths(sub_hg, weight_param='weight'):
    """
    For every combination of source and sink nodes in the sub_hg graph,
    it finds the best scoring path
    """

    # Handle a special case.
    if len(sub_hg.nodes()) == 1:
        best_path = [1, [sub_hg.nodes()[0]], sub_hg.nodes()[0], sub_hg.nodes()[0]]
        return best_path

    sources = [v for v in sub_hg.nodes() if sub_hg.in_degree(v) == 0]
    sinks = [v for v in sub_hg.nodes() if sub_hg.out_degree(v) == 0]

    best_path = None
    for s in sources:
        s_path = []
        for t in sinks:
            try:
                path = nx.shortest_path(sub_hg, s, t, weight=weight_param)
            except nx.exception.NetworkXNoPath:
                path = []
                continue

            path_score = 0
            for v, w in zip(path[:-1], path[1:]):
                path_score += sub_hg.edge[v][w][weight_param]

            s_path.append([path_score, path, s, t])

        s_path.sort(key=lambda x: x[0])

        if len(s_path) == 0:
            continue

        # Best path should have the least weight.
        if best_path == None or s_path[0][0] < best_path[0]:
            best_path = s_path[0]

    return best_path    # Components: [path_score, path, s, t]

def load_aln(sam_path):
    """
    Loads the SAM file and splits the header by '/',
    in case it was a BLASR hack.
    """
    m4 = sam2m4.sam_to_m4(sam_path)
    for aln in m4:
        aln[0] = aln[0].split('/')[0]

    return m4

def generic_nx_to_gfa(graph, fp_out, node_len_dict=None):
    line = 'H\tVN:Z:1.0'
    fp_out.write(line + '\n')

    if node_len_dict != None:
        for v in graph.nodes():
            line = 'S\t%s\t%s\tLN:i:%d' % (v, '*', node_len_dict[v])
            fp_out.write(line + '\n')
    else:
        for v in graph.nodes():
            line = 'S\t%s\t%s\tLN:i:%d' % (v, '*', 1000)
            fp_out.write(line + '\n')

    for v, w in graph.edges():
        line = 'L\t%s\t+\t%s\t+\t0M' % (v, w)
        fp_out.write(line + '\n')

def nx_to_gfa(ctg_id, haplotig_graph, all_haplotig_dict, fp_out):
    line = 'H\tVN:Z:1.0'
    fp_out.write(line + '\n')

    for v in haplotig_graph.nodes():
        line = 'S\t%s\t%s\tLN:i:%d' % (v, '*', len(all_haplotig_dict[v]['seq']))
        fp_out.write(line + '\n')

    for v, w in haplotig_graph.edges():
        line = 'L\t%s\t+\t%s\t+\t0M' % (v, w)
        fp_out.write(line + '\n')

    # Add a virtual source and sink so that both haplotypes appear
    # connected when visualized.
    virtual_source = '%s-source' % (ctg_id)
    line = 'S\t%s\t%s\tLN:i:%d' % (virtual_source, '*', 0)
    fp_out.write(line + '\n')
    virtual_sink = '%s-sink' % (ctg_id)
    line = 'S\t%s\t%s\tLN:i:%d' % (virtual_sink, '*', 0)
    fp_out.write(line + '\n')
    sources = [n for n in haplotig_graph.nodes() if haplotig_graph.in_degree(n) == 0]
    sinks = [n for n in haplotig_graph.nodes() if haplotig_graph.out_degree(n) == 0]
    for v in sources:
        line = 'L\t%s\t+\t%s\t+\t0M' % (virtual_source, v)
        fp_out.write(line + '\n')
    for v in sinks:
        line = 'L\t%s\t+\t%s\t+\t0M' % (v, virtual_sink)
        fp_out.write(line + '\n')


if __name__ == '__main__':  # pragma: no cover
    pass
