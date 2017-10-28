import argparse
import os
import sys

import networkx as nx
from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.FastaReader import FastaReader
from falcon_kit.gfa_graph import GFAGraph
import falcon_kit.mains.gen_gfa_v1 as gen_gfa_v1
from falcon_kit.util.system import find_files


def gfa_from_unzip(fp_out, preads_fasta, p_ctg_fasta, h_ctg_fasta, add_string_graph, unzip_root, write_reads, write_contigs, min_p_len, min_h_len):
    """
    This method produces the GFA-1 formatted output of the
    FALCON assembly.
    The graphical output is produced from either the entire string
    graph (only the non-filtered edges are considered) or from only
    the tiling paths. String graph can show the neighborhood of contig
    breaks, whereas the tiling path output is more sparse.
    Output is written to stdout.
    """

    hasm_dir = os.path.join(unzip_root, '1-hasm')

    # Haplotig tiling paths are not deduplicated.
    # We need the headers of the final haplotigs to filter
    # out the unnecessary tiling paths.
    h_ctg_headers = set()
    f = FastaReader(h_ctg_fasta)
    for r in f:
        h_ctg_headers.add(r.name)

    gfa_graph = GFAGraph()

    # Find and parse all primary contig files.
    for p_ctg_tp_file in find_files(hasm_dir, 'p_ctg_path.*'):
        p_paths, p_edge_to_ctg = gen_gfa_v1.load_tiling_paths(p_ctg_tp_file, 'P')
        _, p_ctg_len = gen_gfa_v1.calc_tiling_paths_len(p_paths)
        p_paths = gen_gfa_v1.filter_tiling_paths_by_len(p_paths, p_ctg_len, min_p_len)
        for ctg_id, path in p_paths.iteritems():
            gfa_graph.add_tiling_path(path, ctg_id)

    # Find and parse all haplotig files.
    for h_ctg_tp_file in find_files(hasm_dir, 'h_ctg_path.*'):
        h_paths_all, h_edge_to_ctg = gen_gfa_v1.load_tiling_paths(h_ctg_tp_file, 'H')
        # Keep only paths which have contigs present in all_h_ctg.fa.
        h_paths = {}
        for k, path in h_paths_all.iteritems():
            if k in h_ctg_headers:
                h_paths[k] = path
        _, h_ctg_len = gen_gfa_v1.calc_tiling_paths_len(h_paths)
        h_paths = gen_gfa_v1.filter_tiling_paths_by_len(h_paths, h_ctg_len, min_h_len)
        for ctg_id, path in h_paths.iteritems():
            gfa_graph.add_tiling_path(path, ctg_id)

    if add_string_graph:
        # Load the unzip graphs.
        for gexf_file in find_files(hasm_dir, 'sg.gexf'):
            sg = nx.read_gexf(gexf_file)
            gfa_graph.add_nx_string_graph(sg)

    gfa_graph.write_gfa_v1(fp_out, preads_fasta, [p_ctg_fasta, h_ctg_fasta], write_reads, write_contigs)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Generates GFA output (on stdout) from FALCON's assembly.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--preads-fasta', type=str, default='preads4falcon.fasta',
        help='path to the preads4falcon.fasta file')
    parser.add_argument(
        '--p-ctg-fasta', type=str, default='all_p_ctg.fa', help='path to the primary contigs file')
    parser.add_argument(
        '--h-ctg-fasta', type=str, default='all_h_ctg.fa', help='path to the haplotig file')
    parser.add_argument(
        '--unzip-root', type=str, default='./', help='path to the 3-unzip directory')
    parser.add_argument(
        '--add-string-graph', action='store_true',
        help="in addition to tiling paths, output other edges and nodes from the final string graph")
    parser.add_argument(
        '--write-reads', '-r', action='store_true', help="output read sequences in S lines")
    parser.add_argument(
        '--write-contigs', '-c', action='store_true', help="output contig sequences as S lines")
    parser.add_argument(
        '--min-p-len', type=int, default=0,
        help='primary contig paths with length smaller than this will not be reported')
    parser.add_argument(
        '--min-h-len', type=int, default=0,
        help='haplotig paths with length smaller than this will not be reported')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    gfa_from_unzip(sys.stdout, **vars(args))


if __name__ == '__main__':
    main()
