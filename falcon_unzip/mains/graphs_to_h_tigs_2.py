"""
I think this has implicit dependencies:
* 3-unzip/2-hasm/p_ctg_tiling_path
* and?
"""
from ..proto import (cigartools, execute, sam2m4, haplotig as Haplotig)
from ..proto.haplotig import Haplotig
from ..proto.graphs_to_h_tigs_2_utils import (
        mkdir,
        extract_unphased_haplotig_paths,
        extract_weakly_unphased_haplotig_paths,
        load_aln, nx_to_gfa,
        load_all_seq, load_sg_seq,
        revcmp_seq, reverse_end,
        write_haplotigs,
)
from falcon_kit import tiling_path
from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.FastaReader import FastaReader
from falcon_kit import io # (de)serialize()
from pypeflow.io import cd
from multiprocessing import Pool
import collections
import copy
import json
import networkx as nx
import os
import re
import intervaltree
import argparse
import sys
import logging

# for shared memory usage
global p_asm_G
global h_asm_G
global all_rid_to_phase
global all_flat_rid_to_phase
global all_haplotigs_for_ctg
global sg_edges
global p_ctg_seqs
global p_ctg_tiling_paths
global LOG
LOG = logging.getLogger() # root, to inherit from sub-loggers

"""
aln = aln_dict[htig_name]
qname, tname = aln[0:2]
score, perc_similarity = aln[2:4]
qstrand, qstart, qend, qlen = aln[4:8]
tstrand, tstart, tend, tlen = aln[8:12]
"""

#######################################################
### The main method for processing a single ctg_id. ###
#######################################################
def run_generate_haplotigs_for_ctg(input_):
    ctg_id, proto_dir, out_dir, base_dir, allow_multiple_primaries = input_
    LOG.info('Entering generate_haplotigs_for_ctg(ctg_id={!r}, out_dir={!r}, base_dir={!r}'.format(
        ctg_id, out_dir, base_dir))
    mkdir(out_dir)

    # Prepare the log output stream for this particular contig.
    logger = logging.getLogger(ctg_id)
    log_fn = os.path.join(out_dir, 'prototype.log')
    LOG.info('New logging FileHandler: {!r}'.format(os.path.abspath(log_fn)))
    hdlr = logging.FileHandler(log_fn, mode='w') #, level=logging.INFO)
    hdlr.setFormatter(logging.Formatter('%(levelname)s: %(message)s')) # Comment this out if you want.
    logger.addHandler(hdlr)
    hdlr.setLevel(logging.DEBUG) # Set to INFO someday?

    try:
        unzip_dir = os.path.join(base_dir, '3-unzip')
        return generate_haplotigs_for_ctg(ctg_id, allow_multiple_primaries, out_dir, unzip_dir, proto_dir, logger)
    except Exception:
        LOG.exception('Failure in generate_haplotigs_for_ctg({!r})'.format(input_))
        raise
    finally:
        logger.removeHandler(hdlr)
        # Note: The logger itself remains registered, but without a handler.
        # Note that the logger was registered in a thread,
        # so the main-thread sees only the top loggers, not this one.
        # Each thread accumulates some portion of the old loggers in disuse.
        # This will not cause any problems, as there is never a O(n) logger search.

def generate_haplotigs_for_ctg(ctg_id, allow_multiple_primaries, out_dir, unzip_dir, proto_dir, logger):
    # proto_dir is specific to this ctg_id.

    global all_haplotigs_for_ctg
    global p_ctg_seqs
    global sg_edges
    global p_ctg_tiling_paths

    fp_proto_log = logger.info

    # min_linear_len = 8

    fp_proto_log('Started processing contig: "{}".'.format(ctg_id))

    #########################################################
    # Load and prepare data.
    #########################################################
    fp_proto_log('Fetching the p_ctg_seq.')
    p_ctg_seq = p_ctg_seqs[ctg_id]

    fp_proto_log('Fetching the p_ctg_tiling_path.')
    p_ctg_tiling_path = p_ctg_tiling_paths[ctg_id]

    # Path to the directory with the output of the main_augment_pb.py.
    #proto_dir = os.path.join(unzip_dir, '0-phasing', ctg_id, 'proto')

    # Load the linear sequences for alignment.
    minced_ctg_path = os.path.join(proto_dir, 'minced.fasta')

    fp_proto_log('Loading minced ctg seqs from {!r} .'.format(minced_ctg_path))
    minced_ctg_seqs = load_all_seq(minced_ctg_path)

    # Load the phase relation graph, and build a dict of the weakly connected components.
    phase_relation_graph_path = os.path.join(proto_dir, "phase_relation_graph.gexf")
    phase_relation_graph = nx.read_gexf(phase_relation_graph_path)
    fp_proto_log('Loading the phase relation graph from {!r} .'.format(phase_relation_graph_path))
    phase_alias_map, alias_to_phase_list = create_phase_alias_map(phase_relation_graph)

    # Load the regions as built previous in the 0-phasing step.
    # The structure of each region in the list is:
    #   type_, edge_start_id, edge_end_id, pos_start, pos_end, htigs = all_regions[i]
    all_regions_path = os.path.join(proto_dir, 'regions.json')
    fp_proto_log('Loading all regions from {!r} .'.format(all_regions_path))
    all_regions = json.load(open(all_regions_path))

    # Create a list of bubble regions, for easier merging with the diploid groups.
    fp_proto_log('Making bubble region list.')
    bubble_region_list = []
    for region in all_regions:
        region_type, first_edge, last_edge, pos_start, pos_end, region_htigs = region
        if region_type != 'linear':
            bubble_region_list.append(region)

    # The json module converts a tuple to list, so revert back here.
    fp_proto_log('Retupling.')
    for region in all_regions:
        region_type, first_edge, last_edge, pos_start, pos_end, region_htigs = region
        for htig_name, htig in region_htigs.iteritems():
            htig['phase'] = tuple(htig['phase'])

    # Assign sequences to all regions.
    fp_proto_log('Assigning sequences to all regions.')
    for region in all_regions:
        region_type, first_edge, last_edge, pos_start, pos_end, region_htigs = region
        # Assigning sequences to the regions. This could have also been stored in the
        # json/msgpack file where the regions are loaded from, the code would be
        # simpler then.
        for htig_name, htig in region_htigs.iteritems():
            htig['seq'] = minced_ctg_seqs[htig_name]

    # Get all the haplotig objects corresponding to this contig only.
    fp_proto_log('Getting snp_haplotigs.')
    snp_haplotigs = all_haplotigs_for_ctg.get(ctg_id, {})

    #########################################################

    #########################################################
    # Debug verbose.
    #########################################################
    logger.debug('phase_alias_map = {}'.format(phase_alias_map))

    logger.debug('Verbose all regions loaded from regions.json:')
    for region_id in xrange(len(all_regions)):
        region = all_regions[region_id]
        region_type, first_edge, last_edge, pos_start, pos_end, region_htigs = region
        logger.debug('[region_id = {}] type = {}, first_edge = {}, last_edge = {}, pos_start = {}, pos_end = {}'.format(region_id, region_type, first_edge, last_edge, pos_start, pos_end))
    logger.debug('') # newline
    #########################################################

    #########################################################
    # Write the haplotigs to disk for alignment.
    #########################################################
    aln_snp_hasm_ctg_path = os.path.join(out_dir, "aln_snp_hasm_ctg.fasta")
    fp_proto_log('Writing haplotigs to disk: {!r}'.format(aln_snp_hasm_ctg_path))
    write_haplotigs(snp_haplotigs, aln_snp_hasm_ctg_path, fp_proto_log, hack_qnames_for_blasr=True)

    #########################################################
    # Align haplotigs.
    #########################################################
    num_threads = 16
    mapping_out_prefix = os.path.join(out_dir, 'aln_snp_hasm_ctg')
    mapping_ref = os.path.join(unzip_dir, 'reads', ctg_id, 'ref.fa')
    sam_path = mapping_out_prefix + '.sam'

    def excomm(cmd):
        execute.execute_command(cmd, logger)

    if len(snp_haplotigs.keys()) > 0:
        # BLASR crashes on empty files, so address that.
        blasr_params = '--minMatch 15 --maxMatch 25 --advanceHalf --advanceExactMatches 10 --bestn 1 --nproc {} --noSplitSubreads'.format(num_threads)
        excomm('blasr {} {} {} --sam --out {}.tmp.sam'.format(
                blasr_params, aln_snp_hasm_ctg_path, mapping_ref, mapping_out_prefix))
        excomm('samtools sort {pre}.tmp.sam -o {pre}.sam'.format(
                pre=mapping_out_prefix))
        excomm('rm -f {pre}.tmp.sam'.format(
                pre=mapping_out_prefix))
    aln_dict = load_and_hash_sam(sam_path, fp_proto_log)

    #########################################################
    # Filter out overlapping haplotigs for each phasing block.
    #########################################################
    # This filters out any overlapping haplotig. We won't have that.
    filtered_snp_haplotigs = snp_haplotigs

    # After haplotigs are aligned, we know which ones are not facing the
    # fwd direction of the corresonding contig. We need to reverse them.
    # Edges need to be looked up from the sg_edges_list for the rev_cmp.
    reorient_haplotigs(filtered_snp_haplotigs, aln_dict, sg_edges, fp_proto_log)

    # Find all of the locations where the haplotigs need to be broken.
    clippoints, bubble_tree = collect_clippoints(all_regions, filtered_snp_haplotigs, aln_dict, fp_proto_log)

    # Fragment the haplotigs on all clippoints.
    fragmented_snp_haplotigs = fragment_haplotigs(filtered_snp_haplotigs, aln_dict, clippoints, bubble_tree, fp_proto_log)

    # Take the fragmented haplotigs, and extract only diploid pairs as regions.
    diploid_region_list = create_diploid_regions(fragmented_snp_haplotigs, fp_proto_log)

    # Collect all bubbles.
    final_phased_regions = diploid_region_list + bubble_region_list

    # Create the linear regions in between bubbles.
    final_linear_regions = get_linear_regions_between_bubbles(ctg_id, final_phased_regions, p_ctg_seq, p_ctg_tiling_path, fp_proto_log)

    # Merge all regions, bubble and linear.
    final_all_regions = final_phased_regions + final_linear_regions
    # final_all_regions = sorted(final_all_regions, key = lambda x: x[3])

    # Convert the linear region list to a graph.
    fp_proto_log('Creating a haplotig graph.')
    haplotig_graph = regions_to_haplotig_graph(ctg_id, final_all_regions, fp_proto_log)
    haplotig_graph = update_haplotig_graph(haplotig_graph, phase_alias_map)

    # Write the nx graph to disk for debugging.
    fp_proto_log('  - Writing the haplotig graph in the gexf format.')
    nx.write_gexf(haplotig_graph, os.path.join(out_dir, "haplotig_graph.gexf"))

    # Write the haplotig graph to disk.
    fp_proto_log('Writing the haplotig_graph.gfa.')
    with open(os.path.join(out_dir, "haplotig_graph.gfa"), 'w') as fp_out:
        nx_to_gfa(ctg_id, haplotig_graph, fp_out)

    # Extract the p_ctg and h_ctg sequences.
    unzipped_p_ctg_seqs, unzipped_p_ctg_edges, unzipped_h_ctg_seqs, unzipped_h_ctg_edges, unzipped_h_ctg_paf = \
                    extract_unzipped_ctgs(ctg_id, haplotig_graph, allow_multiple_primaries, fp_proto_log)

    # Write the results out.
    write_unzipped(out_dir, ctg_id, unzipped_p_ctg_seqs, unzipped_p_ctg_edges,
                    unzipped_h_ctg_seqs, unzipped_h_ctg_edges, unzipped_h_ctg_paf, fp_proto_log)

    #########################################################
    # Debug verbose.
    #########################################################
    for line in sorted(clippoints.iteritems()):
        logger.debug(' clip_point: {}'.format(line))

    logger.debug('Fragmented haplotigs:')
    for hname, htig in fragmented_snp_haplotigs.iteritems():
        logger.debug('  - name = {}, phase = {}, path = {}'.format(htig.name, htig.phase, htig.path))

    logger.debug('Verbose the generated diploid regions:')
    for region_id in xrange(len(final_all_regions)):
        region = final_all_regions[region_id]
        region_type, first_edge, last_edge, pos_start, pos_end, region_htigs = region
        logger.debug('[region_id = {}] type = {}, first_edge = {}, last_edge = {}, pos_start = {}, pos_end = {}'.format(region_id, region_type, first_edge, last_edge, pos_start, pos_end))
    logger.debug('')
    #########################################################

    fp_proto_log('Finished processing contig: %s.' % (ctg_id))

def load_haplotigs(hasm_falcon_path, all_flat_rid_to_phase):
    """
    Returns a dict `haplotigs_for_ctg`, where:
        haplotigs_for_ctg[ctg_id][htig.name] = Haplotig(...)
    """

    haplotigs_for_ctg = {}

    # The haplotigs are named as FALCON contigs, and will have
    # to be renamed. This dict keeps the relation between
    # the new haplotig name and the original p_ctg name from the
    # 3-unzip/1-hasm/asm-falcon stage.
    #
    htig_name_to_original_pctg = {}

    # Load the primary contig tiling paths.
    p_tiling_paths = os.path.join(hasm_falcon_path, 'p_ctg_tiling_path')
    tiling_paths = tiling_path.load_tiling_paths(p_tiling_paths, contig_lens=None, whitelist_seqs=None)

    p_ctg_fasta = os.path.join(hasm_falcon_path, 'p_ctg.fa')
    hasm_p_ctg_seqs = load_all_seq(p_ctg_fasta)

    LOG.info('Loading haplotigs.')

    # Process the tiling paths for each assembled haplotig.
    counter = io.Percenter('tiling_paths', len(tiling_paths))
    num_hctg = 0
    for hctg, hpath in tiling_paths.iteritems():
        num_hctg += 1

        # Sanity check. This should not happen, but if it does,
        # skip the entire haplotig.
        if len(hpath.edges) == 0:
            continue

        # Get the first node of the tiling path.
        vrid, vorient = hpath.edges[0].v.split(':')
        # Sanity check. This should not happen, but if it does,
        # skip the entire haplotig.
        if vrid not in all_flat_rid_to_phase:
            continue
        # By definition, all reads should have the same phase.
        vphase = all_flat_rid_to_phase[vrid]
        phase_ctg, phase_block, phase_id = vphase

        # Sanity check that all the nodes of the tiling path
        # actually of the same phase. Otherwise, there
        # was a bug in the overlap filter.
        for e in hpath.edges:
            wrid, worient = e.w.split(':')
            assert(wrid in all_flat_rid_to_phase)
            wphase = all_flat_rid_to_phase[wrid]
            assert(wphase == vphase)

        # Skip the circular haplotigs, since this makes no sense in the diploid context.
        # This is legal and valid output from FALCON's graph_to_contig.py, but in the
        # context of Unzipping a linear assembly, these are an artifact of the input
        # data (e.g. missing adapters).
        if hpath.edges[0].v == hpath.edges[-1].w:
            LOG.debug('[{}] Skipping haplotig because it\'s circular: {}'.format(num_hctg, hctg))
            continue

        # Sanity check for a case which should absolutely never happen in practice, but better
        # be safe and log + skip it.
        # If the same node is used in multiple places of the same tiling path,
        # there is some degeneracy here and we should skip the entire haplotig.
        occ_count_dict = collections.defaultdict(int)
        do_skip_haplotig = False
        for e in hpath.edges:
            occ_count_dict[e.v] += 1
            occ_count_dict[e.w] += 1
            if occ_count_dict[e.v] > 2 or occ_count_dict[e.w] > 2:
                do_skip_haplotig = True
        if do_skip_haplotig == True:
            LOG.debug('[{}] Skipping an invalid haplotig because contains nodes used multiple times in the path: {}'.format(num_hctg, hctg))
            continue

        ###########################
        # Construct the haplotig. #
        ###########################
        # Info needed to construct a haplotig.
        # Seq is generated from path, to avoid proper/improper confusion.
        h_tig_name = '%s-HAP%s-%s.%s.%s' % (phase_ctg, hctg, phase_ctg, str(phase_block), str(phase_id))
        counter(1, label=h_tig_name)
        complete_phase = vphase
        path = hpath.dump_as_split_lines()
        seq = hasm_p_ctg_seqs[hctg] # The seq should contain the first read (be `proper`).

        new_haplotig = Haplotig(name = h_tig_name, phase = complete_phase, seq = seq, path = path, edges = [])

        # The name relation dicts.
        htig_name_to_original_pctg[h_tig_name] = hctg
        # original_pctg_to_htig_name[hctg] = h_tig_name

        def verbose_haplotig(haplotig):
            return 'phase = %s, h_tig_name = %s, num_edges = %d, len(seq) = %d' % (str(haplotig.phase), haplotig.name, len(haplotig.path), len(haplotig.seq))

        LOG.debug('[{}] Loaded haplotig: {}'.format(num_hctg, verbose_haplotig(new_haplotig)))

        # Append the haplotig to the right place.
        # Secondary haplotigs for any phase will be filtered later.
        haplotigs_for_ctg.setdefault(phase_ctg, {})
        haplotigs_for_ctg[phase_ctg].setdefault(h_tig_name, {})
        haplotigs_for_ctg[phase_ctg][h_tig_name] = new_haplotig

    return haplotigs_for_ctg, htig_name_to_original_pctg

def load_and_hash_sam(sam_path, fp_proto_log):
    """
    Loads the SAM file, converts it to the M4 format, and
    builds a dict of aln_dict[qname] = aln_m4.
    Keeps only one (longest) alignment per qname.
    """
    aln_dict = {}
    if not os.path.exists(sam_path):
        return aln_dict
    fp_proto_log('Loading the alignments.')
    m4 = load_aln(sam_path)
    m4 = sorted(m4, key = lambda x: x[0].split('-')[-1])    # Not required, but simplifies manual debugging.
    for aln in m4:
        if aln[0] not in aln_dict:
            aln_dict[aln[0]] = aln
        else:
            qstart, qend = aln[5], aln[6]
            qdist = qend - qstart
            prev_aln = aln_dict[aln[0]]
            prev_qstart, prev_qend = prev_aln[5], prev_aln[6]
            prev_qdist = prev_qend - prev_qstart
            if qdist > prev_qdist:
                aln_dict[aln[0]] = aln
    return aln_dict

def create_phase_alias_map(phase_relation_graph):
    alias_to_phase_list = []
    phase_alias_map = {}
    for subg in nx.weakly_connected_component_subgraphs(phase_relation_graph):
        # Make all phases point to the same ID.
        num_aliases = len(alias_to_phase_list)
        for v in subg.nodes():
            phase_alias_map[v] = num_aliases
        alias_to_phase_list.append(set(subg.nodes()))
    return phase_alias_map, alias_to_phase_list

def reverse_sg_path(sg_path, sg_edges):
    new_path = []
    for edge in sg_path:
        ctg_id, v, w = edge[0:3]
        rv = reverse_end(w)
        rw = reverse_end(v)
        # sg_edges example:         000000315:E 000000076:E 000000076 7699 7927 7696 99.82 G
        # Tiling path edge example: 000001F 000000315:E 000000076:E 000000076 7699 7927 7696 99.82
        sg_edge = sg_edges[(rv, rw)]
        new_tp_edge = [ctg_id] + sg_edge[0:7] # , sg_edge[0], sg_edge[1], sg_edge[2], int(sg_edge[3]), int(sg_edge[4]), int(sg_edge[5]), float(sg_edge[6])]
        new_path.append(new_tp_edge)
    new_path = new_path[::-1]
    return new_path

def reorient_haplotigs(snp_haplotigs, aln_dict, sg_edges, fp_proto_log):
    """
    Changes the orientation of paths and alignments for haplotigs
    which align toreverse strand.
    Changes are made in-place.
    """
    fp_proto_log('Reorienting haplotigs.')
    for qname, haplotig in snp_haplotigs.iteritems():
        if qname not in aln_dict:
            continue
        aln = aln_dict[qname]
        qstrand, qstart, qend, qlen = aln[4:8]
        tstrand, tstart, tend, tlen = aln[8:12]
        if tstrand == 1:
            fp_proto_log('  - qname = {}'.format(qname))
            # fp_proto_log('  - Before reversing:')
            # for line in haplotig.path:
            #     fp_proto_log('  {}' %.format(line))

            # Reverse the path and the seq.
            haplotig.path = reverse_sg_path(haplotig.path, sg_edges)
            haplotig.seq = revcmp_seq(haplotig.seq)
            fp_proto_log('')
            # Reset the reverse flag.
            # Coordinates in aln are on the fwd strand anyway,
            # so they remain the same.
            aln[8] = 0

            # fp_proto_log('  - After reversing:')
            # for line in haplotig.path:
            #     fp_proto_log('  {}'.format(line))

def collect_clippoints(all_regions, snp_haplotigs, aln_dict, fp_proto_log):
    """
    Creates a dict of all positions where haplotig alignments should be broken.
    This includes:
        - Beginning and ending of bubble regions.
        - Beginning and ending of other haplotigs (e.g. if haplotig for (000000F, 3000001, 0) is longer than (000000F, 3000001, 1).
    """
    clippoints = {}

    # Add the coordinates of the bubble regions.
    bubble_intervals = []
    for region in all_regions:
        region_type, first_edge, last_edge, pos_start, pos_end, region_htigs = region
        if region_type != 'linear':
            clippoints[pos_start] = [region_type + '_start']
            clippoints[pos_end] = [region_type + '_end']
            bubble_intervals.append(intervaltree.Interval(pos_start, pos_end, region))
    bubble_tree = intervaltree.IntervalTree(bubble_intervals)

    for qname, haplotig in snp_haplotigs.iteritems():
        if qname not in aln_dict:
            continue
        aln = aln_dict[qname]
        qstrand, qstart, qend, qlen = aln[4:8]
        tstrand, tstart, tend, tlen = aln[8:12]
        clippoints[tstart] = [qname + '_start']
        clippoints[tend] = [qname + '_end']

    return clippoints, bubble_tree

def fragment_haplotigs(in_haplotigs, aln_dict, clippoints, bubble_tree, fp_proto_log):
    filtered_haplotigs = {}

    for qname, haplotig in in_haplotigs.iteritems():
        if qname not in aln_dict:
            continue
        aln = aln_dict[qname]
        new_haplotigs = fragment_single_haplotig(haplotig, aln, clippoints, bubble_tree, fp_proto_log)
        filtered_haplotigs.update(new_haplotigs)

    return filtered_haplotigs

def fragment_single_haplotig(haplotig, aln, clippoints, bubble_tree, fp_proto_log):
    """
    This method takes a haplotig and a list of coordinates, and breaks
    the haplotig based on the CIGAR string alignment.
    This method expects that the haplotig.seq is already oriented in the fwd
    direction of the target, and that the M4 alignment already reflects that.
    """
    q_name, t_name = aln[0], aln[1]
    q_orient, q_start, q_end, q_len = aln[4], aln[5], aln[6], aln[7]
    t_orient, t_start, t_end, t_len = aln[8], aln[9], aln[10], aln[11]
    sam = aln[13]

    aln_array = cigartools.cigar_to_aln_array(sam.cigar)
    positions = cigartools.find_positions(aln_array, t_start)

    # For each aligned position, check if it's in the clippoints.
    # If it is, store the position in a new set.
    pos_set = set()
    for t_pos, q_pos_array in positions.iteritems():
        if t_pos not in clippoints:
            continue
        # Multiple bases could be stacked at the same t_pos,
        # namely, insertions.
        q_pos = q_pos_array[0]
        pos_set.add((t_pos, q_pos))
    pos_set.add((t_start, q_start))
    pos_set.add((t_end, q_end))

    all_pos_of_interest = sorted(list(pos_set))

    regions_of_interest = []
    if len(all_pos_of_interest) > 1:
        for start, end in zip(all_pos_of_interest[:-1], all_pos_of_interest[1:]):
            fp_proto_log(' start = {}, end = {}'.format(start, end))
            tstart, tend = start[0], end[0]
            found_intervals = bubble_tree.search(tstart, tend)
            if found_intervals:
                continue
            regions_of_interest.append((start, end, q_name, q_len, t_name, t_len, haplotig.phase))

    fp_proto_log('pos_of_interest for q_name: {}'.format(q_name))
    for line in regions_of_interest:
        fp_proto_log('{}'.format(line))
    fp_proto_log('')

    # Here we extract the subsequence and the subpath of the tiling path.
    ret_haplotigs = {}
    tp_ctg_id = None if len(haplotig.path) == 0 else haplotig.path[0][0]    # The original "contig" ID from 3-unzip/1-hasm is still used in the tiling path.
    tp_dict = tiling_path.load_tiling_paths_from_split_lines(haplotig.path, contig_lens={tp_ctg_id: len(haplotig.seq)}, whitelist_seqs=None)
    tp = tp_dict[tp_ctg_id]

    for region_of_interest in regions_of_interest:
        start, end, q_name, q_len, t_name, t_len, q_phase = region_of_interest
        qstart, qend = start[1], end[1]

        new_name = haplotig.name + '-%d' % (len(ret_haplotigs))
        new_seq = haplotig.seq[qstart:qend]
        new_path, new_start_coord, new_end_coord = tp.get_subpath(qstart, qend)
        # new_path = [e.split_line for e in new_path] # Convert the tiling path back to a list of split values.
        new_haplotig = Haplotig(name = new_name, phase = haplotig.phase, seq = new_seq, path = new_path, edges = [])
        # Encode additional attributes.
        new_haplotig.labels['region_of_interest'] = region_of_interest
        new_haplotig.labels['start_in_path'] = new_start_coord
        new_haplotig.labels['end_in_path'] = new_end_coord

        ret_haplotigs[new_name] = new_haplotig

    return ret_haplotigs

def create_diploid_regions(fragmented_snp_haplotigs, fp_proto_log):
    """
    Takes the fragmented haplotigs (based on the alignment coordinates and bubble
    regions), and piles them up based on their coordinates and phasing blocks.
    If there are more than 2 haplotigs per coordinate span, the region is filtered
    out.
    If there are exactly 2 haplotigs per coordinate span, and they are not of
    different phase_id (e.g. when two haplotigs of the same `(ctg_id, phasing_block, phase_id)`
    cover the same region), this is possibly due to a repeat and is also filtered out.
    Any region which is covered by only 1 haplotig is also filtered out.

    This method returns a list of diploid regions, where each haplotig in the region
    perfectly aligns to the start/end coordinates on the p_ctg.fa.
    """

    ret_diploid_regions = []

    groups = {}
    for qname, haplotig in fragmented_snp_haplotigs.iteritems():
        phase = haplotig.phase
        start, end, q_name, q_len, t_name, t_len, q_phase = haplotig.labels['region_of_interest']
        tstart, tend = start[0], end[0]
        qstart, qend = start[1], end[1]

        region = (tstart, tend)
        phase_group = (phase[0], phase[1])
        phase_id = phase[2]
        groups.setdefault(region, {})
        groups[region].setdefault(phase_group, {})
        groups[region][phase_group].setdefault(phase_id, [])
        groups[region][phase_group][phase_id].append(haplotig)

    for region, phase_group_dict in groups.iteritems():
        is_valid = True

        # Skip any region covered by more than 1 phasing block. Possible repeat.
        if len(phase_group_dict.keys()) != 1:
            is_valid = False
            continue

        # Skip any group which doesn't have exactly 2 phases in it.
        # E.g. one haplotig was larger than the other, due to read length.
        for phase_group, phase_id_dict in phase_group_dict.iteritems():
            if len(phase_id_dict.keys()) != 2:
                is_valid = False
                break
            # Skip any phase which has more than 1 haplotig covering the region. Possible repeat.
            for phase_id, haplotigs in phase_id_dict.iteritems():
                if len(haplotigs) != 1:
                    is_valid = False
                    break
            if is_valid == False:
                break

        if is_valid == False:
            continue


        # If we got here, then this is a valid diploid group (exactly 2 phases,
        # and exactly 1 haplotig per phase).

        # Get the placement coordinates on the original collapsed primary contig.
        pos_start, pos_end = region

        # Collect the haplotigs.
        region_htigs = {}
        for phase_group, phase_id_dict in phase_group_dict.iteritems():
            for phase_id, haplotigs in phase_id_dict.iteritems():
                for haplotig in haplotigs:
                    haplotig.cstart = pos_start
                    haplotig.cend = pos_end
                    region_htigs[haplotig.name] = haplotig.__dict__

        # Create a diploid region.
        region_type = 'diploid'
        first_edge = None           # First edge on the p_ctg. Not available right now.
        last_edge = None            # Last edge on the p_ctg. Not available right now.
        new_region = (region_type, first_edge, last_edge, pos_start, pos_end, region_htigs)

        ret_diploid_regions.append(new_region)

    return ret_diploid_regions

def make_linear_region(ctg_id, pos_start, pos_end, p_ctg_seq, p_ctg_tiling_path, fp_proto_log):
    """
    This is a helper method used by get_linear_regions_between_bubbles.
    Formats a new region object for a linear (collapsed region), specified by the
    start and end coordinate on the primary contig.
    This method extracts the subsequence from the primary contig, and also
    the subpath of the tiling path. Everything is encoded in the new region object
    (the sequence and path are part of the Haplotig object).
    """

    fp_proto_log('Entered function: "{}"'.format(__name__))
    if (pos_end - pos_start) <= 0:
        fp_proto_log('Exiting function: "{}"'.format(__name__))
        return None

    region_type = 'linear'
    complete_phase = (ctg_id, '-1', '0')

    htig_name = '%s_%s2_%d:%d_base' % (ctg_id, region_type, pos_start, pos_end)

    fp_proto_log('Info: pos_start = {}, pos_end = {}, seq len: {}'.format(pos_start, pos_end, len(p_ctg_seq)))

    fp_proto_log('Extracting subpath.')
    new_path, new_start_coord, new_end_coord = p_ctg_tiling_path.get_subpath(pos_start, pos_end)

    fp_proto_log('Extracting the seq.')
    new_seq = p_ctg_seq[pos_start:pos_end]

    fp_proto_log('len(new_path) = {}'.format(len(new_path)))
    first_edge = new_path[0]
    last_edge = new_path[-1]

    fp_proto_log('Forming the haplotig.')
    new_haplotig = Haplotig(name = htig_name, phase = complete_phase, seq = new_seq, path = new_path, edges = [])
    new_haplotig.labels['start_in_path'] = new_start_coord
    new_haplotig.labels['end_in_path'] = new_end_coord
    new_haplotig.cstart = pos_start
    new_haplotig.cend = pos_end

    region_htigs = {htig_name: new_haplotig.__dict__}
    new_region = (region_type, first_edge, last_edge, pos_start, pos_end, region_htigs)

    fp_proto_log('Exiting function: "{}"'.format(__name__))
    return new_region

def get_linear_regions_between_bubbles(ctg_id, bubble_regions, p_ctg_seq, p_ctg_tiling_path, fp_proto_log):
    """
    For a given set of non-overlapping bubble regions, creates a set of linear regions
    at the suffix, prefix and in between the bubbles, and returns them as a list.
    """
    ret_linear_regions = []

    sorted_bubble_regions = sorted(bubble_regions, key = lambda x: x[3])

    fp_proto_log('Function: "{}"'.format(__name__))
    fp_proto_log('len(sorted_bubble_regions) = {}'.format(len(sorted_bubble_regions)))

    if len(sorted_bubble_regions) == 0:
        pos_start = 0
        pos_end = len(p_ctg_seq)
        new_region = make_linear_region(ctg_id, pos_start, pos_end, p_ctg_seq, p_ctg_tiling_path, fp_proto_log)
        if new_region != None:
            ret_linear_regions.append(new_region)
        return ret_linear_regions

    # Handle the prefix linear region.
    fp_proto_log('Handling prefix.')
    region1_type, region1_first_edge, region1_last_edge, region1_pos_start, region1_pos_end, region1_htigs = sorted_bubble_regions[0]
    pos_start = 0
    pos_end = region1_pos_start
    new_region = make_linear_region(ctg_id, pos_start, pos_end, p_ctg_seq, p_ctg_tiling_path, fp_proto_log)
    if new_region != None:
        ret_linear_regions.append(new_region)

    fp_proto_log('Handling infix.')
    if len(sorted_bubble_regions) > 1:
        temp_num_regions = 0
        for region1, region2 in zip(sorted_bubble_regions[0:-1], sorted_bubble_regions[1:]):
            temp_num_regions += 1
            fp_proto_log('Checking for in-between region {} / {}.'.format(temp_num_regions, len(sorted_bubble_regions)))
            region1_type, region1_first_edge, region1_last_edge, region1_pos_start, region1_pos_end, region1_htigs = region1
            region2_type, region2_first_edge, region2_last_edge, region2_pos_start, region2_pos_end, region2_htigs = region2
            pos_start = region1_pos_end
            pos_end = region2_pos_start
            new_region = make_linear_region(ctg_id, pos_start, pos_end, p_ctg_seq, p_ctg_tiling_path, fp_proto_log)
            if new_region != None:
                fp_proto_log('Generated the new region.')
                ret_linear_regions.append(new_region)

    # Handle the suffix linear region.
    fp_proto_log('Handling suffix.')
    region2_type, region2_first_edge, region2_last_edge, region2_pos_start, region2_pos_end, region2_htigs = sorted_bubble_regions[-1]
    pos_start = region2_pos_end
    pos_end = len(p_ctg_seq)
    new_region = make_linear_region(ctg_id, pos_start, pos_end, p_ctg_seq, p_ctg_tiling_path, fp_proto_log)
    if new_region != None:
        ret_linear_regions.append(new_region)

    fp_proto_log('Dunn!')

    return ret_linear_regions

def regions_to_haplotig_graph(ctg_id, all_regions, fp_proto_log):
    """
    Takes a list of non-overlapping (but adjacent) regions (basically, a DAG),
    and creates a NetworkX object which describes the relationship between all
    haplotigs.
    """
    all_regions = sorted(all_regions, key = lambda x: x[3])

    haplotig_graph = nx.DiGraph()

    # Add a dummy source node to keep fully-phased front of the contig in place.
    source_node_name = '{ctg_id}-source-node'.format(ctg_id=ctg_id)
    haplotig_graph.add_node(source_node_name, label='source',
                phase='{ctg_id}_-1_0'.format(ctg_id=ctg_id), phase_alias=-1, htig={})

    # Add a dummy sink node to keep fully-phased back of the contig in place.
    sink_node_name = '{ctg_id}-sink-node'.format(ctg_id=ctg_id)
    haplotig_graph.add_node(sink_node_name, label='sink',
                phase='{ctg_id}_-1_0'.format(ctg_id=ctg_id), phase_alias=-1, htig={})

    # Add nodes.
    fp_proto_log('  - Adding nodes.')
    for region_id in xrange(len(all_regions)):
        region = all_regions[region_id]
        region_type, edge_start, edge_end, region_pos_start, region_pos_end, region_haplotigs = region
        fp_proto_log('    - region_id = {}, region_type = {}, region_pos_start = {}, region_pos_end = {}'.format(region_id, region_type, region_pos_start, region_pos_end))
        for key, htig in region_haplotigs.iteritems():
            if region_type == 'complex' and key.endswith('_base') == False:
                # Keep the complex regions as collapsed.
                fp_proto_log('      - [haplotig graph, adding node] key = {} skipped because region_type == "{}" and key does not end in "_base"'.format(key, region_type))
                continue
            fp_proto_log('      - [haplotig graph, adding node] key = {}'.format(key))
            v = htig['name']
            vphase = htig['phase']
            vphase_alias = -1 # This is a placeholder, value `-1` means uninitialized. It will be used in `update_haplotig_graph`.
            haplotig_graph.add_node(v, label='%s_%s_%s_%s' % (ctg_id, region_type, region_pos_start, region_pos_end),
                        phase='_'.join([str(val) for val in vphase]),
                        phase_alias=vphase_alias, htig=htig)

    # Add edges between all regions. Filtering the edges will be performed afterwards.
    # First, add dummy edges.
    if len(all_regions) > 0:
        node_set = set(haplotig_graph.nodes())

        # Connect the first region to the source node.
        region = all_regions[0]
        region_type, edge_start, edge_end, region_pos_start, region_pos_end, region_haplotigs = region
        for htig_name, htig in region_haplotigs.iteritems():
            if htig_name not in node_set:
                continue
            v = source_node_name
            w = htig_name
            haplotig_graph.add_edge(v, w, weight=1)

        # Connect the last region to the sink node.
        region = all_regions[-1]
        region_type, edge_start, edge_end, region_pos_start, region_pos_end, region_haplotigs = region
        for htig_name, htig in region_haplotigs.iteritems():
            if htig_name not in node_set:
                continue
            v = htig_name
            w = sink_node_name
            haplotig_graph.add_edge(v, w, weight=1)

    # Second, add all proper edges.
    fp_proto_log('  - Adding edges.')
    for region_id in xrange(1, len(all_regions)):
        region = all_regions[region_id - 1]
        region_type, edge_start, edge_end, region_pos_start, region_pos_end, region_haplotigs = region
        next_region = all_regions[region_id]
        next_region_type, next_edge_start, next_edge_end, next_region_pos_start, next_region_pos_end, next_region_haplotigs = next_region
        for htig_name, htig in region_haplotigs.iteritems():
            v = htig_name

            if region_type == 'complex' and v.endswith('_base') == False:
                # Keep the complex regions as collapsed.
                fp_proto_log('      - [haplotig graph, adding edge] v = {} skipped because region_type == "{}" and key does not end in "_base"'.format(v, region_type))
                continue

            for next_htig_name, next_htig in next_region_haplotigs.iteritems():
                w = next_htig_name

                if next_region_type == 'complex' and w.endswith('_base') == False:
                    # Keep the complex regions as collapsed.
                    fp_proto_log('      - [haplotig graph, adding edge] w = {} skipped because next_region_type == "{}" and key does not end in "_base"'.format(w, next_region_type))
                    continue

                fp_proto_log('      - [haplotig graph, adding edge] v = {}, w = {}'.format(v, w))
                haplotig_graph.add_edge(v, w, weight=10000)

    return haplotig_graph

def update_haplotig_graph(haplotig_graph, phase_alias_map):
    """
    Iterative updating of a haplotig graph in the form of a
    NetworkX object.
    For every node, it checks the phase_alias_map and updates
    the alias info.
    For any node which has an in-alias edge, removes all other
    edges (cross-phase).
    """
    haplotig_graph2 = haplotig_graph.copy()

    # Update the aliases.
    for v in haplotig_graph2.nodes():
        # v = htig['name']
        vdata = haplotig_graph2.node[v]
        vphase = vdata['phase']
        vphase_alias = phase_alias_map.get(vphase, -1)
        vdata['phase_alias'] = vphase_alias

    edges_to_remove = set()
    for v in haplotig_graph2.nodes():
        # Check if there are any edges which resolve
        # to the same alias in the graph.
        matching_alias = False
        v_alias = haplotig_graph2.node[v]['phase_alias']
        nonmatching_edges = set()
        for vv, ww in haplotig_graph2.out_edges(v):
            ww_alias = haplotig_graph2.node[ww]['phase_alias']
            if ww_alias == v_alias:
                matching_alias = True
            else:
                nonmatching_edges.add((vv, ww))

        # Remove the edges between nodes of different phases,
        # but only if there is at least one edge which keeps
        # the same phase.
        if matching_alias:
            edges_to_remove.update(nonmatching_edges)

    # Actually remove the edges.
    for v, w in edges_to_remove:
        haplotig_graph2.remove_edge(v, w)

    # Update the weights of any remaining edges.
    for v, w in haplotig_graph2.edges():
        v_alias = haplotig_graph2.node[v]['phase_alias']
        w_alias = haplotig_graph2.node[w]['phase_alias']
        if v_alias == -1 and w_alias == -1:
            weight = 50
        elif v_alias == -1 or w_alias == -1:
            weight = 50
        elif v_alias == w_alias:
            weight = 1
        else:
            weight = 10000
        haplotig_graph2[v][w]['weight'] = weight

    return haplotig_graph2

def construct_ctg_seq(haplotig_graph, new_ctg_id, node_path):
    """
    Given a list of nodes which specify the path (e.g. as the output from
    the NetworkX shortest path), this function generates the contig sequence
    and the corresponding list of edges (and the path).
    Nodes here are the diploid or collapsed regions, so concatenating them
    is good by design.
    Inputs:
        haplotig_graph  - A NetworkX object.
        new_ctg_id      - Name of the new contig.
        node_path       - List of node IDs that make up a path in the haplotig_graph.
    Returns:
        new_ctg_seq         - String containing the sequence of the new contig.
        new_ctg_edges       - List of strings, where each string is one of the edges in the contig's tiling path.
        node_start_coords   - Dict where key is the node ID, and value start coordinate of the node in the path.
        node_end_coords     - Dict where key is the node ID, and value end coordinate of the node in the path.
    """

    # Sanity check.
    node_set = set(haplotig_graph.nodes())
    for v in node_path:
        if v not in node_set:
            msg = 'Error while attempting to extract contig sequence. Node {v} does not exist in haplotig graph.'.format(v=v)
            raise Exception(msg)
    edge_set = set(haplotig_graph.edges())
    for v, w in zip(node_path[:-1], node_path[1:]):
        if (v, w) not in edge_set:
            msg = 'Error while attempting to extract contig sequence. Edge ({v}, {w}) does not exist in haplotig graph.'.format(v=v, w=w)
            raise Exception(msg)

    # Create the contig sequence.
    # The sequence is composed of clipped haplotigs, so plain concatenation is valid.
    new_ctg_seq = ''
    for v in node_path:
        node = haplotig_graph.node[v]
        if node['label'] == 'source' or node['label'] == 'sink':
            continue
        new_ctg_seq += node['htig']['seq']

    # Construct the edges for the contig.
    new_ctg_path = []
    new_ctg_edges = []
    node_start_coords = {}
    node_end_coords = {}
    path_seq_len = 0

    for v in node_path:
        node = haplotig_graph.node[v]

        if node['label'] == 'source':
            node_start_coords[v] = 0
            node_end_coords[v] = 0
            continue
        elif node['label'] == 'sink':
            node_start_coords[v] = len(new_ctg_seq)
            node_end_coords[v] = len(new_ctg_seq)
            continue

        # Calculate the position of the node in the new contig.
        node_seq_len = len(node['htig']['seq'])
        node_start_coords[v] = path_seq_len
        node_end_coords[v] = path_seq_len + node_seq_len
        path_seq_len += node_seq_len

        # Get the phase info and the tiling path.
        haplotig = node['htig']
        phase = haplotig['phase']
        tp = haplotig['path']
        cross_phase = 'N'
        source_graph = 'H'

        for edge in tp:
            # Example tiling path line from 2-asm-falcon:
            #   ["000000F", "000000040:E", "000000217:E", "000000217", "14608", "33796", "14608", "99.86"]
            edge_ctg_id, edge_v, edge_w, edge_wid, edge_b, edge_e, edge_score, edge_idt = edge

            # Compose the path line.
            new_path_line = [new_ctg_id, edge_v, edge_w, edge_wid, edge_b, edge_e, edge_score, edge_idt, phase[1], phase[2]]
            new_ctg_path.append(' '.join([str(val) for val in new_path_line]))

            # Compose an edge line. All pread nodes within the same haplotig are phased in the same way.
            new_edge_line = [new_ctg_id, edge_v, edge_w, cross_phase, source_graph, phase[1], phase[2], phase[1], phase[2]]
            new_ctg_edges.append(' '.join([str(val) for val in new_edge_line]))

    return new_ctg_seq, new_ctg_edges, node_start_coords, node_end_coords

def find_haplotig_placement(haplotig_graph,
                            p_ctg_id, p_ctg_seq_len, p_nodes,
                            p_node_start_coords, p_node_end_coords,
                            h_ctg_id, h_ctg_seq_len, htig_node_path):

    """
    Inputs:
        haplotig_graph      - ANetworkX graph object. It doesn't have to be the subgraph
                            for p_ctg_id, it can also be the general graph before subgraph extraction.
        p_ctg_id            - Name (e.g. header) of the currently analyzed contig.
        p_ctg_seq_len       - Length of the primary contig.
        p_nodes             - A list() or a set() of node names that comprise the primary contig called p_ctg_id.
        p_node_start_coords - A dict where key is a node from p_node_set, and value the start
                            coordinate on p_ctg. Every node in p_node_set shold have a key here.
        p_node_end_coords   - A dict where key is a node from p_node_set, and value the end coordinate
                            on p_ctg. Every node in p_node_set shold have a key here.
        h_ctg_id            - The name of the haplotig for which the placement is generated.
        h_ctg_seq_len       - Length of the sequence of the haplotig.
        htig_node_path      - Path of nodes (list of node IDs) comprising the haplotig.
    Returns:
        A PAF-formatted tuple for the mapping of the haplotig to the primary contig.
    """

    p_node_set = set(p_nodes)

    hg_node_set = set(haplotig_graph.nodes())

    for v in p_nodes:
        if v not in hg_node_set:
            msg = 'Not all primary contig nodes can be found in the haplotig graph. p_node_set = "{}", hg_node_set = "{}", missing node = "{}"'.format(str(p_node_set), str(hg_node_set), v)
            raise Exception(msg)

    # Validate the inputs.
    if p_node_set != set(p_node_start_coords.keys()):
        msg = 'The p_node_set and p_node_start_coords contain different keys. p_node_set = "{}", p_node_start_coords = "{}"'.format(str(p_node_set), str(p_node_start_coords))
        raise Exception(msg)

    if p_node_set != set(p_node_end_coords.keys()):
        msg = 'The p_node_set and p_node_end_coords contain different keys. p_node_set = "{}", p_node_end_coords = "{}"'.format(str(p_node_set), str(p_node_end_coords))
        raise Exception(msg)

    if not htig_node_path:
        return None

    # Find the predecessor in the primary contig, if it exists.
    htig_s_node = htig_node_path[0]
    ss_set = set([v for v, w in haplotig_graph.in_edges(htig_s_node)]) & p_node_set
    ss = sorted(list(ss_set), key = lambda x: p_node_end_coords[x])
    # If no predecessors, it's dangling and starts at 0.
    unzipped_start = 0 if len(ss) == 0 else p_node_end_coords[ss[0]]

    # Find the successor in the primary contig, if it exists.
    htig_t_node = htig_node_path[-1]
    tt_set = set([w for v, w in haplotig_graph.out_edges(htig_t_node)]) & p_node_set
    tt = sorted(list(tt_set), key = lambda x: p_node_start_coords[x])
    # If no successors, it's dangling and ends at contig end.
    unzipped_end = p_ctg_seq_len if len(tt) == 0 else p_node_start_coords[tt[-1]]

    # Calculations for the shorthand.
    unzipped_span = unzipped_end - unzipped_start

    # Final PAF formulation of the placement.
    h_ctg_placement = (h_ctg_id, h_ctg_seq_len, 0, h_ctg_seq_len,
                        '+', p_ctg_id, p_ctg_seq_len,
                        unzipped_start, unzipped_end,
                        unzipped_span, unzipped_span, 60)

    return h_ctg_placement

def extract_unzipped_ctgs(ctg_id, haplotig_graph, allow_multiple_primaries, fp_proto_log):
    """
    Finds an arbitrary (shortest) walk down the haplotig DAG, and denotes it
    as "p_ctg.fa".
    "Shortest" is defined in terms of weighted shortest path, where weight is given
    by the notion of "cross-phase" edges.
    It then removes all p_ctg edges from the graph, and all cross-edges,
    and outputs all other weakly connected components as haplotigs in "h_ctg.fa"
    """

    fp_proto_log('Beginning to extract all p_ctg and h_ctg.')

    all_p_seqs = {}
    all_p_edges = {}
    all_h_seqs = {}
    all_h_edges = {}
    all_h_paf = {}

    num_prim_ctg = 0

    for sub_hg_id, sub_hg in enumerate(nx.weakly_connected_component_subgraphs(haplotig_graph)):
        if (not allow_multiple_primaries) and sub_hg_id > 0:
            msg = 'Skipping additional subgraphs of the primary contig: {ctg_id}. The graph has multiple primary components.'.format(ctg_id=ctg_id)
            raise Exception(msg)

        best_path = extract_unphased_haplotig_paths(sub_hg)
        if best_path == None:
            continue

        total_weight, p_node_path, s_node, t_node = best_path
        if len(p_node_path) == 0:
            continue
        p_node_set = set(p_node_path)

        ##################################
        ### Extract the primary contig ###
        ##################################
        # Form the primary contig.
        num_prim_ctg += 1
        p_ctg_id = ctg_id if not allow_multiple_primaries else '%sp%02d' % (ctg_id, num_prim_ctg)
        fp_proto_log('Extracting primary contig: p_ctg_id = {}'.format(p_ctg_id))

        p_ctg_seq, p_ctg_edges, p_node_start_coords, p_node_end_coords = construct_ctg_seq(sub_hg, p_ctg_id, p_node_path)
        all_p_seqs[p_ctg_id] = p_ctg_seq
        all_p_edges[p_ctg_id] = p_ctg_edges

        ##################################
        ### Extract the haplotigs      ###
        ##################################
        fp_proto_log('Extracting the associate haplotigs for p_ctg_id = {}'.format(p_ctg_id))

        # Make a copy where we'll delete stuff.
        sub_hg2 = sub_hg.copy()
        edges_to_remove = set()

        # Mark the primary path for deletion.
        for v, w in zip(p_node_path[:-1], p_node_path[1:]):
            edges_to_remove.add((v, w))

        # Mark any ambiguous edges for deletion.
        for v in sub_hg2.nodes():
            in_edges = set(sub_hg2.in_edges(v))
            if len(in_edges) > 1:
                edges_to_remove.update(in_edges)
            out_edges = set(sub_hg2.out_edges(v))
            if len(out_edges) > 1:
                edges_to_remove.update(out_edges)

        # Actually delete the edges and nodes.
        for v, w in edges_to_remove:
            try:
                sub_hg2.remove_edge(v, w)
            except Exception:
                pass
        for v in p_node_path:
            try:
                sub_hg2.remove_node(v)
            except Exception:
                pass

        # Loop through all associate haplotigs.
        num_hctg = 0
        for vals in extract_weakly_unphased_haplotig_paths(sub_hg2):
            htig_total_weight, htig_node_path, htig_s_node, htig_t_node = vals

            # This should never happen, but let's be sure.
            if not htig_node_path:
                continue

            # Form the haplotig.
            num_hctg += 1
            h_ctg_id = '%s_%03d' % (p_ctg_id, num_hctg)
            h_ctg_seq, h_ctg_edges, _, _ = construct_ctg_seq(sub_hg, h_ctg_id, htig_node_path)

            # Determine placement.
            h_ctg_placement = find_haplotig_placement(haplotig_graph,
                                                      p_ctg_id, len(p_ctg_seq), p_node_set,
                                                      p_node_start_coords, p_node_end_coords,
                                                      h_ctg_id, len(h_ctg_seq), htig_node_path)

            # Store the results.
            all_h_seqs[h_ctg_id] = h_ctg_seq
            all_h_edges[h_ctg_id] = h_ctg_edges
            all_h_paf[h_ctg_id] = h_ctg_placement

    return all_p_seqs, all_p_edges, all_h_seqs, all_h_edges, all_h_paf

def write_unzipped(out_dir, ctg_id, p_ctg_seqs, p_ctg_edges, h_ctg_seqs, h_ctg_edges, h_ctg_paf, fp_proto_log):
    blacklist = set()

    # Any sequence of length 0 should not be output.
    for seq_name, seq in p_ctg_seqs.iteritems():
        if len(seq) == 0:
            fp_proto_log('[write_unzipped] Sequence "{seq_name}" has length zero. Adding to blacklist.'.format(seq_name=seq_name))
            blacklist.add(seq_name)
    # If the length of the edges is 0, something is awry.
    # Rather skip and report, than break execution on the entire assembly.
    for seq_name, edges in p_ctg_edges.iteritems():
        if len(edges) == 0:
            fp_proto_log('[write_unzipped] Sequence "{seq_name}" no edges. Adding to blacklist.'.format(seq_name=seq_name))
            blacklist.add(seq_name)
    # Any haplotig of length 0 should not be written.
    for seq_name, seq in h_ctg_seqs.iteritems():
        if len(seq) == 0:
            fp_proto_log('[write_unzipped] Sequence "{seq_name}" has length zero. Adding to blacklist.'.format(seq_name=seq_name))
            blacklist.add(seq_name)
        # Check if the corresponding primary contig is blacklisted.
        if seq_name.split('_')[0] in blacklist:
            fp_proto_log('[write_unzipped] The primary contig of "{seq_name}" is blacklisted. Adding to blacklist.'.format(seq_name=seq_name))
            blacklist.add(seq_name)
    for seq_name, edges in h_ctg_edges.iteritems():
        if len(edges) == 0:
            fp_proto_log('[write_unzipped] Sequence "{seq_name}" no edges. Adding to blacklist.'.format(seq_name=seq_name))
            blacklist.add(seq_name)
    for seq_name, paf in h_ctg_paf.iteritems():
        if len(paf) == 0:
            fp_proto_log('[write_unzipped] Sequence "{seq_name}" has no placement. This is not critical and will not be added to blacklist, but worth noting.'.format(seq_name=seq_name))
            # blacklist.add(seq_name)

    with open(os.path.join(out_dir, "p_ctg.%s.fa" % ctg_id), "w") as fp_out:
        for seq_name in sorted(p_ctg_seqs.keys()):
            if seq_name in blacklist:
                continue
            seq = p_ctg_seqs[seq_name]
            fp_out.write('>%s\n' % (seq_name))
            fp_out.write(seq)
            fp_out.write('\n')
    with open(os.path.join(out_dir, "p_ctg_edges.%s" % ctg_id), "w") as fp_out:
        for seq_name in sorted(p_ctg_edges.keys()):
            if seq_name in blacklist:
                continue
            edges = p_ctg_edges[seq_name]
            fp_out.write('\n'.join(edges))
            fp_out.write('\n')
    with open(os.path.join(out_dir, "h_ctg.%s.fa" % ctg_id), "w") as fp_out:
        for seq_name in sorted(h_ctg_seqs.keys()):
            if seq_name in blacklist:
                continue
            seq = h_ctg_seqs[seq_name]
            fp_out.write('>%s\n' % (seq_name))
            fp_out.write(seq)
            fp_out.write('\n')
    with open(os.path.join(out_dir, "h_ctg_edges.%s" % ctg_id), "w") as fp_out:
        for seq_name in sorted(h_ctg_edges.keys()):
            if seq_name in blacklist:
                continue
            edges = h_ctg_edges[seq_name]
            fp_out.write('\n'.join(edges))
            fp_out.write('\n')
    with open(os.path.join(out_dir, "h_ctg.%s.paf" % ctg_id), "w") as fp_out:
        for seq_name in sorted(h_ctg_paf.keys()):
            if seq_name in blacklist:
                continue
            paf = h_ctg_paf[seq_name]
            if len(paf) == 0:
                continue
            fp_out.write('\t'.join([str(val) for val in paf]))
            fp_out.write('\n')

def get_rid2proto_dir(gath_fn):
    LOG.info('Reading {!r} to create rid2proto_dir map.'.format(gath_fn))
    gath = io.deserialize(gath_fn)
    result = dict()
    gath_dn = os.path.dirname(gath_fn)
    def fixpath(path):
        if os.path.isabs(path):
            return path
        return os.path.normpath(os.path.join(gath_dn, path))
    for rec in gath:
        rid_to_phase_fn = fixpath(rec['rid_to_phase_out'])
        # That was the output-name for TASK_PHASING_RUN_SCRIPT,
        # specified in phasing_split.py.
        # We now assume that its ctg_id is two dirs up, and
        # the proto/ directory is next to the rid_to_phase file.
        phasing_run_dir = os.path.dirname(rid_to_phase_fn)
        proto_dir = os.path.join(phasing_run_dir, 'proto')
        assert os.path.isdir(proto_dir), 'Missing proto_dir {!r}'.format(proto_dir)
        ctg_id = os.path.basename(os.path.dirname(phasing_run_dir))
        # I guess ctg_id == rid?
        result[ctg_id] = proto_dir
    LOG.info('From {!r}, rid2proto_dir len(dict)={!r}'.format(
        gath_fn, len(result)))
    if result:
        arbitrary_key = next(iter(result))
        LOG.info(' For example: dict[{!r}]->\n  {!r}'.format(
            arbitrary_key, result[arbitrary_key]))
    return result

def define_globals(args):
    # make life easier for now. will refactor it out if possible
    global all_rid_to_phase
    global all_flat_rid_to_phase
    global all_haplotigs_for_ctg
    global p_asm_G
    global h_asm_G
    global p_ctg_seqs
    global sg_edges
    global p_ctg_tiling_paths

    fc_asm_path = args.fc_asm_path
    fc_hasm_path = args.fc_hasm_path
    base_dir = args.base_dir
    fasta_fn = args.fasta

    #hasm_falcon_path = os.path.join(fc_hasm_path, 'asm-falcon')
    hasm_falcon_path = fc_hasm_path # They run in same dir, for now.

    LOG.info('Loading p assembly graph.')
    p_asm_G = AsmGraph(os.path.join(fc_asm_path, "sg_edges_list"),
                       os.path.join(fc_asm_path, "utg_data"),
                       os.path.join(fc_asm_path, "ctg_paths"))
    LOG.info('Loading h assembly graph.')
    h_asm_G = AsmGraph(os.path.join(fc_hasm_path, "sg_edges_list"),
                       os.path.join(fc_hasm_path, "utg_data"),
                       os.path.join(fc_hasm_path, "ctg_paths"))
    assert p_asm_G, 'Empty AsmGraph. Maybe empty inputs?\n{!r}\n{!r}\n{!r}'.format(
        os.path.join(fc_asm_path, "sg_edges_list"),
        os.path.join(fc_asm_path, "utg_data"),
        os.path.join(fc_asm_path, "ctg_paths"),
    )

    LOG.info('Loading phasing info and making the read ID sets.')

    all_rid_to_phase = {}
    all_flat_rid_to_phase = {}
    with io.open_progress(args.rid_phase_map) as f:
        for row in f:
            row = row.strip().split()
            all_rid_to_phase.setdefault(row[1], {})
            all_rid_to_phase[row[1]][row[0]] = (int(row[2]), int(row[3]))
            all_flat_rid_to_phase[row[0]] = (row[1], int(row[2]), int(row[3]))

    # Load the primary contig sequences.
    LOG.info('Loading the 2-asm-falcon primary contigs.')
    p_ctg_seqs = load_all_seq(os.path.join(fc_asm_path, "p_ctg.fa"))
    LOG.info('Done loading 2-asm-falcon primary contigs.')

    LOG.info('Loading tiling paths.')
    # Hash the lengths of the primary contig sequences.
    # Needed to correctly assign node coords when loading tiling paths
    p_ctg_seq_lens = {}
    for p_ctg_id, ctg_seq in p_ctg_seqs.iteritems():
        p_ctg_seq_lens[p_ctg_id] = len(ctg_seq)
    # Load the tiling path of the primary contig, and assign coordinants to nodes.
    p_ctg_tiling_paths = tiling_path.load_tiling_paths(os.path.join(fc_asm_path, "p_ctg_tiling_path"), contig_lens=p_ctg_seq_lens, whitelist_seqs=None)
    LOG.info('Done loading tiling paths.')

    # Load the haplotig sequences (assembled with falcon_kit.mains.graph_to_contig).
    LOG.info('Loading the 1-hasm haplotigs.')
    all_haplotigs_for_ctg, htig_name_to_original_pctg = load_haplotigs(hasm_falcon_path, all_flat_rid_to_phase)
    LOG.info('Done loading haplotigs.')

    # Load all sg_edges_list so that haplotig paths can be reversed if needed.
    LOG.info('Loading sg_edges_list.')
    sg_edges = {}
    sg_edges_list_fn = os.path.join(fc_hasm_path, 'sg_edges_list')
    with io.open_progress(sg_edges_list_fn, 'r') as fp:
        for line in fp:
            sl = line.strip().split()
            sg_edges[(sl[0], sl[1])] = sl
    LOG.info('Done loading sg_edges_list.')

def cmd_apply(args):
    units_of_work_fn = args.units_of_work_fn
    results_fn = args.results_fn

    LOG.info('Loading units-of-work from {!r}'.format(units_of_work_fn))
    units_of_work = io.deserialize(units_of_work_fn)
    if not units_of_work:
        LOG.warning('No units of work in {!r}'.format(units_of_work_fn))
        io.serialize(results_fn, list())
        return
    units_of_work_dn = os.path.dirname(units_of_work_fn)
    def fixpath(path):
        if os.path.isabs(path):
            return path
        return os.path.normpath(os.path.join(units_of_work_dn, path))
    uow = units_of_work[0] # Get the common args from the first UOW.
    sub_args = lambda: None # for attribute storage
    sub_args.fc_asm_path = fixpath(uow['input']['fc_asm_path'])
    sub_args.fc_hasm_path = fixpath(uow['input']['fc_hasm_path'])
    sub_args.base_dir = fixpath(uow['input']['base_dir'])
    sub_args.fasta = fixpath(uow['input']['fasta_fn'])
    sub_args.rid_phase_map = fixpath(uow['input']['rid_phase_map'])

    define_globals(sub_args)

    #LOG.info('Creating the exe list for: {}'.format(str(ctg_id_list)))

    exe_list = list()
    for i, uow in enumerate(units_of_work):
        ctg_id = uow['params']['ctg_id']
        proto_dir = uow['params']['proto_dir']
        out_dir = os.path.join('.', 'uow-{}'.format(ctg_id))
        base_dir = sub_args.base_dir

        exe_list.append((ctg_id, proto_dir, out_dir, base_dir, False))

    LOG.info('Running {} units of work.'.format(len(exe_list)))

    #exec_pool = Pool(args.nproc)  # TODO, make this configurable
    #exec_pool.map(run_generate_haplotigs_for_ctg, exe_list)
    ##map( generate_haplotigs_for_ctg, exe_list)

    results = list()
    for i, exe in enumerate(exe_list):
        out_dir = exe[2] # See prior for-loop. # TODO: with cd(out_dir)

        LOG.info('UOW #{} of {} ...'.format(i, len(exe_list)))
        result = run_generate_haplotigs_for_ctg(exe)

        # We could specify this to 'run_gen', but for now it is an implicit output.
        output_h_ctg_fn = os.path.join(out_dir, 'h_ctg.{}.fa'.format(ctg_id))
        assert os.path.exists(output_h_ctg_fn), 'Missing h_ctg fasta: {!r}'.format(
                os.path.abspath(output_h_ctg_fn))

        result = dict(
                h_ctg=output_h_ctg_fn, # We record the result even if zero-size.
                ctg_id=ctg_id,
        )
        results.append(result)

    io.serialize(results_fn, results)


######
TASK_APPLY_UNITS_OF_WORK = """\
python -m falcon_unzip.mains.graphs_to_h_tigs_2 apply --units-of-work-fn={input.units_of_work} --results-fn={output.results}

#--bash-template-fn= # not needed
"""

def cmd_split(args):
    gath_fn = args.gathered_rid_to_phase
    split_fn = args.split_fn
    bash_template_fn = args.bash_template_fn
    with open(bash_template_fn, 'w') as stream:
        stream.write('echo hi')

    rid2proto_dir = get_rid2proto_dir(gath_fn)

    define_globals(args)

    ctg_id_list = p_asm_G.ctg_data.keys()

    LOG.info('Creating units-of-work for ctg_id_list (though many will be skipped): {}'.format(ctg_id_list))

    uows = []
    for ctg_id in ctg_id_list:
        if ctg_id[-1] != "F":
            continue
        if ctg_id not in all_rid_to_phase:
            continue
        proto_dir = rid2proto_dir[ctg_id]
        uow = dict(
                input=dict(
                    # common inputs:
                    fc_asm_path = args.fc_asm_path,
                    fc_hasm_path = args.fc_hasm_path,
                    base_dir = args.base_dir,
                    fasta_fn = args.fasta,
                    rid_phase_map = args.rid_phase_map,
                ),
                params=dict(
                    ctg_id=ctg_id,
                    proto_dir=proto_dir,
                ),
                wildcards=dict(
                    chunk_id='chunk_{}'.format(ctg_id), # for later substitution
                ),
        )
        uows.append(uow)
    io.serialize(split_fn, uows)

def generate_h_ctg_ids(run_dir, ctg_id):
    cmd = 'grep ">" ./h_ctg.{ctg_id}.fa | sed "s/^>//" >| ./h_ctg_ids.{ctg_id}'.format(
            **locals())
    with cd(run_dir):
        execute.execute_command(cmd, LOG)
def combine(ostream, fns):
    for fn in fns:
        with open(fn) as istream: # IOError would report fn
            ostream.write(istream.read())
def cmd_combine(args):
    results_fn = args.results_fn
    done_fn = args.done_fn
    results = io.deserialize(results_fn)
    combined = collections.defaultdict(list)
    results_dn = os.path.dirname(results_fn)
    def fixpath(path):
        if os.path.isabs(path):
            return os.path.normpath(path)
        return os.path.normpath(os.path.join(results_dn, path))
    for result in results:
        h_ctg_fn = fixpath(result['h_ctg'])
        run_dir = os.path.dirname(h_ctg_fn)
        ctg_id = result['ctg_id']
        combined['h_ctg'].append(h_ctg_fn)
        combined['p_ctg'].append(os.path.join(run_dir, 'p_ctg.{}.fa'.format(ctg_id)))
        combined['h_ctg_ids'].append(os.path.join(run_dir, 'h_ctg_ids.{}'.format(ctg_id)))
        combined['h_ctg_edges'].append(os.path.join(run_dir, 'h_ctg_edges.{}'.format(ctg_id)))
        combined['p_ctg_edges'].append(os.path.join(run_dir, 'p_ctg_edges.{}'.format(ctg_id)))
        combined['h_ctg_placement'].append(os.path.join(run_dir, 'h_ctg.{}.paf'.format(ctg_id)))
        generate_h_ctg_ids(run_dir, ctg_id)
    with open('all_h_ctg.fa', 'w') as stream:
        combine(stream, combined['h_ctg'])
    with open('all_p_ctg.fa', 'w') as stream:
        combine(stream, combined['p_ctg'])
    with open('all_h_ctg_ids', 'w') as stream:
        combine(stream, combined['h_ctg_ids'])
    with open('all_h_ctg_edges', 'w') as stream:
        combine(stream, combined['h_ctg_edges'])
    with open('all_p_ctg_edges', 'w') as stream:
        combine(stream, combined['p_ctg_edges'])
    with open('all_h_ctg.paf', 'w') as stream:
        combine(stream, combined['h_ctg_placement'])
    io.touch(done_fn)
'''
find ./0-phasing -name "phased_reads" | sort | xargs cat >| all_phased_reads
find ./2-htigs -name "h_ctg_ids.*" | sort | xargs cat >| all_h_ctg_ids
find ./2-htigs -name "p_ctg_edges.*" | sort | xargs cat >| all_p_ctg_edges
find ./2-htigs -name "h_ctg_edges.*" | sort | xargs cat >| all_h_ctg_edges
find ./2-htigs -name "p_ctg.*.fa" | sort | xargs cat >| all_p_ctg.fa
find ./2-htigs -name "h_ctg.*.fa" | sort | xargs cat >| all_h_ctg.fa
'''
def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='layout haplotigs from primary assembly graph and phased aseembly graph')
    subparsers = parser.add_subparsers(help='sub-command help')
    description_combine = """
    The run-dirs will have various files. Some of them will be named explicitly in results.json. Others are implicit outputs (to minimize explicits). E.g. p/h ctgs, ids, and edges. We will concatenate these files into single files, excluding any units of work which produced empty fasta. For now, we expect:
  all_h_ctg_edges
  all_h_ctg.fa
  all_h_ctg_ids
  all_p_ctg_edges
  all_p_ctg.fa
  all_phased_reads # ??? This would come from 0-phasing/
    """
    help_combine = 'Combine the results of each "graph_to_h_tigs" application into something useful to the next task.'
    help_split = 'Split input into atomic units-of-work.'
    help_apply = 'Apply "graph_to_h_tigs" to a subset of one or more units-of-work, serially.'
    parser_split = subparsers.add_parser('split',
            description=help_split,
            help=help_split)
    parser_apply = subparsers.add_parser('apply',
            description=help_apply,
            help=help_apply)
    parser_combine = subparsers.add_parser('combine',
            description=help_combine,
            epilog=description_combine,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            help=help_combine)
    parser_combine.set_defaults(func=cmd_combine)

    parser_split.add_argument(
        '--split-fn', required=True,
        help='Output: JSON list of all units of work.')
    parser_split.add_argument(
        '--bash-template-fn', required=True,
        help='Output: bash script to run a unit-of-work, given a record from JSON split-file.')
    parser_split.add_argument(
        '--gathered-rid-to-phase', required=True,
        help='Input. (Possibly 3-unzip/1-hasm/gathered-rid-to-phase/gathered.json) JSON list of dict of "rid_to_phase_out"->filename. This is used to find the proto/ directory for each unit of phasing work.')
    parser_split.add_argument(
        '--fc-asm-path', type=str, required=True,
        help='Common input: path to the primary Falcon assembly output directory')
    parser_split.add_argument(
        '--fc-hasm-path', type=str, required=True,
        help='Common input: path to the phased Falcon assembly output directory')
    parser_split.add_argument(
        '--base-dir', type=str, default='.',
        help='Common input: the output base_dir, default to current working directory')
    parser_split.add_argument(
        '--rid-phase-map', type=str, required=True,
        help="Common input: path to the file that encode the relationship of the read id to phase blocks")
    parser_split.add_argument(
        '--fasta', type=str, required=True,
        help="Common input: sequence file of the p-reads")
    parser_split.set_defaults(func=cmd_split)

    parser_apply.add_argument(
        '--units-of-work-fn', required=True,
        help='Input. JSON list of units of work. This can be on, several, or all of the list from subcommand "split".')
    parser_apply.add_argument(
        '--results-fn', required=True,
        help='Output. JSON list of results, one record per unit-of-work.')
    parser_apply.set_defaults(func=cmd_apply)

    parser_combine.add_argument(
        '--results-fn', required=True,
        help='Input: JSON list of results, one record per unit-of-work. Probably, a "gatherer" has already gathered the "results.json" from each application into one large file, this one.')
    parser_combine.add_argument(
        '--done-fn',
        help='Output: Sentinel.')
    #parser_combine.add_argument(
    #    '--combined-fn',
    #    help='Output: JSON. I think it will be a list of dicts, describing the p/h ctgs, ids, and edges. This might be deduced from the run-dir, to minimize the number of outputs we need to specify explicitly. The list will exclude units-of-work which produced empty h_ctg files. Also, the paths will probably be absolutized here, in case the "gatherer" did not do that already. (There were relative so "apply" could occur in a tmpdir.)')
    parser_combine.set_defaults(func=cmd_combine)

    #parser_run = subparsers.add_parser('run',
    #        help='Actually run "graph_to_h_tigs". This is called for each unit-of-work, in a directory selected by subcommand "apply".',
    #        description='This a single-threaded-program, but it will call multi-threaded blasr.')

    #parser_run.add_argument(
    #    '--ctg-id', type=str, required=True,
    #    help='contig identifier in the bam file')
    #parser_run.add_argument(
    #    '--proto-dir', type=str, required=True,
    #    help='directory of phasing work, e.g. "3-unzip/0-phasing/000000F/uow-00/proto/"')
    #parser_run.add_argument(
    #    '--nproc', type=int, default=8, help="number of processes to use")

    args = parser.parse_args(argv[1:])

    return args


def main(argv=sys.argv):
    global LOG
    args = parse_args(argv)

    # Write all warnings to stderr, including from thread-loggers.
    # (However, the main-thread logger does not currently propagate recs here.)
    hdlr = logging.StreamHandler(sys.stderr)
    #hdlr.setLevel(logging.WARNING - 1) # We want records from execute_command() too.
    hdlr.setLevel(logging.INFO)
    #hdlr.setFormatter(logging.Formatter('[Proto %(asctime)s] %(levelname)s:%(name)s:t%(message)s', '%Y-%m-%d %H:%M:%S'))
    hdlr.setFormatter(logging.Formatter('[%(levelname)s %(asctime)s] %(message)s', '%Y-%m-%d %H:%M:%S'))
    LOG.addHandler(hdlr)
    LOG.setLevel(logging.NOTSET) # Important, as other NOTSET loggers inherit this level.

    # When we were multi-threads, we wanted separate logging per thread.
    #### In main thread, log to a special file.
    #### (Turn this off if too verbose.)
    #### This handler will not see the thread-logger
    #### log-records at all.
    ###LOG = logging.getLogger('mainthread')
    ####hdlr = logging.FileHandler('graphs_to_h_tigs_2.log', 'w')
    ###hdlr = logging.StreamHandler(sys.stderr)
    ###hdlr.setLevel(logging.INFO)
    ###hdlr.setFormatter(logging.Formatter('[Proto %(asctime)s] %(levelname)s:%(message)s', '%Y-%m-%d %H:%M:%S'))
    ###LOG.addHandler(hdlr)
    ###LOG.propagate = False

    #import pdb; pdb.set_trace()
    logging.addLevelName(logging.WARNING-1, 'EXECUTE')

    args.func(args)


if __name__ == '__main__':  # pragma: no cover
    main()
