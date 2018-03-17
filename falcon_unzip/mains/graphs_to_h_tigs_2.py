from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.FastaReader import FastaReader
import os
import networkx as nx
from multiprocessing import Pool
import re
import falcon_unzip.proto.execute as execute
from falcon_unzip.proto.haplotig import Haplotig
import falcon_unzip.proto.tiling_path as tiling_path
import falcon_unzip.proto.sam2m4 as sam2m4
import copy
import pysam
import falcon_unzip.proto.cigartools as cigartools
import json
from graphs_to_h_tigs_2_utils import *
from graphs_to_h_tigs_2_utils import proto_log
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
global seqs
global seq_lens
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
    try:
        ctg_id, out_dir = input_
        sys.stderr.write('[Proto] Entered run_generate_haplotigs_for_ctg for ctg_id = %s, out_dir = %s\n' % (ctg_id, out_dir))
        mkdir(out_dir)

        # Prepare the log output stream for this particular contig.
        fp_proto_log = logging.getLogger(ctg_id)
        fp_proto_log.setLevel(logging.INFO)
        log_fn = os.path.join(out_dir, 'prototype.log')
        LOG.warning('New logging FileHandler: {!r}'.format(os.path.abspath(log_fn)))
        hdlr = logging.FileHandler(log_fn) #, level=logging.INFO)
        fp_proto_log.addHandler(hdlr)

        base_dir = os.path.join(out_dir, '..', '..') # wrong; fix soon
        return generate_haplotigs_for_ctg(ctg_id, out_dir, fp_proto_log)
    except Exception:
        LOG.exception('Failure in generate_haplotigs_for_ctg({!r})'.format(input_))
        raise

def proto_log(message, fp_proto_log):
    fp_proto_log.info(message)

def generate_haplotigs_for_ctg(ctg_id, out_dir, base_dir, fp_proto_log):
    global all_haplotigs_for_ctg
    global seqs
    global seq_lens
    global p_ctg_seqs
    global sg_edges
    global p_ctg_tiling_paths

    # min_linear_len = 8

    #########################################################
    # Load and prepare data.
    #########################################################
    proto_log('Fetching the p_ctg_seq.\n', fp_proto_log)
    p_ctg_seq = p_ctg_seqs[ctg_id]

    proto_log('Fetching the p_ctg_tiling_path.\n', fp_proto_log)
    p_ctg_tiling_path = p_ctg_tiling_paths[ctg_id]

    # Path to the directory with the output of the main_augment_pb.py.
    proto_dir = os.path.join(base_dir, '0-phasing', ctg_id, 'proto')

    # Load the linear sequences for alignment.
    p_ctg_path = os.path.join(base_dir, 'reads', ctg_id, 'ref.fa')
    linear_p_ctg_path = os.path.join(base_dir, '0-phasing', ctg_id, 'proto', 'p_ctg_linear_%s.fasta' % (ctg_id))

    proto_log('Loading linear seqs from %s .\n' % (linear_p_ctg_path), fp_proto_log)
    linear_seqs = load_all_seq(linear_p_ctg_path)

    # Load the phase relation graph, and build a dict of the weakly connected components.
    phase_relation_graph = nx.read_gexf(os.path.join(proto_dir, "phase_relation_graph.gexf"))
    proto_log('Loading the phase relation graph from %s .\n' % (phase_relation_graph), fp_proto_log)
    phase_alias_map, alias_to_phase_list = create_phase_alias_map(phase_relation_graph)

    # Load the regions as built previous in the 0-phasing step.
    # The structure of each region in the list is:
    #   type_, edge_start_id, edge_end_id, pos_start, pos_end, htigs = all_regions[i]
    all_regions_path = os.path.join(proto_dir, 'regions.json')
    proto_log('Loading all regions from %s .\n' % (all_regions_path), fp_proto_log)
    all_regions = json.load(open(all_regions_path))

    # Create a list of bubble regions, for easier merging with the diploid groups.
    proto_log('Making bubble region list.\n', fp_proto_log)
    bubble_region_list = []
    for region in all_regions:
        region_type, first_edge, last_edge, pos_start, pos_end, region_htigs = region
        if region_type != 'linear':
            bubble_region_list.append(region)

    # The json module converts a tuple to list, so revert back here.
    proto_log('Retoupling.\n', fp_proto_log)
    for region in all_regions:
        region_type, first_edge, last_edge, pos_start, pos_end, region_htigs = region
        for htig_name, htig in region_htigs.iteritems():
            htig['phase'] = tuple(htig['phase'])

    # Assign sequences to all regions.
    proto_log('Assigning sequences to all regions.\n', fp_proto_log)
    for region in all_regions:
        region_type, first_edge, last_edge, pos_start, pos_end, region_htigs = region
        if region_type == 'linear':
            # For linear regions, we can take a shortcut by directly
            # extracting them from the
            region_name = region_htigs.keys()[0]
            region_htigs[region_name]['seq'] = p_ctg_seq[pos_start:pos_end]
        else:
            # For any bubble component, generate the sequence here to
            # ensure that the 'proper' approach is not applied.
            # This is important to make good tiling of the sequences.
            for htig_name, htig in region_htigs.iteritems():
                htig['seq'] = path_to_seq(seqs, htig['path'], False)

    # Get all the haplotig objects corresponding to this contig only.
    proto_log('Getting snp_haplotigs.\n', fp_proto_log)
    snp_haplotigs = all_haplotigs_for_ctg.get(ctg_id, {})

    #########################################################

    #########################################################
    # Debug verbose.
    #########################################################
    proto_log('phase_alias_map = %s\n' % (str(phase_alias_map)), fp_proto_log)

    proto_log('Verbose all regions loaded from regions.json:\n', fp_proto_log)
    for region_id in xrange(len(all_regions)):
        region = all_regions[region_id]
        region_type, first_edge, last_edge, pos_start, pos_end, region_htigs = region
        proto_log('[region_id = %d] type = %s, first_edge = %s, last_edge = %s, pos_start = %s, pos_end = %s\n' % (region_id, region_type, str(first_edge), str(last_edge), str(pos_start), str(pos_end)), fp_proto_log)
    proto_log('\n', fp_proto_log)
    #########################################################

    #########################################################
    # Write the haplotigs to disk for alignment.
    #########################################################
    proto_log('[Proto] Writing haplotigs to disk.\n', fp_proto_log)
    aln_snp_hasm_ctg_path = os.path.join(out_dir, "aln_snp_hasm_ctg.fasta")
    write_haplotigs(snp_haplotigs, aln_snp_hasm_ctg_path, fp_proto_log, hack_qnames_for_blasr=True)

    #########################################################
    # Align haplotigs.
    #########################################################
    num_threads = 16
    mapping_out_prefix = os.path.join(out_dir, 'aln_snp_hasm_ctg')
    mapping_ref = linear_p_ctg_path
    mapping_ref = p_ctg_path
    sam_path = mapping_out_prefix + '.sam'
    if len(snp_haplotigs.keys()) > 0:
        # BLASR crashes on empty files, so address that.
        blasr_params = '--minMatch 15 --maxMatch 25 --advanceHalf --advanceExactMatches 10 --bestn 1 --nproc %d --noSplitSubreads' % (num_threads)
        execute.execute_command('blasr %s %s %s --sam --out %s.sam' % (blasr_params, aln_snp_hasm_ctg_path, mapping_ref, mapping_out_prefix), sys.stderr, dry_run = False)
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
    proto_log('[Proto] Creating a haplotig graph.\n', fp_proto_log)
    haplotig_graph = regions_to_haplotig_graph(ctg_id, final_all_regions, phase_alias_map, fp_proto_log)
    haplotig_graph = update_haplotig_graph(haplotig_graph, phase_alias_map)

    # Hash all haplotigs for lookup.
    proto_log('  - Hashing haplotigs.\n', fp_proto_log)
    all_haplotig_dict = {}
    for region in final_all_regions:
        region_type, edge_start, edge_end, region_pos_start, region_pos_end, region_haplotigs = region
        all_haplotig_dict.update({htig_key: htig for htig_key, htig in region_haplotigs.iteritems()})

    # Write the nx graph to disk for debugging.
    proto_log('  - Writing the haplotig graph in the gexf format.\n', fp_proto_log)
    nx.write_gexf(haplotig_graph, os.path.join(out_dir, "haplotig_graph.gexf"))

    # Write the haplotig graph to disk.
    proto_log('[Proto] Writing the haplotig_graph.gfa.\n', fp_proto_log)
    with open(os.path.join(out_dir, "haplotig_graph.gfa"), 'w') as fp_out:
        nx_to_gfa(ctg_id, haplotig_graph, all_haplotig_dict, fp_out)

    # Extract the p_ctg and h_ctg sequences.
    extract_and_write_all_ctg(ctg_id, haplotig_graph, all_haplotig_dict, phase_alias_map, final_all_regions, out_dir, fp_proto_log)

    #########################################################
    # Debug verbose.
    #########################################################
    for line in sorted(clippoints.iteritems()):
        proto_log('clip_point: %s\n' % (str(line)), fp_proto_log)

    proto_log('[Proto] Fragmented haplotigs:\n', fp_proto_log)
    for hname, htig in fragmented_snp_haplotigs.iteritems():
        proto_log('  - name = %s, phase = %s, path = %s\n' % (htig.name, str(htig.phase), str(htig.path)), fp_proto_log)

    proto_log('Verbose the generated diploid regions:\n', fp_proto_log)
    for region_id in xrange(len(final_all_regions)):
        region = final_all_regions[region_id]
        region_type, first_edge, last_edge, pos_start, pos_end, region_htigs = region
        proto_log('[region_id = %d] type = %s, first_edge = %s, last_edge = %s, pos_start = %s, pos_end = %s\n' % (region_id, region_type, str(first_edge), str(last_edge), str(pos_start), str(pos_end)), fp_proto_log)
    proto_log('\n', fp_proto_log)
    #########################################################

    proto_log('[Proto] Finished.\n', fp_proto_log)

    fp_proto_log.close()

def load_haplotigs(hasm_falcon_path, all_flat_rid_to_phase, preads):
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

    sys.stderr.write('[Proto] Loading haplotigs.\n')

    # Load the primary contig tiling paths.
    p_tiling_paths = os.path.join(hasm_falcon_path, 'p_ctg_tiling_path')
    tiling_paths = tiling_path.load_tiling_paths(p_tiling_paths, None, None)

    # Process the tiling paths for each assembled haplotig.
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

        ###########################
        # Construct the haplotig. #
        ###########################
        # Info needed to construct a haplotig.
        # Seq is generated from path, to avoid proper/improper confusion.
        h_tig_name = '%s-HAP%s-%s.%s.%s' % (phase_ctg, hctg, phase_ctg, str(phase_block), str(phase_id))
        complete_phase = vphase
        path = hpath.dump_as_split_lines()
        seq = path_to_seq(preads, path, True)

        new_haplotig = Haplotig(name = h_tig_name, phase = complete_phase, seq = seq, path = path, edges = [])

        # The name relation dicts.
        htig_name_to_original_pctg[h_tig_name] = hctg
        # original_pctg_to_htig_name[hctg] = h_tig_name

        def verbose_haplotig(haplotig):
            return 'phase = %s, h_tig_name = %s, num_edges = %d, len(seq) = %d' % (str(haplotig.phase), haplotig.name, len(haplotig.path), len(haplotig.seq))

        sys.stderr.write('[%d] Loaded haplotig: %s\n' % (num_hctg, verbose_haplotig(new_haplotig)))

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
    proto_log('[Proto] Loading the alignments.\n', fp_proto_log)
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
    proto_log('[Proto] Reorienting haplotigs.\n', fp_proto_log)
    for qname, haplotig in snp_haplotigs.iteritems():
        aln = aln_dict[qname]
        qstrand, qstart, qend, qlen = aln[4:8]
        tstrand, tstart, tend, tlen = aln[8:12]
        if tstrand == 1:
            proto_log('  - qname = %s\n' % (qname), fp_proto_log)
            # proto_log('  - Before reversing:\n', fp_proto_log)
            # for line in haplotig.path:
            #     proto_log('  %s\n' % (str(line)), fp_proto_log)

            # Reverse the path and the seq.
            haplotig.path = reverse_sg_path(haplotig.path, sg_edges)
            haplotig.seq = revcmp_seq(haplotig.seq)
            proto_log('\n', fp_proto_log)
            # Reset the reverse flag.
            # Coordinates in aln are on the fwd strand anyway,
            # so they remain the same.
            aln[8] = 0

            # proto_log('  - After reversing:\n', fp_proto_log)
            # for line in haplotig.path:
            #     proto_log('  %s\n' % (str(line)), fp_proto_log)

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
        aln = aln_dict[qname]
        qstrand, qstart, qend, qlen = aln[4:8]
        tstrand, tstart, tend, tlen = aln[8:12]
        clippoints[tstart] = [qname + '_start']
        clippoints[tend] = [qname + '_end']

    return clippoints, bubble_tree

def fragment_haplotigs(in_haplotigs, aln_dict, clippoints, bubble_tree, fp_proto_log):
    filtered_haplotigs = {}

    for qname, haplotig in in_haplotigs.iteritems():
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
            proto_log('start = %s, end = %s\n' % (str(start), str(end)), fp_proto_log)
            tstart, tend = start[0], end[0]
            found_intervals = bubble_tree.search(tstart, tend)
            if found_intervals:
                continue
            regions_of_interest.append((start, end, q_name, q_len, t_name, t_len, haplotig.phase))

    proto_log('[Proto] pos_of_interest for q_name: %s\n' % (q_name), fp_proto_log)
    for line in regions_of_interest:
        proto_log('%s\n' % (str(line)), fp_proto_log)
    proto_log('\n', fp_proto_log)

    # Here we extract the subsequence and the subpath of the tiling path.
    ret_haplotigs = {}
    tp = tiling_path.convert_split_lines_to_tiling_path(haplotig.path, len(haplotig.seq))

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

        for phase_group, phase_id_dict in phase_group_dict.iteritems():
            # Skip any group which doesn't have exactly 2 phases in it.
            # E.g. one haplotig was larger than the other, due to read length.
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

        # Collect the haplotigs.
        region_htigs = {}
        for phase_group, phase_id_dict in phase_group_dict.iteritems():
            for phase_id, haplotigs in phase_id_dict.iteritems():
                for haplotig in haplotigs:
                    region_htigs[haplotig.name] = haplotig.__dict__

        # Create a diploid region.
        region_type = 'diploid'
        first_edge = None           # First edge on the p_ctg. Not available right now.
        last_edge = None            # Last edge on the p_ctg. Not available right now.
        pos_start, pos_end = region
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

    proto_log('[Proto] Entered function: "%s"\n' % (str(__name__)), fp_proto_log)
    if (pos_end - pos_start) <= 0:
        proto_log('[Proto] Exiting function: "%s"\n' % (str(__name__)), fp_proto_log)
        return None

    region_type = 'linear'
    complete_phase = (ctg_id, '-1', '0')

    htig_name = '%s_%s2_%d:%d_base' % (ctg_id, region_type, pos_start, pos_end)

    proto_log('[Proto] Info: pos_start = %d, pos_end = %d, seq len: %d\n' % (pos_start, pos_end, len(p_ctg_seq)), fp_proto_log)

    proto_log('[Proto] Extracting subpath.\n', fp_proto_log)
    new_path, new_start_coord, new_end_coord = p_ctg_tiling_path.get_subpath(pos_start, pos_end)

    proto_log('[Proto] Extracting the seq.\n', fp_proto_log)
    new_seq = p_ctg_seq[pos_start:pos_end]

    proto_log('[Proto] len(new_path) = %d\n' % (len(new_path)), fp_proto_log)
    first_edge = new_path[0]
    last_edge = new_path[-1]

    proto_log('[Proto] Forming the haplotig.\n', fp_proto_log)
    new_haplotig = Haplotig(name = htig_name, phase = complete_phase, seq = new_seq, path = new_path, edges = [])
    new_haplotig.labels['start_in_path'] = new_start_coord
    new_haplotig.labels['end_in_path'] = new_end_coord

    region_htigs = {htig_name: new_haplotig.__dict__}
    new_region = (region_type, first_edge, last_edge, pos_start, pos_end, region_htigs)

    proto_log('[Proto] Exiting function: "%s"\n' % (str(__name__)), fp_proto_log)
    return new_region

def get_linear_regions_between_bubbles(ctg_id, bubble_regions, p_ctg_seq, p_ctg_tiling_path, fp_proto_log):
    """
    For a given set of non-overlapping bubble regions, creates a set of linear regions
    at the suffix, prefix and in between the bubbles, and returns them as a list.
    """
    ret_linear_regions = []

    sorted_bubble_regions = sorted(bubble_regions, key = lambda x: x[3])

    proto_log('[Proto] Function: "%s"\n' % (str(__name__)), fp_proto_log)
    proto_log('len(sorted_bubble_regions) = %d\n' % (len(sorted_bubble_regions)), fp_proto_log)

    if len(sorted_bubble_regions) == 0:
        pos_start = 0
        pos_end = len(p_ctg_seq)
        new_region = make_linear_region(ctg_id, pos_start, pos_end, p_ctg_seq, p_ctg_tiling_path, fp_proto_log)
        if new_region != None:
            ret_linear_regions.append(new_region)
        return ret_linear_regions

    # Handle the prefix linear region.
    proto_log('Handling prefix.\n', fp_proto_log)
    region1_type, region1_first_edge, region1_last_edge, region1_pos_start, region1_pos_end, region1_htigs = sorted_bubble_regions[0]
    pos_start = 0
    pos_end = region1_pos_start
    new_region = make_linear_region(ctg_id, pos_start, pos_end, p_ctg_seq, p_ctg_tiling_path, fp_proto_log)
    if new_region != None:
        ret_linear_regions.append(new_region)

    proto_log('Handling infix.\n', fp_proto_log)
    if len(sorted_bubble_regions) > 1:
        temp_num_regions = 0
        for region1, region2 in zip(sorted_bubble_regions[0:-1], sorted_bubble_regions[1:]):
            temp_num_regions += 1
            proto_log('Checking for in-between region %d / %d.\n' % (temp_num_regions, len(sorted_bubble_regions)), fp_proto_log)
            region1_type, region1_first_edge, region1_last_edge, region1_pos_start, region1_pos_end, region1_htigs = region1
            region2_type, region2_first_edge, region2_last_edge, region2_pos_start, region2_pos_end, region2_htigs = region2
            pos_start = region1_pos_end
            pos_end = region2_pos_start
            new_region = make_linear_region(ctg_id, pos_start, pos_end, p_ctg_seq, p_ctg_tiling_path, fp_proto_log)
            if new_region != None:
                proto_log('Generated the new region.\n', fp_proto_log)
                ret_linear_regions.append(new_region)

    # Handle the suffix linear region.
    proto_log('Handling suffix.\n', fp_proto_log)
    region2_type, region2_first_edge, region2_last_edge, region2_pos_start, region2_pos_end, region2_htigs = sorted_bubble_regions[-1]
    pos_start = region2_pos_end
    pos_end = len(p_ctg_seq)
    new_region = make_linear_region(ctg_id, pos_start, pos_end, p_ctg_seq, p_ctg_tiling_path, fp_proto_log)
    if new_region != None:
        ret_linear_regions.append(new_region)

    proto_log('Done!\n', fp_proto_log)

    return ret_linear_regions

def regions_to_haplotig_graph(ctg_id, all_regions, phase_alias_map, fp_proto_log):
    """
    Takes a list of non-overlapping (but adjacent) regions (basically, a DAG),
    and creates a NetworkX object which describes the relationship between all
    haplotigs.
    """
    all_regions = sorted(all_regions, key = lambda x: x[3])

    haplotig_graph = nx.DiGraph()
    # Add nodes.
    proto_log('  - Adding nodes.\n', fp_proto_log)
    for region_id in xrange(len(all_regions)):
        region = all_regions[region_id]
        region_type, edge_start, edge_end, region_pos_start, region_pos_end, region_haplotigs = region
        if region_type == 'complex':
            # Complex regions are sketchy, better just break the graph here.
            continue
        proto_log('    - region_id = %d, region_type = %s, region_pos_start = %d, region_pos_end = %d\n' % (region_id, region_type, region_pos_start, region_pos_end), fp_proto_log)
        for key, htig in region_haplotigs.iteritems():
            proto_log('      - key = %s\n' % (key), fp_proto_log)
            v = htig['name']
            vphase = htig['phase']
            vphase_alias = phase_alias_map.get(vphase, -1)
            haplotig_graph.add_node(v, label='%s_%s_%s_%s' % (ctg_id, region_type, region_pos_start, region_pos_end),
                        phase='_'.join([str(val) for val in vphase]),
                        phase_alias=vphase_alias)

    # Add edges between all bubble components. Filtering the edges
    # will be performed afterwards.
    proto_log('  - Adding edges.\n', fp_proto_log)
    for region_id in xrange(1, len(all_regions)):
        region = all_regions[region_id - 1]
        region_type, edge_start, edge_end, region_pos_start, region_pos_end, region_haplotigs = region
        next_region = all_regions[region_id]
        next_region_type, next_edge_start, next_edge_end, next_region_pos_start, next_region_pos_end, next_region_haplotigs = next_region
        if region_type == 'complex' or next_region_type == 'complex':
            # Complex regions are sketchy, better just break the graph here.
            continue
        for htig_name, htig in region_haplotigs.iteritems():
            v = htig_name
            for next_htig_name, next_htig in next_region_haplotigs.iteritems():
                w = next_htig_name
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

def extract_and_write_all_ctg(ctg_id, haplotig_graph, all_haplotig_dict, phase_alias_map, all_regions, out_dir, fp_proto_log):
    """
    Finds an arbitrary (shortest) walk down the haplotig DAG, and denotes it
    as "p_ctg.fa".
    It then removes all p_ctg edges from the graph, and all cross-edges,
    and outputs all other weakly connected components as haplotigs in "h_ctg.fa"
    """
    # Hash all node sequences for easier fetching.
    node_seqs = {}
    for region_id in xrange(len(all_regions)):
        region = all_regions[region_id]
        region_type, edge_start, edge_end, region_pos_start, region_pos_end, region_haplotigs = region
        for key, htig in region_haplotigs.iteritems():
            v = htig['name']
            node_seqs[v] = htig['seq']

    ##################################
    # Write all haplotigs for possible
    # future use.
    ##################################
    proto_log('  - Writing all the haplotigs to disk in haplotigs.fasta.\n', fp_proto_log)
    with open(os.path.join(out_dir, "haplotigs.fasta"), 'w') as fp_out:
        for qname, qseq in node_seqs.iteritems():
            fp_out.write('>%s\n%s\n' % (qname, qseq))

    ##################################
    # Extract the primary collapsed haplotig.
    ##################################
    proto_log('[Proto] Extracting the primary contigs.\n', fp_proto_log)

    # More than one path may be found, if there were multiple connected
    # components. For example, if there was a complex bubble, the graph
    # would be broken.
    best_paths = extract_weakly_unphased_haplotig_paths(haplotig_graph)
    if best_paths == None:
        return
    with open(os.path.join(out_dir, "p_ctg.%s.fa" % ctg_id), "w") as fp_p_tig_fa, \
         open(os.path.join(out_dir, "p_ctg_path.%s" % ctg_id), "w") as fp_p_tig_path, \
         open(os.path.join(out_dir, "p_ctg_edges.%s" % ctg_id), "w") as fp_p_tig_edges:
        for i, best_path in enumerate(best_paths):
            # Write the primary contig.
            total_weight, path, s_node, t_node = best_path
            # Node here is a haplotig, not a pread.
            seq = ''
            new_ctg_id = '%sp%d' % (ctg_id, i)
            for v in path:
                seq += node_seqs[v]
            fp_p_tig_fa.write('>%s\n%s\n' % (new_ctg_id, seq))

            for v in path:
                haplotig = all_haplotig_dict[v]
                phase = haplotig['phase']
                tp = haplotig['path']
                proto_log('[Proto] v = %s\n' % (v), fp_proto_log)

                for edge in tp:
                    # Example tiling path line from 2-asm-falcon:
                    #   ["000000F", "000000040:E", "000000217:E", "000000217", "14608", "33796", "14608", "99.86"]
                    # Example path line from 3-unzip/1-hasm/000000F:
                    #   000000F 002560539:B 002559985:B 002559985 12404 0 8446 97.93 1 0
                    edge_ctg_id, edge_v, edge_w, edge_wid, edge_b, edge_e, edge_score, edge_idt = edge[0:8]

                    # Compose the path line.
                    new_path_line = [new_ctg_id, edge_v, edge_w, edge_wid, edge_b, edge_e, edge_score, edge_idt, phase[1], phase[2]]
                    fp_p_tig_path.write(' '.join([str(val) for val in new_path_line]) + '\n')

                    # Compose an edge line. All pread nodes within the same haplotig
                    # are phased in the same way.
                    cross_phase = 'N'
                    source_graph = 'H'
                    new_edge_line = [new_ctg_id, edge_v, edge_w, cross_phase, source_graph, phase[1], phase[2], phase[1], phase[2]]
                    fp_p_tig_edges.write(' '.join([str(val) for val in new_edge_line]) + '\n')

                proto_log('[Proto] v = %s done.\n' % (v), fp_proto_log)

            # print >>p_tig_path, "%s" % ctg_id, v, w, seq_id, s, t, edge_data[1], edge_data[2], "%d %d" % arid_to_phase.get(seq_id, (-1, 0))
            # print >>f, "%s" % ctg_id, v, w, sg[v][w]["cross_phase"], sg[v][w]["src"], vphase[0], vphase[1], wphase[0], wphase[1]

    #####################################################
    # Create and write the associate contigs.
    #####################################################
    proto_log('[Proto] Extracting the associate haplotigs.\n', fp_proto_log)

    haplotig_graph2 = haplotig_graph.copy()

    # Remove the primary path from the graph.
    for best_path in best_paths:
        total_weight, path, s_node, t_node = best_path
        for v, w in zip(path[:-1], path[1:]):
            haplotig_graph2.remove_edge(v, w)
        for v in path:
            haplotig_graph2.remove_node(v)

    # Mark any non-phased edges and contigs for removal.
    # We want the haplotigs to be only the phased components.
    nodes_to_remove = set()
    edges_to_remove = set()
    for v, w in haplotig_graph2.edges():
        # Get the phase for v, and check if it's phased.
        v_htig = all_haplotig_dict[v]
        vphase = v_htig['phase']
        vphase_str = '_'.join([str(val) for val in vphase])
        vphase_alias = phase_alias_map.get(vphase_str, -1)
        if vphase_alias == -1:
            nodes_to_remove.add(v)
            edges_to_remove.add((v, w))
        # Get the phase for w, and check if it's phased.
        w_htig = all_haplotig_dict[w]
        wphase = w_htig['phase']
        wphase_str = '_'.join([str(val) for val in wphase])
        wphase_alias = phase_alias_map.get(wphase_str, -1)
        if wphase_alias == -1:
            nodes_to_remove.add(v)
            edges_to_remove.add((v, w))
        # Remove edges which do not resolve to the same alias.
        if vphase_alias != wphase_alias:
            edges_to_remove.add((v, w))

    for v, w in edges_to_remove:
        haplotig_graph2.remove_edge(v, w)
    for v in nodes_to_remove:
        haplotig_graph2.remove_node(v)

    proto_log('  - After removing the primary path and all cross-phase edges.\n', fp_proto_log)

    # Each connected component is one haplotig.
    htig_id = 0
    with open(os.path.join(out_dir, "h_ctg.%s.fa" % ctg_id), "w") as fp_h_tig_fa, \
         open(os.path.join(out_dir, "h_ctg_path.%s" % ctg_id), "w") as fp_h_tig_path, \
         open(os.path.join(out_dir, "h_ctg_edges.%s" % ctg_id), "w") as fp_h_tig_edges:
        best_paths = extract_weakly_unphased_haplotig_paths(haplotig_graph2)
        if best_paths == None:
            return
        for htig_id, best_path in enumerate(best_paths):
            total_weight, path, s_node, t_node = best_path
            new_ctg_id = '%s_%d' % (ctg_id, (htig_id + 1))
            seq = ''
            for v in path:
                seq += node_seqs[v]
            fp_h_tig_fa.write('>%s\n%s\n' % (new_ctg_id, seq))

            for v in path:
                haplotig = all_haplotig_dict[v]
                phase = haplotig['phase']
                tp = haplotig['path']
                for edge in tp:
                    # Example tiling path line from 2-asm-falcon:
                    #   ["000000F", "000000040:E", "000000217:E", "000000217", "14608", "33796", "14608", "99.86"]
                    # Example path line from 3-unzip/1-hasm/000000F:
                    #   000000F 002560539:B 002559985:B 002559985 12404 0 8446 97.93 1 0
                    edge_ctg_id, edge_v, edge_w, edge_wid, edge_b, edge_e, edge_score, edge_idt = edge[0:8]

                    # Compose the path line.
                    new_path_line = [new_ctg_id, edge_v, edge_w, edge_wid, edge_b, edge_e, edge_score, edge_idt, phase[1], phase[2]]
                    fp_h_tig_path.write(' '.join([str(val) for val in new_path_line]) + '\n')

                    # Compose an edge line. All pread nodes within the same haplotig
                    # are phased in the same way.
                    cross_phase = 'N'
                    source_graph = 'H'
                    new_edge_line = [new_ctg_id, edge_v, edge_w, cross_phase, source_graph, phase[1], phase[2], phase[1], phase[2]]
                    fp_h_tig_edges.write(' '.join([str(val) for val in new_edge_line]) + '\n')

def run(args):
    # make life easier for now. will refactor it out if possible
    global all_rid_to_phase
    global all_flat_rid_to_phase
    global all_haplotigs_for_ctg
    global p_asm_G
    global h_asm_G
    global seqs
    global seq_lens
    global p_ctg_seqs
    global sg_edges
    global p_ctg_tiling_paths

    fc_asm_path = args.fc_asm_path
    fc_hasm_path = args.fc_hasm_path
    ctg_id = args.ctg_id
    base_dir = args.base_dir
    fasta_fn = args.fasta
    hasm_falcon_path = os.path.join(fc_hasm_path, 'asm-falcon')

    p_asm_G = AsmGraph(os.path.join(fc_asm_path, "sg_edges_list"),
                       os.path.join(fc_asm_path, "utg_data"),
                       os.path.join(fc_asm_path, "ctg_paths"))
    h_asm_G = AsmGraph(os.path.join(fc_hasm_path, "sg_edges_list"),
                       os.path.join(fc_hasm_path, "utg_data"),
                       os.path.join(fc_hasm_path, "ctg_paths"))
    assert p_asm_G, 'Empty AsmGraph. Maybe empty inputs?\n{!r}\n{!r}\n{!r}'.format(
        os.path.join(fc_asm_path, "sg_edges_list"),
        os.path.join(fc_asm_path, "utg_data"),
        os.path.join(fc_asm_path, "ctg_paths"),
    )
    all_rid_to_phase = {}
    all_flat_rid_to_phase = {}
    all_read_ids = set()
    with open(args.rid_phase_map) as f:
        for row in f:
            row = row.strip().split()
            all_rid_to_phase.setdefault(row[1], {})
            all_rid_to_phase[row[1]][row[0]] = (int(row[2]), int(row[3]))
            all_flat_rid_to_phase[row[0]] = (row[1], int(row[2]), int(row[3]))
            all_read_ids.add(row[0])
    assert all_read_ids, 'Empty all_read_ids and all_rid_to_phase. Maybe empty rid_phase_map file? {!r}'.format(
        args.rid_phase_map)
    for v, w in p_asm_G.sg_edges:
        if p_asm_G.sg_edges[(v, w)][-1] != "G":
            continue
        v = v.split(":")[0]
        w = w.split(":")[0]
        all_read_ids.add(v)
        all_read_ids.add(w)

    for v, w in h_asm_G.sg_edges:
        if h_asm_G.sg_edges[(v, w)][-1] != "G":
            continue
        v = v.split(":")[0]
        w = w.split(":")[0]
        all_read_ids.add(v)
        all_read_ids.add(w)

    # Load the preads.
    seqs = load_sg_seq(all_read_ids, fasta_fn)

    # Load the primary contig sequences.
    p_ctg_seqs = load_all_seq(os.path.join(fc_asm_path, "p_ctg.fa"))

    # Hash the lengths of the primary contig sequences.
    # Needed to correctly assign node coords when loading tiling paths
    p_ctg_seq_lens = {}
    for p_ctg_id, ctg_seq in p_ctg_seqs.iteritems():
        p_ctg_seq_lens[p_ctg_id] = len(ctg_seq)

    # Load the tiling path of the primary contig, and assign coordiantes to nodes.
    p_ctg_tiling_paths = tiling_path.load_tiling_paths(os.path.join(fc_asm_path, "p_ctg_tiling_path"), None, p_ctg_seq_lens)

    # Load the haplotig sequences (assembled with falcon_kit.mains.graph_to_contig).
    all_haplotigs_for_ctg, htig_name_to_original_pctg = load_haplotigs(hasm_falcon_path, all_flat_rid_to_phase, seqs)
    # all_asm_haplotig_seqs = load_all_seq(hasm_falcon_path)

    sys.stderr.write('[Proto] Loaded haplotigs.\n')

    # Load all sg_edges_list so that haplotig paths can be reversed if needed.
    sg_edges = {}
    with open(os.path.join(fc_hasm_path, "sg_edges_list"), 'r') as fp:
        for line in fp:
            sl = line.strip().split()
            sg_edges[(sl[0], sl[1])] = sl

    sys.stderr.write('[Proto] Loaded sg_edges_list.\n')

    # Hash the lengths of the preads.
    seq_lens = {}
    for key, val in seqs.iteritems():
        l = len(val)
        seq_lens[key] = l
        seq_lens[key + ':B'] = l
        seq_lens[key + ':E'] = l

    sys.stderr.write('[Proto] Loaded seq lens.\n')

    if ctg_id == "all":
        ctg_id_list = p_asm_G.ctg_data.keys()
    else:
        ctg_id_list = [ctg_id]

    sys.stderr.write('[Proto] Creating the exe list for: %s\n' % (str(ctg_id_list)))

    exe_list = []
    for ctg_id in ctg_id_list:
        if ctg_id[-1] != "F":
            continue
        if ctg_id not in all_rid_to_phase:
            continue
        exe_list.append((ctg_id, os.path.join(".", ctg_id)))

    sys.stderr.write('[Proto] Running jobs.\n')

    exec_pool = Pool(args.nproc)  # TODO, make this configurable
    exec_pool.map(run_generate_haplotigs_for_ctg, exe_list)
    #map( generate_haplotigs_for_ctg, exe_list)


######

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='layout haplotigs from primary assembly graph and phased aseembly graph')

    parser.add_argument(
        '--fc-asm-path', type=str,
        help='path to the primary Falcon assembly output directory', required=True)
    parser.add_argument(
        '--fc-hasm-path', type=str,
        help='path to the phased Falcon assembly output directory', required=True)
    parser.add_argument(
        '--ctg-id', type=str, help='contig identifier in the bam file', default="all", required=True)
    parser.add_argument(
        '--base-dir', type=str, default="./",
        help='the output base_dir, default to current working directory')
    parser.add_argument(
        '--rid-phase-map', type=str,
        help="path to the file that encode the relationship of the read id to phase blocks", required=True)
    parser.add_argument(
        '--fasta', type=str, help="sequence file of the p-reads", required=True)
    parser.add_argument(
        '--nproc', type=int, default=8, help="number of processes to use")

    args = parser.parse_args(argv[1:])

    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig()
    run(args)


if __name__ == '__main__':  # pragma: no cover
    main()
