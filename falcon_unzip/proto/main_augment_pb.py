#! /usr/bin/env python2.7

import os
import sys
import argparse
from falcon_kit import tiling_path
import phasing_block
import intervaltree.intervaltree as intervaltree
#from intervaltree import * # Let's avoid this.
import json
import networkx as nx
import sam2m4
import collections
import copy
import falcon_kit.FastaReader
import haplotig
import logging
import time
LOG = logging.getLogger() # root, to inherit from sub-loggers

def load_seqs(seq_path_list):
    seqs = {}
    seq_lens = {}
    for seq_path in seq_path_list:
        fp_seqs = falcon_kit.FastaReader.FastaReader(seq_path)
        for record in fp_seqs:
            header = record.name.split()[0]
            seq = record.sequence.upper()
            seqs[header] = seq
            seq_lens[header] = len(seq)
    return seqs, seq_lens

def calc_overlap_len(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))

# Find the region with the largest overlap, and assign the read to that region.
# A pread could map to multiple regions if it's really really long.
def find_max_linear_region(all_regions, found_regions, reference_start, reference_end):
    max_linear_region = None # regions[0].data
    max_overlap_len = 0
    for interval in found_regions:
        region_id = interval.data
        type_, edge_start_id, edge_end_id, pos_start, pos_end, htigs = all_regions[region_id]
        overlap_len = calc_overlap_len(pos_start, pos_end, reference_start, reference_end)
        if max_linear_region == None or overlap_len > max_overlap_len:
            max_overlap_len = overlap_len
            max_linear_region = region_id
    return max_linear_region

def create_new_rid2phase(ctg_id, rid2phase, all_regions, m4, a_paths):
    # Hash the regions in an intervaltree
    intervals = []
    for region_id in xrange(len(all_regions)):
        type_, edge_start_id, edge_end_id, pos_start, pos_end, htigs = all_regions[region_id]
        if type_ != 'linear':
            continue
        if pos_end <= pos_start:
            continue # See http://bitbucket.nanofluidics.com:7990/projects/SAT/repos/falcon_unzip_private/pull-requests/62/overview
        intervals.append(intervaltree.Interval(pos_start, pos_end, region_id))
    tree = intervaltree.IntervalTree(intervals)

    """
    When the rid2phase is updated, we need to preserve the
    relation between the old and the new phases, so that
    individual haplotigs can be stitched together later.

    Phases can have multiple aliases, which will be maintained
    as weakly connected components of the phase relation graph.
    Any phase belonging to the same weak component is actually
    the same phase.
    """
    phase_relation_graph = nx.DiGraph()

    # Initialize every read as unphased.
    new_rid2phase = {}
    for rid, phase in rid2phase.iteritems():
        r_ctg_id, phase_block_id, phase_id = phase
        new_rid2phase[rid] = (r_ctg_id, -1, 0)

    # For any read which maps at least partially to a linear region, update the
    # phasing block, so that each linear region is assembled separately.
    for aln in m4:
        rid, r_ctg_id, reference_start, reference_end = aln[0], aln[1], aln[9], aln[10]

        # If there is no phasing info, treat it like unphased.
        if rid not in rid2phase:
            continue

        r_ctg_id, phase_block_id, phase_id = rid2phase[rid]

        if phase_block_id == '-1':
            continue

        found_regions = tree.search(reference_start, reference_end)

        # If the read maps to a linear region at least partially,
        # process it further.
        if found_regions:
            if len(found_regions) > 1:
                LOG.info('[Proto] Read %s has %d regions.' % (rid, len(found_regions)))

            # Change the ID of the phasing block, by adding a large prefix value.
            max_linear_region = find_max_linear_region(all_regions, found_regions, reference_start, reference_end)
            type_, edge_start_id, edge_end_id, pos_start, pos_end, htigs = all_regions[max_linear_region]
            # new_phase = htigs[htigs.keys()[0]]['phase']

            region_phase_block_id = str((max_linear_region + 1) * 1000000 + int(phase_block_id))
            new_phase = (r_ctg_id, region_phase_block_id, phase_id)
            new_rid2phase[rid] = new_phase
            old_phase = rid2phase[rid]

            # Add the phases as nodes in the graph, and add edges.
            old_phase_str = '_'.join([str(val) for val in old_phase])
            new_phase_str = '_'.join([str(val) for val in new_phase])
            phase_relation_graph.add_node(old_phase_str)
            phase_relation_graph.add_node(new_phase_str)
            phase_relation_graph.add_edge(old_phase_str, new_phase_str)
            phase_relation_graph.add_edge(new_phase_str, old_phase_str)

    # Each SV bubble with 2 phases is already phased.
    # Assign new phase to each branch, and add the relations
    # to the alias graph.
    def get_v_phase(edge):
        v, w = edge[1], edge[2]
        vrid = v.split(':')[0]
        # If there is no phasing info, treat it like unphased.
        vphase = rid2phase.get(vrid, (ctg_id, '-1', '0'))
        vphase = (vphase[0], vphase[1], int(vphase[2]))
        return vphase

    new_htig_to_phase = {}
    for region_id in xrange(len(all_regions)):
        type_, edge_start_id, edge_end_id, pos_start, pos_end, htigs = all_regions[region_id]
        if type_ != 'simple':
            continue

        # For each branch of the simple bubble, collect all phases and create
        # a lookup of haplotigs they are related to. If each phase is related to
        # more than 1 haplotig, then there was a problem in the phasing process,
        # and we need to handle this special case.
        phase_to_htig = {}
        for htig_name, htig in htigs.iteritems():
            htig_path = htig['path']
            for edge_id in xrange(1, len(htig_path)):
                edge = htig_path[edge_id]
                vphase = get_v_phase(edge)
                if vphase[1] == '-1':
                    continue
                phase_to_htig.setdefault(vphase, set())
                phase_to_htig[vphase].add(htig_name)

        for htig_name, htig in htigs.iteritems():
            htig_path = htig['path']
            htig_phase = htig['phase']

            # region_phase_block_id = str((region_id + 1) * 1000000)
            # new_phase = (ctg_id, region_phase_block_id, phase_id)
            # htig_name = htigs.keys()[phase_id]
            new_htig_to_phase[htig_name] = htig_phase
            htig_phase_str = '_'.join([str(val) for val in htig_phase])
            phase_relation_graph.add_node(htig_phase_str)

            # The very first and very last edge of each bubble could already
            # be phased in a different phase (the bubble merges to a collapsed
            # region). We need to skip them in order
            # to not create fused graph components.
            #
            # First, check that the phasing of the bubble was sane.
            branch_phases = set()
            is_sane = True
            for edge_id in xrange(1, len(htig_path)):   # Only `v` of each edge is used, so the last `w` is skipped.
                edge = htig_path[edge_id]
                vphase = get_v_phase(edge)
                if vphase[1] == '-1':
                    continue
                other_vphase = (vphase[0], vphase[1], 1 - int(vphase[2]))   # Diploid
                # Check if the other phase is present in the same branch.
                if other_vphase in branch_phases:
                    is_sane = False
                    break
                # Check that the same phase is not present in more than 1 branch.
                # This would mean that the phasing process created phased reads
                # of the same phase but in different branches of the SV bubbles,
                # which should not be possible.
                # But if this case would happen (and it happens in practice),
                # the phase relation graph would become fused.
                if vphase not in phase_to_htig or len(phase_to_htig[vphase]) != 1:
                    is_sane = False
                    break

                branch_phases.add(vphase)

            # Don't add aliases for insanely phased bubbles.
            if is_sane == False:
                continue

            for edge_id in xrange(1, len(htig_path)):
            # for edge in htig_path:
                edge = htig_path[edge_id]
                vphase = get_v_phase(edge)
                if vphase[1] == '-1':
                    continue
                vphase_str = '_'.join([str(val) for val in vphase])
                # Add the phases as nodes in the graph, and add edges.
                phase_relation_graph.add_node(vphase_str)
                phase_relation_graph.add_edge(htig_phase_str, vphase_str)
                phase_relation_graph.add_edge(vphase_str, htig_phase_str)

    return new_rid2phase, phase_relation_graph, new_htig_to_phase




def construct_new_bubble_region(ctg_id, path_type, region_id, p_path, a_paths, path_start, path_end, bubble_branches):

    pos_start = 0 if path_start == 0 else p_path.coords[p_path.edges[path_start - 1].w]
    pos_end = p_path.coords[p_path.edges[path_end - 1].w]

    bubble_a_paths = {}

    phase_id = 0
    for key in bubble_branches:
        region_phase_block_id = str((region_id + 1) * 1000000)
        new_phase = (ctg_id, region_phase_block_id, phase_id)
        path_htig = haplotig.Haplotig(key, phase = new_phase, seq = '', edges = [],
                                        path = [edge.get_split_line() for edge in a_paths[key].edges],
                                        labels = {}, cstart = pos_start, cend = pos_end)
        bubble_a_paths[key] = path_htig.__dict__
        phase_id += 1

    # bubble_a_paths = {key: [edge.split_line for edge in a_paths[key].edges] for key in bubble_branches}

    path_base_name = '%s_%s_%d_base' % (ctg_id, path_type, region_id)
    region_phase_block_id = str((region_id + 1) * 1000000)
    new_phase = (ctg_id, region_phase_block_id, phase_id)
    path_htig = haplotig.Haplotig(path_base_name, phase = new_phase, seq = '', edges = [],
                                        path = [edge.get_split_line() for edge in p_path.edges[path_start:path_end]],
                                        labels = {}, cstart = pos_start, cend = pos_end)
    bubble_a_paths[path_base_name] = path_htig.__dict__

    # bubble_a_paths[path_base_name] = [edge.split_line for edge in p_path.edges[bubble_start:bubble_end]]

    new_region = (path_type, path_start, path_end, pos_start, pos_end, bubble_a_paths)

    return new_region


def delineate_regions(ctg_id, p_path, a_paths, a_placement):
    """
    Performs a linear pass through a contig, and
    finds all bubbles. It classifies subpaths into:
      - linear
      - simple bubbles (multipath bubbles where all
        paths start and end in the same start and end node)
      - complex bubbles (strange furcations in the graph,
        such as interleaved bubbles).
    Usage:
    delineate_regions(p_paths[ctg_id], a_paths, a_placement[ctg_id])
    """
    # For each node of a p_path, keep track
    out_edges = {}
    in_edges = {}

    for a_ctg_id, placement in a_placement.iteritems():
        start, end, p_ctg_id, a_ctg_id, first_node, last_node = placement
        # Add the out edges.
        out_edges.setdefault(first_node, set())
        out_edges[first_node].add(a_ctg_id) # ((a_ctg_id, placement))
        # Add the in edges.
        in_edges.setdefault(last_node, set())
        in_edges[last_node].add(a_ctg_id) # ((a_ctg_id, placement))

    if len(p_path.edges) == 0:
        return

    # The state machine.
    STATE_LINEAR = 0
    STATE_SIMPLE_BUBBLE = 1
    STATE_COMPLEX_BUBBLE = 2
    current_state = STATE_LINEAR
    next_state = current_state
    # Keep trach of where a bubble starts once
    # we hit a fork. Works for both simple and
    # complex bubbles.
    bubble_start = None
    bubble_end = None
    bubble_branches = set()
    num_open_branches = 0
    all_regions = []
    linear_start = 0
    linear_end = 0

    v0 = p_path.edges[0].v
    if v0 in out_edges:
        next_state = STATE_SIMPLE_BUBBLE
        # open_set.union(out_edges[v0])
        num_open_branches += len(out_edges[v0])
        bubble_start = 0
        bubble_branches |= out_edges[v0]
        current_state = next_state

    for i in xrange(len(p_path.edges)):
        w = p_path.edges[i].w
        if current_state == STATE_LINEAR:
            if w in out_edges:
                # Add a linear region.
                path_type = 'linear'
                path_base_name = '%s_%s_%d_base' % (ctg_id, path_type, len(all_regions))

                linear_start = 0 if bubble_end == None else (bubble_end) # + 1)
                linear_end = i + 1
                pos_start = 0 if linear_start == 0 else p_path.coords[p_path.edges[linear_start - 1].w]
                pos_end = p_path.coords[p_path.edges[linear_end - 1].w]

                # linear_paths = {path_base_name: [edge.split_line for edge in p_path.edges[linear_start:linear_end]]}

                path_htig = haplotig.Haplotig(path_base_name, phase = (ctg_id, -1, 0), seq = '', edges = [],
                                                path = [edge.get_split_line() for edge in p_path.edges[linear_start:linear_end]],
                                                labels = {}, cstart = pos_start, cend = pos_end)

                linear_paths = {path_base_name: path_htig.__dict__}

                new_region = (path_type, linear_start, linear_end, pos_start, pos_end, linear_paths)

                if pos_end > pos_start: # Avoid zero-length regions.
                    all_regions.append(new_region)

                # Change the state.
                num_open_branches += len(out_edges[w])
                bubble_start = i + 1
                bubble_branches |= out_edges[w]
                next_state = STATE_SIMPLE_BUBBLE if num_open_branches == 1 else STATE_COMPLEX_BUBBLE

            if w in in_edges:
                msg = "Cycles in the graph? w={!r}, len(in_edges)={}. Not a fatal error - skipping this contig because it's too sketchy.".format(w, len(in_edges))
                LOG.warning(msg)
                all_regions = []
                return all_regions

        elif current_state == STATE_SIMPLE_BUBBLE:

            # First close the bubbles, in case new bubbles start branching
            # from this very node. In this case, we do not want to call a
            # complex bubble.

            if w in in_edges:
                num_open_branches -= len(in_edges[w])

                if num_open_branches == 0:
                    bubble_end = i + 1  # Non-inclusive end edge ID.

                    path_type = 'simple'
                    new_region = construct_new_bubble_region(ctg_id, path_type, len(all_regions), p_path, a_paths, bubble_start, bubble_end, bubble_branches)

                    _, _, _, pos_start, pos_end, _ = new_region
                    if pos_end > pos_start: # Avoid zero-length regions.
                        all_regions.append(new_region)

                    bubble_branches = set()
                    next_state = STATE_LINEAR

                else:
                    next_state = STATE_COMPLEX_BUBBLE

            if w in out_edges:
                if num_open_branches == 0:
                    num_open_branches += len(out_edges[w])
                    bubble_start = i + 1
                    next_state = STATE_SIMPLE_BUBBLE if num_open_branches == 1 else STATE_COMPLEX_BUBBLE

                else:
                    next_state = STATE_COMPLEX_BUBBLE
                    num_open_branches += len(out_edges[w])

                bubble_branches |= out_edges[w]

        elif current_state == STATE_COMPLEX_BUBBLE:
            if w in out_edges:
                next_state = STATE_COMPLEX_BUBBLE
                num_open_branches += len(out_edges[w])
                bubble_branches |= out_edges[w]

            if w in in_edges:
                num_open_branches -= len(in_edges[w])

                if num_open_branches == 0:
                    bubble_end = i + 1  # Non-inclusive end edge ID.

                    path_type = 'complex'
                    new_region = construct_new_bubble_region(ctg_id, path_type, len(all_regions), p_path, a_paths, bubble_start, bubble_end, bubble_branches)

                    _, _, _, pos_start, pos_end, _ = new_region
                    if pos_end > pos_start: # Avoid zero-length regions.
                        all_regions.append(new_region)

                    bubble_branches = set()
                    next_state = STATE_LINEAR

                else:
                    next_state = STATE_COMPLEX_BUBBLE

        current_state = next_state

    # Handle the leftover path.
    if current_state != STATE_LINEAR:
        bubble_end = len(p_path.edges)

        path_type = 'simple' if (current_state == STATE_SIMPLE_BUBBLE) else 'complex'
        new_region = construct_new_bubble_region(ctg_id, path_type, len(all_regions), p_path, a_paths, bubble_start, bubble_end, bubble_branches)

        _, _, _, pos_start, pos_end, _ = new_region
        if pos_end > pos_start: # Avoid zero-length regions.
            all_regions.append(new_region)

    else:
        path_type = 'linear'
        path_base_name = '%s_%s_%d_base' % (ctg_id, path_type, len(all_regions))

        linear_start = 0 if bubble_end == None else (bubble_end) # + 1)
        linear_end = len(p_path.edges)
        pos_start = 0 if linear_start == 0 else p_path.coords[p_path.edges[linear_start - 1].w]
        pos_end = p_path.coords[p_path.edges[linear_end - 1].w]

        path_htig = haplotig.Haplotig(path_base_name, phase = (ctg_id, -1, 0), seq = '', edges = [],
                                        path = [edge.get_split_line() for edge in p_path.edges[linear_start:linear_end]],
                                        labels = {}, cstart = pos_start, cend = pos_end)

        linear_paths = {path_base_name: path_htig.__dict__}

        new_region = (path_type, linear_start, linear_end, pos_start, pos_end, linear_paths)

        _, _, _, pos_start, pos_end, _ = new_region
        if pos_end > pos_start: # Avoid zero-length regions.
            all_regions.append(new_region)

    return all_regions

def load_preads_alignments(sam_path):
    # Load alignments, convert to M4 for easier coordinate lookup,
    # and filter multiple mappings.
    qname_dict = {}
    for aln in sam2m4.sam_to_m4(sam_path):
        qname, tname = aln[0], aln[1]
        q_orient, q_start, q_end, q_len = aln[4], aln[5], aln[6], aln[7]
        t_orient, t_start, t_end, t_len = aln[8], aln[9], aln[10], aln[11]
        sam = aln[13]

        # Take only one alignment per qname, the longest spanning one.
        if qname in qname_dict:
            prev_aln = qname_dict[qname]
            prev_t_start, prev_t_end = prev_aln[9], prev_aln[10]
            if (t_end - t_start) < (prev_t_end - prev_t_start):
                continue
        qname_dict[qname] = aln
        # qname_set.add(qname)
        # m4.append(aln)
    m4 = [val for key, val in qname_dict.iteritems()]
    m4 = sorted(m4, key = lambda x: x[9])

    return m4

def run(wd, ctg_id, extracted_ctg_fasta, p_ctg, p_ctg_tiling_path, a_ctg, a_ctg_tiling_path, p_variant_fn, rid_phase_map, preads_sam, rawread_bam, read_to_contig_map, out_updated_rid_phase_map, threads):
    num_threads = threads

    ###################################################
    # Load all contigs and their lengths.
    ###################################################
    seqs, seq_lens = load_seqs([p_ctg, a_ctg])

    seqs_for_ctg = {key: val for key, val in seqs.iteritems() if key.startswith(ctg_id)}

    ###################################################
    # Load tiling paths with coords, and filter the
    # deduplicated tiling paths (based on their presence
    # in the a_ctg.fa file).
    ###################################################
    p_paths = tiling_path.load_tiling_paths(p_ctg_tiling_path, whitelist_seqs=seqs_for_ctg, contig_lens=seq_lens)
    a_paths = tiling_path.load_tiling_paths(a_ctg_tiling_path, whitelist_seqs=seqs_for_ctg, contig_lens=None)
    # Find the associate contig placement. `a_placement` is a dict:
    #   placement[p_ctg_id][a_ctg_id] = (start, end, p_ctg_id, a_ctg_id, first_node, last_node)
    a_placement = tiling_path.find_a_ctg_placement(p_paths, a_paths)

    ###################################################
    # Classify regions in the p_ctg tiling path as either
    # linear or bubbles (simple and complex).
    # Each region has the following values:
    #   (type_, edge_start_id, edge_end_id, pos_start, pos_end, paths)
    # where type_ is one of 'linear', 'simple' or 'complex,
    # and paths is a dict with a mandatory 'base' key describing
    # the path belonging to the primary contig, and additional
    # keys corresponding th a_ctg paths.
    ###################################################
    a_placement_for_ctg = a_placement.get(ctg_id, {})
    all_regions = delineate_regions(ctg_id, p_paths[ctg_id], a_paths, a_placement_for_ctg)
    regions_json = os.path.join(wd, 'regions.json')
    with open(regions_json, 'w') as fp_out:
        fp_out.write(json.dumps(all_regions))

    ###################################################
    # Load the phasing information for every raw read.
    ###################################################
    rid2phase = phasing_block.load_rid_to_phase(rid_phase_map)

    ###################################################
    # Extract preads for this contig and map them to
    # the contig, to get the coordinates of each pread.
    # This is needed to filter preads which map
    # to bubbles.
    ###################################################
    m4 = load_preads_alignments(preads_sam)

    ###################################################
    # Modify the phases.
    # This is the most important part to ensure that
    # the linear regions will get assembled into haplotigs.
    ###################################################
    updated_rid2phase, phase_relation_graph, htig_to_phase = create_new_rid2phase(ctg_id, rid2phase, all_regions, m4, a_paths)

    # Write the phase relation graph to disk.
    nx.write_gexf(phase_relation_graph, os.path.join(wd, "phase_relation_graph.gexf"))

    ###################################################
    # Write the updated_rid_phase_map to be used
    # to construct the haplotigs.
    ###################################################
    with open(out_updated_rid_phase_map, 'w') as fp_out:
        for key, val in updated_rid2phase.iteritems():
            fp_out.write('%s %s %s %s\n' % (key, val[0], val[1], val[2]))

    ###################################################
    # Write the FASTA sequences for the linear
    # and the bubble regions.
    ###################################################
    minced_ctg_path = os.path.join(wd, 'minced.fasta')
    with open(minced_ctg_path, 'w') as fp_out:
        for region_id in xrange(len(all_regions)):
            type_, edge_start_id, edge_end_id, pos_start, pos_end, htigs = all_regions[region_id]
            if type_ == 'linear':
                region_name = htigs.keys()[0]
                region_seq = seqs[ctg_id][pos_start:pos_end]
                fp_out.write('>%s\n%s\n' % (region_name, region_seq))
            else:
                for htig_name, htig in htigs.iteritems():
                    region_name = htig_name
                    if htig_name.endswith('_base'):
                        # This is the 'base' part of any bubble. Literally, a region
                        # of the primary contig.
                        # The primary contigs should ideally be `proper`.
                        region_seq = seqs[ctg_id][pos_start:pos_end]
                        fp_out.write('>%s\n%s\n' % (region_name, region_seq))
                    else:
                        # Branches other than 'base' are actually the a_ctg sequences.
                        # All p_ctg and a_ctg are loaded in the seqs dict.
                        # If the region_name is not here, let the process fail.
                        # The a_ctg here should not contain the first read
                        # (they should be `improper`).
                        region_seq = seqs[region_name]
                        fp_out.write('>%s\n%s\n' % (region_name, region_seq))

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='This script prepares all info needed for prototyping the Unzip workflow.',  # better description?
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--wd', required=False,
        default='./',
        help='Working directory to store the output files.'
    )
    parser.add_argument(
        '--ctg-id', required=True,
        help='The contig ID to process.'
    )
    parser.add_argument(
        '--extracted-ctg-fasta', required=True,
        help='Path to the FASTA file containing extracted contig for analysis.'
    )
    parser.add_argument(
        '--rawread-bam', required=True,
        help='BAM file containing raw read to current ctg alignments.'
    )
    parser.add_argument(
        '--p-ctg', required=True,
        help='The primary contigs FASTA file.'
    )
    parser.add_argument(
        '--a-ctg', required=True,
        help='The associate contigs FASTA file.'
    )
    parser.add_argument(
        '--p-ctg-tiling-path', required=True,
        help='The primary contig tiling path from the 2-asm-falcon step'
    )
    parser.add_argument(
        '--a-ctg-tiling-path', required=True,
        help='The associate contig tiling path from the 2-asm-falcon step'
    )
    parser.add_argument(
        '--p-variant-fn', required=True,
        help='File containing a list of reads and assigned phasing blocks and phases.'
    )
    parser.add_argument(
        '--rid-phase-map', type=str,
        help='The file that encode the relationship of the read id to phase blocks.', required=True
    )
    parser.add_argument(
        '--preads-sam', type=str,
        help='SAM alignments of the preads to ref, for this particular ctg_id.', required=True
    )
    parser.add_argument(
        '--out-updated-rid-phase_map', type=str,
        help='Path where the augmented rid to phase map will be written.', required=True
    )
    parser.add_argument(
        # ./3-unzip/reads/get_read_ctg_map/read_to_contig_map
        '--read-to-contig-map', type=str, default="../../../3-unzip/reads/get_read_ctg_map/read_to_contig_map",
        help='read_to_contig_map, from fc_get_read_hctg_map'
    )
    parser.add_argument(
        '--threads', type=int,
        default=4,
        help='(Ignored. Single-threaded.) Number of threads to run the alignment with.'
    )

    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    logging.Formatter.converter = time.gmtime
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s %(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    args = parse_args(argv)
    run(**vars(args))

if __name__ == '__main__':  # pragma: no cover
    main()
