import falcon_unzip.mains.graphs_to_h_tigs_2 as mod
import helpers
import pytest
import os
import sys
import networkx as nx
from falcon_unzip.proto.haplotig import Haplotig
from falcon_kit.tiling_path import TilingPathEdge

"""
The haplotig_graph nodes are defined as:
            haplotig_graph.add_node(v, label='%s_%s_%s_%s' % (ctg_id, region_type, region_pos_start, region_pos_end),
                        phase='_'.join([str(val) for val in vphase]),
                        phase_alias=vphase_alias, htig=htig)
"""

#######################
### Utility methods ###
#######################
def mock_fp_proto_log(line):
    sys.stderr.write(line + '\n')
    pass

def make_dummy_linear_region(ctg_id, seq, pos_start, pos_end):
    region_type = 'linear'
    first_edge = None
    last_edge = None

    htig_name = '%s_linear2_%d:%d_base' % (ctg_id, pos_start, pos_end)
    complete_phase = (ctg_id, '-1', '0')
    path = [[ctg_id, '000000001:B', '000000002:B', '000000002', '4', '10', '1986', '99.95']]
    new_haplotig = Haplotig(name = htig_name, phase = complete_phase, seq = seq, path = path, edges = [])
    new_haplotig.cstart = pos_start
    new_haplotig.cend = pos_end

    new_region_htigs = {htig_name: new_haplotig.__dict__}

    new_region = (region_type, first_edge, last_edge, pos_start, pos_end, new_region_htigs)

    return new_region, htig_name

def make_dummy_diploid_region(ctg_id, seq1, seq2, pos_start, pos_end, phasing_block):
    region_type = 'diploid'
    first_edge = None
    last_edge = None

    htig_name_1 = '%s_diploid_%d:%d_phase0' % (ctg_id, pos_start, pos_end)
    complete_phase = (ctg_id, str(phasing_block), '0')
    path_1 = [[ctg_id, '000000003:B', '000000004:B', '000000004', '30', '6', '1986', '99.99']]
    new_haplotig_1 = Haplotig(name = htig_name_1, phase = complete_phase, seq = seq1, path = path_1, edges = [])
    # "Collapsed-start/end" are required to correctly determine placement coordinates.
    # All haplotigs from the same bubble have the same cstart and cend coordinates
    # (this is the coordinate on the collapsed primary contig (2-asm-falcon/p_ctg).
    new_haplotig_1.cstart = pos_start
    new_haplotig_1.cend = pos_end

    htig_name_2 = '%s_diploid_%d:%d_phase1' % (ctg_id, pos_start, pos_end)
    path_2 = [[ctg_id, '000000005:B', '000000006:B', '000000006', '8', '6', '2019', '100.00']]
    complete_phase = (ctg_id, str(phasing_block), '1')
    new_haplotig_2 = Haplotig(name = htig_name_2, phase = complete_phase, seq = seq2, path = path_2, edges = [])
    # "Collapsed-start/end" are required to correctly determine placement coordinates.
    # All haplotigs from the same bubble have the same cstart and cend coordinates
    # (this is the coordinate on the collapsed primary contig (2-asm-falcon/p_ctg).
    new_haplotig_2.cstart = pos_start
    new_haplotig_2.cend = pos_end

    new_region_htigs = {htig_name_1: new_haplotig_1.__dict__, htig_name_2: new_haplotig_2.__dict__}

    new_region = (region_type, first_edge, last_edge, pos_start, pos_end, new_region_htigs)

    return new_region, htig_name_1, htig_name_2

def make_source_and_sink_nodes(ctg_id):
    nodes = ['{ctg_id}-source-node'.format(ctg_id=ctg_id), '{ctg_id}-sink-node'.format(ctg_id=ctg_id)]
    return nodes

def evaluate_update_haplotig_graph(test_haplotig_graph, expected_nodes, expected_edges):
    expected_node_set = set(expected_nodes)
    expected_edge_set = set(expected_edges)

    # Evaluate nodes.
    assert len(test_haplotig_graph.nodes()) == len(expected_node_set)
    for v in test_haplotig_graph.nodes():
        assert v in expected_node_set

    # Evaluate edges.
    assert len(test_haplotig_graph.edges()) == len(expected_edge_set)
    for e in test_haplotig_graph.edges():
        assert e in expected_edge_set

def evaluate_extract_and_write_all_ctg(tmpdir, expected):
    # Check the number of generated files.
    assert len(tmpdir.listdir()) == 7

    # Compare generated files with expectations.
    fns = {os.path.basename(str(val)):str(val) for val in tmpdir.listdir()}
    for fn, expected_val in expected.iteritems():
        # Check if the file actually exists.
        assert fn in fns
        # Load the data and compare.
        result_val = open(fns[fn], 'r').read()
        assert result_val == expected_val, fn

#######################

def test_regions_to_haplotig_graph_1():
    """
    Test on empty input.
    """

    def create_test():
        ctg_id = '000000F'

        ### Inputs.
        regions = []

        ### Expected results.
        source_node, sink_node = make_source_and_sink_nodes(ctg_id)
        expected_nodes = [source_node, sink_node]
        expected_edges = []
        # Return.
        return ctg_id, regions, expected_nodes, expected_edges

    # Get inputs and outputs.
    ctg_id, regions, expected_nodes, expected_edges = create_test()

    # Run the unit under test.
    haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

    # Evaluate.
    evaluate_update_haplotig_graph(haplotig_graph, expected_nodes, expected_edges)

def test_regions_to_haplotig_graph_2():
    """
    Test on a single linear (collapsed) region.
    """

    def create_test():
        ctg_id = '000000F'

        ### Inputs.
        regions = []
        # Create a linear region.
        new_region, new_region_header = make_dummy_linear_region(ctg_id, 'ACTG', 0, 1000)
        regions.append(new_region)

        ### Expected results.
        source_node, sink_node = make_source_and_sink_nodes(ctg_id)
        expected_nodes = [new_region_header] + [source_node, sink_node]
        expected_edges = [  (source_node, new_region_header),
                            (new_region_header, sink_node)]
        # Return.
        return ctg_id, regions, expected_nodes, expected_edges

    # Get inputs and outputs.
    ctg_id, regions, expected_nodes, expected_edges = create_test()

    # Run the unit under test.
    haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

    # Evaluate.
    evaluate_update_haplotig_graph(haplotig_graph, expected_nodes, expected_edges)

def test_regions_to_haplotig_graph_3():
    """
    Test on two consecutive linear (collapsed) regions.
    """

    def create_test():
        ctg_id = '000000F'

        ### Inputs.
        regions = []
        # Create the first linear region.
        new_region, region_1_header = make_dummy_linear_region(ctg_id, 'ACTG', 0, 1000)
        regions.append(new_region)
        # Create the second linear region.
        new_region, region_2_header = make_dummy_linear_region(ctg_id, 'ACTG', 1000, 2000)
        regions.append(new_region)

        ### Expected results.
        source_node, sink_node = make_source_and_sink_nodes(ctg_id)
        expected_nodes = [region_1_header, region_2_header] + [source_node, sink_node]
        expected_edges = [  (source_node, region_1_header),
                            (region_1_header, region_2_header),
                            (region_2_header, sink_node)
                        ]
        # Return.
        return ctg_id, regions, expected_nodes, expected_edges

    # Get inputs and outputs.
    ctg_id, regions, expected_nodes, expected_edges = create_test()

    # Run the unit under test.
    haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

    # Evaluate.
    evaluate_update_haplotig_graph(haplotig_graph, expected_nodes, expected_edges)

def test_regions_to_haplotig_graph_4():
    """
    Test on a slightly more complicated case of:
        - linear region
        - diploid region
        - linear region
    """

    def create_test():
        ctg_id = '000000F'

        ### Inputs.
        regions = []
        # Create the first linear region.
        new_region, region_1_header = make_dummy_linear_region(ctg_id, 'ACTG', 0, 1000)
        regions.append(new_region)
        # Create the diploid region.
        new_region, region_2_header_1, region_2_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 1000, 2000, 0)
        regions.append(new_region)
        # Create the second linear region.
        new_region, region_3_header = make_dummy_linear_region(ctg_id, 'ACTG', 1000, 2000)
        regions.append(new_region)

        ### Expected results.
        source_node, sink_node = make_source_and_sink_nodes(ctg_id)
        expected_nodes = [region_1_header, region_2_header_1, region_2_header_2, region_3_header] + [source_node, sink_node]
        expected_edges = [  (source_node, region_1_header),
                            (region_1_header, region_2_header_1),
                            (region_1_header, region_2_header_2),
                            (region_2_header_1, region_3_header),
                            (region_2_header_2, region_3_header),
                            (region_3_header, sink_node)
                        ]
        # Return.
        return ctg_id, regions, expected_nodes, expected_edges

    # Get inputs and outputs.
    ctg_id, regions, expected_nodes, expected_edges = create_test()

    # Run the unit under test.
    haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

    # Evaluate.
    evaluate_update_haplotig_graph(haplotig_graph, expected_nodes, expected_edges)

def test_regions_to_haplotig_graph_5():
    """
    Test on a more complicated case of:
        - linear region
        - diploid region
        - diploid region
        - diploid region
        - linear region
    """

    def create_test():
        ctg_id = '000000F'

        ### Inputs.
        regions = []
        # Create the first linear region.
        new_region, region_1_header = make_dummy_linear_region(ctg_id, 'ACTG', 0, 1000)
        regions.append(new_region)
        # Create a diploid region.
        new_region, region_2_header_1, region_2_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 1000, 2000, 0)
        regions.append(new_region)
        # Create a diploid region.
        new_region, region_3_header_1, region_3_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 2000, 3000, 1)
        regions.append(new_region)
        # Create a diploid region.
        new_region, region_4_header_1, region_4_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 3000, 4000, 2)
        regions.append(new_region)
        # Create the second linear region.
        new_region, region_5_header = make_dummy_linear_region(ctg_id, 'ACTG', 4000, 5000)
        regions.append(new_region)

        ### Expected results.
        source_node, sink_node = make_source_and_sink_nodes(ctg_id)
        expected_nodes = [  region_1_header,
                            region_2_header_1, region_2_header_2,
                            region_3_header_1, region_3_header_2,
                            region_4_header_1, region_4_header_2,
                            region_5_header] + [source_node, sink_node]
        expected_edges = [  # Source to linear 1.
                            (source_node, region_1_header),
                            # Linear 1 to diploid 1.
                            (region_1_header, region_2_header_1),
                            (region_1_header, region_2_header_2),
                            # Diploid 1 to diploid 2.
                            (region_2_header_1, region_3_header_1),
                            (region_2_header_1, region_3_header_2),
                            (region_2_header_2, region_3_header_1),
                            (region_2_header_2, region_3_header_2),
                            # Diploid 2 to diploid 3.
                            (region_3_header_1, region_4_header_1),
                            (region_3_header_1, region_4_header_2),
                            (region_3_header_2, region_4_header_1),
                            (region_3_header_2, region_4_header_2),
                            # Diploid 3 to linear 2.
                            (region_4_header_1, region_5_header),
                            (region_4_header_2, region_5_header),
                            # Linear 2 to sink.
                            (region_5_header, sink_node)
                        ]
        # Return.
        return ctg_id, regions, expected_nodes, expected_edges

    # Get inputs and outputs.
    ctg_id, regions, expected_nodes, expected_edges = create_test()

    # Run the unit under test.
    haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

    # Evaluate.
    evaluate_update_haplotig_graph(haplotig_graph, expected_nodes, expected_edges)

def test_regions_to_haplotig_graph_6():
    """
    Test on a single diploid region.
    """

    # regions = []
    # ctg_id = '000000F'
    # # Create the diploid region.
    # new_region, new_header_1, new_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 1000, 2000, 0)
    # regions.append(new_region)

    # haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

    # assert len(haplotig_graph.nodes()) == 2
    # assert len(haplotig_graph.edges()) == 0

    def create_test():
        ctg_id = '000000F'

        ### Inputs.
        regions = []
        # Create a diploid region.
        new_region, region_1_header_1, region_1_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 1000, 2000, 0)
        regions.append(new_region)

        ### Expected results.
        source_node, sink_node = make_source_and_sink_nodes(ctg_id)
        expected_nodes = [
                            region_1_header_1, region_1_header_2,
                        ] + [source_node, sink_node]
        expected_edges = [  # Source to diploid.
                            (source_node, region_1_header_1),
                            (source_node, region_1_header_2),
                            # Diploid to sink.
                            (region_1_header_1, sink_node),
                            (region_1_header_2, sink_node),
                        ]
        # Return.
        return ctg_id, regions, expected_nodes, expected_edges

    # Get inputs and outputs.
    ctg_id, regions, expected_nodes, expected_edges = create_test()

    # Run the unit under test.
    haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

    # Evaluate.
    evaluate_update_haplotig_graph(haplotig_graph, expected_nodes, expected_edges)

def test_update_haplotig_graph_1():
    """
    This method should only remove edges, and retain all nodes.
    Test on a case with a single linear region and no additional edges.
    """

    def create_test():
        ctg_id = '000000F'
        phase_alias_map = {}

        # Compile a list of regions.
        regions = []

        # Convert regions to graph.
        haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

        # Expected results.
        # No nodes or edges should be removed.
        expected_node_set = set(haplotig_graph.nodes())
        expected_edge_set = set(haplotig_graph.edges())

        return haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set

    # Create test data.
    haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set = create_test()

    # Run the unit under test.
    updated_haplotig_graph = mod.update_haplotig_graph(haplotig_graph, phase_alias_map)

    # Evaluate
    evaluate_update_haplotig_graph(updated_haplotig_graph, expected_node_set, expected_edge_set)

def test_update_haplotig_graph_2():
    """
    This method should only remove edges, and retain all nodes.
    Test on a case with a single linear region and no additional edges.
    """

    def create_test():
        ctg_id = '000000F'
        phase_alias_map = {}

        # Compile a list of regions.
        regions = []
        new_region, new_region_header = make_dummy_linear_region(ctg_id, 'ACTG', 0, 1000)
        regions.append(new_region)

        # Convert regions to graph.
        haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

        # Expected results.
        # No nodes or edges should be removed.
        expected_node_set = set(haplotig_graph.nodes())
        expected_edge_set = set(haplotig_graph.edges())

        return haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set

    # Create test data.
    haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set = create_test()

    # Run the unit under test.
    updated_haplotig_graph = mod.update_haplotig_graph(haplotig_graph, phase_alias_map)

    # Evaluate
    evaluate_update_haplotig_graph(updated_haplotig_graph, expected_node_set, expected_edge_set)

def test_update_haplotig_graph_3():
    """
    This method should only remove edges, and retain all nodes.
    Test on a case with a:
        1. Linear region
        2. Diploid region
        3. Linear region
    The phase alias map in this case is empty, so no edges should be removed.
    """

    def create_test():
        ctg_id = '000000F'
        phase_alias_map = {}

        # Compile a list of regions.
        regions = []
        # Create the first linear region.
        new_region, new_region_header = make_dummy_linear_region(ctg_id, 'ACTG', 0, 1000)
        regions.append(new_region)
        # Create the diploid region.
        new_region, new_header_1, new_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 1000, 2000, 0)
        regions.append(new_region)
        # Create the second linear region.
        new_region, new_region_header = make_dummy_linear_region(ctg_id, 'ACTG', 2000, 3000)
        regions.append(new_region)

        # Convert regions to graph.
        haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

        # Expected results.
        # No nodes or edges should be removed.
        expected_node_set = set(haplotig_graph.nodes())
        expected_edge_set = set(haplotig_graph.edges())

        return haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set

    # Create test data.
    haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set = create_test()

    # Run the unit under test.
    updated_haplotig_graph = mod.update_haplotig_graph(haplotig_graph, phase_alias_map)

    # Evaluate
    evaluate_update_haplotig_graph(updated_haplotig_graph, expected_node_set, expected_edge_set)

def test_update_haplotig_graph_4():
    """
    This method should only remove edges, and retain all nodes.
    Test on a case with a:
        1. Linear region
        2. Diploid region
        3. Linear region
    The phase alias map in this case contains phases not present in the graph. It should have no impact.
    """

    def create_test():
        ctg_id = '000000F'
        phase_alias_map = { '000123F_0_0': 0,
                            '000123F_0_1': 1,}

        # Compile a list of regions.
        regions = []
        # Create the first linear region.
        new_region, new_region_header = make_dummy_linear_region(ctg_id, 'ACTG', 0, 1000)
        regions.append(new_region)
        # Create the diploid region.
        new_region, new_header_1, new_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 1000, 2000, 0)
        regions.append(new_region)
        # Create the second linear region.
        new_region, new_region_header = make_dummy_linear_region(ctg_id, 'ACTG', 2000, 3000)
        regions.append(new_region)

        # Convert regions to graph.
        haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

        # Expected results.
        # No nodes or edges should be removed.
        expected_node_set = set(haplotig_graph.nodes())
        expected_edge_set = set(haplotig_graph.edges())

        return haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set

    # Create test data.
    haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set = create_test()

    # Run the unit under test.
    updated_haplotig_graph = mod.update_haplotig_graph(haplotig_graph, phase_alias_map)

    # Evaluate
    evaluate_update_haplotig_graph(updated_haplotig_graph, expected_node_set, expected_edge_set)


def test_update_haplotig_graph_5():
    """
    This method should only remove edges, and retain all nodes.
    Test on a case with a:
        1. Linear region (unphased, phasing block is -1).
        2. Diploid region (phasing block 0)
        3. Diploid region (phasing block 1)
        4. Diploid region (phasing block 2)
        5. Linear region  (unphased, phasing block is -1).
    Non-empty phase_alias_map. Edges between nodes that belong to different phase aliases
    should be removed.
    """

    def create_test():
        ctg_id = '000000F'

        # Compile a list of regions.
        regions = []
        # Create the first linear region.
        new_region, new_region_header = make_dummy_linear_region(ctg_id, 'ACTG', 0, 1000)
        regions.append(new_region)
        # Create the diploid region.
        dip_0_region, dip_0_header_1, dip_0_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 1000, 2000, 0)
        regions.append(dip_0_region)
        # Create the diploid region.
        dip_1_region, dip_1_header_1, dip_1_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 2000, 3000, 1)
        regions.append(dip_1_region)
        # Create the diploid region.
        dip_2_region, dip_2_header_1, dip_2_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 3000, 4000, 2)
        regions.append(dip_2_region)
        # Create the second linear region.
        new_region, new_region_header = make_dummy_linear_region(ctg_id, 'ACTG', 4000, 5000)
        regions.append(new_region)

        # Convert regions to graph.
        haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

        # Declare that the phases '{ctg_id}_0_0' and '{ctg_id}_1_1' are the same (they have the same alias).
        # Also, '{ctg_id}_0_1' is the same as '{ctg_id}_1_0'.
        phase_alias_map = { '{ctg_id}_0_0'.format(ctg_id=ctg_id): 0,    # dip_0_header_1
                            '{ctg_id}_0_1'.format(ctg_id=ctg_id): 1,    # dip_0_header_2
                            '{ctg_id}_1_0'.format(ctg_id=ctg_id): 1,    # dip_1_header_1
                            '{ctg_id}_1_1'.format(ctg_id=ctg_id): 0,    # dip_1_header_2
                        }

        # Expected results.
        # No nodes should be removed.
        expected_node_set = set(haplotig_graph.nodes())

        # Two edges should be removed.
        expected_edge_set = set(haplotig_graph.edges())
        expected_edge_set.remove((dip_0_header_1, dip_1_header_1))
        expected_edge_set.remove((dip_0_header_2, dip_1_header_2))

        return haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set

    # Create test data.
    haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set = create_test()

    # Run the unit under test.
    updated_haplotig_graph = mod.update_haplotig_graph(haplotig_graph, phase_alias_map)

    # Evaluate
    evaluate_update_haplotig_graph(updated_haplotig_graph, expected_node_set, expected_edge_set)

def test_update_haplotig_graph_6():
    """
    This method should only remove edges, and retain all nodes.
    Test on a case with a:
        1. Linear region (unphased, phasing block is -1).
        2. Diploid region (phasing block 0)
        3. Diploid region (phasing block 1)
        4. Diploid region (phasing block 2)
        5. Linear region  (unphased, phasing block is -1).
    Non-empty phase_alias_map, but the aliases are for distant bubbles, so no neighboring edges.
    No edges should be removed here.
    """

    def create_test():
        ctg_id = '000000F'

        # Compile a list of regions.
        regions = []
        # Create the first linear region.
        new_region, new_region_header = make_dummy_linear_region(ctg_id, 'ACTG', 0, 1000)
        regions.append(new_region)
        # Create the diploid region.
        dip_0_region, dip_0_header_1, dip_0_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 1000, 2000, 0)
        regions.append(dip_0_region)
        # Create the diploid region.
        dip_1_region, dip_1_header_1, dip_1_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 2000, 3000, 1)
        regions.append(dip_1_region)
        # Create the diploid region.
        dip_2_region, dip_2_header_1, dip_2_header_2 = make_dummy_diploid_region(ctg_id, 'ACTG', 'TT', 3000, 4000, 2)
        regions.append(dip_2_region)
        # Create the second linear region.
        new_region, new_region_header = make_dummy_linear_region(ctg_id, 'ACTG', 4000, 5000)
        regions.append(new_region)

        # Convert regions to graph.
        haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

        # Declare that the phases '{ctg_id}_0_0' and '{ctg_id}_1_1' are the same (they have the same alias).
        # Also, '{ctg_id}_0_1' is the same as '{ctg_id}_1_0'.
        phase_alias_map = { '{ctg_id}_0_0'.format(ctg_id=ctg_id): 0,    # dip_0_header_1
                            '{ctg_id}_0_1'.format(ctg_id=ctg_id): 1,    # dip_0_header_2
                            '{ctg_id}_2_0'.format(ctg_id=ctg_id): 1,    # dip_2_header_1
                            '{ctg_id}_2_1'.format(ctg_id=ctg_id): 0,    # dip_2_header_2
                        }

        # Expected results.
        # No nodes or edges should be removed.
        expected_node_set = set(haplotig_graph.nodes())

        expected_edge_set = set(haplotig_graph.edges())

        return haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set

    # Create test data.
    haplotig_graph, phase_alias_map, expected_node_set, expected_edge_set = create_test()

    # Run the unit under test.
    updated_haplotig_graph = mod.update_haplotig_graph(haplotig_graph, phase_alias_map)

    # Evaluate.
    evaluate_update_haplotig_graph(updated_haplotig_graph, expected_node_set, expected_edge_set)

def test_extract_and_write_all_ctg_1(tmpdir):
    """
    A test case of only one linear region without any additional haplotigs.
    """

    def create_test(ctg_id):
        # Make a list of regions (DAG).
        regions = []

        # Linear region. Nothing is unzipped.
        region_1_seq = 'ACTG'
        region_1, region_1_name = make_dummy_linear_region(ctg_id, region_1_seq, 1000, 2000)
        regions.append(region_1)

        # Create a haplotig graph (NetworkX object) from the regions.
        haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

        # Expected results.
        p_ctg_fasta = '>{ctg_id}\n{seq}\n'.format(ctg_id=ctg_id, seq=region_1_seq)
        p_ctg_edges = '\n'.join([' '.join(edge[0:3]) + ' N H -1 0 -1 0' for edge in region_1[5][region_1_name]['path']]) + '\n'

        expected = {}
        expected['p_ctg.{ctg_id}.fa'.format(ctg_id=ctg_id)] = p_ctg_fasta
        expected['p_ctg_edges.{ctg_id}'.format(ctg_id=ctg_id)] = p_ctg_edges
        expected['h_ctg.{ctg_id}.fa'.format(ctg_id=ctg_id)] = ''
        expected['h_ctg_edges.{ctg_id}'.format(ctg_id=ctg_id)] = ''
        expected['h_ctg.{ctg_id}.paf'.format(ctg_id=ctg_id)] = ''
        # expected['regions.fasta'] = '>{header}\n{seq}\n'.format(header=region_1_name, seq=region_1_seq)
        # expected['regions.paf'] = regions_paf

        return haplotig_graph, expected

    ctg_id = '000000F'
    allow_multiple_primaries = False

    # Make the test inputs and expected outputs.
    haplotig_graph, expected = create_test(ctg_id)

    # Run unit under test.
    mod.extract_and_write_all_ctg(ctg_id, haplotig_graph, str(tmpdir), allow_multiple_primaries, mock_fp_proto_log)

    # Evaluate.
    evaluate_extract_and_write_all_ctg(tmpdir, expected)

def test_extract_and_write_all_ctg_2(tmpdir):
    """
    A test case of only one *diploid* region.
    This test should successfully generate one primary and one fully phased haplotig.
    """

    def create_test(ctg_id):
        # Make a list of regions (DAG).
        regions = []

        # Diploid region. The entire contig is unzipped.
        region_1_seq_1 = 'ACTG'
        region_1_seq_2 = 'AT'
        region_1, region_1_name_1, region_1_name_2 = make_dummy_diploid_region(ctg_id, region_1_seq_1, region_1_seq_2, 1000, 2000, 0)
        regions.append(region_1)

        # Create a haplotig graph (NetworkX object) from the regions.
        haplotig_graph = mod.regions_to_haplotig_graph(ctg_id, regions, mock_fp_proto_log)

        # Expected results.
        expected = {}
        # Expected results.
        ctg_seq = region_1_seq_1
        p_ctg_fasta = '>{ctg_id}\n{seq}\n'.format(ctg_id=ctg_id, seq=ctg_seq)
        p_ctg_edges = '\n'.join([' '.join(edge[0:3]) + ' N H 0 0 0 0' for edge in region_1[5][region_1_name_1]['path']]) + '\n'

        h_ctg_id = '{ctg_id}_001'.format(ctg_id=ctg_id)
        h_ctg_seq = region_1_seq_2
        h_ctg_fasta = ''
        h_ctg_paf = ''
        h_ctg_fasta += '>{ctg_id}\n{seq}\n'.format(ctg_id=h_ctg_id, seq=h_ctg_seq)
        h_ctg_paf += '{h_ctg_id}\t{qlen}\t{qstart}\t{qend}\t+\t{ctg_id}\t{tlen}\t{tstart}\t{tend}\t{tspan}\t{tspan}\t{mapq}\n'.format(
                        h_ctg_id=h_ctg_id, qlen=len(h_ctg_seq), qstart=0, qend=len(h_ctg_seq),
                        ctg_id=ctg_id, tlen=len(ctg_seq), tstart=0, tend=len(ctg_seq), tspan=len(ctg_seq), mapq=60)
        h_ctg_edges = '\n'.join([h_ctg_id + ' ' + ' '.join(edge[1:3]) + ' N H 0 1 0 1' for edge in region_1[5][region_1_name_2]['path']]) + '\n'

        expected = {}
        expected['p_ctg.{ctg_id}.fa'.format(ctg_id=ctg_id)] = p_ctg_fasta
        expected['p_ctg_edges.{ctg_id}'.format(ctg_id=ctg_id)] = p_ctg_edges
        expected['h_ctg.{ctg_id}.fa'.format(ctg_id=ctg_id)] = h_ctg_fasta
        expected['h_ctg_edges.{ctg_id}'.format(ctg_id=ctg_id)] = h_ctg_edges
        expected['h_ctg.{ctg_id}.paf'.format(ctg_id=ctg_id)] = h_ctg_paf

        return haplotig_graph, expected

    ctg_id = '000000F'
    allow_multiple_primaries = False

    # Make the test inputs and expected outputs.
    haplotig_graph, expected = create_test(ctg_id)

    # Run unit under test.
    mod.extract_and_write_all_ctg(ctg_id, haplotig_graph, str(tmpdir), allow_multiple_primaries, mock_fp_proto_log)

    # Evaluate.
    evaluate_extract_and_write_all_ctg(tmpdir, expected)

def test_extract_and_write_all_ctg_3(tmpdir):
    """
    A degenerate case which should throw.
    This test case should not occur in practice,
    but we need to check for it still.
    This tests for two graph components in the same haplotig graph, which are not connected with
    a source and a sink node.
    """

    def create_test(ctg_id):
        # Make a list of regions (DAG).
        # Diploid region. The entire contig is unzipped.
        regions_1 = []
        region_1_seq_1 = 'ACTG'
        region_1_seq_2 = 'AT'
        region_1, region_1_name_1, region_1_name_2 = make_dummy_diploid_region(ctg_id + 'p0', region_1_seq_1, region_1_seq_2, 1000, 2000, 0)
        regions_1.append(region_1)

        # Another diploid region. Could be linear, or any kind.
        regions_2 = []
        region_2_seq_1 = 'ACTG'
        region_2_seq_2 = 'AT'
        region_2, region_2_name_1, region_2_name_2 = make_dummy_diploid_region(ctg_id + 'p1', region_2_seq_1, region_2_seq_2, 4000, 5000, 0)
        regions_2.append(region_2)

        # Create haplotig graphs (NetworkX object) from the regions.
        haplotig_graph_1 = mod.regions_to_haplotig_graph(ctg_id + 'p0', regions_1, mock_fp_proto_log)
        haplotig_graph_2 = mod.regions_to_haplotig_graph(ctg_id + 'p1', regions_2, mock_fp_proto_log)

        # Make a degenerate graph.
        haplotig_graph = nx.compose(haplotig_graph_1, haplotig_graph_2)

        return haplotig_graph

    ctg_id = '000000F'
    allow_multiple_primaries = False

    # Make the test inputs and expected outputs.
    haplotig_graph = create_test(ctg_id)

    # Run unit under test.
    with pytest.raises(Exception, match=r'Skipping additional subgraphs of the primary contig.*'):
        mod.extract_and_write_all_ctg(ctg_id, haplotig_graph, str(tmpdir), allow_multiple_primaries, mock_fp_proto_log)

def test_extract_and_write_all_ctg_4(tmpdir):
    """
    A degenerate case which should throw.
    This test should allow multiple primary componets to be extracted from the same graph.
    """

    def create_test(ctg_id):
        # Linear region.
        regions = []
        region_1_seq = 'AAAA'
        region_1, region_1_name = make_dummy_linear_region(ctg_id + 'a', region_1_seq, 1000, 2000)
        regions.append(region_1)
        haplotig_graph_1 = mod.regions_to_haplotig_graph(ctg_id + 'a', regions, mock_fp_proto_log)

        # Another separate linear region.
        regions = []
        region_2_seq = 'TTTT'
        region_2, region_2_name = make_dummy_linear_region(ctg_id + 'b', region_2_seq, 4000, 5000)
        regions.append(region_2)
        haplotig_graph_2 = mod.regions_to_haplotig_graph(ctg_id + 'b', regions, mock_fp_proto_log)

        # Make a multi-component graph.
        haplotig_graph = nx.compose(haplotig_graph_1, haplotig_graph_2)

        # Expected results.
        expected = {}
        # Expected results.
        p_ctg_fasta =   '>{ctg_id}p01\n{seq}\n'.format(ctg_id=ctg_id, seq=region_2_seq) + \
                        '>{ctg_id}p02\n{seq}\n'.format(ctg_id=ctg_id, seq=region_1_seq)
        p_ctg_edges =   '\n'.join(['{ctg_id}p01 '.format(ctg_id=ctg_id) + ' '.join(edge[1:3]) + ' N H -1 0 -1 0' for edge in haplotig_graph.node[region_2_name]['htig']['path']]) + '\n' + \
                        '\n'.join(['{ctg_id}p02 '.format(ctg_id=ctg_id) + ' '.join(edge[1:3]) + ' N H -1 0 -1 0' for edge in haplotig_graph.node[region_1_name]['htig']['path']]) + '\n'

        h_ctg_fasta =   ''
        h_ctg_edges =   ''
        h_ctg_paf = ''

        expected = {}
        expected['p_ctg.{ctg_id}.fa'.format(ctg_id=ctg_id)] = p_ctg_fasta
        expected['p_ctg_edges.{ctg_id}'.format(ctg_id=ctg_id)] = p_ctg_edges
        expected['h_ctg.{ctg_id}.fa'.format(ctg_id=ctg_id)] = h_ctg_fasta
        expected['h_ctg_edges.{ctg_id}'.format(ctg_id=ctg_id)] = h_ctg_edges
        expected['h_ctg.{ctg_id}.paf'.format(ctg_id=ctg_id)] = h_ctg_paf

        return haplotig_graph, expected

    ctg_id = '000000F'

    # Make the test inputs and expected outputs.
    haplotig_graph, expected = create_test(ctg_id)

    # Run unit under test. Test for failure, because there are multiple components.
    allow_multiple_primaries = False
    with pytest.raises(Exception, match=r'Skipping additional subgraphs of the primary contig.*'):
        mod.extract_and_write_all_ctg(ctg_id, haplotig_graph, str(tmpdir), allow_multiple_primaries, mock_fp_proto_log)

    # Run unit under test. Allow multiple components, and test if it will actually generate valid output.
    allow_multiple_primaries = True
    mod.extract_and_write_all_ctg(ctg_id, haplotig_graph, str(tmpdir), allow_multiple_primaries, mock_fp_proto_log)

    # Evaluate.
    evaluate_extract_and_write_all_ctg(tmpdir, expected)

def test_extract_and_write_all_ctg_5(tmpdir):
    """
    A degenerate case which should throw.
    This test should allow multiple primary componets to be extracted from the same graph.
    """

    def create_test(ctg_id):
        # Linear region.
        regions = []
        region_1_seq = 'AAAA'
        region_1, region_1_name = make_dummy_linear_region(ctg_id + 'a', region_1_seq, 1000, 2000)
        regions.append(region_1)
        haplotig_graph_1 = mod.regions_to_haplotig_graph(ctg_id + 'a', regions, mock_fp_proto_log)

        # Another separate diploid region.
        regions = []
        region_2_seq_1 = 'TTTT'
        region_2_seq_2 = 'GG'
        region_2, region_2_name_1, region_2_name_2 = make_dummy_diploid_region(ctg_id + 'b', region_2_seq_1, region_2_seq_2, 4000, 5000, 1)
        regions.append(region_2)
        haplotig_graph_2 = mod.regions_to_haplotig_graph(ctg_id + 'b', regions, mock_fp_proto_log)

        # Make a multi-component graph.
        haplotig_graph = nx.compose(haplotig_graph_1, haplotig_graph_2)

        # Expected results.
        expected = {}
        # Expected results.
        # The order in which primaries are listed is abritrary and depends on NetworkX.
        # The selection of path through diploid bubbles is arbitrary and depends on NetworkX.
        p_ctg_fasta =   '>{ctg_id}p01\n{seq}\n'.format(ctg_id=ctg_id, seq=region_2_seq_2) + \
                        '>{ctg_id}p02\n{seq}\n'.format(ctg_id=ctg_id, seq=region_1_seq)
        p_ctg_edges =   '\n'.join(['{ctg_id}p01 '.format(ctg_id=ctg_id) + ' '.join(edge[1:3]) + ' N H 1 1 1 1' for edge in haplotig_graph.node[region_2_name_2]['htig']['path']]) + '\n' + \
                        '\n'.join(['{ctg_id}p02 '.format(ctg_id=ctg_id) + ' '.join(edge[1:3]) + ' N H -1 0 -1 0' for edge in haplotig_graph.node[region_1_name]['htig']['path']]) + '\n'

        h_ctg_fasta =   '>{ctg_id}p01_001\n{seq}\n'.format(ctg_id=ctg_id, seq=region_2_seq_1)
        h_ctg_edges =   '\n'.join(['{ctg_id}p01_001 '.format(ctg_id=ctg_id) + ' '.join(edge[1:3]) + ' N H 1 0 1 0' for edge in haplotig_graph.node[region_2_name_1]['htig']['path']]) + '\n'
        h_ctg_paf =     '{ctg_id}p01_001\t{qlen}\t{qstart}\t{qend}\t+\t{ctg_id}p01\t{tlen}\t{tstart}\t{tend}\t{tspan}\t{tspan}\t{mapq}\n'.format(
                        ctg_id=ctg_id,
                        qlen=len(region_2_seq_1), qstart=0, qend=len(region_2_seq_1),
                        tlen=len(region_2_seq_2), tstart=0, tend=len(region_2_seq_2), tspan=len(region_2_seq_2), mapq=60)

        expected = {}
        expected['p_ctg.{ctg_id}.fa'.format(ctg_id=ctg_id)] = p_ctg_fasta
        expected['p_ctg_edges.{ctg_id}'.format(ctg_id=ctg_id)] = p_ctg_edges
        expected['h_ctg.{ctg_id}.fa'.format(ctg_id=ctg_id)] = h_ctg_fasta
        expected['h_ctg_edges.{ctg_id}'.format(ctg_id=ctg_id)] = h_ctg_edges
        expected['h_ctg.{ctg_id}.paf'.format(ctg_id=ctg_id)] = h_ctg_paf

        return haplotig_graph, expected

    ctg_id = '000000F'

    # Make the test inputs and expected outputs.
    haplotig_graph, expected = create_test(ctg_id)

    # Run unit under test. Test for failure, because there are multiple components.
    allow_multiple_primaries = False
    with pytest.raises(Exception, match=r'Skipping additional subgraphs of the primary contig.*'):
        mod.extract_and_write_all_ctg(ctg_id, haplotig_graph, str(tmpdir), allow_multiple_primaries, mock_fp_proto_log)

    # Run unit under test. Allow multiple components, and test if it will actually generate valid output.
    allow_multiple_primaries = True
    mod.extract_and_write_all_ctg(ctg_id, haplotig_graph, str(tmpdir), allow_multiple_primaries, mock_fp_proto_log)

    # Evaluate.
    evaluate_extract_and_write_all_ctg(tmpdir, expected)
