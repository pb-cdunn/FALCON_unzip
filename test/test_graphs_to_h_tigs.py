import falcon_unzip.mains.graphs_to_h_tigs as mod
import helpers
import pytest
import os
import networkx as nx

def add_node(sg, v, vphase):
    if v not in sg.nodes():
        sg.add_node(v, label="%d_%d" % vphase, phase="%d_%d" % vphase, src="P")

def add_edge(sg, p_asm_graph, h_asm_graph, v, vphase, w, wphase, src, cross_phase, overlap_len):
    add_node(sg, v, vphase)
    add_node(sg, w, wphase)
    sg.add_edge(v, w, src=src, cross_phase=cross_phase)
    rv = mod.reverse_end(v)
    rw = mod.reverse_end(w)
    add_node(sg, rv, vphase)
    add_node(sg, rw, wphase)
    sg.add_edge(rw, rv, src=src, cross_phase=cross_phase)
    # The structure of an sg_edges object looks like this:
    #   p_asm_graph.sg_edges[(v, w)] = ((seq_id, b, e), score, idt, type_)
    # We don't care about the actual values, just that abs(b - e) == overlap_len.
    # Mirror all edges in both graphs for simpler testing.
    p_asm_graph.sg_edges[(v, w)] = ((0, 0, overlap_len), 10000, 100.0, 'G')
    h_asm_graph.sg_edges[(v, w)] = ((0, 0, overlap_len), 10000, 100.0, 'G')

class MockAsmGraph:
    def __init__(self):
        self.sg_edges = {}

def node_name(node_id):
    return '%04d:B' % (node_id)

def test_find_phased_neighbors_1():
    """
    Test on an empty input.
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()
   
    s_node, s_phase = None, None
    t_node, t_phase = None, None
    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    expected = (None, None)

    assert(result == expected)

def test_find_phased_neighbors_2():
    """
    Test on a case where there is nothing to be done.
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()

    # Just a linear contig, all nodes in phase
    for i in xrange(1, 10):
        add_edge(sg, p_asm_graph, h_asm_graph, node_name(i), (1, 0), node_name(i + 1), (1, 0), 'OP', 'N', 1000)

    # Setup our test source and sink nodes.
    s_node = node_name(1)
    t_node = node_name(10)
    
    # This is how phases are extracted in Unzip.
    s_phase = sg.node[s_node]["phase"].split("_")
    t_phase = sg.node[t_node]["phase"].split("_")

    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    expected = (s_node, t_node)

    assert(result == expected)

def test_find_phased_neighbors_3():
    """
    Test a simple cross-phase case where the first
    node is not in phase with the rest of the contig.
    This tests the lookup on the out edge of the wrong-phase
    source node.
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()

    # Just a linear contig, all nodes in phase
    for i in xrange(1, 10):
        add_edge(sg, p_asm_graph, h_asm_graph, node_name(i), (1, 0), node_name(i + 1), (1, 0), 'OP', 'N', 1000)
    
    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(11), (1, 1), node_name(1), (1, 0), 'OP', 'Y', 1000)

    # Setup our test source and sink nodes.
    s_node = '%04d:B' % (11)
    t_node = '%04d:B' % (10)

    # This is how phases are extracted in Unzip.
    s_phase = sg.node[s_node]["phase"].split("_")
    t_phase = sg.node[t_node]["phase"].split("_")

    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    # The function call should move the source node down the out-edge to node_name(1).
    expected = (node_name(1), t_node)

    assert(result == expected)

def test_find_phased_neighbors_4():
    """
    Test a cross-phase case where the first
    node is out of phase with the rest of the contig,
    and there are two alternative out-edges to choose from
    (the one with larger overlap len should prevail, and
    the correct result should wind up in node_name(1)).
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()

    # Just a linear contig, all nodes in phase
    for i in xrange(1, 10):
        add_edge(sg, p_asm_graph, h_asm_graph, node_name(i), (1, 0), node_name(i + 1), (1, 0), 'OP', 'N', 1000)

    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(11), (1, 1), node_name(2), (1, 0), 'OP', 'Y', 700)
    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(11), (1, 1), node_name(1), (1, 0), 'OP', 'Y', 900)

    # Setup our test source and sink nodes.
    s_node = node_name(11)
    t_node = node_name(10)

    # This is how phases are extracted in Unzip.
    s_phase = sg.node[s_node]["phase"].split("_")
    t_phase = sg.node[t_node]["phase"].split("_")

    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    expected = (node_name(1), t_node)

    assert(result == expected)

def test_find_phased_neighbors_5():
    """
    This should look at the in-edges
    of the source node to find the candidate.
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()

    # Just a linear contig, all nodes in phase
    for i in xrange(1, 10):
        add_edge(sg, p_asm_graph, h_asm_graph, node_name(i), (1, 0), node_name(i + 1), (1, 0), 'OP', 'N', 1000)

    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(1), (1, 0), node_name(11), (1, 1), 'OP', 'Y', 1000)

    # Setup our test source and sink nodes.
    s_node = node_name(11)
    t_node = node_name(10)

    # This is how phases are extracted in Unzip.
    s_phase = sg.node[s_node]["phase"].split("_")
    t_phase = sg.node[t_node]["phase"].split("_")

    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    expected = (node_name(1), t_node)

    assert(result == expected)

def test_find_phased_neighbors_6():
    """
    This should look at the in-edges
    of the source node to find the candidate.
    There are two edges to choose from, and
    the better one points to node_name(2).
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()

    # Just a linear contig, all nodes in phase
    for i in xrange(1, 10):
        add_edge(sg, p_asm_graph, h_asm_graph, node_name(i), (1, 0), node_name(i + 1), (1, 0), 'OP', 'N', 1000)

    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(1), (1, 0), node_name(11), (1, 1), 'OP', 'Y', 700)
    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(2), (1, 0), node_name(11), (1, 1), 'OP', 'Y', 900)

    # Setup our test source and sink nodes.
    s_node = node_name(11)
    t_node = node_name(10)

    # This is how phases are extracted in Unzip.
    s_phase = sg.node[s_node]["phase"].split("_")
    t_phase = sg.node[t_node]["phase"].split("_")

    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    expected = (node_name(2), t_node)

    assert(result == expected)

def test_find_phased_neighbors_7():
    """
    Test a simple cross-phase case where the last
    node is not in phase with the rest of the contig.
    This tests the lookup on the out edge of the wrong-phase
    node.
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()

    # Just a linear contig, all nodes in phase
    for i in xrange(1, 10):
        add_edge(sg, p_asm_graph, h_asm_graph, node_name(i), (1, 0), node_name(i + 1), (1, 0), 'OP', 'N', 1000)

    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(11), (1, 1), node_name(10), (1, 0), 'OP', 'Y', 1000)

    # Setup our test source and sink nodes.
    s_node = '%04d:B' % (1)
    t_node = '%04d:B' % (11)

    # This is how phases are extracted in Unzip.
    s_phase = sg.node[s_node]["phase"].split("_")
    t_phase = sg.node[t_node]["phase"].split("_")

    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    # The function call should move the source node down the out-edge to node_name(1).
    expected = (s_node, node_name(10))

    assert(result == expected)

def test_find_phased_neighbors_8():
    """
    Test a cross-phase case where the last
    node is out of phase with the rest of the contig,
    and there are two alternative out-edges to choose from
    (the one with larger overlap len should prevail, and
    the correct result should wind up in node_name(10)).
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()

    # Just a linear contig, all nodes in phase
    for i in xrange(1, 10):
        add_edge(sg, p_asm_graph, h_asm_graph, node_name(i), (1, 0), node_name(i + 1), (1, 0), 'OP', 'N', 1000)

    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(11), (1, 1), node_name(9), (1, 0), 'OP', 'Y', 700)
    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(11), (1, 1), node_name(10), (1, 0), 'OP', 'Y', 900)

    # Setup our test source and sink nodes.
    s_node = node_name(1)
    t_node = node_name(11)

    # This is how phases are extracted in Unzip.
    s_phase = sg.node[s_node]["phase"].split("_")
    t_phase = sg.node[t_node]["phase"].split("_")

    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    expected = (s_node, node_name(10))

    assert(result == expected)

def test_find_phased_neighbors_9():
    """
    This should look at the in-edges
    of the sink node to find the candidate.
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()

    # Just a linear contig, all nodes in phase
    for i in xrange(1, 10):
        add_edge(sg, p_asm_graph, h_asm_graph, node_name(i), (1, 0), node_name(i + 1), (1, 0), 'OP', 'N', 1000)

    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(10), (1, 0), node_name(11), (1, 1), 'OP', 'Y', 1000)

    # Setup our test source and sink nodes.
    s_node = node_name(1)
    t_node = node_name(11)

    # This is how phases are extracted in Unzip.
    s_phase = sg.node[s_node]["phase"].split("_")
    t_phase = sg.node[t_node]["phase"].split("_")

    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    expected = (s_node, node_name(10))

    assert(result == expected)

def test_find_phased_neighbors_10():
    """
    This should look at the in-edges
    of the sink node to find the candidate.
    There are two edges to choose from, and
    the better one points to node_name(2).
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()

    # Just a linear contig, all nodes in phase
    for i in xrange(1, 10):
        add_edge(sg, p_asm_graph, h_asm_graph, node_name(i), (1, 0), node_name(i + 1), (1, 0), 'OP', 'N', 1000)

    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(10), (1, 0), node_name(11), (1, 1), 'OP', 'Y', 700)
    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(9), (1, 0), node_name(11), (1, 1), 'OP', 'Y', 900)

    # Setup our test source and sink nodes.
    s_node = node_name(1)
    t_node = node_name(11)

    # This is how phases are extracted in Unzip.
    s_phase = sg.node[s_node]["phase"].split("_")
    t_phase = sg.node[t_node]["phase"].split("_")

    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    expected = (s_node, node_name(9))

    assert(result == expected)

def test_find_phased_neighbors_11():
    """
    This is identical to test 10, but it builds
    everything in the 'H' graph. This test
    checks if that branch of calc_overlap_len works
    correctly.
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()

    # Just a linear contig, all nodes in phase
    for i in xrange(1, 10):
        add_edge(sg, p_asm_graph, h_asm_graph, node_name(i), (1, 0), node_name(i + 1), (1, 0), 'H', 'N', 1000)

    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(10), (1, 0), node_name(11), (1, 1), 'H', 'Y', 700)
    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(9), (1, 0), node_name(11), (1, 1), 'H', 'Y', 900)

    # Setup our test source and sink nodes.
    s_node = node_name(1)
    t_node = node_name(11)

    # This is how phases are extracted in Unzip.
    s_phase = sg.node[s_node]["phase"].split("_")
    t_phase = sg.node[t_node]["phase"].split("_")

    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    expected = (s_node, node_name(9))

    assert(result == expected)

def test_find_phased_neighbors_12():
    """
    Tests multiple equal branches.
    """
    sg = nx.DiGraph()
    p_asm_graph = MockAsmGraph()
    h_asm_graph = MockAsmGraph()

    # Just a linear contig, all nodes in phase
    for i in xrange(1, 10):
        add_edge(sg, p_asm_graph, h_asm_graph, node_name(i), (1, 0), node_name(i + 1), (1, 0), 'OP', 'N', 1000)

    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(10), (1, 0), node_name(11), (1, 1), 'OP', 'Y', 0)
    # Add a cross-phase edge.
    add_edge(sg, p_asm_graph, h_asm_graph, node_name(9), (1, 0), node_name(11), (1, 1), 'OP', 'Y', 0)

    # Setup our test source and sink nodes.
    s_node = node_name(1)
    t_node = node_name(11)

    # This is how phases are extracted in Unzip.
    s_phase = sg.node[s_node]["phase"].split("_")
    t_phase = sg.node[t_node]["phase"].split("_")

    result = mod.find_phased_neighbors(sg, p_asm_graph, h_asm_graph, s_node, t_node, s_phase, t_phase)
    expected1 = (s_node, node_name(9))
    expected2 = (s_node, node_name(10))

    assert(result == expected1 or result == expected2)

