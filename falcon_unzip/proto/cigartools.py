import pysam
import copy

CIGAR_OPS = 'MIDNSHP=X'
CIGAR_IS_MATCH = [True, False, False, False, False, False, False, True, True]
CONSUMES_REF = [True, False, True, True, False, False, False, True, True]
CONSUMES_QUERY = [True, True, False, False, True, False, False, True, True]
CIGAR_OP_M = 0
CIGAR_OP_I = 1
CIGAR_OP_D = 2
CIGAR_OP_N = 3
CIGAR_OP_S = 4
CIGAR_OP_H = 5
CIGAR_OP_P = 6
CIGAR_OP_EQ = 7
CIGAR_OP_X = 8
CIGAR_CHAR_TO_OP = {'M': CIGAR_OP_M, 'I': CIGAR_OP_I, 'D': CIGAR_OP_D, 'N': CIGAR_OP_N, 'S': CIGAR_OP_S, 'H': CIGAR_OP_H, 'P': CIGAR_OP_P, '=': CIGAR_OP_EQ, 'X': CIGAR_OP_X}

def cigar_to_aln_array(cigar):
    """
    Converts a Pysam CIGAR list into an array of CIGAR ops, where each element
    is a single CIGAR operation of count 1.
    """
    ret = []
    for op, count in cigar:
        ret += [op for i in xrange(count)]
    return ret

def aln_array_to_cigar(aln_array):
    """
    Converts the alignment array back to a Pysam CIGAR list.
    """
    cigar = []
    if len(aln_array) == 0:
        return cigar
    count = 0
    prev_op = aln_array[0]
    for i in xrange(len(aln_array)):
        op = aln_array[i]
        if op == prev_op:
            count += 1
        else:
            cigar.append((prev_op, count))
            count = 1
        prev_op = op
    if count > 0:
        if len(cigar) > 0 and cigar[-1][0] == prev_op:
            cigar[-1] = (cigar[-1][0], cigar[-1][1] + count)    # pragma: no cover
        else:
            cigar.append((prev_op, count))
    return cigar

def find_positions(aln_array, reference_start):
    """
    Finds the query bases which align to the specified
    reference positions.
    """
    ret = {}

    ref_pos = reference_start + 0
    query_pos = 0

    for op_id in xrange(len(aln_array)):
        op = aln_array[op_id]
        if op != CIGAR_OP_S:
            op_start = op_id
            break

    # op_start = 0
    # while aln_array[op_start] == CIGAR_OP_H:
    #     op_start += 1
    # while aln_array[op_start] == CIGAR_OP_S:
    #     op_start += 1
    #     q_start += 1

    # op_end = len(aln_array) - 1
    # while aln_array[op_end] == CIGAR_OP_S or aln_array[op_end] == CIGAR_OP_H:
    #     op_end -= 1

    for op_id in xrange(len(aln_array)):
    # for op_id in xrange(op_start, op_end + 1):
        op = aln_array[op_id]
        cq = CONSUMES_QUERY[op]
        cr = CONSUMES_REF[op]

        if op == CIGAR_OP_S:
            query_pos += 1 # + (1 if cq else 0)
            continue

        ret.setdefault(ref_pos, [])
        ret[ref_pos].append(query_pos)

        if cr:
            ref_pos += 1
        if cq:
            query_pos += 1

    # ret.setdefault(ref_pos, [])
    # ret[ref_pos].append(query_pos)
    # import sys
    # sys.stderr.write('positions: %s\n' % (str(ret)))

    return ret

def find_split_positions(aln_array, reference_start, reference_split_pos):
    """
    Finds the position of the rightmost match operation which is located
    left of the linearized contig end, and the leftmost match operation
    which is located right of linearized contig end (or at the very end).
    These coordinates can then be used to split the alignment array
    into two alignments: the left part corresponds to the alignment at
    the sequence end, whereas the right part corresponds to the alignment
    at the reference beginning.
    """
    last_left_M_op = None
    op_id = 0
    for op_id in xrange(len(aln_array)):
        op = aln_array[op_id]
        cq = CONSUMES_QUERY[op]
        cr = CONSUMES_REF[op]
        if reference_start >= reference_split_pos:
            break
        if CIGAR_IS_MATCH[aln_array[op_id]]:
            last_left_M_op = op_id
        reference_start = reference_start + (1 if cr else 0)
        # query_pos = query_pos + (1 if cq else 0)

    split_op_id = op_id

    first_right_M_op = None
    for op_id in xrange(split_op_id, len(aln_array)):
        op = aln_array[op_id]
        if CIGAR_IS_MATCH[aln_array[op_id]]:
            first_right_M_op = op_id
            break;

    return [last_left_M_op, first_right_M_op]

def split_aln_array(aln_array, reference_start, last_left_M_op, first_right_M_op):
    """
    Takes the alignment array (each array element is a single CIGAR op),
    the position of the alignment on the reference (reference_start), the
    position of the rightmost match operation which is left of linearized
    contig end, and the leftmost match operation which is right of linearized
    contig end, and produces two alignment arrays corresponding to the left and
    the right parts.
    Both arrays are padded with soft clipping operations to maintain query
    length.
    If the original alignment had hard clippings, those will be preserved.
    """
    # Sanity checks.
    if None in [last_left_M_op, first_right_M_op] or last_left_M_op >= first_right_M_op:
        return [None, None, None]

    # Initialize.
    left_aln_array = aln_array[0:(last_left_M_op + 1)]
    right_aln_array = []

    # Process bases which need clipping.
    for op_id in xrange(0, len(aln_array)):
        op = aln_array[op_id]
        cq = CONSUMES_QUERY[op]
        cr = CONSUMES_REF[op]

        if op_id > last_left_M_op:
            # If op consumes query, push S op instead so that query length is preserved.
            # The only alternative is if there was a H CIGAR op. Push that as well.
            if cq:
                left_aln_array.append(CIGAR_OP_S)
            elif op == CIGAR_OP_H:
                left_aln_array.append(CIGAR_OP_H)

        # Find where the right alignment part begins. We need this
        # to update the reference mapping position.
        if op_id == first_right_M_op:
            right_ref_start = reference_start

        if op_id < first_right_M_op:
            # If op consumes query, push S op instead so that query length is preserved.
            # The only alternative is if there was a H CIGAR op. Push that as well.
            if cq:
                right_aln_array.append(CIGAR_OP_S)
            elif op == CIGAR_OP_H:
                right_aln_array.append(CIGAR_OP_H)

        reference_start = reference_start + (1 if cr else 0)

    # Fill out the rest of the right alignment.
    right_aln_array += aln_array[first_right_M_op:]

    return [left_aln_array, right_aln_array, right_ref_start]

def split_cigar(cigar, reference_start, reference_split_pos):
    # Convert the CIGAR string to a longer alignment array, where each
    # element of the array is a single CIGAR op. This is less efficient
    # but more readable.
    # Efficiency will still not be impacted greatly, because only a small
    # fraction of reads need to be split like this.
    aln_array = cigar_to_aln_array(cigar)

    # Find the position of the last match/mismatch base right before the
    # reference end (linearized reference end), and also right
    # after (or right at) the reference end.
    [last_left_M_op, first_right_M_op] = find_split_positions(aln_array, reference_start, reference_split_pos)

    # Sanity checks.
    if last_left_M_op == None or first_right_M_op == None:
        return [None, None, None]
    if last_left_M_op >= first_right_M_op:
        return [None, None, None]           # pragma: no cover

    # This creates two alignment arrays: one for the part left of the linear
    # reference end with soft clipping CIGAR ops added as padding after
    # the last_left_M_op, and one for the part which starts after the linear
    # reference end, with soft clipping CIGAR ops prepended to the actual
    # alignment for padding. Hard clipping CIGAR ops are preserved, if there
    # were any.
    left_aln_array, right_aln_array, right_ref_start = split_aln_array(aln_array, reference_start, last_left_M_op, first_right_M_op)

    if left_aln_array == None or right_aln_array == None or right_ref_start == None: return [None, None, None]

    # The alignment arrays can now be converted to CIGAR strings.
    cigar_left = aln_array_to_cigar(left_aln_array)
    cigar_right = aln_array_to_cigar(right_aln_array)

    return [cigar_left, cigar_right, right_ref_start]

def split_sam(sam, linear_ref_len):
    reference_start = sam.reference_start + 0           # Pysam reference_start is 0-based.

    # Split the CIGAR string into two pieces.
    cigar_left, cigar_right, right_ref_start = split_cigar(sam.cigar, reference_start, linear_ref_len)

    # Create a SAM entry for the part left of linear reference end.
    sam_left = copy.copy(sam)
    sam_left.cigar = cigar_left

    # Create a SAM entry for the part right of linear reference end.
    # The coordinate needs to be modified as well.
    sam_right = copy.copy(sam)
    sam_right.reference_start = right_ref_start % linear_ref_len
    sam_right.cigar = cigar_right

    # Set the smaller part as supplementary.
    if (linear_ref_len - sam.reference_start) >= (sam.reference_length) / 2:
        sam_right.flag |= 2048
    else:
        sam_left.flag |= 2048

    return [sam_left, sam_right]
