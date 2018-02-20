#! /usr/bin/env python2.7

import re
import sys
import pysam
# import argparse

################################
######### Utility tools ########
import sys;
import os;
from time import gmtime, strftime
import traceback;
import fnmatch;
### Logs messages to STDERR and an output log file if provided (opened elsewhere).
def log(message, fp_log=sys.stderr):
  timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
  if (fp_log != None):
      fp_log.write('[%s] %s\n' % (timestamp, message))
      fp_log.flush();
################################

CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_H = 5
CIGAR_P = 6
CIGAR_EQ = 7
CIGAR_X = 8
CIGAR_B = 9

def is_primary(sam):
    return not (sam.is_secondary or sam.is_supplementary)

def open_sam_bam_for_reading(file_path):
    if (file_path.endswith('bam')):
        fp = pysam.AlignmentFile(file_path, 'rb', check_sq=False)
        return fp, True
    fp = pysam.AlignmentFile(file_path, 'r', check_sq=False)
    return fp, False

def open_sam_bam_for_writing(file_path, header):
    if (file_path.endswith('bam')):
        fp = pysam.AlignmentFile(file_path, 'wb', header=header)
    else:
        fp = pysam.AlignmentFile(file_path, 'wh', header=header)
    return fp

def pysam_to_m4(aln, ref_lens = None, skip_supplementary=True, skip_secondary=True):
    ret = None

    if aln.is_unmapped == True:
        return ret
    if aln.is_supplementary == True and skip_supplementary == True:
        return ret
    if aln.is_secondary == True and skip_secondary == True:
        return ret
    if len(aln.cigar) == 0:
        return ret
    if aln.query_sequence == None or len(aln.query_sequence) == 0:
        return ret

    ref_len = (ref_lens[aln.reference_name]) if ref_lens != None else 0

    m = 0
    matches = 0
    mismatches = 0
    ins = 0
    dels = 0

    for op, count in aln.cigar:

        if op == CIGAR_EQ:
            matches += count
            m += count
        elif op == CIGAR_X:
            mismatches += count
            m += count
        elif op == CIGAR_M:
            # assert(op != CIGAR_M and "Cannot calculate match rate from 'M' operations.'")
            matches += count
            mismatches += count
            m += count
        elif op == CIGAR_I:
            ins += count
        elif op == CIGAR_D:
            dels += count

    score = 0 - (mismatches + ins + dels)
    identity = 100.0 * float(matches) / float(aln.query_alignment_end - aln.query_alignment_start)
    ret = [aln.qname, aln.reference_name, score, identity,
                0, aln.query_alignment_start, aln.query_alignment_end, aln.query_length,
                (1 if aln.is_reverse == True else 0), aln.reference_start, aln.reference_end, ref_len, 254, aln]

    return ret

def sam_to_m4(in_path):
    ret = []

    fp_in, is_bam = open_sam_bam_for_reading(in_path)
    it_fp_in = fp_in.fetch() if not is_bam else fp_in

    ref_lens = {}
    for val in fp_in.header['SQ']:
        ref_lens[val['SN']] = val['LN']

    num_lines = 0
    for sam in it_fp_in:
        num_lines += 1
        if (num_lines % 10000 == 0):
            log('Loaded %d lines.' % ((num_lines - 1)));

        m4 = pysam_to_m4(sam, ref_lens)
        if m4 == None:
            continue
        ret.append(m4)

        # if sam.is_unmapped == True: continue

        # m = 0
        # matches = 0
        # mismatches = 0
        # ins = 0
        # dels = 0
        # for op, count in sam.cigar:
        #     if op == CIGAR_EQ:
        #         matches += count
        #         m += count
        #     elif op == CIGAR_X:
        #         mismatches += count
        #         m += count
        #     elif op == CIGAR_M:
        #         # assert(op != CIGAR_M and "Cannot calculate match rate from 'M' operations.'")
        #         matches += count
        #         mismatches += count
        #         m += count
        #     elif op == CIGAR_I:
        #         ins += count
        #     elif op == CIGAR_D:
        #         dels += count

        # score = 0 - (mismatches + ins + dels)
        # identity = 100.0 * float(matches) / float(sam.query_alignment_end - sam.query_alignment_start)
        # ret.append([sam.qname, sam.reference_name, score, identity,
        #             0, sam.query_alignment_start, sam.query_alignment_end, sam.query_length,
        #             (1 if sam.is_reverse == True else 0), sam.reference_start, sam.reference_end, ref_lens[sam.reference_name], 254, sam])

    return ret
