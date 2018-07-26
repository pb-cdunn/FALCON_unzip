import falcon_unzip.proto.sam2m4 as M
import pytest

class Aln(object):
    def __init__(self,
            is_unmapped=False,
            is_supplementary=False,
            is_secondary=False,
            is_reverse=False,
            cigar=None,
            query_sequence=list(),
            query_length=0,
            query_alignment_start=0,
            query_alignment_end=0,
            reference_start=0,
            reference_end=0,
            reference_name='',
            qname='',
        ):
        self.is_unmapped = is_unmapped
        self.is_supplementary = is_supplementary
        self.is_secondary = is_secondary
        self.is_reverse = is_reverse
        self.cigar = cigar
        self.query_sequence = query_sequence
        self.query_length = query_length
        self.query_alignment_start = query_alignment_start
        self.query_alignment_end = query_alignment_end
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.reference_name = reference_name
        self.qname = qname

def test_pysam_to_m4_unmapped():
    aln = Aln(is_unmapped=True)
    assert M.pysam_to_m4(aln) == None
def test_pysam_to_m4_supplementary():
    aln = Aln(is_supplementary=True)
    assert M.pysam_to_m4(aln) == None
def test_pysam_to_m4_secondary():
    aln = Aln(is_secondary=True)
    assert M.pysam_to_m4(aln) == None
def test_pysam_to_m4_no_cigar():
    aln = Aln(cigar=[])
    assert M.pysam_to_m4(aln) == None
def test_pysam_to_m4_no_query():
    aln = Aln(cigar=[(7,0)], query_sequence='')
    assert M.pysam_to_m4(aln) == None
    aln = Aln(cigar=[(7,0)], query_sequence=None)
    assert M.pysam_to_m4(aln) == None
def test_pysam_to_m4_no_aln():
    aln = Aln(cigar=[(7,0)], query_sequence='A', query_alignment_start=0, query_alignment_end=0)
    assert M.pysam_to_m4(aln) == None
def test_pysam_to_m4():
    aln = Aln(
            cigar=[(7,0)], query_sequence='A', query_length=9,
            query_alignment_start=0, query_alignment_end=42,
            reference_start=3, reference_end=5,
            qname='Q', reference_name='REF',
    )
    identity = 0.0
    assert M.pysam_to_m4(aln) == ['Q', 'REF', 0, identity, 0, 0, 42, 9, 0, 3, 5, 0, 254, aln]
def test_pysam_to_m4_identity():
    aln = Aln(
            cigar=[(7,1)], query_sequence='A', query_length=9,
            query_alignment_start=0, query_alignment_end=2,
            reference_start=3, reference_end=5,
            qname='Q', reference_name='REF',
    )
    identity = 50.0
    assert M.pysam_to_m4(aln) == ['Q', 'REF', 0, identity, 0, 0, 2, 9, 0, 3, 5, 0, 254, aln]
