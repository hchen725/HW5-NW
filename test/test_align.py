# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat",
                         gap_open = -10,
                         gap_extend = -1)
    
    score, seqA, seqB = NW.align(seq1, seq2)
    assert score == 4
    assert seqA == "MYQR"
    assert seqB == "M-QR"

    i = -np.inf
    m = np.array([[0,   i,    i,   i],
                  [i,   5,  -11,   13],
                  [i, -12,    4,   -8],
                  [i, -12,   -1,    5],
                  [i, -14,   -6,    4]])
    assert m.all() == NW._align_matrix.all()
    
    known_A = np.array([[-10,   i,    i,    i],
                        [-11, -12,   -6,   -7],
                        [-12, -13,  -14,   -7],
                        [-13, -14,  -15,  -12],
                        [-14, -15,  -16,  -17]])
    assert known_A.all() == NW._gapA_matrix.all()

    known_B = np.array([[-10, -11,  -12, -13],
                        [i,   -12,  -13, -14],
                        [i,    -6,  -14, -15],
                        [i,    -7,   -7, -16],
                        [i,    -8,   -8,  -6]])
    assert known_B.all() == NW._gapB_matrix.all()
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat",
                         gap_open = -10,
                         gap_extend = -1)
    
    assert NW.align(seq3, seq4) == (17, 
                                    'MAVHQLIRRP', 
                                    'M---QLIRHP') 




