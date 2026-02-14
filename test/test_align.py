# Importing Dependencies
import pytest
from align.align import NeedlemanWunsch, read_fasta
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
    #init
    aligning = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    
    #results
    score, seq1_aligned, seq2_aligned = aligning.align(seq1, seq2)
    

    assert aligning._align_matrix.shape == (5, 4), "Wrong dim"
    assert aligning._back.shape == (5, 4), "Backtrace wrong dim"
    

    expected_first_col = np.array([0, -10, -20, -30, -40])
    assert np.allclose(aligning._align_matrix[:, 0], expected_first_col), "First column wrong init"
    
    expected_first_row = np.array([0, -10, -20, -30])
    assert np.allclose(aligning._align_matrix[0, :], expected_first_row), "First row wrong"
    
    # Googled this as edge case
    assert not np.isnan(score), "Alignment score is NaN"
    assert score != 0, "Alignment score should not be zero"
    

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

    aligning = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    
    score, seq3_aligned, seq4_aligned = aligning.align(seq3, seq4)
    
    # Test if 17
    assert score == 17, f"Expected alignment score of 17, but got {score}"
    
    # MAVHQLIRRP
    # M---QLIRHP
    expected_seq3_aligned = "MAVHQLIRRP"
    expected_seq4_aligned = "M---QLIRHP"
    
    assert seq3_aligned == expected_seq3_aligned, f"seq3 alignment incorrect."
    assert seq4_aligned == expected_seq4_aligned, f"seq4 alignment incorrect."
    









    





