"""
Util. functions to perform misc. things
"""

import numpy as np
from config import ALN_GAP, SEQ_GAP


def make_score_matrix(seq_a, seq_b, score_funct):
    """
    Compute a score matrix between two sequences

    :param seq_a: an Iterable object
    :param seq_b: an Iterable object
    :param score_funct: a function giving a score between 2 elements of seq_a and seq_b
    :return: a score matrix of seq_a and seq_b
    """
    matrix = np.zeros(shape=(len(seq_a), len(seq_b)), dtype=np.dtype('f8'))
    for i, a in enumerate(seq_a):
        for j, b in enumerate(seq_b):
            matrix[i, j] = score_funct(a, b)
    return matrix


def get_common_pos(aln_a, aln_b):
    """
    Return the common positions (wihout gap) between 2 alignments

    :param aln_a: alignment a of size x
    :param aln_b: alignment b of size x
    :return: the common positions
    """
    cmn_a = [aln_a[i] for i in range(len(aln_a)) if aln_a[i] != ALN_GAP and aln_b[i] != ALN_GAP]
    cmn_b = [aln_b[i] for i in range(len(aln_b)) if aln_a[i] != ALN_GAP and aln_b[i] != ALN_GAP]
    return cmn_a, cmn_b


def get_aln_idx(aligned_seq):
    """
    Retrieve the alignment indices from an aligned sequence

    :param aligned_seq: an aligned sequence
    :return: the aligned index list
    """
    aln, idx = [], 0
    for r in aligned_seq:
        if r != SEQ_GAP:
            aln.append(idx)
            idx += 1
        else:
            aln.append(ALN_GAP)
    return aln


def get_alignment(seq, aln):
    """
    Get an aligned sequence

    :param seq: an Iterable object
    :param aln: an alignment from dp_align
    """
    return ''.join([SEQ_GAP if i == ALN_GAP else str(seq[i]) for i in aln])


def get_ratio_gapless(alignments):
    """
    Compute the ratio of positions in an alignments that don't
    contain any gap

    :param alignments: the result of a msa
    :return: the ratio of gap-less positions
    """
    gap_less, aln_len = 0, len(list(alignments.values())[0])
    for i in range(aln_len):
        if __is_gap_less(alignments, i):
            gap_less += 1
    return np.around(gap_less / aln_len, 3)


def get_equivalence_score(ref_alignments, alignments):
    """
    Return the percentage of gap-less positions which are present in the reference alignments

    :param ref_alignments: the reference alignments
    :param alignments: the alignments to compare
    :return: the percentage of gap-less positions compared to the ref alignments
    """
    equivalence, gap_less, aln_len = 0, 0, len(list(alignments.values())[0])
    for i in range(aln_len):
        if __is_gap_less(alignments, i):
            gap_less += 1
            if __is_gap_less(ref_alignments, i):
                equivalence += 1
    return np.around(equivalence / gap_less, 3)


def __is_gap_less(alignments, idx):
    """ Return True if the alignment is gap-less at the current index """
    non_gap_count = 0
    for k in alignments.keys():
        if idx < len(alignments[k]) and alignments[k][idx] != SEQ_GAP:
            non_gap_count += 1
    return non_gap_count == len(alignments)
