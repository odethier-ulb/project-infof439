"""
Functions to perform secondary and signal superposition
"""

import numpy as np
from utils.dtw import dp_align
from utils.helper import get_common_pos, make_score_matrix
from utils.kabsch import kabsch_superpose


def caretta_score(sup_c_a, sup_c_b, gamma=.03):
    """ Measure the similarity score between two sets of superposed coordinates """
    score = 0
    for i in range(len(sup_c_a)):
        score += caretta_score_func(sup_c_a[i], sup_c_b[i], gamma)
    return score


def caretta_score_func(coord_a, coord_b, gamma=.03):
    """ Measure the similarity score between two coordinates """
    return np.exp(-gamma * np.sum((coord_a - coord_b) ** 2))


def secondary_superpose(coord_a, coord_b, ss_a, ss_b, gap_open_penalty=1, gap_extend_penalty=0):
    """
    Align two proteins according to their secondary structure elements

    :param coord_a: the (x, y, z) coordinates of alpha-carbon atoms of protein A
    :param coord_b: the (x, y, z) coordinates of alpha-carbon atoms of protein B
    :param ss_a: the secondary structure codes for each residue of protein A
    :param ss_b: the secondary structure codes for each residue of protein B
    :param gap_open_penalty: gap open penalty for the dp alignment
    :param gap_extend_penalty: gap extend penalty for the dp alignment
    :return: the score of the superposition (score_sec),
        aligned positions of protein A (pos_sec_a) and protein B (pos_sec_b)
    """
    score_matrix = np.zeros((len(ss_a), len(ss_b)), np.short)

    for i in range(len(ss_a)):
        for j in range(len(ss_b)):
            if ss_a[i] == ss_b[j]:
                if ss_a[i] != '-':
                    score_matrix[i, j] = 1
            else:
                score_matrix[i, j] = -1

    _, aln_a, aln_b = dp_align(score_matrix, gap_open_penalty, gap_extend_penalty)
    pos_sec_a, pos_sec_b = get_common_pos(aln_a, aln_b)
    sup_c_a, sup_c_b = kabsch_superpose(coord_a, coord_b, coord_a[pos_sec_a], coord_b[pos_sec_b])
    score_sec = caretta_score(sup_c_a, sup_c_b)
    return score_sec, pos_sec_a, pos_sec_b


def signal_superpose(coord_a, coord_b, from_first=True, size=30, gap_open_penalty=1, gap_extend_penalty=.01):
    """
    Align two proteins on one-dimensional signals of distance from all residues to the first
    or last residue in a segment

    :param coord_a: the (x, y, z) coordinates of alpha-carbon atoms of protein A
    :param coord_b: the (x, y, z) coordinates of alpha-carbon atoms of protein B
    :param from_first: True if distance from the first residue, False for the last one
    :param size: segment size
    :param gap_open_penalty: gap open penalty for the dp alignment
    :param gap_extend_penalty: gap extend penalty for the dp alignment
    :return: the score of superposition (score_signal), aligned positions of
        protein A (pos_signal_a) and protein B (pos_signal_b)
    """
    sup_c_a, sup_c_b = __signal_superpose(coord_a, coord_b, from_first, size)
    score_matrix = make_score_matrix(sup_c_a, sup_c_b, caretta_score_func)
    _, aln_a, aln_b = dp_align(score_matrix, gap_open_penalty, gap_extend_penalty)
    pos_a, pos_b = get_common_pos(aln_a, aln_b)
    sup_c_a, sup_c_b = kabsch_superpose(sup_c_a, sup_c_b, sup_c_a[pos_a], sup_c_b[pos_b])
    return caretta_score(sup_c_a, sup_c_b), pos_a, pos_b


def __signal_superpose(coord_a, coord_b, from_first, size):
    """ Initial superposition of the whole structure based on the signal """
    signal_a = np.zeros((coord_a.shape[0] - size, size))
    signal_b = np.zeros((coord_b.shape[0] - size, size))
    # ref is the reference coordinates in each segment
    ref_a = np.zeros((signal_a.shape[0], coord_a.shape[1]))
    ref_b = np.zeros((signal_b.shape[0], coord_b.shape[1]))
    offset = 0 if from_first else size - 1

    for i in range(len(coord_a) - size):
        for j in range(size):
            signal_a[i, j] = np.linalg.norm(coord_a[i + j] - coord_a[i + offset])
        ref_a[i] = coord_a[i + offset]

    for i in range(len(coord_b) - size):
        for j in range(size):
            signal_b[i, j] = np.linalg.norm(coord_b[i + j] - coord_b[i + offset])
        ref_b[i] = coord_b[i + offset]

    score_matrix = make_score_matrix(signal_a, signal_b, lambda x, y: np.median(np.exp(-((x - y) ** 2) / 10)))
    _, aln_a, aln_b = dp_align(score_matrix, 0, 0)
    cmn_a, cmn_b = get_common_pos(aln_a, aln_b)
    aln_coord_a = np.array([ref_a[idx] for idx in cmn_a])
    aln_coord_b = np.array([ref_b[idx] for idx in cmn_b])
    sup_c_a, sup_c_b = kabsch_superpose(coord_a, coord_b, aln_coord_a, aln_coord_b, False)
    return sup_c_a, sup_c_b
