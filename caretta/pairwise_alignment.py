"""
Functions to perform a pairwise alignment
"""

from superposition import *
from utils.kabsch import kabsch_superpose
from utils.helper import make_score_matrix


def pairwise_alignment(protein_a, protein_b, gap_open_penalty=1, gap_extend_penalty=.01):
    """
    Perform pairwise alignment of two proteins using secondary_superpose and signal_superpose
    to make an initial superposition

    :param protein_a: a pdb.Protein
    :param protein_b: a pdb.Protein
    :param gap_open_penalty: gap open penalty for the final dp alignment
    :param gap_extend_penalty: gap extend penalty for final the dp alignment
    :return: an alignment score (aln_score) and the alignment of the first (aln_a) and
        second (aln_b) protein
    """
    coord_a, coord_b = protein_a.get_ca_coordinates(), protein_b.get_ca_coordinates()
    ss_a, ss_b = protein_a.get_dssp(), protein_b.get_dssp()
    score_sig_1, pos_sig_a_1, pos_sig_b_1 = signal_superpose(coord_a, coord_b, True)
    score_sig_2, pos_sig_a_2, pos_sig_b_2 = signal_superpose(coord_a, coord_b, False)
    score_sec, pos_sec_a, pos_sec_b = secondary_superpose(coord_a, coord_b, ss_a, ss_b)

    if score_sig_1 > score_sig_2:
        if score_sig_1 > score_sec:
            pos_a, pos_b = pos_sig_a_1, pos_sig_b_1
        else:
            pos_a, pos_b = pos_sec_a, pos_sec_b
    else:
        if score_sig_2 > score_sec:
            pos_a, pos_b = pos_sig_a_2, pos_sig_b_2
        else:
            pos_a, pos_b = pos_sec_a, pos_sec_b

    sup_ca, sup_cb = kabsch_superpose(coord_a, coord_b, coord_a[pos_a], coord_b[pos_b], False)
    # concatenate the consensus to the coordinates before computing Score_c
    sup_ca_cw = np.array([np.append(sup_ca[i], protein_a.get_consensus()[i]) for i in range(len(sup_ca))])
    sup_cb_cw = np.array([np.append(sup_cb[i], protein_b.get_consensus()[i]) for i in range(len(sup_cb))])
    score_matrix = make_score_matrix(sup_ca_cw, sup_cb_cw, caretta_score_func)
    aln_score, aln_a, aln_b = dp_align(score_matrix, gap_open_penalty, gap_extend_penalty)
    return aln_score, aln_a, aln_b
