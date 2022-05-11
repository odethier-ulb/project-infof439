"""
Functions to perform multiple alignments
"""

import copy
import logging
import numpy as np
from config import ALN_GAP
from pairwise_alignment import pairwise_alignment
from utils.neighbor_joining import neighbor_joining
from utils.helper import get_common_pos, get_alignment
from utils.kabsch import kabsch_superpose
from utils.pdb import Protein


def multiple_alignment(proteins, cw=1, gap_open_penalty=1, gap_extend_penalty=.01):
    """
    Align all input proteins using the pairwise_alignement and make_intermediate_node function

    :param proteins: a list of pdb.Protein
    :param cw: the initial value for the consensus row
    :param gap_open_penalty: gap open penalty for the pairwise alignment
    :param gap_extend_penalty: gap extend penalty for the pairwise alignment
    :return: aln, the final alignment of all input proteins and protein, a protein containing the mean
        coordinates of all the input protein superposed
    """
    proteins = copy.deepcopy(proteins)

    for p in proteins:
        p.set_consensus([cw for _ in range(len(p))])

    logging.info('Guide tree construction...')

    score_matrix = np.zeros((len(proteins), len(proteins)))

    for i in range(score_matrix.shape[0] - 1):
        for j in range(i + 1, score_matrix.shape[1]):
            score, _, _ = pairwise_alignment(proteins[i], proteins[j], gap_open_penalty, gap_extend_penalty)
            score_matrix[i, j] = -score

    score_matrix += score_matrix.T

    tree = neighbor_joining(score_matrix)

    logging.info('Protein alignment...')

    alignments = {p.name: p.get_sequence() for p in proteins}

    for i in range(0, len(tree) - 1, 2):
        protein_a, protein_b = proteins[tree[i]], proteins[tree[i + 1]]
        protein_ab, aln_a, aln_b = make_intermediate_node(protein_a, protein_b, gap_open_penalty, gap_extend_penalty)
        __add_alignment(protein_a, aln_a, alignments)
        __add_alignment(protein_b, aln_b, alignments)
        proteins.append(protein_ab)

    protein_a, protein_b = proteins[-1], proteins[tree[-1]]
    final_protein, aln_a, aln_b = make_intermediate_node(proteins[-1], proteins[tree[-1]],
                                                         gap_open_penalty, gap_extend_penalty)
    __add_alignment(protein_a, aln_a, alignments)
    __add_alignment(protein_b, aln_b, alignments)
    return alignments, final_protein


def make_intermediate_node(protein_a, protein_b, gap_open_penalty, gap_extend_penalty):
    """
    Combines two input proteins into an intermediate node

    :param protein_a: a pdb.Protein
    :param protein_b: a pdb.Protein
    :param gap_open_penalty: gap open penalty for the pairwise alignment
    :param gap_extend_penalty: gap extend penalty for the pairwise alignment
    :return: a pdb.Protein from a pairwise alignment between protein_a and protein_b
    """
    protein_a.set_consensus([cw * protein_b.get_len_aln() / 2 for cw in protein_a.get_consensus()])
    protein_b.set_consensus([cw * protein_a.get_len_aln() / 2 for cw in protein_b.get_consensus()])
    _, aln_a, aln_b = pairwise_alignment(protein_a, protein_b, gap_open_penalty, gap_extend_penalty)
    protein_a.set_consensus([cw * 2 / protein_b.get_len_aln() for cw in protein_a.get_consensus()])
    protein_b.set_consensus([cw * 2 / protein_a.get_len_aln() for cw in protein_b.get_consensus()])
    coord_ab = __combine_coordinates(protein_a, protein_b, aln_a, aln_b)
    sec_ab = __combine_sec_struct(protein_a, protein_b, aln_a, aln_b)
    con_ab = [0 for _ in range(len(aln_a))]
    for i in range(len(con_ab)):
        if aln_a[i] != ALN_GAP:
            con_ab[i] += protein_a.get_consensus()[aln_a[i]]
        if aln_b[i] != ALN_GAP:
            con_ab[i] += protein_b.get_consensus()[aln_b[i]]
    name = ';'.join([protein_a.name, protein_b.name])
    return Protein.from_intermediate_node(name, coord_ab, sec_ab, con_ab), aln_a, aln_b


def __combine_coordinates(protein_a, protein_b, aln_a, aln_b):
    """
    Return the mean of the coordinates in protein_a and protein_b after aligning them according
    to aln_a and aln_b
    """
    # Superpose
    coord_a, coord_b, coord_ab = protein_a.get_ca_coordinates(), protein_b.get_ca_coordinates(), []
    cmn_pos_a, cmn_pos_b = get_common_pos(aln_a, aln_b)
    sup_a, sup_b = kabsch_superpose(coord_a, coord_b, coord_a[cmn_pos_a], coord_b[cmn_pos_b], False)
    # Combine
    for i in range(len(aln_a)):
        if aln_a[i] != -1 and aln_b[i] == ALN_GAP:
            coord_ab.append(sup_a[aln_a[i]])
        elif aln_a[i] == -1 and aln_b[i] != ALN_GAP:
            coord_ab.append(sup_b[aln_b[i]])
        else:
            coord_ab.append((sup_a[aln_a[i]] + sup_b[aln_b[i]]) / 2)
    return np.array(coord_ab)


def __combine_sec_struct(protein_a, protein_b, aln_a, aln_b):
    """
    Combine secondary structures by taking the non-gap secondary structure value from the
    two proteins after alignment
    """
    ss_a, ss_b, ss_ab = protein_a.get_dssp(), protein_b.get_dssp(), []
    for i in range(len(aln_a)):
        if aln_a[i] == ALN_GAP and aln_b[i] != ALN_GAP:
            ss_ab.append(ss_b[aln_b[i]])
        else:
            ss_ab.append(ss_a[aln_a[i]])
    return ss_ab


def __add_alignment(protein, aln, alignments):
    """ Add alignment of the given protein to the current ones """
    for prot_name in protein.name.split(';'):
        alignments[prot_name] = get_alignment(alignments[prot_name], aln)
