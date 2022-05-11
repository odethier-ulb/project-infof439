"""
Implements the Kabsch algorithm in order to calculate the optimal rotation matrix that minimizes the
RMSD between two paired sets of points
"""

import numpy as np
from utils.helper import get_aln_idx, get_common_pos


def kabsch_superpose(coord_a, coord_b, cmn_coord_a, cmn_coord_b, superpose_common=True):
    """
    Performs superposition of the pos_a positions of coord_a coordinates onto the pos_b
    positions of coord_b coordinates using the Kabsch algorithm

    :param coord_a: an array of (x, y, z) coordinates
    :param coord_b: an array of (x, y, z) coordinates
    :param cmn_coord_a: the common coordinates in a
    :param cmn_coord_b: the common coordinates in b
    :param superpose_common: True if the superposition is applied only on given position
    :return: the superposed coordinates sup_ca and sup_cb
    """
    a_centroid, b_centroid = np.average(cmn_coord_a, axis=0), np.average(cmn_coord_b, axis=0)
    p, q = cmn_coord_a - a_centroid, cmn_coord_b - b_centroid
    h = p.T @ q
    u, _, vt = np.linalg.svd(h)
    d = np.linalg.det(vt.T @ u.T)
    r = vt.T @ np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1 if d < 0 else 1]]) @ u.T
    t = a_centroid - (b_centroid @ r)
    if superpose_common:
        return cmn_coord_a, (cmn_coord_b @ r) + t
    else:
        return coord_a - a_centroid, (coord_b - b_centroid) @ r


def get_rmsd(coord_a, coord_b):
    """
    Get the Root-mean-square deviation of two set of points

    :param coord_a: an array of (x, y, z) coordinates
    :param coord_b: an array of (x, y, z) coordinates
    :return: the rsmd between coord_a and coord_b
    """
    return np.sqrt(np.sum((coord_a - coord_b)**2) / coord_a.shape[0])


def get_aligned_coordinates(final_protein, proteins, alignments):
    """
    Align the coordinates of each protein to the final one created by the msa

    :param final_protein: the final protein resulting of the msa
    :param proteins: a List of pdb.Protein
    :param alignments: the msa result
    :return: the aligned coordinates
    """
    aligned_coordinates = {}
    coord_a = final_protein.get_ca_coordinates()
    aln_a = [i for i in range(len(coord_a))]

    for p in proteins:
        coord_b = p.get_coordinates()
        aln_b = get_aln_idx(alignments[p.name])
        cmn_a, cmn_b = get_common_pos(aln_a, aln_b)
        cmn_pos_a = coord_a[cmn_a]
        cmn_pos_b = p.get_ca_coordinates()[cmn_b]
        _, coordinates = kabsch_superpose(coord_a, coord_b, cmn_pos_a, cmn_pos_b, False)
        aligned_coordinates[p.name] = coordinates

    return aligned_coordinates


def get_final_rmsd(alignments, proteins):
    """
    Get the Root-mean-square deviation of an alignment. It is calculated
    for every pair of structures in the msa after superposing all structures to
    one reference structure, the longest protein

    :param alignments: the alignments of the proteins
    :param proteins: a List of pdb.Protein

    :return: the rsmd computed as the average of each pair of protein
    """
    # 1. superpose to ref
    ref_idx, superposed_coord = np.argmax([len(p) for p in proteins]), []
    coord_a = proteins[ref_idx].get_ca_coordinates() - np.average(proteins[ref_idx].get_ca_coordinates(), axis=0)
    aln_a = get_aln_idx(alignments[proteins[ref_idx].name])

    for i in range(len(proteins)):
        if i != ref_idx:
            coord_b = proteins[i].get_ca_coordinates()
            aln_b = get_aln_idx(alignments[proteins[i].name])
            cmn_a, cmn_b = get_common_pos(aln_a, aln_b)
            cmn_pos_b = proteins[i].get_ca_coordinates()[cmn_b]
            _, coordinates = kabsch_superpose(coord_a, coord_b, coord_a[cmn_a], cmn_pos_b, False)
            superposed_coord.append(coordinates)
        else:
            superposed_coord.append(coord_a)

    # 2. compute rmsd for each pair
    pairwise_rmsd = np.zeros((len(proteins), len(proteins)))
    pairwise_rmsd[:] = np.nan
    for i in range(len(proteins)):
        for j in range(i + 1, len(proteins)):
            aln_a = get_aln_idx(alignments[proteins[i].name])
            aln_b = get_aln_idx(alignments[proteins[j].name])
            cmn_a, cmn_b = get_common_pos(aln_a, aln_b)
            pairwise_rmsd[i, j] = pairwise_rmsd[j, i] = \
                get_rmsd(superposed_coord[i][cmn_a], superposed_coord[j][cmn_b])

    return np.nanmean(pairwise_rmsd)
