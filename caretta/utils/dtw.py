"""
Util. functions to perform dynamic programming alignment

Inspired from https://www.youtube.com/watch?v=NqYY0PJbD3s for the traceback
"""

import numpy as np


def dp_align(score_matrix, gap_open_penalty=0, gap_extend_penalty=0):
    """
    Perform a generic affine gap cost dynamic programming alignment

    :param score_matrix: a 2D score matrix as a numpy ndarray
    :param gap_open_penalty: the gap open penalty
    :param gap_extend_penalty: the gap extension penalty
    :return: aln_score, aln_a, aln_b : the alignment score, the alignment of the first and second
        proteins as index arrays
    """
    scores, backtrack = __compute_matrix(score_matrix, gap_open_penalty, gap_extend_penalty)
    max_score = max(scores)
    aln_a, aln_b = __get_alignment(backtrack, scores[scores.index(max_score)] + 1)
    return max_score, aln_a, aln_b


def __compute_matrix(score_matrix, gap_open_penalty, gap_extend_penalty):
    n, m = score_matrix.shape[0] + 1, score_matrix.shape[1] + 1
    I, E = abs(gap_open_penalty), abs(gap_extend_penalty)
    S = np.zeros((n, m), np.dtype('f8'))
    V = np.zeros((n, m), np.dtype('f8'))
    W = np.zeros((n, m), np.dtype('f8'))
    backtrack = np.zeros((n, m, 3), np.ushort)
    V[0, :], W[:, 0] = -np.Infinity, -np.Infinity
    V[0, 0], W[0, 0] = 0, 0

    for i in range(1, n):
        for j in range(1, m):
            V[i, j] = max(S[i - 1, j] - I, V[i - 1, j] - E)
            W[i, j] = max(S[i, j - 1] - I, W[i, j - 1] - E)
            S[i, j] = max(S[i - 1, j - 1] + score_matrix[i - 1, j - 1], V[i, j], W[i, j])
            backtrack[i, j, 0] = 1 if V[i, j] == V[i - 1, j] - E else 3
            backtrack[i, j, 1] = 2 if W[i, j] == W[i, j - 1] - E else 3
            backtrack[i, j, 2] = 1 if S[i, j] == V[i, j] else 2 if S[i, j] == W[i, j] else 3

    return (S[n - 1, m - 1], V[n - 1, m - 1], W[n - 1, m - 1]), backtrack


def __get_alignment(backtrack, trace):
    aln_a, aln_b = [], []
    i = backtrack.shape[0] - 1
    j = backtrack.shape[1] - 1

    while not (i == 0 and j == 0):
        if i == 0:
            j -= 1
            aln_a.insert(0, -1)
            aln_b.insert(0, max(-1, j))
        elif j == 0:
            i -= 1
            aln_a.insert(0, max(-1, i))
            aln_b.insert(0, -1)
        else:
            if trace == 1:
                trace = backtrack[i, j, 0]
                i -= 1
                aln_a.insert(0, max(-1, i))
                aln_b.insert(0, -1)
            elif trace == 2:
                trace = backtrack[i, j, 1]
                j -= 1
                aln_a.insert(0, -1)
                aln_b.insert(0, max(-1, j))
            else:
                trace = backtrack[i, j, 2]
                if trace == 3:
                    i -= 1
                    j -= 1
                    aln_a.insert(0, max(-1, i))
                    aln_b.insert(0, max(-1, j))

    return aln_a, aln_b
