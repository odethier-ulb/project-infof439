"""
Functions to do the Neighbor Joining

From https://en.wikipedia.org/wiki/Neighbor_joining
"""

import numpy as np


def neighbor_joining(score_matrix):
    """
    Uses the neighbor joining algorithm with a maximum linkage scheme to compute a guide tree for
    the progression of multiple alignment

    :param score_matrix: a score (distance) matrix
    :return: a List of indexes, looping through it 2 by 2 gives the 2 nodes to join
    """
    n = score_matrix.shape[0]
    tree, new_node = [], n
    indices = [i for i in range(n)]
    dist_matrix = score_matrix

    while n > 3:
        min_i, min_j = __get_closest(dist_matrix)
        tree.append(indices[min_i])
        tree.append(indices[min_j])
        dist_matrix, tmp_indices = __get_new_distance_matrix(dist_matrix, min_i, min_j, indices)
        indices = [new_node] + [indices[i] for i in tmp_indices]
        new_node += 1
        n -= 1

    tree.append(indices[1])
    tree.append(indices[2])
    tree.append(indices[0])
    return tree


def __get_closest(q_matrix):
    """ Return the pair of closest indexes """
    min_dist, min_i, min_j = np.inf, -1, -1
    q_matrix = __get_q_matrix(q_matrix)
    for i in range(q_matrix.shape[0]):
        for j in range(q_matrix.shape[1]):
            if i != j and q_matrix[i, j] < min_dist:
                min_dist, min_i, min_j = q_matrix[i, j], i, j
    return min_i, min_j


def __get_q_matrix(dist_matrix):
    """ Get the distance matrix Q based on the given distance matrix """
    n = dist_matrix.shape[0]
    q_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                q_matrix[i, j] = (n - 2) * dist_matrix[i, j] - np.sum(dist_matrix[i, :]) - np.sum(dist_matrix[j, :])
    return q_matrix


def __get_new_distance_matrix(dist_matrix, min_i, min_j, indices):
    """ Compute a new distance matrix by creating a new node from min i and min j """
    new_matrix = np.zeros((dist_matrix.shape[0] - 1, dist_matrix.shape[0] - 1))
    new_indices = [i for i in range(len(indices)) if i != min_i and i != min_j]
    new_matrix[1:, 1:] = dist_matrix[new_indices, :][:, new_indices]
    for i in range(len(new_indices)):
        new_matrix[0, i + 1] = new_matrix[i + 1, 0] = \
            0.5 * (dist_matrix[min_i, new_indices[i]] + dist_matrix[min_j, new_indices[i]] - dist_matrix[min_i, min_j])
    return new_matrix, new_indices
