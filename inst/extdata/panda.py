import math
import time
import pandas as pd
from scipy.stats import zscore
import numpy as np
import os
import math

#
# GPU Calculation functions:
# These functions were defined in Panda, but are
# also shared by Puma and possibly others.
# We are going to start putting these here so they can be
# shared by all classes as they are independent from the class


def gt_function(x, y=None):
    """
    Description:
        Continuous Tanimoto similarity function computed on the GPU
    Inputs:
        x: First object to measure the distance from. If only this matrix is provided, then the distance is meausred between the columns of x.
        y: Second object to measure the distance to
    Ouputs:
        a_matrix: Matrix containing the pairwsie distances.
    """
    if y is None:
        a_matrix = cp.dot(x, x.T)
        s = cp.square(x).sum(axis=1)
        a_matrix /= cp.sqrt(s + s.reshape(-1, 1) - cp.abs(a_matrix))
    else:
        a_matrix = cp.dot(x, y)
        a_matrix /= cp.sqrt(
            cp.square(y).sum(axis=0)
            + cp.square(x).sum(axis=1).reshape(-1, 1)
            - cp.abs(a_matrix)
        )
    return a_matrix


def gupdate_diagonal(diagonal_matrix, num, alpha, step):
    """
    Description:
        Updates the diagonal of the input matrix in the message passing computed on the GPU
    Inputs:
        diagonal_matrix: Input diagonal matrix.
        num            : Number of rows/columns.
        alpha          : Learning rate.
        step           : The current step in the algorithm.
    """
    cp.fill_diagonal(diagonal_matrix, cp.nan)
    diagonal_std = cp.nanstd(diagonal_matrix, 1)
    diagonal_fill = diagonal_std * num * math.exp(2 * alpha * step)
    cp.fill_diagonal(diagonal_matrix, diagonal_fill)
    return diagonal_matrix


def compute_panda_gpu(
    correlation_matrix,
    ppi_matrix,
    motif_matrix,
    computing="cpu",
    threshold=0.001,
    alpha=0.1,
):

    """
    Panda network optimization
    Args:
        correlation_matrix (numpy float): coexpression matrix
        ppi_matrix (numpy float): PPI network matrix
        motif_matrix (numpy float): motif matrix
        computing (str) : either cpu or gpu. Defaults to 'cpu'
        threshold (float, optional): hamming distance threshold for stop. Defaults to 0.001.
        alpha (float, optional): learning rate. Defaults to 0.1
    """
    num_tfs, num_genes = motif_matrix.shape
    step = 0
    hamming = 1

    ppi_matrix = cp.array(ppi_matrix)
    motif_matrix = cp.array(motif_matrix)
    correlation_matrix = cp.array(correlation_matrix)

    if hamming > 0.001:

        W = 0.5 * (
            gt_function(ppi_matrix, motif_matrix)
            + gt_function(motif_matrix, correlation_matrix)
        )  # W = (R + A) / 2
        hamming = cp.abs(motif_matrix - W).mean()
        # update motif matrix
        motif_matrix *= 1 - alpha
        motif_matrix += alpha * W
        # Update ppi_matrix
        ppi = gt_function(motif_matrix)  # t_func(X, X.T)
        ppi = gupdate_diagonal(ppi, num_tfs, alpha, step)
        ppi_matrix *= 1 - alpha
        ppi_matrix += alpha * ppi

        # Update correlation_matrix
        motif = gt_function(motif_matrix.T)
        motif = gupdate_diagonal(motif, num_genes, alpha, step)
        correlation_matrix *= 1 - alpha
        correlation_matrix += alpha * motif
        # del W, ppi, motif  # release memory for next step
        print("step: {}, hamming: {}".format(step, hamming))
        step = step + 1

    motif_matrix = cp.asnumpy(motif_matrix)
    return motif_matrix

# this scirpt is derived from https://github.com/netZoo/netZooPy/blob/master/netZooPy/panda/panda.py
# and https://github.com/netZoo/netZooPy/blob/master/netZooPy/panda/timer.py

class Timer(object):
    def __init__(self, name=None):
        if name:
            print(name)

    def __enter__(self):
        self.tic = time.time()

    def __exit__(self, type, value, traceback):
        print('  Elapsed time: %.2f sec.' % (time.time() - self.tic))



#
# Calculation functions:
# These functions were defined in Panda, but are
# also shared by Puma and possibly others.
# We are going to start putting these here so they can be
# shared by all classes as they are independent from the class


def t_function(x, y=None):
    """
    Description:
        Continuous Tanimoto similarity function computed on the CPU.
    Inputs:
                x: First object to measure the distance from. If only this matrix is provided, then the distance is meausred between the columns of x.
        y: Second object to measure the distance to.
    Ouputs:
        a_matrix: Matrix containing the pairwise distances.
    """
    if y is None:
        a_matrix = np.dot(x, x.T)
        s = np.square(x).sum(axis=1)
        a_matrix /= np.sqrt(s + s.reshape(-1, 1) - np.abs(a_matrix))
    else:
        a_matrix = np.dot(x, y)
        a_matrix /= np.sqrt(
            np.square(y).sum(axis=0)
            + np.square(x).sum(axis=1).reshape(-1, 1)
            - np.abs(a_matrix)
        )
    return a_matrix


def update_diagonal(diagonal_matrix, num, alpha, step):
    """
    Description:
        Updates the diagonal of the input matrix in the message passing computed on the CPU
    Inputs:
        diagonal_matrix: Input diagonal matrix.
        num            : Number of rows/columns.
        alpha          : Learning rate.
        step           : The current step in the algorithm.
    """
    np.fill_diagonal(diagonal_matrix, np.nan)
    diagonal_std = np.nanstd(diagonal_matrix, 1)
    diagonal_fill = diagonal_std * num * math.exp(2 * alpha * step)
    np.fill_diagonal(diagonal_matrix, diagonal_fill)
    return diagonal_matrix


def compute_panda_cpu(
    correlation_matrix,
    ppi_matrix,
    motif_matrix,
    threshold=0.001,
    alpha=0.1,
):
    """Panda network optimization
    Args:
        correlation_matrix (numpy float): coexpression matrix
        ppi_matrix (numpy float): PPI network matrix
        motif_matrix (numpy float): motif matrix
        threshold (float, optional): hamming distance threshold for stop. Defaults to 0.001.
        alpha (float, optional): learning rate. Defaults to 0.1
    """
    motif_matrix = motif_matrix.copy()
    ppi_matrix = ppi_matrix.copy()
    correlation_matrix = correlation_matrix.copy()

    num_tfs, num_genes = motif_matrix.shape
    step = 0
    hamming = 1

    while hamming > threshold:
        W = 0.5 * (
            t_function(ppi_matrix, motif_matrix)
            + t_function(motif_matrix, correlation_matrix)
        )  # W = (R + A) / 2
        hamming = np.abs(motif_matrix - W).mean()
        # Update motif matrix
        motif_matrix *= 1 - alpha
        motif_matrix += alpha * W
        # Update ppi_matrix
        ppi = t_function(motif_matrix)  # t_func(X, X.T)
        ppi = update_diagonal(ppi, num_tfs, alpha, step)
        ppi_matrix *= 1 - alpha
        ppi_matrix += alpha * ppi
        # Update correlation_matrix
        motif = t_function(motif_matrix.T)
        motif = update_diagonal(motif, num_genes, alpha, step)
        correlation_matrix *= 1 - alpha
        correlation_matrix += alpha * motif
        # del W, ppi, motif  # release memory for next step
        print("step: {}, hamming: {}".format(step, hamming))
        step = step + 1

    return motif_matrix


def compute_panda(
    correlation_matrix,
    ppi_matrix,
    motif_matrix,
    computing="cpu",
    threshold=0.001,
    alpha=0.1,
):
    """Panda network optimization
    Args:
        correlation_matrix (numpy float): coexpression matrix
        ppi_matrix (numpy float): PPI network matrix
        motif_matrix (numpy float): motif matrix
        computing (str) : either cpu or gpu. Defaults to 'cpu'
        threshold (float, optional): hamming distance threshold for stop. Defaults to 0.001.
        alpha (float, optional): learning rate. Defaults to 0.1
    """

    if computing == "cpu":
        # Initialise W and hamming
        motif_matrix = compute_panda_cpu(
            correlation_matrix,
            ppi_matrix,
            motif_matrix,
            alpha=alpha,
        )

    elif computing == "gpu":
        import cupy as cp
        from netZooPy.panda.calculations_gpu import compute_panda_gpu

        motif_matrix = compute_panda_gpu(
            correlation_matrix,
            ppi_matrix,
            motif_matrix,
            alpha=alpha,
        )
    else:
        sys.error("ERROR: %s is not an existing computing device" % str(computing))
    return motif_matrix


def normalize_network(x):
    """
    Description:
        Standardizes the input data matrices.
    Inputs:
        x     : Input adjacency matrix.
    Outputs:
        normalized_matrix: Standardized adjacency matrix.
    """
    norm_col = zscore(x, axis=0)
    if x.shape[0] == x.shape[1]:
        norm_row = norm_col.T
    else:
        norm_row = zscore(x, axis=1)
    # Alessandro: replace nan values
    normalized_matrix = (norm_col + norm_row) / math.sqrt(2)
    norm_total = (x - np.mean(x)) / np.std(x)  # NB zscore(x) is not the same
    nan_col = np.isnan(norm_col)
    nan_row = np.isnan(norm_row)
    normalized_matrix[nan_col] = (norm_row[nan_col] + norm_total[nan_col]) / math.sqrt(
        2
    )
    normalized_matrix[nan_row] = (norm_col[nan_row] + norm_total[nan_row]) / math.sqrt(
        2
    )
    normalized_matrix[nan_col & nan_row] = (
        2 * norm_total[nan_col & nan_row] / math.sqrt(2)
    )
    return normalized_matrix
