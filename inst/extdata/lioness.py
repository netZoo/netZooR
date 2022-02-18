import os, os.path, sys
import numpy as np
import pandas as pd

from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects

import math
from scipy.stats import zscore

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
    
class Lioness(Panda):
    """
    Description:
       Using LIONESS to infer single-sample gene regulatory networks.

    Usage:
       1. Reading in PANDA network and preprocessed middle data
       2. Computing coexpression network
       3. Normalizing coexpression network
       4. Running PANDA algorithm
       5. Writing out LIONESS networks

    Inputs:
       obj: PANDA object, generated with keep_expression_matrix=True.
       obj.motif_matrix: TF DNA motif binding data in tf-by-gene format.
                         If set to None, Lioness will be performed on gene coexpression network.
       computing  : 'cpu' uses Central Processing Unit (CPU) to run PANDA
                    'gpu' use the Graphical Processing Unit (GPU) to run PANDA

    Authors: 
       cychen, davidvi
    """

    def __init__(self, obj, computing='cpu', start=1, end=None, save_single_network=False, save_dir='lioness_output', save_fmt='npy'):
        # Load data
        with Timer("Loading input data ..."):
            self.export_panda_results = obj.export_panda_results
            self.expression_matrix = obj.expression_matrix
            self.motif_matrix = obj.motif_matrix
            self.ppi_matrix = obj.ppi_matrix
            self.computing=computing
            self.save_single_network=save_single_network
            if hasattr(obj,'panda_network'):
                self.network = obj.panda_network
            elif hasattr(obj,'puma_network'):
                self.network = obj.puma_network
            else:
                print('Cannot find panda or puma network in object')
                raise AttributeError('Cannot find panda or puma network in object')
            del obj

        # Get sample range to iterate
        self.n_conditions = self.expression_matrix.shape[1]
        self.indexes = range(self.n_conditions)[start-1:end]  # sample indexes to include

        if self.save_single_network:
            # Create the output folder if not exists
            self.save_dir = save_dir
            self.save_fmt = save_fmt
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)

        # Run LIONESS
        self.total_lioness_network = self.__lioness_loop()

        # create result data frame
        self.export_lioness_results = pd.DataFrame(self.total_lioness_network)

    def __lioness_loop(self):
        for i in self.indexes:
            print("Running LIONESS for sample %d:" % (i+1))
            idx = [x for x in range(self.n_conditions) if x != i]  # all samples except i

            with Timer("Computing coexpression network:"):
                if self.computing=='gpu':
                    import cupy as cp
                    correlation_matrix = cp.corrcoef(self.expression_matrix[:, idx])
                    if cp.isnan(correlation_matrix).any():
                        cp.fill_diagonal(correlation_matrix, 1)
                        correlation_matrix = cp.nan_to_num(correlation_matrix)
                    correlation_matrix=cp.asnumpy(correlation_matrix)
                else:
                    correlation_matrix = np.corrcoef(self.expression_matrix[:, idx])
                    if np.isnan(correlation_matrix).any():
                        np.fill_diagonal(correlation_matrix, 1)
                        correlation_matrix = np.nan_to_num(correlation_matrix)

            with Timer("Normalizing networks:"):
                correlation_matrix_orig = correlation_matrix # save matrix before normalization
                correlation_matrix = self._normalize_network(correlation_matrix)

            with Timer("Inferring LIONESS network:"):
                if self.motif_matrix is not None:
                    del correlation_matrix_orig
                    subset_panda_network = self.panda_loop(correlation_matrix, np.copy(self.motif_matrix), np.copy(self.ppi_matrix),self.computing)
                else:
                    del correlation_matrix
                    subset_panda_network = correlation_matrix_orig

            lioness_network = self.n_conditions * (self.network - subset_panda_network) + subset_panda_network
            if self.save_single_network:
                with Timer("Saving LIONESS network %d to %s using %s format:" % (i+1, self.save_dir, self.save_fmt)):
                    path = os.path.join(self.save_dir, "lioness.%d.%s" % (i+1, self.save_fmt))
                    if self.save_fmt == 'txt':
                        np.savetxt(path, lioness_network)
                    elif self.save_fmt == 'npy':
                        np.save(path, lioness_network)
                    elif self.save_fmt == 'mat':
                        from scipy.io import savemat
                        savemat(path, {'PredNet': lioness_network})
                    else:
                        print("Unknown format %s! Use npy format instead." % self.save_fmt)
                        np.save(path, lioness_network)
            if i == 0:
                self.total_lioness_network = np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)
            else:
                self.total_lioness_network=np.column_stack((self.total_lioness_network ,np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)))

        return self.total_lioness_network

    def save_lioness_results(self, file='lioness.txt'):
        '''Write lioness results to file.'''
        #self.lioness_network.to_csv(file, index=False, header=False, sep="\t")
        np.savetxt(file, self.total_lioness_network, delimiter="\t",header="")
        return None
