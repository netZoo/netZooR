import os, os.path, sys
import numpy as np
import pandas as pd

from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects

import math
import numpy as np
from scipy.stats import zscore
import sys

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
    """ Using LIONESS to infer single-sample gene regulatory networks.
        1. Reading in PANDA network and preprocessed middle data
        2. Computing coexpression network
        3. Normalizing coexpression network
        4. Running PANDA algorithm
        5. Writing out LIONESS networks

    Parameters
    ----------

            obj             : object
                PANDA object, generated with keep_expression_matrix=True.
            computing       : str
                'cpu' uses Central Processing Unit (CPU) to run PANDA
                'gpu' use the Graphical Processing Unit (GPU) to run PANDA
            precision       : str
                'double' computes the regulatory network in double precision (15 decimal digits).
                'single' computes the regulatory network in single precision (7 decimal digits) which is fastaer, requires half the memory but less accurate.
            subset_numbers  : list
                List of sample index onto which lioness should be run. ([1,10,20])
            subset_names  : list
                List of sample names onto which lioness should be run. (['s1','s2','s3'])
            start           : int
                Index of first sample to compute the network. If subset_numbers or subset_names is passed, 
                this is ignored
            end             : int
                Index of last sample to compute the network.If subset_numbers or subset_names is passed, 
                this is ignored
            all_background : bool
                Pass the flag if you want to keep the whole samples as background
            save_dir        : str
                Directory to save the networks.
            save_fmt        : str
                Save format.
                - '.npy': (Default) Numpy file of the network.
                - '.txt': Text file, only values are saved, no tf or gene names. Will be deprecated.
                - '.csv': text file with index (tf) and column (gene) names
                - '.h5': hdf file, fastest way to save the lioness dataframe with index/column names
                - '.mat': MATLAB file.
            output          : str
                - 'network' returns all networks in a single edge-by-sample matrix (lioness_obj.total_lioness_network is the unlabeled variable and lioness_obj.export_lioness_results is the row-labeled variable). For large sample sizes, this variable requires large RAM memory.
                - 'gene_targeting' returns gene targeting scores for all networks in a single gene-by-sample matrix (lioness_obj.total_lioness_network).
                - 'tf_targeting' returns tf targeting scores for all networks in a single gene-by-sample matrix (lioness_obj.total_lioness_network).
            alpha            : float
                learning rate, set to 0.1 by default but has to be changed manually to match the learning rate of the PANDA object.
            save_single: bool
                when set to True it will save each lioness network with its sample name inside the lioness output folder
            export_filename: str
                if passed, the final lioness table will be saved with all tf-gene edges as dataframe index and 
                samples as column name
            ignore_final: bool
                if True, no lioness network is kept in memory. This requires saving single networks at each step
    Returns
    --------
    export_lioness_results : _
        Depeding on the output argument, this can be either all the lioness 
        networks or their gene/tf targeting scores.

    Notes
    -------
    Example on how to use Lioness and plot the network

        >>> from netZooPy.lioness.lioness import Lioness
        >>> #To run the Lioness algorithm for single sample networks, first run PANDA using the keep_expression_matrix flag, then use Lioness as follows:
        >>> panda_obj = Panda('../../tests/ToyData/ToyExpressionData.txt', '../../tests/ToyData/ToyMotifData.txt', '../../tests/ToyData/ToyPPIData.txt', remove_missing=False, keep_expression_matrix=True)
        >>> lioness_obj = Lioness(panda_obj)

        >>> #Save Lioness results:
        >>> lioness_obj.save_lioness_results('Toy_Lioness.txt')
        >>> #Return a network plot for one of the Lioness single sample networks:
        >>> plot = AnalyzeLioness(lioness_obj)
        >>> plot.top_network_plot(column= 0, top=100, file='top_100_genes.png')

    Example lioness output:
        TF, Gene and Motif order is identical to the panda output file.

        - Sample1 Sample2 Sample3 Sample4\n
        - -------------------------------\n
        - -0.667452814003	-1.70433776179	-0.158129613892	-0.655795512803\n
        - -0.843366539284	-0.733709815256	-0.84849895139	-0.915217389738\n
        - 3.23445386464	2.68888472802	3.35809757371	3.05297381396\n
        - 2.39500370135	1.84608635425	2.80179804094	2.67540878165\n
        - -0.117475863987	0.494923925853	0.0518448588965	-0.0584810456421\n

        

    References
    -----------

    .. [1]__ Kuijjer, Marieke Lydia, et al. "Estimating sample-specific regulatory networks." 
        Iscience 14 (2019): 226-240.
    
    Authors: Cho-Yi Chen, David Vi, Daniel Morgan
    """

    def __init__(
        self,
        obj,
        computing="cpu",
        precision="double",
        ncores=1,
        start=1,
        end=None,
        subset_numbers=None,
        subset_names=None,
        save_dir="lioness_output",
        save_fmt="npy",
        output="network",
        alpha=0.1,
        save_single = False,
        export_filename = None,
        ignore_final=False, 
    ):
        """ Initialize instance of Lioness class and load data.
        """
        # Load data
        with Timer("Loading input data ..."):
            self.export_panda_results = obj.export_panda_results
            self.expression_samples = obj.expression_samples
            if precision == "single":
                self.expression_matrix = np.float32(obj.expression_matrix)
                self.correlation_matrix = np.float32(obj.correlation_matrix)
                self.motif_matrix = np.float32(obj.motif_matrix)
                self.ppi_matrix = np.float32(obj.ppi_matrix)
                self.alpha = np.float32(alpha)
            else:
                self.expression_matrix = obj.expression_matrix
                self.motif_matrix = obj.motif_matrix
                self.ppi_matrix = obj.ppi_matrix
                self.correlation_matrix = obj.correlation_matrix
                self.alpha = alpha
                
            self.computing = computing
            self.n_cores = int(ncores)
            self.save_single = save_single
            self.precision = precision
            if precision == "single":
                self.np_dtype = np.float32
            else:
                self.np_dtype = np.float64
            if hasattr(obj, "panda_network"):
                self.network = obj.panda_network.to_numpy()
            elif hasattr(obj, "puma_network"):
                self.network = obj.puma_network
            else:
                print("Cannot find panda or puma network in object")
                raise AttributeError("Cannot find panda or puma network in object")
            gene_names = obj.gene_names
            tf_names = obj.unique_tfs
            origmotif = obj.motif_data  # save state of original motif matrix
            self.gene_names = gene_names
            self.tf_names = tf_names
            del obj

        # Get sample range to iterate
        # the number of conditions is the N parameter used for the number of samples in the whole background
        self.n_conditions = self.expression_matrix.shape[1]
        self.n_lio_samples = self.n_conditions
        if (subset_numbers!=None or subset_names!=None):
            if (subset_numbers!=None and subset_names!=None):
                sys.exit('Pass only one between subset_numbers and subset_names')
            elif (subset_numbers!=None and subset_names==None):
                # select using indexes
                assert isinstance(subset_numbers, list)
                assert isinstance(subset_numbers[0], int)
                self.indexes = [int(i) for i in subset_numbers]
            else:
                #select using sample names
                assert isinstance(subset_names, list)
                assert isinstance(subset_names[0], int)
                self.indexes = [self.expression_samples.index(int(i)) for i in subset_names]
            #self.expression_samples = self.expression_samples[self.indexes]
            # number of lioness networks to be computed
            self.n_lio_samples = len(self.indexes)
        else:
            # if no subset is selected, we just use the start and end numbers to decide
            # which samples need to be analyses. The background is always what is used for PANDA
            # and stays the same
            
            self.indexes = range(self.n_conditions)[
                start - 1 : end
            ]  # sample indexes to include
            #self.expression_samples = self.expression_samples[start-1:end]
            self.n_lio_samples = len(self.indexes)
            
        print("Number of total samples:", self.n_conditions)
        print("Number of computed samples:", len(self.indexes))
        print("Number of parallel cores:", self.n_cores)

        # Create the output folder if not exists
        self.save_dir = save_dir
        self.save_fmt = save_fmt
        # if true, no final large matrix is saved
        self.ignore_final=ignore_final
        
        # Here we make sure that at least one between the complete lioness (edgex x samples) 
        # dataframe or each single network is saved
        if ((self.save_single==False) and (self.ignore_final==True)):
            sys.exit("ERROR: you are passing save_single=False, and ignore_final=True. This way no output is saved, change one of the two to have either the single networks or the final large table")
        # We create the folder
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        #############
        # Run LIONESS
        #############
        if int(self.n_conditions) >= int(self.n_cores) and self.computing == "cpu":
            # the total_lioness_network here is a list of 1d
            # arrays (network:(tfXgene,1),gene_targeting:(gene,1),tf_targeting:(tf,1))
            self.total_lioness_network = Parallel(n_jobs=self.n_cores)(
                self.__par_lioness_loop(i, output) for i in (self.indexes)
            ) 

        elif self.computing == "gpu":
            for i in self.indexes:
                self.total_lioness_network = self.__lioness_loop(i)
            self.total_lioness_network = self.total_lioness_network.T
        # create result data frame
        if self.ignore_final:
            print('WARNING: we do not keep all lionesses in memory. They have been saved singularly.')
        else:
            if output == "network":
                if isinstance(origmotif, pd.DataFrame):
                    # get row of all TFs
                    total_tfs = tf_names * len(gene_names)
                    # get row of all genes
                    total_genes = [i for i in gene_names for _ in range(len(tf_names))]
                    # first dataframe is made of tf and gene names
                    indDF = pd.DataFrame([total_tfs, total_genes], index=["tf", "gene"])
                    # concatenate with dataframe of data, rows are samples, columns the edges
                    indDF = indDF.append(
                        pd.DataFrame(self.total_lioness_network, index = self.expression_samples[self.indexes])
                    ).transpose()
                else:  # if equal to None to be specific
                    total_genes1 = gene_names * len(gene_names)
                    total_genes2 = [i for i in gene_names for _ in range(len(gene_names))]
                    indDF = pd.DataFrame(
                        [total_genes1, total_genes2], index=["gene1", "gene2"]
                    )
                    indDF = indDF.append(
                        pd.DataFrame(self.total_lioness_network, index = self.expression_samples[self.indexes])
                    ).transpose()
                
                # keep the df as the export results
                self.export_lioness_results = indDF
                del indDF
            elif output == "gene_targeting":
                self.export_lioness_results = pd.DataFrame(
                    self.total_lioness_network, columns=gene_names, index =self.expression_samples[self.indexes]
                ).transpose()
            elif output == "tf_targeting":
                self.export_lioness_results = pd.DataFrame(
                    self.total_lioness_network, columns=tf_names, index = self.expression_samples[self.indexes]
                ).transpose()
            
            # if export filename is passed, the full lioness table is saved
            if export_filename:
                self.export_lioness_table(output_filename = export_filename)
            else:
                self.save_lioness_results()

    def __lioness_loop(self, i):
        #TODO: this is now for GPU only in practice
        """ Initialize instance of Lioness class and load data.

        Returns
        --------
            self.total_lioness_network: array
                An edge-by-sample matrix containing sample-specific networks.
        """
        # for i in self.indexes:
        print("Running LIONESS for sample %d/%d:" %((i),(self.n_conditions)))
        idx = [x for x in range(self.n_conditions) if x != i]  # all samples except i
        with Timer("Computing coexpression network:"):
            if self.computing == "gpu":
                import cupy as cp
                
                correlation_matrix_cp = cp.corrcoef(self.expression_matrix[:, idx].astype(self.np_dtype)).astype(self.np_dtype)
                if cp.isnan(correlation_matrix_cp).any():
                    cp.fill_diagonal(correlation_matrix_cp, 1)
                    correlation_matrix_cp = cp.nan_to_num(correlation_matrix_cp)
                correlation_matrix = cp.asnumpy(correlation_matrix_cp)
                del correlation_matrix_cp
                cp._default_memory_pool.free_all_blocks()
            else:
                correlation_matrix = np.corrcoef(self.expression_matrix[:, idx])
                if np.isnan(correlation_matrix).any():
                    np.fill_diagonal(correlation_matrix, 1)
                    correlation_matrix = np.nan_to_num(correlation_matrix)

        with Timer("Normalizing networks:"):
            correlation_matrix_orig = (
                correlation_matrix  # save matrix before normalization
            )
            correlation_matrix = self._normalize_network(correlation_matrix)

        with Timer("Inferring LIONESS network:"):
            if self.motif_matrix is not None:
                del correlation_matrix_orig
                subset_panda_network = compute_panda(
                    correlation_matrix,
                    np.copy(self.ppi_matrix),
                    np.copy(self.motif_matrix),
                    computing = self.computing,
                    alpha = self.alpha,
                )
            else:
                del correlation_matrix
                subset_panda_network = correlation_matrix_orig

        # For consistency with R, we are using the N panda_all - (N-1) panda_all_but_q
        lioness_network = (self.n_conditions * self.network) - (
            (self.n_conditions - 1) * subset_panda_network
        )
        # old
        #lioness_network = self.n_conditions * (self.network - subset_panda_network) + subset_panda_network

        if self.save_single:
            with Timer(
                "Saving LIONESS network %d (%s) to %s using %s format:"
                % (i, self.expression_samples[i], self.save_dir, self.save_fmt)
            ):
                
                path = os.path.join(self.save_dir, "lioness.%s.%s.%s" % (self.expression_samples[i],str(i),self.save_fmt))
                self.__lioness_to_disk(lioness_network, path)

        if self.ignore_final:
            return(np.array([0]))
        else:
            if self.computing == "gpu" and i == 0:
                self.total_lioness_network = np.fromstring(
                    np.transpose(lioness_network).tostring(), dtype=lioness_network.dtype
                )[:,np.newaxis]
                
            elif self.computing == "gpu" and i != 0:
                self.total_lioness_network = np.column_stack(
                    (
                        self.total_lioness_network,
                        np.fromstring(
                            np.transpose(lioness_network).tostring(),
                            dtype=lioness_network.dtype,
                        ),
                    )
                )
            
            return self.total_lioness_network

    @delayed
    @wrap_non_picklable_objects
    def __par_lioness_loop(self, i, output):
        """ Initialize instance of Lioness class and load data.

        Returns
        ---------
            self.total_lioness_network: array
                An edge-by-sample matrix containing sample-specific networks.
        """
        # for i in self.indexes:
        print("Running LIONESS for sample %d:" % (i + 1))
        idx = [x for x in range(self.n_conditions) if x != i]  # all samples except i
        with Timer("Computing coexpression network:"):
            if self.computing == "gpu":
                import cupy as cp

                correlation_matrix = cp.corrcoef(self.expression_matrix[:, idx])
                if cp.isnan(correlation_matrix).any():
                    cp.fill_diagonal(correlation_matrix, 1)
                    correlation_matrix = cp.nan_to_num(correlation_matrix)
                correlation_matrix = cp.asnumpy(correlation_matrix)
            else:
                # run on CPU with numpy
                correlation_matrix = np.corrcoef(self.expression_matrix[:, idx])
                if np.isnan(correlation_matrix).any():
                    np.fill_diagonal(correlation_matrix, 1)
                    correlation_matrix = np.nan_to_num(correlation_matrix)

        with Timer("Normalizing networks:"):
            correlation_matrix_orig = (
                correlation_matrix  # save matrix before normalization
            )
            correlation_matrix = self._normalize_network(correlation_matrix)

        with Timer("Inferring LIONESS network:"):
            # TODO: fix this correlation matrix+delete
            if self.motif_matrix is not None:
                del correlation_matrix_orig
                subset_panda_network = compute_panda(
                    correlation_matrix,
                    np.copy(self.ppi_matrix),
                    np.copy(self.motif_matrix),
                    computing = self.computing,
                    alpha = self.alpha,
                )
            else:
                del correlation_matrix
                subset_panda_network = correlation_matrix_orig

        # For consistency with R, we are using the N panda_all - (N-1) panda_all_but_q
        #lioness_network = self.n_conditions * (self.network - subset_panda_network) + subset_panda_network

        lioness_network = (self.n_conditions * self.network) - (
            (self.n_conditions - 1) * subset_panda_network
        )
        # the lioness network here is a TFxGENE np array
        
        # if save_single flag is passed, save each single lioness sample
        if self.save_single:
            # TODO: here we need to decide whether to add the tf and gene name
            with Timer(
                "Saving LIONESS network %d (%s) to %s using %s format:"
                % (i, self.expression_samples[i], self.save_dir, self.save_fmt)
            ):
                
                path = os.path.join(self.save_dir, "lioness.%s.%s.%s" % (self.expression_samples[i],str(i),self.save_fmt))
                print(path)
                self.__lioness_to_disk(lioness_network, path)
        
        # TODO: why this? Should we remove it?
        # if i == 0:
        # self.total_lioness_network = np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)
        # else:
        #    self.total_lioness_network=np.column_stack((self.total_lioness_network ,np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)))
        if self.ignore_final:
            return(np.array([0]))
        else:
            if output == "network":
                self.total_lioness_network = np.transpose(lioness_network).flatten()
            elif output == "gene_targeting":
                self.total_lioness_network = np.sum(lioness_network, axis=0)
            elif output == "tf_targeting":
                self.total_lioness_network = np.sum(lioness_network, axis=1)
            return self.total_lioness_network

    def save_lioness_results(self, lioness_output_filename = None):
        """ Saves LIONESS network.
            Uses self.save_fmt, self.save_dir to save the data
            into self.total_lioness_network
        """
        # self.lioness_network.to_csv(file, index=False, header=False, sep="\t")
        if lioness_output_filename:
            fullpath = lioness_output_filename
        else:
            fullpath = os.path.join(self.save_dir, "lioness.%s" % (self.save_fmt))

        if fullpath.endswith("txt"):
            np.savetxt(
                fullpath,
                np.transpose(self.total_lioness_network),
                delimiter="\t",
            )
        elif fullpath.endswith("csv"):
            np.savetxt(
                fullpath,
                np.transpose(self.total_lioness_network),
                delimiter=",",
            )
        elif fullpath.endswith("npy"):
            np.save(fullpath, np.transpose(self.total_lioness_network))
        elif fullpath.endswith("h5"):
            print('Warning: saving as txt for now. Need to rewrite this')
            np.savetxt(fullpath,
                            np.transpose(self.total_lioness_network),
                            delimiter="\t",
                        )
        elif fullpath.endswith("mat"):
            from scipy.io import savemat
            mdic = {"results": np.transpose(self.total_lioness_network), "label": "lioness"}
            savemat(fullpath, mdic)
        else:
            print('Trying to save lioness output. Format %s not recognised' %str(fullpath))
        return None

    def __lioness_to_disk(self, net, path):
        if self.save_fmt == "txt":
            np.savetxt(path, net)
        elif self.save_fmt == "h5":
            pd.DataFrame(data = net, columns=self.gene_names, index = self.tf_names).to_hdf(path, key = 'lioness' ,mode='w')
        elif self.save_fmt == 'csv':
            pd.DataFrame(data = net, columns=self.gene_names, index = self.tf_names).to_csv(path)
        elif self.save_fmt == "npy":
            np.save(path, net)
        elif self.save_fmt == "mat":
            from scipy.io import savemat
            savemat(path, {"PredNet": net})
        else:
            print("Unknown format %s! Use npy format instead." % self.save_fmt)
            np.save(path, net)

    def export_lioness_table(self, output_filename="lioness_table.txt", header=False, output = 'network'):
        """ 
            Saves LIONESS network with edge names. This saves a dataframe with the corresponding
            header and indexes.

        Parameters
        ------------
            output_filename: str
                Path to save the network. Specify relative path
                and format. Choose between .csv, .tsv and .txt. 
                (Defaults to .lioness_table.txt))
        """
        # TODO: add case where there is tf_targeting or gene_targeting
        if (output=='network'):
            # we first get the names of first two columns (tf,gene) or (gene1,gene2)
            sort_cols = self.export_lioness_results.columns.tolist()[:2]
            if output_filename.endswith("txt"):
                self.export_lioness_results.sort_values(by=sort_cols).to_csv(output_filename, sep=" ", index = False)
            elif output_filename.endswith("csv"):
                self.export_lioness_results.sort_values(by=sort_cols).to_csv(output_filename, sep=",", index = False)
            elif output_filename.endswith("tsv"):
                self.export_lioness_results.sort_values(by=sort_cols).to_csv(output_filename, sep="\t", index = False)
            else:
                sys.exit('Export output unknown: use txt/csv/tsv')

        else:
            # otherwise we only need one column and we sort by index
            if output_filename.endswith("txt"):
                self.export_lioness_results.sort_index().to_csv(output_filename, sep=" ")
            elif output_filename.endswith("csv"):
                self.export_lioness_results.sort_index().to_csv(output_filename, sep=",")
            elif output_filename.endswith("tsv"):
                self.export_lioness_results.sort_index().to_csv(output_filename, sep="\t")
            else:
                sys.exit('Export output unknown: use txt/csv/tsv')
