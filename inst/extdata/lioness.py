import os, os.path, sys
import numpy as np
import pandas as pd
from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects

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
    
class Lioness(Panda):
    """
    Description:
        Using LIONESS to infer single-sample gene regulatory networks.
        1. Reading in PANDA network and preprocessed middle data
        2. Computing coexpression network
        3. Normalizing coexpression network
        4. Running PANDA algorithm
        5. Writing out LIONESS networks

    Inputs:
        Panda: Panda object.

    Methods:
        __init__              : Intialize instance of Lioness class and load data.
        __lioness_loop        : The LIONESS algorithm.
        save_lioness_results  : Saves LIONESS network.

    Example:
        from netZooPy.lioness.lioness import Lioness
        To run the Lioness algorithm for single sample networks, first run PANDA using the keep_expression_matrix flag, then use Lioness as follows:
        panda_obj = Panda('../../tests/ToyData/ToyExpressionData.txt', '../../tests/ToyData/ToyMotifData.txt', '../../tests/ToyData/ToyPPIData.txt', remove_missing=False, keep_expression_matrix=True)
        lioness_obj = Lioness(panda_obj)

        Save Lioness results:
        lioness_obj.save_lioness_results('Toy_Lioness.txt')
        Return a network plot for one of the Lioness single sample networks:
        plot = AnalyzeLioness(lioness_obj)
        plot.top_network_plot(column= 0, top=100, file='top_100_genes.png')

        Example lioness output:
        Sample1 Sample2 Sample3 Sample4
        -------------------------------
        -0.667452814003	-1.70433776179	-0.158129613892	-0.655795512803
        -0.843366539284	-0.733709815256	-0.84849895139	-0.915217389738
        3.23445386464	2.68888472802	3.35809757371	3.05297381396
        2.39500370135	1.84608635425	2.80179804094	2.67540878165
        -0.117475863987	0.494923925853	0.0518448588965	-0.0584810456421

        TF, Gene and Motif order is identical to the panda output file.

    Authors:
        Cho-Yi Chen, David Vi, Daniel Morgan

    Reference:
        Kuijjer, Marieke Lydia, et al. "Estimating sample-specific regulatory networks." Iscience 14 (2019): 226-240.
    """

    def __init__(
        self,
        obj,
        computing="cpu",
        precision="double",
        ncores=1,
        start=1,
        end=None,
        save_dir="lioness_output",
        save_fmt="npy",
        output="network",
        alpha=0.1,
    ):
        """
        Description:
            Initialize instance of Lioness class and load data.

        Inputs:
            obj             : PANDA object, generated with keep_expression_matrix=True.
            obj.motif_matrix: TF DNA motif binding data in tf-by-gene format.
                              If set to None, Lioness will be performed on gene coexpression network.
            computing       : 'cpu' uses Central Processing Unit (CPU) to run PANDA
                              'gpu' use the Graphical Processing Unit (GPU) to run PANDA
            precision       : 'double' computes the regulatory network in double precision (15 decimal digits).
                              'single' computes the regulatory network in single precision (7 decimal digits) which is fastaer, requires half the memory but less accurate.
            start           : Index of first sample to compute the network.
            end             : Index of last sample to compute the network.
            save_dir        : Directory to save the networks.
            save_fmt        : Save format.
                              '.npy': (Default) Numpy file.
                              '.txt': Text file.
                              '.mat': MATLAB file.
            output          : 'network' returns all networks in a single edge-by-sample matrix (lioness_obj.total_lioness_network). For large sample sizes, this variable requires large RAM memory.
                              'gene_targeting' returns gene targeting scores for all networks in a single gene-by-sample matrix (lioness_obj.total_lioness_network).
                              'tf_targeting' returns tf targeting scores for all networks in a single gene-by-sample matrix (lioness_obj.total_lioness_network).
            alpha            : learning rate, set to 0.1 by default but has to be changed manually to match the learning rate of the PANDA object.

        Output:
            export_lioness_results : Depeding on the output argument, this can be either all the lioness networks or their gene/tf targeting scores.
        """
        # Load data
        with Timer("Loading input data ..."):
            self.export_panda_results = obj.export_panda_results
            self.expression_matrix = obj.expression_matrix
            self.motif_matrix = obj.motif_matrix
            self.ppi_matrix = obj.ppi_matrix
            self.correlation_matrix = obj.correlation_matrix
            if precision == "single":
                self.correlation_matrix = np.float32(self.correlation_matrix)
                self.motif_matrix = np.float32(self.motif_matrix)
                self.ppi_matrix = np.float32(self.ppi_matrix)
            self.computing = computing
            self.alpha = alpha
            self.n_cores = int(ncores)
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
            del obj

        # Get sample range to iterate
        self.n_conditions = self.expression_matrix.shape[1]
        self.indexes = range(self.n_conditions)[
            start - 1 : end
        ]  # sample indexes to include
        print("Number of total samples:", self.n_conditions)
        print("Number of computed samples:", len(self.indexes))
        print("Number of parallel cores:", self.n_cores)

        # Create the output folder if not exists
        self.save_dir = save_dir
        self.save_fmt = save_fmt
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        # Run LIONESS
        if int(self.n_conditions) >= int(self.n_cores) and self.computing == "cpu":
            self.total_lioness_network = Parallel(n_jobs=self.n_cores)(
                self.__par_lioness_loop(i, output) for i in (self.indexes)
            )

        elif self.computing == "gpu":
            for i in self.indexes:
                self.total_lioness_network = self.__lioness_loop(i)
        #        # self.export_lioness_results = pd.DataFrame(self.total_lioness_network)

        # create result data frame
        if output == "network":
            if isinstance(origmotif, pd.DataFrame):
                total_tfs = tf_names * len(gene_names)
                total_genes = [i for i in gene_names for _ in range(len(tf_names))]
                indDF = pd.DataFrame([total_tfs, total_genes], index=["tf", "gene"])
                indDF = indDF.append(
                    pd.DataFrame(self.total_lioness_network)
                ).transpose()
            else:  # if equal to None to be specific
                total_genes1 = gene_names * len(gene_names)
                total_genes2 = [i for i in gene_names for _ in range(len(gene_names))]
                indDF = pd.DataFrame(
                    [total_genes1, total_genes2], index=["gene1", "gene2"]
                )
                indDF = indDF.append(
                    pd.DataFrame(self.total_lioness_network)
                ).transpose()
            self.export_lioness_results = indDF
        elif output == "gene_targeting":
            self.export_lioness_results = pd.DataFrame(
                self.total_lioness_network, columns=gene_names
            ).transpose()
        elif output == "tf_targeting":
            self.export_lioness_results = pd.DataFrame(
                self.total_lioness_network, columns=tf_names
            ).transpose()
        self.save_lioness_results()

    def __lioness_loop(self, i):
        """
        Description:
            Initialize instance of Lioness class and load data.

        Outputs:
            self.total_lioness_network: An edge-by-sample matrix containing sample-specific networks.
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

        with Timer(
            "Saving LIONESS network %d to %s using %s format:"
            % (i + 1, self.save_dir, self.save_fmt)
        ):
            path = os.path.join(self.save_dir, "lioness.%d.%s" % (i + 1, self.save_fmt))
            if self.save_fmt == "txt":
                np.savetxt(path, lioness_network)
            elif self.save_fmt == "npy":
                np.save(path, lioness_network)
            elif self.save_fmt == "mat":
                from scipy.io import savemat

                savemat(path, {"PredNet": lioness_network})
            else:
                print("Unknown format %s! Use npy format instead." % self.save_fmt)
                np.save(path, lioness_network)
        if self.computing == "gpu" and i == 0:
            self.total_lioness_network = np.fromstring(
                np.transpose(lioness_network).tostring(), dtype=lioness_network.dtype
            )
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
        """
        Description:
            Initialize instance of Lioness class and load data.

        Outputs:
            self.total_lioness_network: An edge-by-sample matrix containing sample-specific networks.
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
        #lioness_network = self.n_conditions * (self.network - subset_panda_network) + subset_panda_network

        lioness_network = (self.n_conditions * self.network) - (
            (self.n_conditions - 1) * subset_panda_network
        )

        with Timer(
            "Saving LIONESS network %d to %s using %s format:"
            % (i + 1, self.save_dir, self.save_fmt)
        ):
            path = os.path.join(self.save_dir, "lioness.%d.%s" % (i + 1, self.save_fmt))
            if self.save_fmt == "txt":
                np.savetxt(path, lioness_network)
            elif self.save_fmt == "npy":
                np.save(path, lioness_network)
            elif self.save_fmt == "mat":
                from scipy.io import savemat

                savemat(path, {"PredNet": lioness_network})
            else:
                print("Unknown format %s! File will not be saved." % self.save_fmt)
                # np.save(path, lioness_network)
        # if i == 0:
        # self.total_lioness_network = np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)
        # else:
        #    self.total_lioness_network=np.column_stack((self.total_lioness_network ,np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)))
        if output == "network":
            self.total_lioness_network = np.transpose(lioness_network).flatten()
        elif output == "gene_targeting":
            self.total_lioness_network = np.sum(lioness_network, axis=0)
        elif output == "tf_targeting":
            self.total_lioness_network = np.sum(lioness_network, axis=1)
        return self.total_lioness_network

    def save_lioness_results(self):
        """
        Description:
            Saves LIONESS network.

        Outputs:
            file: Path to save the network.
        """
        # self.lioness_network.to_csv(file, index=False, header=False, sep="\t")
        fullpath = os.path.join(self.save_dir, "lioness.%s" % (self.save_fmt))
        if self.save_fmt == "txt":
            np.savetxt(
                fullpath,
                np.transpose(self.total_lioness_network),
                delimiter="\t",
                header="",
            )
        elif self.save_fmt == "npy":
            np.save(fullpath, np.transpose(self.total_lioness_network))
        elif self.save_fmt == "mat":
            from scipy.io import savemat

            savemat(fullpath, np.transpose(self.total_lioness_network))
        return None

    def export_lioness_table(self, output_filename="lioness_table.txt", header=False):
        """
        Description:
            Saves LIONESS network with edge names. This saves a dataframe with the corresponding
            header and indexes.
            So far we
        Outputs:
            output_filename: Path to save the network. Specify relative path
            and format. Choose between .csv, .tsv and .txt. (Defaults to .lioness_table.txt))
        """
        df  = self.export_lioness_results
        df = df.sort_values(by=['tf','gene'])
        if output_filename.endswith("txt"):
            df.to_csv(output_filename, sep=" ", header=False, index = False)
        elif output_filename.endswith("csv"):
            df.to_csv(output_filename, sep=",", header=False, index = False)
        elif output_filename.endswith("tsv"):
            df.to_csv(output_filename, sep="\t", header=False, index = False)
