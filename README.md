[![master](https://github.com/netZoo/netZooR/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/netZoo/netZooR/actions/workflows/main.yml)
[![devel](https://github.com/netZoo/netZooR/actions/workflows/main.yml/badge.svg?branch=devel)](https://github.com/netZoo/netZooR/actions/workflows/main.yml)
[![codecov](https://codecov.io/gh/netZoo/netZooR/branch/devel/graph/badge.svg)](https://codecov.io/gh/netZoo/netZooR)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<a href="https://netzoo.github.io/netZooR/"><img src="https://img.shields.io/badge/docs-passing-green"></a>
[![tutorials](https://img.shields.io/badge/netZooR-tutorials-9cf)](https://github.com/netZoo/netZooR/tree/master/vignettes)
[![Netbooks](https://img.shields.io/badge/netZooR-Netbooks-ff69b4)](http://netbooks.networkmedicine.org/user/marouenbg/notebooks/Welcome_to_netBooks.ipynb?)
[![discussions](https://img.shields.io/badge/netZooR-discussions-orange)](https://github.com/netZoo/netZooR/discussions)
[![DOI](https://zenodo.org/badge/190646802.svg)](https://zenodo.org/badge/latestdoi/190646802)

netZooR is tested on: (OS: Ubuntu + Macos) X (Language: R v4.2)

## Description
netZooR is an R package to reconstruct, analyse, and plot biological networks.

## Features

netZooR currently integrates:

<details>
<summary>PANDA</summary>
<b>PANDA</b> (Passing Attributes between Networks for Data Assimilation) <a href="http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832">Glass et al. 2013</a>: PANDA is a method for estimating bipartite gene regulatory networks (GRNs) consisting of two types of nodes: transcription factors (TFs) and genes. An edge between TF $i$ and gene $j$ indicates that gene $j$ is regulated by TF $i$. The edge weight represents the strength of evidence for this regulatory relationship obtained by integrating three types of biological data: gene expression data, protein-protein interaction (PPI) data, and transcription factor binding motif (TFBM) data. PANDA is an iterative approach that begins with a seed GRN estimated from TFBMs and uses message passing between data types to refine the seed network to a final GRN that is consistent with the information contained in gene expression, PPI, and TFBM data. 
</details>

<details>
<summary>CONDOR</summary>
<b>CONDOR</b> (COmplex Network Description Of Regulators) <a href="http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005033">Platig et al. 2016</a>: CONDOR is a tool for community detection in bipartite networks. Many community detection methods for unipartite networks are based on the concept of maximizing a modularity metric that compares the weight of edges within communities to the weight of edges between communities, prioritizing community assignments with higher values of the former relative to the latter. CONDOR extends this concept to bipartite networks by optimizing a bipartite version of modularity defined by <a href="https://pubmed.ncbi.nlm.nih.gov/18233893/">Barber (2007)</a>. To enable bipartite community detection on large networks such gene regulatory networks, CONDOR uses a fast unipartite modularity maximization method on one of the two unipartite projections of the bipartite network.  In Platig et al. (2016), CONDOR is applied to bipartite networks of single nucleotide polymorphisms (SNPs) and gene expression, where a network edge from a SNP node to a gene node is indicative of an association between the SNP and the gene expression level, commonly known as an expression quantitative trait locus (eQTL). Communities detected with CONDOR contained local hub nodes ("core SNPs") enriched for association with disease, suggesting that functional eQTL relationships are encoded at the community level.
</details>

<details>
<summary>LIONESS</summary>
<b>LIONESS</b> (Linear Interpolation to Obtain Network Estimates for Single Samples) <a href="https://doi.org/10.1016/j.isci.2019.03.021">Kuijjer et al. 2019</a>: LIONESS is a flexible method for single-sample network integration. The machinery behind LIONESS is a leave-one-out approach. To construct a single-sample network for sample $i$, a first network is estimated on the full dataset and a second network is estimated on the dataset with sample $i$ withheld. The single-sample network is then estimated based on the difference between these two networks. Any method that can be used to estimate a network can be used with LIONESS to estimate single-sample networks. Two common use cases are the use of LIONESS to generate single-sample GRNs based on PANDA and the use of LIONESS to generate single-sample Pearson correlation networks.
</details>

<details>
<summary>ALPACA</summary>
<b>ALPACA</b> (ALtered Partitions Across Community Architectures) <a href="https://www.nature.com/articles/s41540-018-0052-5">Padi and Quackenbush 2018</a>: ALPACA is a method for differential network analysis that is based on a novel approach to comparison of network community structures. Comparisons of community structure have typically been accomplished by assessing which nodes switch community membership between networks ("community comparison") or by computing the edge weight differences by subtracting the adjacency matrices of two networks and then performing community detection on the resulting differential network ("edge subtraction"). Both these approaches have important limitations. Community comparison is subject to a resolution limit and cannot detect differences smaller than the average community size in a network. Edge subtraction transfers noise from both of the original networks to the differential network, leading to an imprecise estimator. Moreover, positive and negative edge differences cannot be distinguished in the subsequent community detection performed on the differential network. 

  In contrast to community comparison and edge subtraction, ALPACA compares the community structure of two networks by optimizing a new metric: "differential modularity". In the ALPACA algorithm, one network is defined as the reference network and the second is defined as the perturbed network. The differential modularity metric measures the extent to which edges in a community in the perturbed network differ from those that would be expected by random chance according to a null distribution based on the reference network. Community structure of the perturbed network is determined by maximizing this differential modularity. The resulting communities are "differential modules" that show how the perturbed network differs from the reference network at the community level.
</details>

<details>
<summary>SAMBAR</summary>
<b>SAMBAR</b> (Subtyping Agglomerated Mutations By Annotation Relations) <a href="https://www.nature.com/articles/s41416-018-0109-7">Kuijjer et al.</a>: SAMBAR is a tool for studying cancer subtypes based on patterns of somatic mutations in curated biological pathways. Rather than characterize cancer according to mutations at the gene level, SAMBAR agglomerates mutations within pathways to define a pathway mutation score. To avoid bias based on pathway representation, these pathway mutation scores correct for the number of genes in each pathway as well as the number of times each gene is represented in the universe of pathways. By taking a pathway rather than gene-by-gene lens, SAMBAR both de-sparsifies somatic mutation data and incorporates important prior biological knowledge. Kuijjer et al. (2018) demonstrate that SAMBAR is capable of outperforming other methods for cancer subtyping, producing subtypes with greater between-subtype distances; the authors use SAMBAR for a pan-cancer subtyping analysis that identifies four diverse pan-cancer subtypes linked to distinct molecular processes. 
</details>

<details>
<summary>MONSTER</summary>
<b>MONSTER</b> (Modeling Network State Transitions from Expression and Regulatory data) <a href="https://doi.org/10.1186/s12918-017-0517-y">Schlauch et al.</a>: MONSTER is a method for estimating transitions between network states by modeling the adjacency matrix of one state as a linear transformation of the adjacency matrix of another. Like LIONESS, MONSTER is a flexible method that does not require a particular type of network structure. MONSTER models the perturbation of an initial network A into a perturbed network B according to a matrix product B = AT. T is a transition matrix encoding the changes that map A to B. When A and B are gene regulatory networks, i.e., bipartite networks between TFs and genes, the MONSTER framework leads naturally to the definition of TF involvement as the sum of the off-diagonal weights for a transcription factor $i$ in the transition matrix T. This perspective enables MONSTER to identify differentially involved TFs that contribute to network transitions differently between different conditions. This dimension cannot be captured from a traditional differential expression analysis of TFs, which will not detect TFs that have the same concentration between conditions.
</details>

<details>
<summary>OTTER</summary>
<b>OTTER</b> (Optimization to Estimate Regulation) <a href="https://www.biorxiv.org/content/10.1101/2020.06.23.167999v2.abstract">Weighill et al.</a>: OTTER is a GRN inference method based on the idea that observed biological data (PPI data and gene co-expression data) are projections of a bipartite GRN between TFs and genes. Specifically, PPI data represent the projection of the GRN onto the TF-TF space and gene co-expression data represent the projection of the GRN onto the gene-gene space. OTTER reframes the problem of GRN inference as a problem of relaxed graph matching and finds a GRN that has optimal agreement with the observed PPI and coexpression data. The OTTER objective function is tunable in two ways: first, one can prioritize matching the PPI data or the coexpression data more heavily depending on one's confidence in the data source; second, there is a regularization parameter that can be applied to induce sparsity on the estimated GRN. The OTTER objective function can be solved using spectral decomposition techniques and gradient descent; the latter is shown to be closely related to the PANDA message-passing approach (Glass et al. 2013).
</details>

<details>
<summary>CRANE</summary>
<b>CRANE</b> (Constrained Random Alteration of Network Edges) <a href="https://doi.org/10.3389/fgene.2020.603264">Lim et al.</a>: CRANE is a method for determining statistical significance of structural differences between networks.  Analysis with CRANE is a four-phase process. The first step of CRANE is to estimate two networks: a reference network and a perturbed network. In the same spirit as LIONESS, CRANE is flexible: any network inference method (e.g., correlation, partial correlation, PANDA) can be used at this stage. In the second step, differential features are determined by comparing the reference and perturbed networks. Here, CRANE is again flexible: such differential features could arise from simple measures such as a comparison of node degree or centrality, or from more nuanced techniques such as differential module detection with ALPACA. Third, a large number of constrained random networks are developed based on the network structure of the reference network. By comparing each random network with the original reference network, a set of null differential measures is obtained. Fourth, the observed differential features from step two can be compared with the null distribution from step three to generate empirical p-values. A typical workflow for applying CRANE in NetZooR would involve fitting PANDA networks in step one and using ALPACA to estimate differential modules in step two. 
</details>

<details>
<summary>EGRET</summary>
<b>EGRET</b> (Estimating the Genetic Regulatory effects on TFs) <a href="https://www.genome.org/cgi/doi/10.1101/gr.275107.120">Weighill et al.</a>: EGRET incorporates genetic variants as a fourth data type in the PANDA message-passing framework, enabling the estimation of genotype-specific GRNs. Genetic variants can alter transcription factor binding by affecting the composition of motif sites on the DNA. Not every genetic variant has such an affect; EGRET incorporates only genetic variants which have (1) been shown to be associated with gene expression (expression quantitative trait loci, or eQTL), and (2) are predicted to affect transcription factor binding based on a tool called QBiC (Martin et al. 2019). This information is used in combination with TFBM predictions as input to the PANDA message-passing framework. The resulting EGRET network is a genotype-specific bipartite GRN that is similar to a PANDA network but incorporates the information contained by individual genetic variation.
</details>

<details>
<summary>YARN</summary>
<b>YARN</b> (Yet Another RNa-seq package) <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1847-x">Paulsson et al.</a>: YARN is a package that combines quality control, gene filtering, and normalization steps to streamline the preprocessing of large-scale, multi-tissue gene expression data from resources such as the Genotype-Tissue Expression (GTEx) project.  Among other steps, YARN uses principal coordinate analysis (PCoA) to determine if samples collected from different sites on the same tissue (for example, transverse and sigmoid colon) can be treated as "transcriptionally indistinguishable" and grouped together to increase power for downstream analyses. Paulsson et al. (2017) demonstrate the use of YARN to develop a pan-cancer RNA-seq dataset for 30,333 genes from 9435 samples across 38 tissues from the GTEx dataset.
</details>

<details>
<summary>PUMA</summary>
<b>PUMA</b> (PANDA Using MicroRNA Associations) <a href="https://www.sciencedirect.com/science/article/pii/S2589004219300872">Kuijjer et al.</a> extends the PANDA framework to model how microRNAs (miRNAs) participate in gene regulatory networks. PUMA networks are bipartite networks that consist of a regulatory layer and a layer of genes being regulated, similar to PANDA networks. While the regulatory layer of PANDA networks consists only of transcription factors (TFs), the regulatory layer of PUMA networks consists of both TFs and miRNAs. A PUMA network is seeded using a combination of input data sources such as motif scans or ChIP-seq data (for TF-gene edges) and an miRNA target prediction tool such as TargetScan or miRanda (for miRNA-gene edges). PUMA uses a message passing framework similar to PANDA to integrate this prior information with gene-gene coexpression and protein-protein interactions to estimate a final regulatory network incorporating miRNAs. Kuijjer and colleagues [7] apply PUMA to 38 GTEx tissues and demonstrate that PUMA can identify important patterns in tissue-specific regulation of genes by miRNA.
</details>

<details>
<summary>SPIDER</summary>
<b>SPIDER</b> (Seeding PANDA Interactions to Derive Epigenetic Regulation) <a href="https://www.nature.com/articles/s41540-021-00208-3">Sonawane et al.</a> extends the PANDA framework by incorporating DNase-Seq data to account for chromatin state for the prediction of TF binding sites. The method consists of processing DNase-Seq data to find open chromatin regions and build a “mask” matrix that is then overlaid on the TF-gene motif network to select binding sites that are available fro TF binding. This method can be applied for various biological contexts such as cell lines and tissues. Sonawane and colleagues have employed their method to model cell- type specific GRNs using DNase-Seq data from ENCODE and showed that integrating epigenetic data in SPIDER networks allows building more accurate networks.
</details>

<details>
<summary>DRAGON</summary>
<b>DRAGON</b> (Determining Regulatory Associations using Graphical models on Omics Networks) <a href="https://arxiv.org/abs/2104.01690">Shutta, Weighill et al.</a> is a method for estimating multiomic Gaussian graphical models (GGMs, also known as partial correlation networks) that incorporate two different omics data types. DRAGON builds off of the popular covariance shrinkage method of Ledoit and Wolf with an optimization approach that explicitly accounts for the differences in two separate omics "layers" in the shrinkage estimator. The resulting sparse covariance matrix is then inverted to obtain a precision matrix estimate and a corresponding GGM.  Although GGMs assume normally distributed data, DRAGON can be used on any type of continuous data by transforming data to approximate normality prior to network estimation. Currently, DRAGON can be applied to estimate networks with two different types of omics data. Investigators interested in applying DRAGON to more than two types of omics data can consider estimating pairwise networks and "chaining" them together.
</details>

<details>
<summary>TIGER</summary>
<b>TIGER</b> (Transcription Inference using Gene Expression and Regulatory Data) <a href="https://www.biorxiv.org/content/10.1101/2022.12.12.520141v1">Chen et al.</a> is a Bayesian matrix factorization framework that combines prior TF binding knowledge, such as from the DoRothEA database, with gene expression data from experiments. It estimates individual-level TF activities (TFA) and context-specific gene regulatory networks (GRN). Unlike other methods, TIGER can flexibly model activation and inhibition events, prioritize essential edges, shrink irrelevant edges towards zero using a sparse Bayesian prior, and simultaneously estimate TF activity levels and the underlying regulatory network. It is important to note that TIGER works most appropriately with large sample size datasets like TCGA to include a wide range of TFs due to its lower rank constraint.
</details>

<details>
<summary>COBRA</summary>
<b>COBRA</b> (Co-expression Batch Reduction Adjustment) Micheletti, Schlauch et al. is method for correction of high-order batch effects such as those that persist in co-expression networks. Batch effects and other covariates are known to induce spurious associations in co-expression networks and confound differential gene expression analyses. These effects are corrected for using various methods prior to downstream analyses such as the inference of co-expression networks and computing differences between them. In differential co-expression analysis, the pairwise joint distribution of genes is considered rather than independently analyzing the distribution of expression levels for each individual gene. Computing co-expression matrices after standard batch correction on gene expression data is not sufficient to account for the possibility of batch-induced changes in the correlation between genes as existing batch correction methods act solely on the marginal distribution of each gene. Consequently, uncorrected, artifactual differential co-expression can skew the correlation structure such that network-based methods that use gene co-expression can produce false, nonbiological associations even using data corrected using standard batch correction. Co-expression Batch Reduction Adjustment (COBRA) addresses this question by computing a batch-corrected gene co-expression matrix based on estimating a conditional covariance matrix. COBRA estimates a reduced set of parameters that express the co-expression matrix as a function of the sample covariates and can be used to control for continuous and categorical covariates. The method is computationally fast and makes use of the inherently modular structure of genomic data to estimate accurate gene regulatory associations and enable functional analysis for high-dimensional genomic data.
</details>

* Source protein-protein interaction network from [STRINGdb](https://string-db.org/) based on a list of protein of interest.

* Plot one PANDA network in [Cytoscape](https://cytoscape.org/).

* Plot two differential PANDA networks in Cytoscape.

## Requirements, installation and basic configuration.

- netZooR is compatible with R (>= 4.1),  click [here](https://www.r-project.org/) for more installation details.

- To use PANDA and LIONESS, there are two options: 

  1. Use `panda.py()` and `lioness.py()` by invoking the respective Python implementations in [netZooPy]((https://github.com/netZoo/netZooPy/tree/netZoo)). Because the native R linear algebra libraries can be slow, this way is recommended for faster analysis. However, optimized parallel libraries can give reasonable run times (option ii). To invoke Python scripts, there are some requirements to meet before using netZooR:

     a) [**Python**](https://www.python.org/downloads/) (>= 3.5.0) installed;

     b) Python libraries [pandas](https://pandas.pydata.org/), [numpy](https://numpy.org/), and [scipy](https://www.scipy.org/) installed;

     c) Internet access as package `reticulate` will link the R wrapper to the Python scripts located [here](https://github.com/netZoo/netZooPy/tree/netZoo) for those two methods.

  2. Use `panda()` and `lioness()` for the pure R implementations of PANDA and LIONESS. To speed up the run time, it is highly recommended to install an optimized linear algebra library, particularly for Ubunutu. Macos generally comes with optimized linear algebra libraries. You can check the `BLAS/LAPACK` fields in `sessionInfo()` in your R console. Detailed instructions can be found [here](https://csantill.github.io/RPerformanceWBLAS/).

     :warning: However, we found that Intel MKL linear algebra library with R 4.0.3 on Ubuntu 18.04 gave inconsistent results for the multiplication of large matrices and the results of PANDA were inconsistent. Therefore, Intel MKL is not currently recommended. 


- Most of plotting function can be realized by functions in [igraph](https://igraph.org/redirect.html), which will be loaded with netZooR through `library(netZooR)`. Some plotting functions like `vis.panda.in.cytoscape()` and `vis.diff.panda.in.cytoscape()` are able to plot interactive PANDA networks in [Cytoscape](https://cytoscape.org/), but installation of Cytoscape is required before using these plotting functions. Also, please make sure that Cytoscape is open when these functions are called.

### Installation

#### Using devtools/remotes

```r
# install.packages("devtools") 
# install netZooR pkg with vignettes, otherwise remove the "build_vignettes = TRUE" argument.
devtools::install_github("netZoo/netZooR", build_vignettes = TRUE)
library(netZooR)
```
You can use `remotes` instead of `devtools` because it is faster to install and run. The synatx is the following:

```r
# install.packages("remotes") 
# install netZooR pkg with vignettes, otherwise remove the "build_vignettes = TRUE" argument.
remotes::install_github("netZoo/netZooR", build_vignettes = TRUE)
library(netZooR)
```

#### Using Bioconductor

netZooR is also available through [Bioconductor](https://bioconductor.org/packages/netZooR)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("netZooR")
```

For more details please refer to the [documentation website](https://netzoo.github.io/netZooR/).

#### Using bioconda

netZooR is also available through [Bioconda](https://bioconda.github.io/recipes/bioconductor-netzoor/README.html#package-bioconductor-netzoor)

```bash
conda install -c bioconda bioconductor-netzoor
```

### Python binding

This package will invoke Python programming language in R environment through [reticulate](https://rstudio.github.io/reticulate/) package, by default setting there is no additional configuration needed.
Configuring which version of Python to use , here in netZooR, Python 3.X is required. More details can be found [here](https://cran.r-project.org/web/packages/reticulate/vignettes/versions.html).

```r
#check your Python configuration and the specific version of Python in use currently
py_config()

# reset to Python 3.X if necessary, like below:
# use_python("/usr/local/bin/python3")

```

## Help

If you need help or if you have any question about netZoo, feel free to start with [discussions](https://github.com/netZoo/netZooR/discussions).
To report a bug, please open a new [issue](https://github.com/netZoo/netZooR/issues).


## Tutorials
For more details please refer to the [documentation website](https://netzoo.github.io/netZooR/). Tutorials are available at the top navigation bar **Articles/** for basic usage and application cases.
Or use `browseVignettes("netZooR")` after installing the package. [Netbooks](http://netbooks.networkmedicine.org) deploys tutorials on a Jupyter notebook cloud server to get you running without any installation.

## Contribution and Development
Contributions are welcome! The contribution guide to netZooR can be found [here](https://netzoo.github.io/contribute/contribute/). 

We follow the [Bioconductor code guidelines](https://bioconductor.org/packages/devel/bioc/vignettes/BiocCheck/inst/doc/BiocCheck.html). Before pushing a contribution, please run

```r
library(BiocCheck)
BiocCheck("packageDirOrTarball")
```
And resolve any warnings, notes, and errors before committing the code.

After adding new features or optimizing a function in the package, please re-build the package and run `R CMD check .` in the terminal or follow the instructions below before doing the pull request to the devel branch.
To run only the tests:
```r
# install.packages('rcmdcheck')
# setwd('path/to/netZooR/root') # Set the working directory to netZooR root
rcmdcheck::rcmdcheck(args = c("--no-manual","--ignore-vignettes"), error_on = "error", build_args="--no-build-vignettes")
```
To rebuild vignettes, documentation, and tests:
```r
# document the description of function
# setwd('path/to/netZooR/root') # Set the working directory to netZooR root
devtools::document()
# build vignettes
devtools::build_vignettes() # You can skip building the vignettes if you are not contributing a vignette
# build documentation website
pkgdown::build_site(examples=FALSE)

# Install and build the package using devtools
devtools::install() # To install the dependencies
devtools::build() # To build the package
#devtools::build(vignettes = FALSE) # You can skip building the vignettes if you are not contributing a vignette

# CMD check, if passed all tests here, it means this package is ready to pull request to the devel branch. Otherwise, fix the bug before pulling request.
devtools::check()
#devtools::check(vignettes = FALSE) #You can skip building the vignettes if you are not contributing a vignette
```

The master branch on github should always be in good shape, so please to pull request to the **devel** branch.
If the contribution is specific to pandaR, please contribute to its separate GitHub page by [pull request](https://github.com/jnpaulson/pandaR). 

## License
The software is free and is licensed under the GNU General License v3.0, see the file [LICENSE](LICENSE) for details.

## Code of conduct
Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

