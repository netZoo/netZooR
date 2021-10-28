[![master](https://github.com/netZoo/netZooR/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/netZoo/netZooR/actions/workflows/main.yml)
[![devel](https://github.com/netZoo/netZooR/actions/workflows/main.yml/badge.svg?branch=devel)](https://github.com/netZoo/netZooR/actions/workflows/main.yml)
[![codecov](https://codecov.io/gh/netZoo/netZooR/branch/devel/graph/badge.svg)](https://codecov.io/gh/netZoo/netZooR)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<a href="https://netzoo.github.io/netZooR/"><img src="https://img.shields.io/badge/docs-passing-green"></a>
[![tutorials](https://img.shields.io/badge/netZooR-tutorials-9cf)](https://github.com/netZoo/netZooR/tree/master/vignettes)
[![netBooks](https://img.shields.io/badge/netZooR-netBooks-ff69b4)](http://netbooks.networkmedicine.org/user/marouenbg/notebooks/Welcome_to_netBooks.ipynb?)
[![discussions](https://img.shields.io/badge/netZooR-discussions-orange)](https://github.com/netZoo/netZooR/discussions)


netZooR is tested on: (OS: Ubuntu + Macos) X (Language: R v3.6 + R v4.0)

## Description
netZooR is an R package to reconstruct, analyse, and plot biological networks.

## Features

netZooR currently integrates:
* **PANDA** (Passing Attributes between Networks for Data Assimilation) [[Glass et al. 2013]](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832): constructs gene regulatory network from gene expression data, protein-protein interaction data, and transcription factor binding motifs (TFBMs) data.

* **CONDOR** (COmplex Network Description Of Regulators) [[Platig et al. 2016]](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005033): analyses bipartite community structure analysis of biological networks.

* **LIONESS** (Linear Interpolation to Obtain Network Estimates for Single Samples) [[Kuijjer et al. 2019]](https://doi.org/10.1016/j.isci.2019.03.021): reconstructs single-sample gene regulatory networks.

* **ALPACA** (ALtered Partitions Across Community Architectures) [[Padi and Quackenbush 2018]](https://www.nature.com/articles/s41540-018-0052-5): compares two networks and identify changes in modular structure.

* **SAMBAR** (Subtyping Agglomerated Mutations By Annotation Relations) [[Kuijjer et al.]](https://www.nature.com/articles/s41416-018-0109-7): identifies subtypes based on somatic mutation data.

* **MONSTER** (Modeling Network State Transitions from Expression and Regulatory data) [[Schlauch et al.]](https://doi.org/10.1186/s12918-017-0517-y): infers transcription factor drivers of cell state conditions at the gene regulatory network level.

* **OTTER** (Optimization to Estimate Regulation) [[Weighill et al.]](https://www.biorxiv.org/content/10.1101/2020.06.23.167999v2.abstract): models gene regulation estimation as a graph matching problem.

* **CRANE** (Constrained Random Alteration of Network Edges) [[Lim et al.]](https://www.biorxiv.org/content/10.1101/2020.07.12.198747v1): generates ensembles of gene regulatory networks to identify disease modules.

* **EGRET** (Estimating the Genetic Regulatory effects on TFs) [[Weighill et al.]]() In preparation: models individual-specific gene regulatory networks using their genetic variants.

* **YARN** (Yet Another RNa-seq package) [[Paulsson et al.]](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1847-x): combines quality control and preprocessing, gene filtering, and normalization steps to streamline the processing of large-scale gene expression data such as GTEx.

The package also integrates additional functions to:
* Source protein-protein interaction network from [STRINGdb](https://string-db.org/) based on a list of protein of interest.

* Plot one PANDA network in [Cytoscape](https://cytoscape.org/).

* Plot two differential PANDA networks in Cytoscape.

## Requirements, installation and basic configuration.

- netZooR is compatible with R (>= 3.3.3) including R (>= 4.0),  click [here](https://www.r-project.org/) for more installation details.

- To use PANDA and LIONESS, there are two options: 

  1. Use `panda.py()` and `lioness.py()` by invoking the respective Python implementations in [netZooPy]((https://github.com/netZoo/netZooPy/tree/netZoo)). Because the native R linear algebra libraries can be slow, this way is recommended for faster analysis. However, optimized parallel libraries can give reasonable run times (option ii). To invoke Python scripts, there are some requirements to meet before using netZooR:

     a) [**Python**](https://www.python.org/downloads/) (>= 3.5.0) installed;

     b) Python libraries [pandas](https://pandas.pydata.org/), [numpy](https://numpy.org/), and [scipy](https://www.scipy.org/) installed;

     c) Internet access as package `reticulate` will link the R wrapper to the Python scripts located [here](https://github.com/netZoo/netZooPy/tree/netZoo) for those two methods.

  2. Use `panda()` and `lioness()` for the pure R implementations of PANDA and LIONESS. To speed up the run time, it is highly recommended to install an optimized linear algebra library, particularly for Ubunutu. Macos generally comes with optimized linear algebra libraries. You can check the `BLAS/LAPACK` fields in `sessionInfo()` in your R console. Detailed instructions can be found [here](https://csantill.github.io/RPerformanceWBLAS/).

     :warning: However, we found that Intel MKL linear algebra library with R 4.0.3 on Ubuntu 18.04 gave inconsistent results for the multiplication of large matrices and the results of PANDA were inconsistent. Therefore, Intel MKL is not currently recommended. 


- Most of plotting function can be realized by functions in [igraph](https://igraph.org/redirect.html), which will be loaded with netZooR through `library(netZooR)`. Some plotting functions like `vis.panda.in.cytoscape()` and `vis.diff.panda.in.cytoscape()` are able to plot interactive PANDA networks in [Cytoscape](https://cytoscape.org/), but installation of Cytoscape is required before using these plotting functions. Also, please make sure that Cytoscape is open when these functions are called.

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

For more details please refer to the [documentation website](https://netzoo.github.io/netZooR/).

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
Or use `browseVignettes("netZooR")` after installing the package. Also check [netbooks](http://netbooks.networkmedicine.org) to go through the  tutorials on a Jupyter notebook cloud server.

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
devtools::install() # To install the dependecies
devtools::build() # To build the package
#devtools::build(vignettes = FALSE) # You can skip building the vignettes if you are not contributing a vignette

# CMD check, if passed all tests here, it means this package is ready to pull request to the devel branch. Otherwise, fix the bug before pulling request.
devtools::check()
#devtools::check(vignettes = FALSE) #You can skip building the vignettes if you are not contributing a vignette
```

The master branch on github should always be in good shape, so please to pull request to the **devel** branch.
If the contribution is specific to pandaR, please contribute to its seperate GitHub page by [pull request](https://github.com/jnpaulson/pandaR). 

## License
The software is free and is licensed under the GNU General License v3.0, see the file [LICENSE](LICENSE) for details.

## Code of conduct
Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

