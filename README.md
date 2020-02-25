
[![Build Status](https://travis-ci.org/netZoo/netZooR.svg?branch=devel)](https://travis-ci.org/netZoo/netZooR)
[![codecov](https://codecov.io/gh/netZoo/netZooR/branch/devel/graph/badge.svg)](https://codecov.io/gh/netZoo/netZooR)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Linux[![LINUX](https://travis-ci-job-status.herokuapp.com/badge/netZoo/netZooR/devel/linux)](https://travis-ci.org/netZoo/netZooR)

Macos[![MAC](https://travis-ci-job-status.herokuapp.com/badge/netZoo/netZooR/devel/macos)](https://travis-ci.org/netZoo/netZooR)

## Description
netZooR is an R package to construct, analyse and plot gene regulatory networks.

## Features

netZooR currently integrates with:
* **PANDA**(Passing Attributes between Networks for Data Assimilation)[[Glass et al. 2013]](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832): construct gene regulatory network from gene expression data, protein-protein interaction data, and transcription factor binding motifs (TFBMs) data.

* **CONDOR**(COmplex Network Description Of Regulators)[[Platig et al. 2016]](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005033): analyse bipartite community structure analysis of biological networks.

* **LIONESS**(Linear Interpolation to Obtain Network Estimates for Single Samples)[[Kuijjer et al. 2019]](https://doi.org/10.1016/j.isci.2019.03.021): reconstruct single-sample gene regulatory networks.

* **ALPACA**(ALtered Partitions Across Community Architectures)[[Padi and Quackenbush 2018]](https://www.nature.com/articles/s41540-018-0052-5): compare two networks and identify changes in modular structure.

* **SAMBAR**(Subtyping Agglomerated Mutations By Annotation Relations)[[Kuijjer et al.]](https://www.nature.com/articles/s41416-018-0109-7): identify subtypes based on somatic mutation data.


* **MONSTER**(Modeling Network State Transitions from Expression and Regulatory data)[[Schlauch et al.]](https://doi.org/10.1186/s12918-017-0517-y): infer transcription factor drivers of cell state conditions at the gene regulatory network level.



* Source protein-protein interaction network from [STRINGdb](https://string-db.org/) based on a list of protein of interest.

* Plot one PANDA network in [Cytoscape](https://cytoscape.org/).

* Plot two differential PANDA networks in Cytoscape.

## Requirement, installation and basic configuration.
Using this pacakage requires [**Python**](https://www.python.org/downloads/) (3.X) and some [Python libraries](#required-python-libraries), [**R**](https://cran.r-project.org/) (>= 3.3.3), and stable **Internet access**.

Some plotting functions will require the [**Cytoscape**](https://cytoscape.org/) installed.

```r
# install.packages("devtools") 
library(devtools)
# install netZooR pkg with vignettes, otherwise remove the "build_vignettes = TRUE" argument.
devtools::install_github("netZoo/netZooR", build_vignettes = TRUE)
library(netZooR)
```

For more details please refer to the [documentation website](https://netzoo.github.io/netZooR/).

This package will invoke the Python in R environment through reticulate package.
Configure which version of Python to use if necessary, here in netZooR, Python 3.X is required. More details can be found [here](https://cran.r-project.org/web/packages/reticulate/vignettes/versions.html).

```r
#check your Python configuration and the specific version of Python in use currently
py_config()

# reset to Python 3.X if necessary, like below:
# use_python("/usr/local/bin/python3")

```

## Issues

For **data.table** installation issue please refer to [issue #40](https://github.com/netZoo/netZooR/issues/40).

Please report any further issue to the [issues page](https://github.com/netZoo/netZooR/issues).


## Tutorials
Please refer to the top navigation bar **Articles/** for basic usage and application cases.
Or use `browseVignettes("netZooR")` after installing package.

## Contribution and Development
Contributions are welcome! Instruction of how to contribute netZooR repository can be found [here](https://netzoo.github.io/contribute/contribute/). 
After adding new features or debugging, please re-build package and run CDM check by using below codes before pulling request to the devel branch.
```r
library(devtools)
library(pkgdown)
# document the description of function
devtools::document()
# build vignettes
devtools::build_vignettes()
# build documentation website
pkgdown::build_site()
# devtools::build()

# CMD check, if passed all tests here, it means this package is ready to pull request to the devel branch. Otherwise, fix the bug before pulling request.
devtools::check()
```

The master branch on github should always be in good shape, so please to pull request to the **devel** branch.
If the contribution is specific to pandaR, please contribute to its seperate GitHub page by [pull request](https://github.com/jnpaulson/pandaR). 

## License
The software is free and is licensed under the GNU General License v3.0, see the file [LICENSE](LICENSE) for details.

## Code of conduct
Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
