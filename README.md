[![Build Status](https://travis-ci.org/netZoo/netZooR.svg?branch=devel)](https://travis-ci.org/netZoo/netZooR)
[![codecov](https://codecov.io/gh/netZoo/netZooR/branch/devel/graph/badge.svg)](https://codecov.io/gh/netZoo/netZooR)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Linux[![LINUX](https://travis-ci-job-status.herokuapp.com/badge/netZoo/netZooR/devel/linux)](https://travis-ci.org/netZoo/netZooR)

Macos[![MAC](https://travis-ci-job-status.herokuapp.com/badge/netZoo/netZooR/devel/macos)](https://travis-ci.org/netZoo/netZooR)

## Description
netZooR is an R package to construct,analyse and plot gene regulatory networks.

## Features

netZooR currently integrates with:
* **PANDA**(Passing Attributes between Networks for Data Assimilation)[[Glass et al. 2013]](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832): construct gene regulatory network from gene expression data, protein-protein interaction data, and transcription factor binding motifs (TFBMs) data.

* **CONDOR**(COmplex Network Description Of Regulators)[[Platig et al. 2016]](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005033): analyse bipartite community structure analysis of biological networks.

* **LIONESS**(Linear Interpolation to Obtain Network Estimates for Single Samples)[[Kuijjer et al. 2019]](https://doi.org/10.1016/j.isci.2019.03.021): reconstruct single-sample gene regulatory networks.

* **ALPACA**(ALtered Partitions Across Community Architectures)[[Padi and Quackenbush 2018]](https://www.nature.com/articles/s41540-018-0052-5): compare two networks and identify changes in modular structure.

* **SAMBAR**(Subtyping Agglomerated Mutations By Annotation Relations)[[Kuijjer et al.]](https://www.nature.com/articles/s41416-018-0109-7): identify subtypes based on somatic mutation data.

* Source protein-protein interaction network from [STRINGdb](https://string-db.org/) based on a list of protein of interest.

* Plot one PANDA network in [Cytoscape](https://cytoscape.org/).

* Plot two differential PANDA networks in Cytoscape.

## Installation

```r
library(devtools)
devtools::install_github("netZoo/netZooR")`
library(netZooR)
```

## User guide & tutorials
Please refer to the top navigation bar **Get started** for basic usage and application cases.

## Development
Contributions are welcome! Instruction of how to contribute netZooR repository can be found [here](https://netzoo.github.io/contribute/contribute/). 
The master branch on github should always be in good shape, so please to pull request to the **devel** branch.

## Feedback/Issues
Please report any issue to the [issues page](https://github.com/netZoo/netZooR/issues).

## License
The software is free and is licensed under the GNU General License v3.0, see the file [LICENSE](LICENSE) for details.

## Code of conduct
Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
