# netZoo

An R package to integrate [pypanda](https://github.com/davidvi/pypanda)--Python implementation of PANDA and LIONESS, R implementation of [CONDOR](https://github.com/jplatig/condor), R implementation of [ALPACA](https://github.com/meghapadi/ALPACA), and realated downsteaming analysis.

**PANDA**(Passing Attributes between Networks for Data Assimilation) is a message-passing model to gene regulatory network reconstruction. It integrates multiple sources of biological data, including protein-protein interaction, gene expression, and sequence motif information, in order to reconstruct genome-wide, condition-specific regulatory networks.[[Glass et al. 2013]](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832)

**LIONESS**(Linear Interpolation to Obtain Network Estimates for Single Samples) is a method to estimate sample-specific regulatory networks by applying linear interpolation to the predictions made by existing aggregate network inference 		approaches.[[LIONESS arxiv paper]](https://arxiv.org/abs/1505.06440)

**CONDOR** (COmplex Network Description Of Regulators) implements methods for clustering biapartite networks
and estimatiing the contribution of each node to its community's modularity.[[Platig et al. 2016]](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005033)

**ALPACA**(ALtered Partitions Across Community Architectures) is a method for comparing two genome-scale networks derived from different phenotypic states to identify condition-specific modules. [[Padi and Quackenbush 2018]](https://www.nature.com/articles/s41540-018-0052-5)



## Table of Contents
* [Getting Started](#getting-started) 
  * [Prerequisites](#prerequisites)
  * [Installing](#installing)
* [Data Resources](#Data-Resources)
  * [PANDA-ready Motif mapping data](#PANDA-ready-Motif-mapping-data)
  * [PPI](#PPI)
* [Running the sample datasets](#running-the-sample-datasets)
  * [PANDA and plot PANDA network](#panda-and-plot-panda-network)
  * [Cytoscape Plotting](#Cytoscape-Plotting)
  * [LIONESS and plot LIONESS network](#lioness-and-plot-lioness-network)
  * [CONDOR](#condor)
  * [ALPACA](#ALPACA)
* [Downstream Analysis](#Downstream-Analysis)
  * [Pre-ranked GSEA](#Pre-ranked-GSEA)
* [Further information](#further-information)
  * [Future Work](#future-work)
  * [Note](#note)
 


## Getting Started

### Prerequisites
Using this pacakage requires [**Python**](https://www.python.org/downloads/) 2.7, [**R**](https://cran.r-project.org/) (>= 3.3.3), [**Cytoscape**](https://cytoscape.org/) (>=3.6.1), [**GSEA**](http://software.broadinstitute.org/gsea/index.jsp), and **Internet access**.
There are also some Python packages required to apply Python implementation of PANDA and LIONESS.

How to install packages in different platforms could be find [here](https://packaging.python.org/tutorials/installing-packages/). 

The required Python packages are:
[pandas](https://pandas.pydata.org/), [numpy](http://www.numpy.org/), [networkx](https://networkx.github.io/), [matplotlib.pyplot](https://matplotlib.org/api/pyplot_api.html).

### Installing
This package could be downloaded via `install_github()` function from `devtools` package.

```R
install.packages("devtools")
library(devtools)
devtools::install_github("twangxxx/netZoo")

```
## Data Resources

### PANDA-ready Motif mappings data
Here is some [PANDA-ready Motif mappings data](https://sites.google.com/a/channing.harvard.edu/kimberlyglass/tools/resources) by different species.
### PPI
This package is included a function `sourcePPI` will create a Protein-Protein Interactions(PPI) througt STRING database among the proteins of interest. The [STRINGdb](http://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html) is already loaded with loading netZOO for further use.

```R
# TF is a data frame with single column filled with TFs of Mycobacterium tuberculosis H37Rv.
PPI <- sourcePPI(TF,species.index=83332, score_streshold=0)

```

## Running the sample datasets

Package [CONDOR](https://github.com/jplatig/condor), [igprah](http://igraph.org/r/), [viridisLite](https://cran.r-project.org/web/packages/viridisLite/index.html), and [ALPACA](https://github.com/meghapadi/ALPACA) are loaded with this package for further downstream analyses.

Use search() to check all loaded package currently.
```R
library(netZoo)
search()
```
Access help pages for usage of six core functions.
```
?runPanda
?plotPanda
?runLioness
?plotLioness
?runCondor
?runAlpaca
```
Use example datasets within package to test this package.

Refer to four input datasets files: one TB control expression dataset, one TB treated expression dataset, one motif sequence dataset, and one protein-protein interaction datasets in inst/extdat. All datasets are public data.

```R
treated_expression_file_path <- system.file("extdata", "expr4.txt", package = "netZoo", mustWork = TRUE)
control_expression_file_path <- system.file("extdata", "expr10.txt", package = "netZoo", mustWork = TRUE)
motif_file_path <- system.file("extdata", "chip.txt", package = "netZoo", mustWork = TRUE)
ppi_file_path <- system.file("extdata", "ppi.txt", package = "netZoo", mustWork = TRUE)
```

### PANDA and plot PANDA network

Assign the paths of treated expression dataset, motif dataset, and ppi dataset above to flag `e`(refers to "expression dataset"), `m`(refers to "=motif dataset"), and `ppi`(refers to "PPI" dataset) respectively. Then set option `rm_missing` to `TRUE` to run **PANDA** to generate an aggregate network for treated.
Repeat but alter the paths of treated expression dataset to control expression datasets to generate an aggregate network for control. 
vector `treated_all_panda_result` and vector `control_all_panda_result` below are large lists with three elements: the entire PANDA network, indegree ("to" nodes) nodes and score, outdegree ("from" nodes) nodes and score. Use `$panda`,`$indegree` and '$outdegree' to access each item resepctively.

```R
treated_all_panda_result <- runPanda(e = treated_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )
control_all_panda_result <- runPanda(e = control_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )
```
Use `$panda`to access the entire PANDA network.
```R
treated_net <- treated_all_panda_result$panda
```
Plot the 100 edge with the largest weight of two PANDA network. Besides, one message will be returned to indicate the location of output .png plot.
```R
plotPanda(top =100, file="treated_panda_100.png")
```
Repeat with networl of control. 
```R
control_net <- control_all_panda_result$panda
plotPanda(top =100, file="control_panda_100.png")
```

### Cytoscape Plotting
Cytoscape is an interactivity network visualization tool highly recommanded to explore the PANDA/LIONESS network. Before using this function `runCytoscapePlot`, please install and launch Cytoscape (3.6.1 or greater) and keep it running whenever using.

```R
runCytoscapePlot(control_net, top = 200, network.name="TB_control")
```

### LIONESS and plot LIONESS network
The method how to run LIONESS is mostly idential with method how to run PANDA in this package, unless the return values of `runLioness` is a data frame where first two columns represent TFs (regulators) and Genes (targets) while the rest columns represent each sample. each cell filled with estimated score calculated by LIONESS.

```R
treated_lioness <- runLioness(e = treated_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )
```
Plot LIONESS network should clarify which sample by 0-based index.

```R
plotLioness(col = 0, top = 100, file = "treat_lioness_sample1_100.png")
```

Repeat with control.
```R
control_lioness <- runLioness(e = control_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )
plotLioness(col = 0, top = 100, file = "control_lioness_sample1_100.png")
```

### CONDOR 

run CONDOR with a threshold to select edges. 
Defaults to `threshold` is the average of [median weight of non-prior edges] and [median weight of prior edges], all weights mentioned previous are transformationed with formula `w'=ln(e^w+1)` before calculating the median and average. But all the edges selected will remain the orginal weights calculated by PANDA before applying CONDOR.

```R
treated_condor_object <- runCondor(treated_net, threshold = 0)
control_condor_object <- runCondor(control_net, threshold = 0)
```

plot communities. package igraph and package viridisLite (a color map package) are already loaded with this package.

*treated network*:
```R
treated_color_num <- max(treated_condor_object$red.memb$com)
treated_color <- viridis(treated_color_num, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
condor.plot.communities(treated_condor_object, color_list=treated_color, point.size=0.04, xlab="Target", ylab="Regulator")
```

*control network*:
```R
control_color_num <- max(control_condor_object$red.memb$com)
control_color <- viridis(control_color_num, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
condor.plot.communities(control_condor_object, color_list=control_color , point.size=0.04, xlab="Target", ylab="Regulator")
```
### ALPACA

run LIONESS with two PANDA network above as first two arguments.

```R
alpaca_result<- runAlpaca(treated_net, control_net, "~/Desktop/TB", verbose=T)
```



## Downstream Analysis
### Pre-ranked GSEA

## Further information

### Future Work
Use `vignette("condor")` to access the vignette page of `condor` package and `vignette("ALPACA")` to access the vignette page of `ALPACA` package for any downstream analyses.

### Note
If there is an error like `Error in fetch(key) : lazy-load database.rdb' is corrupt` when accessing the help pages of functions in this package after being loaded. It's [a limitation of base R](https://github.com/r-lib/devtools/issues/1660) and has not been solved yet. Restart R session and re-load this package will help.

