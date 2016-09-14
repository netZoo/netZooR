# library(bereR)
library(pandaR)
# library(bptools)
library(reshape2)
library(penalized)
library(Biobase)
library(org.Hs.eg.db)
library(foreach)
library(doParallel)
library(limma)
library(igraph)
library(ggrepel)

analysisCode <- sample(100000,1)

# Prepares analysis including data loading, cleaning, algorithm generation...
# source("./load_TM_data.R")

# Runs the analysis- Includes network inference, transition analysis, etc.  Main computation.
# source("./run_TM.R")

# Process Results- Generates tables, images, etc...
# source("./process_TM.R")

# Odyssey2
