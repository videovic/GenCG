## arch--copula--graphon--v1
## Here to develop a network to a target level of assortativity
## gives different network structures in the copula graphon/
    ## copula density graphon / tensor product framework

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(igraph)
library(transport)
library(tidyverse)
library(kableExtra)
library(viridis)
library(ggplot2)
library(xtable)
library(cubature)
library(gridExtra)
library(cowplot)
library(latex2exp)
library(RColorBrewer)

set.seed(100)

###---------------------------
outputs_folder <- "./outputs/"
# Check if "outputs" folder exists
if (!dir.exists(outputs_folder)) {
  # If not, create it
  dir.create(outputs_folder)
  print("Folder 'outputs' created.")
} else {
  print("Folder 'outputs' already exists.")
}

###---------------------------
#### FUNCTIONS 
source("functions/graphon.R") # gives graphon simulations
source("functions/homfun.R") # gives functions for homorphism counts / densities
source("functions/copula.R") # reads in the generator function/ copulas for some Archimedean copulas # used to be called gen
source("functions/fun.R") # gives functions for network properties

#### ANALYSIS
source("analysis.R") # analysis on graphons, graphon tensor - produces graphs in paper
source("analysis3b.R") 

### EXAMPLE
source("examples.R") # produces additional graphs and checks for the thesis

