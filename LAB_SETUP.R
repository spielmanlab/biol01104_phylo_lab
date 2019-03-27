all_packages <- c("shiny", 
                  "msaR", 
                  "ape", 
                  "phangorn", 
                  "dplyr", 
                  "tidyr", 
                  "ggplot2", 
                  "alignfigR")
                  "BiocManager")
install.packages(all_packages)
library(BiocManager)
BiocManager::install("ggtree", version = "3.8")
