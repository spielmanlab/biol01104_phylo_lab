library(phangorn)
library(dplyr)
library(tidyr)
library(alignfigR)


primate_msa_file     <- "data/primate-mtDNA.fasta"
primates_phydat   <- phangorn::read.phyDat(primate_msa_file, format="fasta")
primates_ntaxa   <- length(primates_phydat)
primates_names   <- names(primates_phydat)
primates_seqlen  <- 898 ### hardcoded wompwomp


