library(phangorn)
library(dplyr)
library(tidyr)
library(alignfigR)


primate_msa_file     <- "data/primate-mtDNA.fasta"
pars_choice <- sample(1:2,1)  ### There are two most parsimonious trees. Students will randomly get one or the other (for some reason phangorn only ever gives one, need to understand phangorn better and then can estimate parsimony here.)
primates_pars_file   <- paste0("data/primates_parsimony_tree_rooted_v", pars_choice, ".tre")

primates_phydat   <- phangorn::read.phyDat(primate_msa_file, format="fasta")

primates_parstree <- ape::read.tree(primates_pars_file)
primates_parstree_rooted <- ape::root(primates_parstree, outgroup = "Ring-tailed_lemur", resolve.root=TRUE)
primates_pscore          <- phangorn::fitch(primates_parstree, primates_phydat)

primates_ntaxa   <- length(primates_phydat)
primates_names   <- names(primates_phydat)
primates_seqlen  <- 898 ### hardcoded wompwomp




##################### Ancestral data #################
primates_pars_file_anc   <- "data/primates_parsimony_tree_labeled.tre"
primates_parstree_anc    <- ape::read.tree(primates_pars_file_anc)
primates_parstree_rooted <- ape::root(primates_parstree_anc, outgroup = "Ring-tailed_lemur", resolve.root=TRUE)

primate_msa_file_anc              <- "data/primates_ancestors_alignment.fasta"
primates_alignment_with_ancestors <- read_alignment(primate_msa_file_anc)
as_tibble(primates_alignment_with_ancestors) %>% 
  mutate(column=1:n()) %>% 
  gather(taxa, character, `Ring-tailed_lemur`:N11) %>%
  select(taxa, column, character) -> primates_ancestors_sequences
