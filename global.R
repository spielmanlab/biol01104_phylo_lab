#library(phangorn)


read_alignment <- function(file){
  raw_data <- readLines( file, warn = FALSE )
  seq_vector <- c()
  seq_name <- ""
  for (line in raw_data){
    # New sequence record? Reset numbering
    if ( grepl("^>", line) ){
      seq_name <- sub("^>", "", line)
      seq_vector[seq_name] <- ""
    }
    else {
      temp_seq <- gsub(" ","",line)
      temp_seq <- gsub("\n","",temp_seq)
      seq_vector[seq_name] <- paste( seq_vector[seq_name], temp_seq, sep="" )
    }
  }
  # Is this an alignment?
  seq_list <- strsplit(seq_vector, split = "")
  lengths <- sapply(seq_list, length)
  if ( sum(lengths != lengths[1]) != 0 )
    stop("Your provided file is not an alignment. Please provide an alignment file in FASTA format to use alignfigR.")
  # Return sequence data parsed into named list
  seq_list 
}


primate_msa_file     <- "data/primate-mtDNA.fasta"
primates_phydat   <- phangorn::read.phyDat(primate_msa_file, format="fasta")
primates_ntaxa   <- length(primates_phydat)
primates_names   <- names(primates_phydat)
primates_seqlen  <- 898 ### hardcoded wompwomp


