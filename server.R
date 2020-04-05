library(shiny)
library(msaR)
library(ape)
library(phangorn)
library(ggtree)
library(dplyr)
library(tidyr)
library(ggplot2)
library(alignfigR)

dna_colors <- c("A"="#A2FA8C", "C"="#FCCE8A", "G"="#F38D8A", "T"="#8AB8F5", "-" = "grey85")
outgroup <- "Ring-tailed_lemur"

######################### Parsimony for Phylogeny panel ############################
#pars_choice <- sample(1:2,1)  ### There are two most parsimonious trees. Students will randomly get one or the other (for some reason phangorn only ever gives one, need to understand phangorn better and then can estimate parsimony here.)
primates_pars_file1   <- paste0("data/primates_parsimony_tree_rooted_v1.tre")
primates_pars_file2   <- paste0("data/primates_parsimony_tree_rooted_v2.tre")

primates_parstree1 <- ape::read.tree(primates_pars_file1)
primates_parstree1_rooted <- ape::root(primates_parstree1, outgroup = outgroup, resolve.root=TRUE)
primates_parstree2 <- ape::read.tree(primates_pars_file2)
primates_parstree2_rooted <- ape::root(primates_parstree2, outgroup = outgroup, resolve.root=TRUE)

primates_pscore          <- phangorn::fitch(primates_parstree1, primates_phydat)
##################################################################################


################################# Ancestral panel ####################################
primates_pars_file_anc   <- "data/primates_parsimony_tree_labeled.tre"
primates_parstree_anc    <- ape::read.tree(primates_pars_file_anc)
primates_parstree_rooted <- ape::root(primates_parstree_anc, outgroup = outgroup, resolve.root=TRUE)

primate_msa_file_anc              <- "data/primates_ancestors_alignment.fasta"
primates_alignment_with_ancestors <- alignfigR::read_alignment(primate_msa_file_anc)
as_tibble(primates_alignment_with_ancestors) %>% 
  mutate(column=1:n()) %>% 
  gather(taxa, character, `Ring-tailed_lemur`:N11) %>%
  select(taxa, column, character) -> primates_ancestors_sequences
##################################################################################




shinyServer(function(input, output) {
   
  ##################################################################################
  ##################################### MSA ########################################
  ##################################################################################
  output$primate_msa1 <- renderMsaR({
    msaR(primate_msa_file, colorscheme = "nucleotide", labelNameLength = 150, labelid = FALSE, menu = FALSE)
  })
  output$primate_msa2 <- renderMsaR({
    msaR(primate_msa_file, colorscheme = "nucleotide", overviewbox = FALSE, labelNameLength = 150, labelid = FALSE, menu = FALSE)
  })
  ##################################################################################
  
  
  ##################################################################################
  ################################# TREE SEARCH ####################################
  ##################################################################################
  best_pscore <- reactiveVal(100000)  
  best_tree   <- reactiveVal(NULL)
  last_pscore  <- reactiveVal(100000)
  last_tree  <- reactiveVal(NULL)

  
  current_tree <- eventReactive(input$update_tree, {
    search_type <- isolate(input$search_type)
    if (input$update_tree == 1){
        new_tree <- ape::as.phylo(ape::rtree(primates_ntaxa, primates_names, rooted = TRUE))  
    } else
    { 
        if (search_type == "best"){
            all_nni <- phangorn::nni(best_tree())
        } else {
            all_nni <- phangorn::nni(last_tree())
        }
        new_tree <- ape::as.phylo( all_nni[[ sample(1:length(all_nni), 1) ]] )
    }

    new_pscore <- phangorn::fitch(new_tree, primates_phydat)   
    

    outmain <- paste0("Search #", input$update_tree,":",new_pscore, "\n")
    if (new_pscore < best_pscore()) {
         outmain <- paste0(outmain,"New best tree!")
    } else if (new_pscore > best_pscore()){
         outmain <- paste0(outmain,"Worse than best.")
    }  else {
        outmain <- paste0(outmain,"Same as best.")
    }
    
    if (is.null(best_tree()) | new_pscore < best_pscore())
    {
        best_pscore(new_pscore)
        best_tree(new_tree)
    }

    last_tree(new_tree)
    last_pscore(new_pscore) 
    out <- list("tree" = new_tree, "pscore" = new_pscore, "main" = outmain)

    if (input$update_tree == 1){
        out$main <- paste0("Score:",new_pscore, "\nInitial Tree")
    }
    out
  })
  
  
  
  output$tree_finding <- renderPlot({
    
    par(mfrow=c(1,2))
    
    tr <- current_tree()$tree
    tr <- ape::chronos( ape::root(tr, outgroup = outgroup, resolve.root=TRUE) )
    plot(tr, font=2, cex=1.25, edge.width=1.75, main = current_tree()$main, cex.main=1.5)

    best_tr <- ape::chronos( ape::root(best_tree(), outgroup = outgroup, resolve.root=TRUE))
    plot(best_tr, font=2, cex=1.25, edge.width=1.75, main = paste("Score:", best_pscore(), "\nBest Tree Found"), cex.main=1.5)
    
  }) 
  
  #observeEvent(input$reveal_parsimony, {
    output$parsimony_tree <- renderPlot({
      
      par(mfrow=c(1,2))      
      plot(primates_parstree1_rooted, cex=1.5, font=2, edge.width=2, main = paste("Most parsimonious tree. Score:", primates_pscore), cex.main=1.5)
      plot(primates_parstree2_rooted, cex=1.5, font=2, edge.width=2, main = paste("Most parsimonious tree. Score:", primates_pscore), cex.main=1.5)
  
    }) 
 # })
  ##################################################################################
  

  
  ##################################################################################
  ################################# ANCESTRAL STATES ###############################
  ##################################################################################
  
  output$treesite <- renderPlot({
    
    primates_ancestors_sequences %>% filter(column==input$site) -> primates_column
    ggtree(primates_parstree_anc, size=1.5) %<+% primates_column +
          geom_tiplab(color="black", offset=0.2, size=10) +
          geom_tippoint(aes(fill = character), color = "black", size=10, shape=22) + 
          scale_fill_manual(values=dna_colors) + 
          guides(fill=FALSE) + 
          xlim_tree(15) -> primate_ggtree
    
    if (input$color_branches) {
        primate_ggtree <- primate_ggtree + 
                            aes(color = character) + 
                            scale_color_manual(values=dna_colors) + 
                            guides(color = FALSE)
    } 
    primate_ggtree
  })
  
  
  ##################################################################################
  
})