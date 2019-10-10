library(shiny)
library(msaR)
library(ape)
library(phangorn)
library(ggtree)
library(dplyr)
library(tidyr)
library(ggplot2)
library(alignfigR)

dna_colors <- c("A"="dodgerblue", "C"="red", "G"="forestgreen", "T"="gold", "-" = "grey85")


######################### Parsimony for Phylogeny panel ############################
pars_choice <- sample(1:2,1)  ### There are two most parsimonious trees. Students will randomly get one or the other (for some reason phangorn only ever gives one, need to understand phangorn better and then can estimate parsimony here.)
primates_pars_file   <- paste0("data/primates_parsimony_tree_rooted_v", pars_choice, ".tre")

primates_parstree <- ape::read.tree(primates_pars_file)
primates_parstree_rooted <- ape::root(primates_parstree, outgroup = "Ring-tailed_lemur", resolve.root=TRUE)
primates_pscore          <- phangorn::fitch(primates_parstree, primates_phydat)
##################################################################################


################################# Ancestral panel ####################################
primates_pars_file_anc   <- "data/primates_parsimony_tree_labeled.tre"
primates_parstree_anc    <- ape::read.tree(primates_pars_file_anc)
primates_parstree_rooted <- ape::root(primates_parstree_anc, outgroup = "Ring-tailed_lemur", resolve.root=TRUE)

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
  output$primate_msa <- renderMsaR({
    msaR(primate_msa_file, colorscheme = "nucleotide", labelNameLength = 150, labelid = FALSE, menu = FALSE)
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
    update_tree_type <- isolate(input$update_tree_type)
    force_improvement <- isolate(input$force_improvement)

    if (update_tree_type == "New random tree")      new_tree <- ape::as.phylo(ape::rtree(primates_ntaxa, primates_names, rooted = TRUE))  
    if (update_tree_type == "Move one branch only") 
    {
        all_nni <- phangorn::nni(last_tree())
        new_tree <- ape::as.phylo( all_nni[[ sample(1:length(all_nni), 1) ]] )
    }
    new_pscore <- phangorn::fitch(new_tree, primates_phydat)   
    
    first_search_title         <- paste0("Initial Tree.\nScore:",new_pscore, "\n")
    improved_title             <- paste0("Search #", input$update_tree,": Tree improved from last search!\nScore:",new_pscore, "\n")
    not_improved_title_force   <- paste0("Search #", input$update_tree,": Tree not improved from last search.\nScore:",last_pscore(), "\n")
    not_improved_title_noforce <- paste0("Search #", input$update_tree,": Tree not improved from last search.\nScore:",new_pscore, "\n")

    
    if (is.null(best_tree()) | new_pscore < best_pscore())
    {
        best_pscore(new_pscore)
        best_tree(new_tree)
    }

    if (force_improvement)
    {   
        if (new_pscore < last_pscore()) {
            ### Force: Tree improved
            out <- list("tree" = new_tree, "pscore" = new_pscore, "main" = improved_title)
            last_tree(new_tree)
            last_pscore(new_pscore)
        } else 
        {   ### Force: Tree not improved
            out <- list("tree" = last_tree(), "pscore" = last_pscore(), "main" = not_improved_title_force)
        }
    } else
    {
        if (new_pscore < last_pscore()) 
        {
            outmain <- improved_title
        } else {
            outmain <- not_improved_title_noforce
        }  
        last_tree(new_tree)
        last_pscore(new_pscore) 
        out <- list("tree" = new_tree, "pscore" = new_pscore, "main" = outmain)
    }
    if (input$update_tree == 1){
        out$main <- first_search_title
    }
    out
  })
  
  
  
  output$tree_finding <- renderPlot({
    
    par(mfrow=c(1,2))
    
    tr <- current_tree()$tree
    tr <- ape::chronos( ape::root(tr, outgroup = input$outgroup_random, resolve.root=TRUE) )
    plot(tr, font=2, cex=1.25, edge.width=1.75, main = current_tree()$main, cex.main=1.5)

    best_tr <- ape::chronos( ape::root(best_tree(), outgroup = input$outgroup_random, resolve.root=TRUE))
    plot(best_tr, font=2, cex=1.25, edge.width=1.75, main = paste("Best Tree Found.\nScore:", best_pscore(), "\n"), cex.main=1.5)
    
  }) 
  
  observeEvent(input$reveal_parsimony, {
    output$parsimony_tree <- renderPlot({
      
      primates_parstree_rooted <- ape::root(primates_parstree, outgroup = input$outgroup_parsimony, resolve.root=TRUE)
      
      plot(primates_parstree_rooted, cex=1.5, font=2, edge.width=2, main = paste("Most parsimonious tree. Score:", primates_pscore), cex.main=1.5)
  
    }) 
  })
  ##################################################################################
  

  
  ##################################################################################
  ################################# ANCESTRAL STATES ###############################
  ##################################################################################
  
  output$treesite <- renderPlot({
    
    primates_ancestors_sequences %>% filter(column==input$site) -> primates_column
    
    if (input$color_branches) {
        ggtree(primates_parstree_anc, aes(color = character), size=1.25) %<+% primates_column +
          geom_tiplab(color="black", offset=0.2, size=6) +
          geom_tippoint(aes(color = character), size=4, shape=15) + 
          scale_color_manual(values=dna_colors) + xlim_tree(12)
    } else{
        ggtree(primates_parstree_anc, size=1.25) %<+% primates_column +
          geom_tiplab(color="black", offset=0.2, size=6) +
          geom_tippoint(aes(color = character), size=4, shape=15) + 
          scale_color_manual(values=dna_colors) + xlim_tree(12)    
    }
  })
  
  
  ##################################################################################
  
})