library(shiny)
library(msaR)
library(ape)
library(phangorn)
library(ggtree)
library(dplyr)

dna_colors <- c("A"="dodgerblue", "C"="red", "G"="forestgreen", "T"="gold", "-" = "grey85")

shinyServer(function(input, output) {
   
  ##################################################################################
  ##################################### MSA ########################################
  ##################################################################################
  output$primate_msa <- renderMsaR({
    msaR(primate_msa_file, colorscheme = "nucleotide", labelNameLength = 150, labelid = FALSE, menu = FALSE, overviewbox = FALSE)
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
    improved_title             <- paste0("Search #", input$update_tree,": Tree improved!\nScore:",new_pscore, "\n")
    not_improved_title_force   <- paste0("Search #", input$update_tree,": Tree not improved.\nScore:",last_pscore(), "\n")
    not_improved_title_noforce <- paste0("Search #", input$update_tree,": Tree not improved.\nScore:",new_pscore, "\n")

    
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
    
    ggtree(primates_parstree_anc, aes(color = character), size=1.25) %<+% primates_column +
      geom_tiplab(color="black", offset=0.3) +
      geom_tippoint(aes(color = character), size=4, shape=15) + 
      scale_color_manual(values=dna_colors) + xlim_tree(12)
  })
  
  
  ##################################################################################
  
})
#' 
#' 
#' output_vector <- generate_tree_vis(sample_df = sample_df,
#'                                    #'                                    alignment = aln_path, tree = tree,
#'                                    #'                                    phy_mat = bears, pscore = TRUE,
#'                                    #'                                    lscore = TRUE)
#'                                    #' }
#'                                    
#'                                    generate_tree_vis <- function(sample_df, alignment, tree, phy_mat,
#'                                                                  pscore = FALSE, lscore = FALSE, random_tree = FALSE){
#'                                      
#'                                      vis_vec <- list()
#'                                      phy_mat <- phangorn::phyDat(phy_mat, levels = c(0, 1), type = "USER")
#'                                      sample_df <- check_subs(sample_df = sample_df, phy_mat = phy_mat)
#'                                      for (i in 1:nrow(sample_df)){
#'                                        if (random_tree == FALSE){
#'                                          chars <- c(sample_df$starting_val[i], sample_df$stop_val[i])
#'                                          tr <- ape::as.phylo(generate_tree_vec(phy_mat, sample_df$starting_val[i],
#'                                                                                sample_df$stop_val[i], tree))
#'                                          pl <- ggtree::msaplot(p=ggtree::ggtree(tr), fasta=alignment, window = chars,                                    width = .1, offset = 9 ) + ggtree::geom_tiplab() +
#'                                            ggplot2::ggtitle(paste0(chars[1],"\n",chars[2]))
#' 
#' 
#' 


