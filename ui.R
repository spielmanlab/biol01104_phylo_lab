library(shiny)
library(msaR)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("BIOL01104 Phylogenetics Lab"),
  br(),

    # sidebarLayout(
    # sidebarPanel(
    #   helpText("This program analyzes a character matrix ('sequence alignment') of mitochondrial DNA from twelve primate species.")
    # ),
    

    fixedPanel(width = "75%", right = "15%", left="15%",
      tabsetPanel(
        tabPanel("Multiple sequence alignment",
            msaROutput("primate_msa1")
        ),
      
        tabPanel("Search for the most parsimonious phylogeny",
                 br(),br(),
                 
                 fluidRow(
                    style = "max-height:80vh; overflow-y: scroll;" ,
                    actionButton("update_tree", "Search for another tree"),
                    radioButtons("search_type", "How to search for a new tree?",
                                        c("Search from last tree found" = "last",
                                          "Search from best tree found" = "best")),
                    br(),br(),               
                 plotOutput("tree_finding", height = "500px")
                )
            ),
        tabPanel("Reveal the most parsimonious phylogeny!",
                 br(),
                 plotOutput("parsimony_tree")
        ), 
        tabPanel("Ancestral state reconstruction",
        fluidRow(
                    style = "max-height:80vh; overflow-y: scroll;" ,
                msaROutput("primate_msa2"),
                numericInput("site", label = "Column in character matrix", value = 1, min=1, max = primates_seqlen),
                checkboxInput("color_branches", label = "Color branches by inferred ancestry?", value = FALSE),
                plotOutput("treesite", width="60%", height="400px"),
                img(src="./dna_legend.png", style="width:200px;", align="left") 
        ) 
                              

        )
      )
    )
))
