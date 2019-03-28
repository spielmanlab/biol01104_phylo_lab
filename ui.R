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
        tabPanel("Character Matrix", 
                 br(),
                 div(style = "margin:1em;padding-top:0;font-size:12px;",
                     msaROutput("primate_msa", height="200%")
                 )
        ),
        tabPanel("Phylogeny",
                 br(),br(),
                 
                 fluidRow(
                        style = "max-height: 75vh; overflow-y: scroll;" ,
                 div(style="display: inline-block;vertical-align:top; width: 150px;",
                     actionButton("update_tree", "Update tree!")
                  ),
                 div(style="display: inline-block;vertical-align:top; width: 200px;",
                     radioButtons("update_tree_type", "How to update tree?", choices=c("New random tree", "Move one branch only"))
                  ),
                 div(style="display: inline-block;vertical-align:top; width: 225px;",
                     checkboxInput("force_improvement", label = "Force tree improvement", value = FALSE)
                  ),
                 div(style="display: inline-block;vertical-align:top; width: 300px;",
                     selectInput("outgroup_random", "Select outgroup:", primates_names, selected = "Ring-tailed_lemur")
                  ),     
                br(),br(),               
                 plotOutput("tree_finding"),
                 
                 div(style="display: inline-block;vertical-align:top; width: 300px;",
                     actionButton("reveal_parsimony", "Reveal the most parsimonious tree")
                  ),
                 div(style="display: inline-block;vertical-align:top; width: 200px;",
                     selectInput("outgroup_parsimony", "Select outgroup:", primates_names, selected = "Ring-tailed_lemur")
                  ),  
                 br(),
                 plotOutput("parsimony_tree")
                )
        ), 
        tabPanel("Character History",
                 br(),br(),
                numericInput("site", label = "Column in character matrix", value = 1, min=1, max = primates_seqlen),
                checkboxInput("color_branches", label = "Color branches by inferred ancestry?", value = TRUE),
                plotOutput("treesite", width="90%", height="300px"),
                img(src="./dna_legend.png", style="width:150px;", align="center")
                 #div(style = "width:100px",
        )
      )
    )
))
