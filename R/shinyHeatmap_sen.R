#'Run an interactive ComplexHeatmap session
#'
#'Dependencies are:
#' library(shiny)
#' library(ComplexHeatmap)
#' library(circlize)
#' library(dplyr)
#' library(openxlsx)
#' library(edgeR)
#' library(stringr)
#' 
#'date 14/02/2023
#'
#'
#'@return Opens a shiny session
#'@export


shinyHeatmap_sen <- function(){
library(shiny)
options(shiny.maxRequestSize = 200*1024^2)
ui <- fluidPage(
  splitLayout(cellWidths = c('45%','65%'),
    mainPanel(
      splitLayout(cellWidths = 300,
        mainPanel(
          actionButton('click','Generate Heatmap'),
          fileInput('expression_table','Expression file\n(xlsx format)', accept = '.xlsx'),
      selectInput('normalization','Normalize the values with', 
                  choices = c('None','CPM','LogCPM','Z-score','Log only')),
      selectInput('color_scale','Choose your color scale', choices = c('Blue-White-Red',
                                                                       'Red-White-Blue',
                                                                       'Green-Black-Purple',
                                                                       'Purple-Black-Green')
      ),
      fileInput('metadata','Metadata file (xlsx format)', accept = '.xlsx'),
      textInput('annotcol', 'Annotation column\nof metadata file', value = NULL),
      selectInput('annot_row_or_col','Annotate rows or columns?', choices = c('column','row')),
      textAreaInput('annotcolors',HTML(paste0('Annotation colors (1 color per line)<br/># of colors must equal<br/> # of annotation categories')), value = NULL),

      textAreaInput('gene_list','List of genes to plot(1 gene per line)', value = NULL)

      
      
      ),
      mainPanel(
      textInput('heatmap_title','Heatmap title', value = ''),
      textInput('legend_title','Legend title', value = NULL),
      textInput('legend_annotation_title','Annotation title', value = 'Annotation'),
      sliderInput('colsize','Column labels size', min=0, max=20, step=1, value=8 ),
      sliderInput('colkm','Number of Column clusters', min=1, max=10, step=1, value=1 ),
      selectInput('coldend','Show column dendrogram', choices = c(TRUE,FALSE)),

      sliderInput('rowsize','Row labels size', min=0, max=20, step=1, value=8 ),
      sliderInput('rowkm','Number of Row clusters', min=1, max=10, step=1, value=1 ),
      selectInput('rowdend','Show row dendrogram', choices = c(TRUE,FALSE))
    )
    )
    ),
    
    mainPanel(
     plotOutput('heatmap', height= '600px')
    )
  )
  
)

server <- function(input, output, session) {
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(openxlsx)
  library(edgeR)
  library(stringr)
  
  #debugging
  # output$text <- reactive(paste0(!is.null(input$metadata), 
  #                                str(input$annotcol)))
  
  
  #output
  observeEvent(input$click, {
    output$heatmap <-  renderPlot( expr = {
      
      ####expression file####
      expression <-   read.xlsx(input$expression_table$datapath, rowNames = TRUE, colNames=TRUE)
      
      
      ####normalization type####
        if(input$normalization == 'CPM'){
          expression_normalized <- as.data.frame(cpm(expression))
        } else if(input$normalization == 'LogCPM'){
          expression_normalized <- as.data.frame(log(cpm(expression)+1))
        } else if(input$normalization == 'Z-score'){
          expression_normalized <- as.data.frame(t(scale(t(cpm(expression)))))
        } else if(input$normalization == 'None'){
          expression_normalized <- expression
        } else if(input$normalization == 'Log only'){
          expression_normalized <- as.data.frame(log(expression + 1))
        }
      
      ####legend title####
      if(input$normalization == 'CPM' & input$legend_title == ''){
        legend_title <- 'CPM\nvalues'
      } else if(input$normalization == 'LogCPM' & input$legend_title == ''){
        legend_title <- 'Log Normalized\nCPM values'
      } else if(input$normalization == 'Z-score' & input$legend_title == ''){
        legend_title <- 'Z-score\nvalues'
      } else if(input$normalization == 'None' & input$legend_title == ''){
        legend_title <- 'Raw\nvalues'
      } else if(input$normalization == 'Log only' & input$legend_title == ''){
        legend_title <- 'Log\nNormalized\nValues'
      } else {
        legend_title <- input$legend_title
      }
      
      
      ####metadata file####
        if(!is.null(input$metadata)){
        metadata <- read.xlsx(input$metadata$datapath, rowNames = TRUE, colNames=TRUE)
        } 
      
      ####Gene list####
        if(input$gene_list!=''){
          genes_list <- unlist(strsplit(input$gene_list, split = '\n'))
          expression_isolated <- expression_normalized %>% dplyr::filter(rownames(.) %in% genes_list)
        } else {
          expression_isolated <- expression_normalized
        }
      
      ####Expression color scale####
      
        if(input$color_scale == 'Blue-White-Red'){
          CS <-colorRamp2(c(min(expression_isolated),0,max(expression_isolated)),
                     c("steelblue3","white","firebrick3"))
        } else if(input$color_scale == 'Red-White-Blue'){
          CS <-colorRamp2(c(min(expression_isolated),0,max(expression_isolated)),
                     c("firebrick3","white","steelblue3"))
        } else if(input$color_scale == 'Green-Black-Purple'){
          CS <-colorRamp2(c(min(expression_isolated),0,max(expression_isolated)),
                     c("olivedrab3","black","mediumorchid3"))
        } else if(input$color_scale == 'Purple-Black-Green'){
          CS <-colorRamp2(c(min(expression_isolated),0,max(expression_isolated)),
                     c("mediumorchid3","black","olivedrab3"))
        }
      
    
      
      ####Plot heatmap####
    if(input$annotcol!=''){
    
      if(input$annotcolors!=''){
      annotcolors <- unlist(strsplit(input$annotcolors, split = '\n'))
      
        if(length(annotcolors)==length(levels(as.factor(metadata[,input$annotcol])))){
        ####Annotation colors####
        names(annotcolors)<-levels(as.factor(metadata[,input$annotcol]))
        ####Annotation list####
        annotlist <- list(annotcolors)
        names(annotlist) <- input$legend_annotation_title
        
        ####Annotation data.frame####
        annotdf <- data.frame(metadata[,input$annotcol]) 
        names(annotdf) <- input$legend_annotation_title
        
        Heatmap(as.matrix(expression_isolated),
            col = CS,
            row_names_gp = gpar(fontsize = input$rowsize),
            column_names_gp = gpar(fontsize = input$colsize),
            column_km = input$colkm,
            row_km= input$rowkm,
            show_row_dend = as.logical(input$rowdend),
            show_column_dend = as.logical(input$coldend),
            column_title = input$heatmap_title,
            heatmap_legend_param = list(title = legend_title),
            top_annotation = HeatmapAnnotation(df = annotdf, which = input$annot_row_or_col,
                                               col = annotlist )
            
        )
        }
      } else {
        
        ####Annotation data.frame####
        annotdf <- data.frame(metadata[,input$annotcol]) 
        names(annotdf) <- input$legend_annotation_title
        
        
        Heatmap(as.matrix(expression_isolated),
                col = CS,
                
                row_names_gp = gpar(fontsize = input$rowsize),
                column_names_gp = gpar(fontsize = input$colsize),
                column_km = input$colkm,
                row_km= input$rowkm,
                show_row_dend = as.logical(input$rowdend),
                show_column_dend = as.logical(input$coldend),
                column_title = input$heatmap_title,
                heatmap_legend_param = list(title = legend_title),
                top_annotation = HeatmapAnnotation(df = annotdf, which = input$annot_row_or_col)
        )
              
            }
      
    } else {
      Heatmap(as.matrix(expression_isolated), 
              col = CS, 
              row_names_gp = gpar(fontsize = input$rowsize),
              column_names_gp = gpar(fontsize = input$colsize),
              column_km = input$colkm,
              row_km= input$rowkm,
              show_row_dend = as.logical(input$rowdend),
              show_column_dend = as.logical(input$coldend),
              column_title = input$heatmap_title,
              heatmap_legend_param = list(title = legend_title)
              )
      
    }
    
  })
  })
    
}

shinyApp(ui, server)
}