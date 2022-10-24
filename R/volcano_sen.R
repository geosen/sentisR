#'Create a volcano plot
#'
#'Creates a volcano plot based on a gene list with logFC and (fdr or PValue) columns
#'Works best with edgeR results. 
#'
#'
#'
#'date 24/10/2022
#'@param table Input table with logFC and (fdr or PValue) columns
#'@param title Title of the plot
#'@param significanceMeasure the significance Measure to be used for statistical significance. Must match the corresponding column' name
#'@param thresSig The statistical Measure significance value
#'@param thresLogFC The Log Fold change significance threshold
#'@param colorUp The upregulated genes colors
#'@param colorDown The downregulated genes colors
#'@param colorNeutral The color of non-statistically significant or non-logFC significant genes
#'
#'
#'
#'@return ggplot object
#'@export


volcano_sen <- function(table, title = "Cond_Up vs Cond_Down", significanceMeasure = c("fdr","PValue"), thresSig = 0.05, thresLogFC = 0.58, colorUp="darkgreen",colorDown = "steelblue", colorNeutral = "lightgrey" )
{
  if (length(significanceMeasure) >1) {
    significanceMeasure <- 'fdr'
  }
  require(ggplot2)
  require(dplyr)
  ##Volcano code
  table_sigup <- table %>%
    filter(get(significanceMeasure) < thresSig, logFC > thresLogFC)
  table_sigdown <- table %>%
    filter(get(significanceMeasure) < thresSig, logFC < -thresLogFC)
  
  ggplot(table, aes(logFC, -log10(get(significanceMeasure)))) + 
    geom_point(color =colorNeutral) +
    geom_point(data = table_sigup, color = colorUp)+
    geom_point(data = table_sigdown, color = colorDown) + 
    geom_hline(yintercept= -log10(thresSig))+
    geom_vline(xintercept = -thresLogFC)+
    geom_vline(xintercept = thresLogFC)+
    scale_x_continuous(breaks = c(-3,-thresLogFC,0,thresLogFC,3,6,9), labels = c(-3,-thresLogFC,0,thresLogFC,3,6,9))+
    scale_y_continuous(breaks = c(0,-log10(thresSig),2.5,5,7.5,10), labels = c(0,paste0(toupper(significanceMeasure)," = ", thresSig),2.5,5,7.5,10))+
    labs(title = title, x = "LogFC", y = paste0('-Log10(',toupper(significanceMeasure),')'))+
    theme_minimal()
  
}