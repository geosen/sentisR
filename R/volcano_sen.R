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


volcano_sen <- function(table, title = "Cond_Up vs Cond_Down", significanceMeasure = c("fdr","PValue"), thresSig = 0.05, thresLogFC = 0.58, colorUp="darkgreen",colorDown = "steelblue", colorNeutral = "lightgrey" , xlim = c(-10,10), ylim = c(0,15))
{
  if (length(significanceMeasure) >1) {
    significanceMeasure <- 'fdr'
  }
  ##Volcano code
  upregulated <- table %>%
    filter(get(significanceMeasure) < thresSig, logFC > thresLogFC) %>%
    mutate(
      change = 'Upregulated'
    )
  
  downregulated <- table %>%
    filter(get(significanceMeasure) < thresSig, logFC < -thresLogFC) %>%
    mutate(
      change = 'Downregulated'
    )
  
  insignificant <- table %>%
    filter(get(significanceMeasure) > thresSig | get(significanceMeasure) < thresSig & logFC <= thresLogFC & logFC >= -thresLogFC ) %>%
    mutate(
      change = 'Insignificant'
    )
  
  table <- rbind(upregulated,downregulated,insignificant)
  
  ggplot(table, aes(logFC, -log10(get(significanceMeasure)),color = as.factor(change))) + 
    geom_point() +
    geom_hline(yintercept= -log10(thresSig))+
    geom_vline(xintercept = -thresLogFC)+
    geom_vline(xintercept = thresLogFC)+
    scale_x_continuous(breaks = c(seq(xlim[1],xlim[2],2),-thresLogFC,0,thresLogFC), labels = c(seq(xlim[1],xlim[2],2),-thresLogFC,0,thresLogFC), limits = xlim)+
    scale_y_continuous(breaks = c(seq(ylim[1],ylim[2],2),-log10(thresSig)), labels = c(seq(ylim[1],ylim[2],2),paste0(toupper(significanceMeasure)," = ", thresSig)), limits = ylim)+
    scale_color_manual(values = c('Upregulated' = colorUp, 'Downregulated' = colorDown, 'Insignificant' = colorNeutral))+
    labs(title = title, x = 'LogFC', y = paste0('-Log10(',toupper(significanceMeasure),')'), color = 'Differential \nExpression \nResult')+
    theme_minimal()
  
}