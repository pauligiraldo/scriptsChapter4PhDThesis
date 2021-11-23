# Build boxplot for this specific taxon
buildBoxPlotOTULevel <- function(data, factorsToTest, factorLevels, otuID, outputFile ){
  
  # data = otuData
  # factorsToTest = factor
  # factorLevels = factorLevels
  # otuID
  
  dataForBoxPlot <- data %>% 
    dplyr::select(SampleID, abundance, one_of(factorsToTest)) %>%
    dplyr::mutate(taxLevel = otuID) %>%
    tidyr::gather( factor, subFactor, -SampleID, -abundance, -taxLevel ) %>%
    dplyr::mutate(subFactor = factor(subFactor, 
                                     levels = rev(factorLevels)) )
  
  pBoxPlot <- ggplot( dataForBoxPlot, 
                      aes(x = subFactor, y = abundance, color = subFactor) ) +
    geom_boxplot() +
    coord_flip() +
    # Define labels of X and Y axes
    labs( x = "", y = "\nAbundance" ) +
    # Define the title
    ggtitle( otuID ) +
    # Remove background color and stripes
    theme( legend.position = "none",
           panel.background = element_blank(),
           legend.key=element_blank(),
           strip.background =element_rect(fill="white"),
           plot.title = element_text(hjust = 0.5) )
  
  # Save graph as PDF
  pdf(file = outputFile, width = 4, height = 4)
  print(pBoxPlot)
  dev.off()
}

