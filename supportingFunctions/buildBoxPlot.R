# Build boxplot for this specific taxon
buildBoxPlot <- function(dataNumberOfTaxa, factorsToTest, factorLevels, title, outputFile ){
  
  dataForBoxPlot <- dataNumberOfTaxa %>% 
    dplyr::select(SampleID, numberOfTaxa, one_of(factorsToTest)) %>%
    dplyr::mutate(taxLevel = title) %>%
    tidyr::gather( factor, subFactor, -SampleID, -numberOfTaxa, -taxLevel ) %>%
    dplyr::mutate(subFactor = factor(subFactor, 
                                     levels = rev(factorLevels)) )
  
  pBoxPlot <- ggplot( dataForBoxPlot, 
                      aes(x = subFactor, y = numberOfTaxa, color = subFactor) ) +
    geom_boxplot() +
    coord_flip() +
    # Define labels of X and Y axes
    labs( x = "", y = "\nNumber of taxa" ) +
    # Define the title
    ggtitle( title) +
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

