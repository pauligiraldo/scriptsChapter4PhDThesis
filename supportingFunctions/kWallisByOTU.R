kWallisByOTU <- function(data, factor, factorLevels, taxLevel, subfactors = NA){
# 
  # data = dataByTaxaLevel2
  # factor = factorsToTest
  # factorLevels = factorLevels
  # taxLevel = taxonomicLevels[i]
  # subfactors = paste(subfactorsToFilter, collapse = "_")
  
  resultsDataPairwise <- data.frame()
  
  resultsDataInteraction <- data.frame()
  
  # Find unique OTUs
  uniqueOTUs <- unique(as.character(data$otuNames))
  
  # Get data by single OTU
  for(otuID in uniqueOTUs){
    
    print(otuID)
    
    # Filter Data
    otuData <- data %>%
      dplyr::filter( otuNames == otuID )
    
    #if(nrow(otuData) > 0){
    
        # Build boxplot for this specific taxon
        source("supportingFunctions/buildBoxPlotOTULevel.R")
        if(is.na(subfactors)){
        
          buildBoxPlotOTULevel( otuData, factor, factorLevels, otuID, outputFile = paste0(outputFolder, "/analysis_top_level_boxplot_", taxLevel, "_", otuID, ".pdf" ))
        
        } else {
          
          buildBoxPlotOTULevel( otuData, factor, factorLevels, otuID, outputFile = paste0(outputFolder, "/analysis_top_level_boxplot_", taxLevel, "_", otuID, "_", subfactors, ".pdf" ))
          
        }
        
        # Generate pairwise comparisons for all the factors we want to test
        source("supportingFunctions/kWallisPairWise.R")
        pairWiseComps <- kWallisPairWise(data = otuData, nameColumnWithNumberOfTaxa = "abundance", factor = factor ) %>%
          dplyr::mutate( taxon = taxLevel,
                         otuName = otuID,
                         subfactorTested = NA )
        
        resultsDataPairwise <- resultsDataPairwise %>%
          dplyr::bind_rows( pairWiseComps )
        
        # Test for interactions
        source("supportingFunctions/kWallisInteractions.R")
        interactionsComps <- kWallisInteractions(data = otuData, nameColumnWithNumberOfTaxa = "abundance", factor = factor) %>%
          dplyr::mutate( taxon = taxLevel,
                         otuName = otuID,
                         subfactorTested = NA )
        
        resultsDataInteraction <- resultsDataInteraction %>%
          dplyr::bind_rows( interactionsComps )
        
    #}
    
  }
  
  return(list(pairwise = resultsDataPairwise, interaction = resultsDataInteraction))
  
}