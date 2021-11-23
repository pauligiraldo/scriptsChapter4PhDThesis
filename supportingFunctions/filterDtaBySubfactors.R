# Filter data by subfactors and repeat test for interactions and pairwises
filterDataBySubfactors <- function(data, dataLongFormat, subfactorToFilter){
  
  filteredData <- dataLongFormat %>%
    dplyr::filter( subfactor %in% subfactorToFilter ) %>%
    dplyr::group_by( SampleID ) %>%
    dplyr::mutate( n = n() ) %>%
    dplyr::filter( n == length(subfactorToFilter) )
  
  results <- data %>%
    dplyr::filter( SampleID %in% filteredData$SampleID )
    
  return(results)

}

