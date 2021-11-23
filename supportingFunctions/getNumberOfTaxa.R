getNumberOfTaxa <- function( biomDataObject, taxaLevel ){
  
  dataNumberOfTaxa <- 
    # Get and groom OTU table
    otu_table(biomDataObject) %>%
    data.frame(check.names = F) %>%
    rownames_to_column( var = "taxCode" ) %>%
    tidyr::gather( SampleID, isPresent, -taxCode ) %>%
    # Change data to binary or presence and absence
    dplyr::mutate( isPresent = ifelse(isPresent > 0, 1, isPresent) ) %>%
    # Group by sample and count number of taxa
    dplyr::group_by( SampleID ) %>%
    dplyr::summarise( numberOfTaxa = sum(isPresent, na.rm = T) ) %>%
    data.frame() %>%
    # Add the columns with factor values
    dplyr::left_join( sample_data(biomDataObject) %>%
                        data.frame(check.names = F),
                      by = "SampleID" ) %>%
    dplyr::mutate( taxon = taxaLevel)
  
  return(dataNumberOfTaxa)
  
}