getDataByTaxa <- function( biomDataObject, taxaLevel ){
  
  #biomDataObject = taxData
  #taxaLevel = taxonomicLevels[i]
  
  dataByTaxa <- 
    otu_table(biomDataObject) %>%
    data.frame(check.names = F, stringsAsFactors = F) %>%
    rownames_to_column( var = "taxCode" ) %>%
    tidyr::gather( SampleID, abundance, -taxCode ) %>%
    # Add the columns with factor values
    dplyr::left_join( tax_table(biomDataObject)[,taxaLevel] %>%
                        data.frame(check.names = F, stringsAsFactors = F) %>%
                        dplyr::rename( "otuNames" = names(.) ) %>%
                        rownames_to_column( var = "taxCode" ) ,
                      by = "taxCode" ) %>%
    # Add the columns with factor values
    dplyr::left_join( sample_data(biomDataObject) %>%
                        data.frame(check.names = F, stringsAsFactors = F),
                      by = "SampleID" )
  
  return(dataByTaxa)
  
}