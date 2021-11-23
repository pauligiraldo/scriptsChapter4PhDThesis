

#
# Kruskal-Wallis with interaction testing per taxonomic level
#

# The taxonomic levels that will be tested are defined through the object taxonomicLevels,
# which is defined in the script analyses/fourYears/analysisFourYearsPipeline.R

# The factors that will be tested are defined through the object factorsToTest,
# which is defined in the script analyses/fourYears/analysisFourYearsPipeline.R

# Start loop by taxonomic levels and interactions
print(taxonomicLevels)

for(i in 1:length(taxonomicLevels)) {
  print(i)
  message(paste("Testing taxonomic level:", taxonomicLevels[i]))
  
  # Check if the analysis should be done by OTUs only or by a taxonomic level
  if(taxonomicLevels[i] == "OTUs"){
    # No need to filter data if analysis for OTUs only
    taxData <- biomDataNormalised
  } else {
    # Filter the data for the taxonomic level we need
    taxData <- tax_glom(biomDataNormalised, taxrank = taxonomicLevels[i], NArm=F)
  }
  
  taxData <- prune_taxa(taxa_sums(taxData) > 0, taxData)
  
  #
  # Analysis by OTU
  #
  source("supportingFunctions/getDataByTaxa.R")
  dataByTaxa <- getDataByTaxa(taxData, taxonomicLevels[i]) %>%
    dplyr::filter( !(is.na(otuNames)) )
  
  # Kruskal Wallis by OTU
  source("supportingFunctions/kWallisByOTU.R")
  message("KWallis top level....")
  message(paste0("Factor to test: ", factorsToTest))
  message(paste0("Factor levels: ", factorLevels))
  message(paste0("TaxLevel: ", taxonomicLevels[i]))

  kWallisResults <- kWallisByOTU(dataByTaxa, factor = factorsToTest, factorLevels = factorLevels, taxLevel = taxonomicLevels[i])
  
  
  if(file.exists(paste0(outputFolder, "/analysis_kwallis_by_organism.csv"))){
    write.table(kWallisResults$pairwise, file = paste0(outputFolder, "/analysis_kwallis_by_organism.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_kwallis_by_organism.csv")), append = T, row.names = F)
  } else {
    write.csv(kWallisResults$pairwise, file = paste0(outputFolder, "/analysis_kwallis_by_organism.csv"), row.names = F )
  }
  
  if(file.exists(paste0(outputFolder, "/analysis_kwallis_interaction_by_organism.csv"))){
    write.table(kWallisResults$interaction, file = paste0(outputFolder, "/analysis_kwallis_interaction_by_organism.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_kwallis_interaction_by_organism.csv")), append = T, row.names = F)
  } else {
    write.csv(kWallisResults$interaction, file = paste0(outputFolder, "/analysis_kwallis_interaction_by_organism.csv"), row.names = F )
  }
  
  
  # Build list for comparison to do at level 2
  factorsToTestLevel2 <- factorsToTest #[factorsToTest != "Management"]
  
  rawData <- sample_data(taxData) %>%
    data.frame(stringsAsFactors = F)
    
  compList <- list()
  
  for(j in 1:length(factorsToTestLevel2)){
    
    compList[factorsToTestLevel2[j]] <- list(unique(as.character(rawData[, factorsToTestLevel2[j]])))
    
  }
  
  if( length(factorsToTestLevel2) == 1 ){
    
    comparisonsToDo <- list(
      expand.grid(compList[[1]])
    )
    
  }
  
  if( length(factorsToTestLevel2) == 2 ){
    
    comparisonsToDo <- list(
      expand.grid(compList[[1]]),
      expand.grid(compList[[2]]),
      expand.grid(compList[[1]], compList[[2]])
    )
    
  }
  
  if( length(factorsToTestLevel2) == 3 ){
    
    comparisonsToDo <- list(
      expand.grid(compList[[1]]),
      expand.grid(compList[[2]]),
      expand.grid(compList[[3]]),
      expand.grid(compList[[1]], compList[[2]]),
      expand.grid(compList[[1]], compList[[3]]),
      expand.grid(compList[[2]], compList[[3]]),
      expand.grid(compList[[1]], compList[[2]], compList[[3]])
    )
    
  }
  
  if( length(factorsToTestLevel2) == 4 ){
    
    comparisonsToDo <- list(
      expand.grid(compList[[1]]),
      expand.grid(compList[[2]]),
      expand.grid(compList[[3]]),
      expand.grid(compList[[4]]),
      expand.grid(compList[[1]], compList[[2]]),
      expand.grid(compList[[1]], compList[[3]]),
      expand.grid(compList[[1]], compList[[4]]),
      expand.grid(compList[[2]], compList[[3]]),
      expand.grid(compList[[2]], compList[[4]]),
      expand.grid(compList[[3]], compList[[4]]),
      expand.grid(compList[[1]], compList[[2]], compList[[3]]),
      expand.grid(compList[[1]], compList[[2]], compList[[4]]),
      expand.grid(compList[[1]], compList[[3]], compList[[4]]),
      expand.grid(compList[[2]], compList[[3]], compList[[4]]) 
    )
    
  }
  
  comparisonsToDo <- comparisonsToDo[c(1, 4, 5, 7)]
  
  # Convert data to long format
  dataNumberOfTaxaLong <- rawData %>%
    dplyr::select( SampleID, one_of(factorsToTest) ) %>%
    tidyr::gather( factor, subfactor, -SampleID )
  
  dataNumberOfTaxaSelected <- rawData %>%
    dplyr::select( SampleID, one_of(factorsToTest) )
  
  # Filter data by subfactors and repeat test for interactions and pairwises
  source("supportingFunctions/filterDtaBySubfactors.R")
  
  for(j in 1:length(comparisonsToDo)){
    
    for(k in 1:nrow(comparisonsToDo[[j]])){
      
      print(j)
      print(k)
      
      if(ncol(comparisonsToDo[[j]]) == 1){
        
        subfactorsToFilter <- as.character(comparisonsToDo[[j]][k,])
        
      } else {
        
        subfactorsToFilter <- (comparisonsToDo[[j]][k,] %>% 
                                 gather( var, value ))$value
        
      }
      
      filteredData <- filterDataBySubfactors(dataNumberOfTaxaSelected, dataNumberOfTaxaLong, subfactorsToFilter )
      
      #
      # Analysis by OTU
      #
      dataByTaxaLevel2 <- dataByTaxa %>%
        dplyr::filter( SampleID %in% unique(filteredData$SampleID) )
      
      # Kruskal Wallis by OTU
      source("supportingFunctions/kWallisByOTU.R")
      message("KWallis lower levels....")
      #message(paste0("Factor to test: ", factorsToTest))
      #message(paste0("Factor levels: ", factorLevels))
      #message(paste0("TaxLevel: ", taxonomicLevels[i]))
      #message(paste0("Subfactors: ", paste(subfactorsToFilter, collapse = "_")))
      
      kWallisResults <- kWallisByOTU(dataByTaxaLevel2, factor = factorsToTest, factorLevels = factorLevels, taxLevel = taxonomicLevels[i], subfactors = paste(subfactorsToFilter, collapse = "_"))
      
      
      if(file.exists(paste0(outputFolder, "/analysis_kwallis_by_organism.csv"))){
        write.table(kWallisResults$pairwise, file = paste0(outputFolder, "/analysis_kwallis_by_organism.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_kwallis_by_organism.csv")), append = T, row.names = F)
      } else {
        write.csv(kWallisResults$pairwise, file = paste0(outputFolder, "/analysis_kwallis_by_organism.csv"), row.names = F )
      }
      
      if(file.exists(paste0(outputFolder, "/analysis_kwallis_interaction_by_organism.csv"))){
        write.table(kWallisResults$interaction, file = paste0(outputFolder, "/analysis_kwallis_interaction_by_organism.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_kwallis_interaction_by_organism.csv")), append = T, row.names = F)
      } else {
        write.csv(kWallisResults$interaction, file = paste0(outputFolder, "/analysis_kwallis_interaction_by_organism.csv"), row.names = F )
      }
      
    }
  }
}







