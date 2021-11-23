# Purpose: run Kruskal Wallis testing for all potential interactions for Shannon and Simpson
# Date: 25/04/2020
# Author: Paulina Giraldo Perez





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
  
  richnessMeasureToUse <- c("Shannon", "Simpson", "Chao1")
  
  shannonAndSimpsonPre <- estimate_richness(taxData, measures = richnessMeasureToUse) %>%
    tibble::rownames_to_column( "SampleID" ) 
  
  message("Estimated richness. Merging with sample_data")
  
  shannonAndSimpson <- shannonAndSimpsonPre %>%
    dplyr::mutate( SampleID = gsub("^X", "", SampleID) ) %>%
    dplyr::left_join( data.frame(sample_data(taxData)), by = "SampleID" )
  
  # Check if all samples matched with sample data
  if(nrow(shannonAndSimpsonPre) == nrow(shannonAndSimpson)){
    message("Merge OK!")
  } else {
    stop("Sample names are not matching... you will need to check it.")
  }
  
  # Save data Shannon and Simpson to use in the future, e.g. to build the box plots later
  save(shannonAndSimpson, file = paste0(intermediateFolder, "/datashannonAndSimpsonLevelTaxon_", taxonomicLevels[i], ".RData"))
  
  # Build boxplot for Shannon
  source("supportingFunctions/buildBoxPlotGeneric.R")
  
  for(measureToUse in richnessMeasureToUse){
    
    buildBoxPlotGeneric( data = shannonAndSimpson, variable = measureToUse, factorsToTest = factorsToTest, 
                  factorLevels = factorLevels, title = taxonomicLevels[i], 
                  yTitle = measureToUse,
                  outputFile = paste0(outputFolder, "/analysis_numbers_top_level_boxplot_measure_", measureToUse, "_taxLevel_", taxonomicLevels[i], ".pdf" ))
  
    # Generate pairwise comparisons for all the factors we want to test
    source("supportingFunctions/kWallisPairWise.R")
    pairWiseComps <- kWallisPairWise(data = shannonAndSimpson, nameColumnWithNumberOfTaxa = measureToUse, factor = factorsToTest) %>%
      dplyr::mutate( taxon = taxonomicLevels[i],
                    subfactorTested = NA )
  
    if(file.exists(paste0(outputFolder, "/analysis_numbers_pairwise.csv"))){
      write.table(pairWiseComps, file = paste0(outputFolder, "/analysis_numbers_pairwise.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_numbers_pairwise.csv")), append = T, row.names = F, )
    } else {
      write.csv(pairWiseComps, file = paste0(outputFolder, "/analysis_numbers_pairwise.csv"), row.names = F )
    }

    # Test for interactions
    source("supportingFunctions/kWallisInteractions.R")
    interactionsComps <- kWallisInteractions(data = shannonAndSimpson, nameColumnWithNumberOfTaxa = measureToUse, factor = factorsToTest) %>%
      dplyr::mutate( taxon = taxonomicLevels[i],
                    subfactorTested = NA )
  
  
    if(file.exists(paste0(outputFolder, "/analysis_numbers_interaction.csv"))){
      write.table(interactionsComps, file = paste0(outputFolder, "/analysis_numbers_interaction.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_numbers_interaction.csv")), append = T, row.names = F)
    } else {
      write.csv(interactionsComps, file = paste0(outputFolder, "/analysis_numbers_interaction.csv"), row.names = F )
    }
  
  }
  
  # Build list for comparison to do at level 2
  factorsToTestLevel2 <- factorsToTest #[factorsToTest != "Management"]
    
  compList <- list()
  
  for(j in 1:length(factorsToTestLevel2)){
    
    compList[factorsToTestLevel2[j]] <- list(unique(as.character(shannonAndSimpson[, factorsToTestLevel2[j]])))
    
  }
  
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
    expand.grid(compList[[2]], compList[[3]], compList[[4]]) )
  
  
  
  for(measureToUse in richnessMeasureToUse){
  
    # Convert data to long format
    dataLong <- shannonAndSimpson %>%
      dplyr::select( SampleID, one_of(measureToUse), one_of(factorsToTest) ) %>%
      tidyr::gather( factor, subfactor, -SampleID, -one_of(measureToUse) )
  
    dataSelected <- shannonAndSimpson %>%
      dplyr::select( SampleID, one_of(measureToUse), one_of(factorsToTest) )
  
    # Filter data by subfactors and repeat test for interactions and pairwises
    source("supportingFunctions/filterDtaBySubfactors.R")
  
    for(j in 1:length(comparisonsToDo)){
    
      for(k in 1:nrow(comparisonsToDo[[j]])){
      
        if(ncol(comparisonsToDo[[j]]) == 1){
       
          subfactorsToFilter <- as.character(comparisonsToDo[[j]][k,])
        
        } else {
        
          subfactorsToFilter <- (comparisonsToDo[[j]][k,] %>% 
            gather( var, value ))$value
        
        }
      
        filteredData <- filterDataBySubfactors(dataSelected, dataLong, subfactorsToFilter )
      
        # Save data per taxon to use in the future, e.g. to build the box plots later
        #save(filteredData, file = paste0(intermediateFolder, "/dataNumberOfTaxaLevelTaxon_", taxonomicLevels[i], "_", paste(subfactorsToFilter, collapse = "_"), ".RData"))
      
        source("supportingFunctions/buildBoxPlotGeneric.R")
        buildBoxPlotGeneric( data = filteredData, variable = measureToUse, factorsToTest = factorsToTest, 
                             factorLevels = factorLevels, title = taxonomicLevels[i], 
                             yTitle = measureToUse, 
                             outputFile = paste0(outputFolder, "/analysis_boxplot_measure_", measureToUse, "_taxLevel_", taxonomicLevels[i], "_factor_", paste(subfactorsToFilter, collapse = "_"), ".pdf" ) )
        
        # Generate pairwise comparisons for all the factors we want to test
        source("supportingFunctions/kWallisPairWise.R")
        factorsToBeRemoved <- dataLong %>%
          dplyr::filter( subfactor %in% subfactorsToFilter )
      
        pairWiseComps <- kWallisPairWise(data = filteredData, nameColumnWithNumberOfTaxa = measureToUse, factor = factorsToTest[!(factorsToTest %in% unique(factorsToBeRemoved$factor))]) %>%
          dplyr::mutate( taxon = taxonomicLevels[i],
                        subfactorTested = paste(subfactorsToFilter, collapse = "_") ) 
      
        if(file.exists(paste0(outputFolder, "/analysis_numbers_pairwise.csv"))){
          write.table(pairWiseComps, file = paste0(outputFolder, "/analysis_numbers_pairwise.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_numbers_pairwise.csv")), append = T, row.names = F)
        } else {
          write.csv(pairWiseComps, file = paste0(outputFolder, "/analysis_numbers_pairwise.csv"), row.names = F )
        }
      
        # Test for interactions
        source("supportingFunctions/kWallisInteractions.R")
        factorsToTestForInteraction <- factorsToTest[!(factorsToTest %in% unique(factorsToBeRemoved$factor))]
      
        if(length(factorsToTestForInteraction) > 1){
          interactionsComps <- kWallisInteractions(data = filteredData, nameColumnWithNumberOfTaxa = measureToUse, factor = factorsToTest[!(factorsToTest %in% unique(factorsToBeRemoved$factor))]) %>%
            dplyr::mutate( taxon = taxonomicLevels[i],
                          subfactorTested = paste(subfactorsToFilter, collapse = "_") ) 
      
        
          if(file.exists(paste0(outputFolder, "/analysis_numbers_interaction.csv"))){
            write.table(interactionsComps, file = paste0(outputFolder, "/analysis_numbers_interaction.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_numbers_interaction.csv")), append = T, row.names = F)
          } else {
            write.csv(interactionsComps, file = paste0(outputFolder, "/analysis_numbers_interaction.csv"), row.names = F )
          }
        
        }
      }
    }
    
  }
}
  
  
  
  
  


