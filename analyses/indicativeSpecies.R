#
##
###
####
##### Indicative species
####
###
##
#

rm(list = ls())  

# The following script will: 
# * load the required libraries
# * load the input data "data/itsASVTable.csv" | "data/metadata_fungi_coding2.csv" | "data/itsTaxonomy.csv"
# * merge the input data into a phyloseq object

# Define the organisms to analyse: "bacteria", "fungi", "metazoa", "all"
targetOrganisms <- c("bacteria", "fungi", "metazoa", "all")
#targetOrganisms <- c("bacteria")

# Define years to use
yearsToUse <- 1:3

# Permutation value to use in the multipatti function
permutationValue <- 100

#
# Choose the normalisation method to use: "rarefy", "bc" [bias correction] or "notNormalised"
#
normaliseMethod <- "bc"

# If the normalisation method chosen is "rarefy", use the options below for further
# configure how rarefying will be applied:

# Test sample size for rarefication:
# The variable or objects below depthRangeFrom,
# depthRangeTo and depthRangeBy can be used to define the range of reads or sample depths that should be tested.
# The default values of depthRangeFrom = 500, depthRangeTo = 10000 and depthRangeBy = 500. These values would result
# in testing from 500 reads to 10000 reads by 500 reads, or the ranges below:

# 500  1000  1500  2000  2500  3000  3500  4000  4500  5000  5500  6000  6500  7000  7500  8000  8500  9000  9500 10000

# If you don't want to rarefy, just sets depthRangeFrom, depthRangeTo and depthRangeBy to 0, or zero.

# Your turn:
depthRangeFrom = 500
depthRangeTo = 10000
depthRangeBy = 100

#### Finishes choosing normalisation method ##############

# Define minimum proportion os samples per data class
# For example, after rarefying, my final dataset cannot contain a data class, like Management:Future, represented by
# less than the minimumPropToUse:
minimumPropToUse <- 0.25

shouldGeneratePlotNReads <- F

# The object dataSizeLog is used to store the profile of samples per factor level for each organism before and after rarefying.
dataSizeLog <- data.frame()

# Define objects to save
toKeep <- ls()
toKeep <- ls()
toKeep <- c(toKeep, "biomData", "biomDataPre", "organism", "shouldTestValuesForRarefying", "readsValueForRarefying",
            "outputFolder", "intermediateFolder", "shouldGeneratePlotNReads", "dataSizeLog",
            "mappingFile", "taxonomyDataFile", "normaliseMethod", "taxLevels", "taxonomicLevelsToTest")

for(organism in targetOrganisms){
  
  message(paste0("Analysing data for ", organism, "...."))
  
  switch(organism,
         bacteria = {
           dataFile <- "data/bacteriaASVTable.csv"
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomicLevelsToTest <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "OTUs")
           taxonomyDataFile <- "data/bacteriaTaxonomy.csv"
           readsValueForRarefying <- 2600
         },
         fungi = {
           dataFile <- "data/itsASVTable.csv"
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomicLevelsToTest <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "OTUs")
           taxonomyDataFile <- "data/itsTaxonomy.csv"
           readsValueForRarefying <- 2000
         },
         metazoa = {
           dataFile <- "data/metazoaASVTable.csv"
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomicLevelsToTest <- c("OTUs")
           taxonomyDataFile <- "data/metazoaTaxonomy.csv"
           readsValueForRarefying <- 2200
         },
         all = {
           dataFile <- c("data/bacteriaASVTable.csv", "data/itsASVTable.csv", "data/metazoaASVTable.csv") 
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomicLevelsToTest <- c("OTUs")
           taxonomyDataFile <- c("data/bacteriaTaxonomy.csv", "data/itsTaxonomy.csv", "data/metazoaTaxonomy.csv")
           readsValueForRarefying <- 2000
         } 
  )
  
  source("supportingFunctions/loadSamples.R")
  
  # Select samples from the desired years
  biomData <- subset_samples(biomDataPre, Year_of_program %in% yearsToUse)
  print(biomData)
  
  # Loaded samples
  message("Loaded samples.....")
  
  dataProfile <- sample_data(biomData) %>%
    data.frame() %>%
    dplyr::select( SampleID, Year_of_program, Season, Management, Region) %>%
    tidyr::gather( Variable, Value, -SampleID  ) %>%
    dplyr::group_by( Variable, Value ) %>%
    dplyr::summarise( Samples = n(),
                      organism = organism,
                      normaliseMethod = normaliseMethod,
                      time = "before")
  
  if(!file.exists(paste0("output/", organism, "/", normaliseMethod, "/a_indicativeSpecies" ))){
    dir.create(paste0("output/", organism, "/", normaliseMethod, "/a_indicativeSpecies"), recursive = T)
  }
  
  outputFolder <- paste0("output/", organism, "/", normaliseMethod, "/a_indicativeSpecies")
  
  if(!file.exists(paste0("intermediateData/", organism, "/", normaliseMethod, "/a_indicativeSpecies"))){
    dir.create(paste0("intermediateData/", organism, "/", normaliseMethod, "/a_indicativeSpecies"), recursive = T)
  }
  
  intermediateFolder <- paste0("intermediateData/", organism, "/", normaliseMethod, "/a_indicativeSpecies")
  
  # Outputs: 
  # * biomData | This is the phyloseq object built from the input files
  
  # Check the content of the phyloseq data
  (biomData)
  tax_table(biomData) # taxonomic table
  sample_data(biomData) # metadata
  otu_table(biomData) # OTU table
  
  ########
  ################ Apply normalisation method
  ########
  
  
  ##############
  ## Rarefy if normalise method is "rarefy"
  ##############
  if(normaliseMethod == "notNormalised"){
    
    biomDataNormalised <- biomData
    
    save(biomData, file = paste0(intermediateFolder, "/biomData.RData"))
    
  }
  
  ##############
  ## Rarefy if normalise method is "rarefy"
  ##############
  if(normaliseMethod == "rarefy"){
    
    source("supportingFunctions/normaliseRarefy.R", local = environment())
    
  }
  
  ##########################
  ## Bias correction if normalise method is "bc"
  ##########################
  if(normaliseMethod == "bc"){
    
    source("supportingFunctions/normaliseBC.R", local = environment())
    
  }
  
  # Add kResults as an object we would like to keep
  
  #
  # Numbers - Kruskal Wallis on the number of taxa per taxonomic level and factor (i.e. "Management", "Region", "Season" and "Year_of_program").
  #
  
  ### Approach
  # 
  
  # List the factors we would like to test for interaction. Factors will also be tested individually. 
  factorsToTest <- c("Management", "Region", "Season", "Year_of_program")
  factorLevels <- c("C", "F", "HB", "MB", "Budburst", "Veraison", "Harvest", "1", "2", "3")
  # Add kResults as an object we would like to keep
  toKeep <- c(toKeep, "biomDataNormalised", "factorsToTest", "factorLevels", "taxonomicLevels")
  
  
  
  
  #
  # Kruskal-Wallis with interaction testing per taxonomic level
  #

  
  # Start loop by taxonomic levels and interactions
  message("The following taxonomic levels will be tested:")
  print(taxonomicLevelsToTest)
  
  taxonomicLevels <- taxonomicLevelsToTest
  
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
    
    # Convert data to long format
    dataTaxaLong <- sample_data(taxData) %>%
      data.frame( stringsAsFactors = F ) %>%
      dplyr::select( SampleID, one_of(factorsToTest) ) %>%
      tidyr::gather( factor, subfactor, -SampleID )
    
    dataTaxaSelected <- sample_data(taxData) %>%
      data.frame( stringsAsFactors = F ) %>%
      dplyr::select( SampleID, one_of(factorsToTest) )
    
    #
    # Build comparisons at top level for each factor to test
    #
    for(xfactor in factorsToTest){
      
      message(paste0("Testing factor: ", xfactor))
      
      factorBySample <- dataTaxaSelected[, xfactor ]
      
      abundanceData <- otu_table(taxData) %>%
        data.frame(stringsAsFactors = F) %>%
        t()
      
      indicativeResults = multipatt(abundanceData, factorBySample, func = "r.g", control = how(nperm = permutationValue))
      
      if(taxonomicLevels[i] == "OTUs"){
        
        taxDataToMerge <- indicativeResults$sign %>%
          data.frame(stringsAsFactors = F) %>%
          dplyr::add_rownames( "taxaID" ) %>%
          dplyr::mutate( taxaName = paste0("OTU_", taxaID) ) %>%
          dplyr::select( taxaID, taxaName )
        
      } else {
        
        taxDataToMerge <- tax_table(taxData) %>%
          data.frame(stringsAsFactors = F) %>%
          dplyr::select( taxonomicLevels[i] ) %>%
          dplyr::add_rownames( "taxaID" ) %>%
          dplyr::rename( taxaName = taxonomicLevels[i] )
        
      }
      
      significativeTaxa <- indicativeResults$sign %>%
        dplyr::add_rownames( "taxaID" ) %>%
        dplyr::filter( p.value <= 0.05 )
      
      indicativeResultsProcessed <- significativeTaxa %>%
        tidyr::gather( subfactorFromMultipatt, presenceFromMultipatt, -taxaID, -index, -stat, -p.value ) %>%
        dplyr::left_join( taxDataToMerge,
                          by = "taxaID" ) %>%
        dplyr::mutate( taxLevel = taxonomicLevels[i],
                       factor = xfactor,
                       subfactor = NA )
      
      # Save results in separated files
      write.csv(indicativeResultsProcessed, file = paste0(outputFolder, "/analysis_indicativeSpecies_", taxonomicLevels[i], "_", xfactor, ".csv"), row.names = F )
      
      # if(file.exists(paste0(outputFolder, "/analysis_indicativeSpecies.csv"))){
      #   write.table(indicativeResultsProcessed, file = paste0(outputFolder, "/analysis_indicativeSpecies.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_indicativeSpecies.csv")), append = T, row.names = F, )
      # } else {
      #   write.csv(indicativeResultsProcessed, file = paste0(outputFolder, "/analysis_indicativeSpecies.csv"), row.names = F )
      # }
      
    }
    
    #
    # Build comparisons at sublevel for each factor to test
    #
    factorsToTestLevel2 <- factorsToTest #[factorsToTest != "Management"]
    
    compList <- list()
    
    sampleData <- data.frame(sample_data(taxData))
    
    for(j in 1:length(factorsToTestLevel2)){
      
      compList[factorsToTestLevel2[j]] <- list(unique(as.character(sampleData[, factorsToTestLevel2[j]])))
      
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
    
    for(j in 1:length(comparisonsToDo)){
      
      for(k in 1:nrow(comparisonsToDo[[j]])){
    
        if(ncol(comparisonsToDo[[j]]) == 1){
          
          subfactorsToFilter <- as.character(comparisonsToDo[[j]][k,])
          
        } else {
          
          subfactorsToFilter <- (comparisonsToDo[[j]][k,] %>% 
                                   gather( var, value ))$value
          
        }
        
        subfactorSamples <- dataTaxaLong %>%
          dplyr::filter( subfactor %in% subfactorsToFilter ) %>%
          dplyr::group_by( SampleID ) %>%
          dplyr::mutate( n = n() ) %>%
          dplyr::filter( n == length(subfactorsToFilter) )
        
        filteredData <- subset_samples(taxData, SampleID %in% unique(subfactorSamples$SampleID))
        
        # For each factor to test, check how many of them has 2 or more subfactors so the multipatti can be applied
        
        dataTaxaSelected <- sample_data(filteredData) %>%
          data.frame( stringsAsFactors = F ) %>%
          dplyr::select( SampleID, one_of(factorsToTest) )
        
        for(xfactor in factorsToTest){
          
          message(paste0("Checking if factor ", xfactor, " has 2 or more levels to run multipatti..."))
          
          factorBySample <- dataTaxaSelected[, xfactor ]
          
          if(length(unique(factorBySample)) > 1){
            
            message(paste0("Yes! It has ", length(unique(factorBySample)), " levels:"))
            message(paste( unique(factorBySample), collapse = " | "))
            
            abundanceData <- otu_table(filteredData) %>%
              data.frame(stringsAsFactors = F) %>%
              t()
            
            indicativeResults = multipatt(abundanceData, factorBySample, func = "r.g", control = how(nperm = permutationValue))
            
            if(taxonomicLevels[i] == "OTUs"){
              
              taxDataToMerge <- indicativeResults$sign %>%
                data.frame(stringsAsFactors = F) %>%
                dplyr::add_rownames( "taxaID" ) %>%
                dplyr::mutate( taxaName = paste0("OTU_", taxaID) ) %>%
                dplyr::select( taxaID, taxaName )
              
            } else {
              
              taxDataToMerge <- tax_table(taxData) %>%
                data.frame(stringsAsFactors = F) %>%
                dplyr::select( taxonomicLevels[i] ) %>%
                dplyr::add_rownames( "taxaID" ) %>%
                dplyr::rename( taxaName = taxonomicLevels[i] )
              
            }
            
            significativeTaxa <- indicativeResults$sign %>%
              dplyr::add_rownames( "taxaID" ) %>%
              dplyr::filter( p.value <= 0.05 )
            
            indicativeResultsProcessed <- significativeTaxa %>%
              tidyr::gather( subfactorFromMultipatt, presenceFromMultipatt, -taxaID, -index, -stat, -p.value ) %>%
              dplyr::left_join( taxDataToMerge,
                                by = "taxaID" ) %>%
              dplyr::mutate( taxLevel = taxonomicLevels[i],
                             factor = xfactor,
                             subfactor = paste(subfactorsToFilter, collapse = "_") )
            
            # Save results in separated files
            write.csv(indicativeResultsProcessed, file = paste0(outputFolder, "/analysis_indicativeSpecies_", taxonomicLevels[i], "_", xfactor, "_", paste(subfactorsToFilter, collapse = "_"), ".csv"), row.names = F )
            
            # if(file.exists(paste0(outputFolder, "/analysis_indicativeSpecies.csv"))){
            #   write.table(indicativeResultsProcessed, file = paste0(outputFolder, "/analysis_indicativeSpecies.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_indicativeSpecies.csv")), append = T, row.names = F, )
            # } else {
            #   write.csv(indicativeResultsProcessed, file = paste0(outputFolder, "/analysis_indicativeSpecies.csv"), row.names = F )
            # }
            
          } else {
            
            message(paste0("No! It has only ", length(unique(factorBySample)), " levels:"))
            message(paste( unique(factorBySample), collapse = " | "))
            
          }
          
        }
        
      }
      
    }
    
  }
  
}
  