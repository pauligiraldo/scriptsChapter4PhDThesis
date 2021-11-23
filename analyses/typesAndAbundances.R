#
##
###
####
##### Types and abundances
####
###
##
#

rm(list = ls())  

# The following script will: 
# * load the input data "data/itsASVTable.csv" | "data/metadata_fungi_coding2.csv" | "data/itsTaxonomy.csv"
# * merge the input data into a phyloseq object

# Define the organisms to analyse: "bacteria", "fungi", "metazoa", "all"
targetOrganisms <- c("bacteria", "fungi", "metazoa", "all")

# Define years to use
yearsToUse <- 1:3


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
            "mappingFile", "taxonomyDataFile", "normaliseMethod", "taxLevels")

for(organism in targetOrganisms){
  
  message(paste0("Analysing data for ", organism, "...."))
  
  switch(organism,
         bacteria = {
           dataFile <- "data/bacteriaASVTable.csv"
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomyDataFile <- "data/bacteriaTaxonomy.csv"
           readsValueForRarefying <- 2600
         },
         fungi = {
           dataFile <- "data/itsASVTable.csv"
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomyDataFile <- "data/itsTaxonomy.csv"
           readsValueForRarefying <- 2000
         },
         metazoa = {
           dataFile <- "data/metazoaASVTable.csv"
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomyDataFile <- "data/metazoaTaxonomy.csv"
           readsValueForRarefying <- 2200
         },
         all = {
           dataFile <- c("data/bacteriaASVTable.csv", "data/itsASVTable.csv", "data/metazoaASVTable.csv") 
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
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
  
  if(!file.exists(paste0("output/", organism, "/", normaliseMethod))){
    dir.create(paste0("output/", organism, "/", normaliseMethod), recursive = T)
  }
  
  outputFolder <- paste0("output/", organism, "/", normaliseMethod)
  
  if(!file.exists(paste0("intermediateData/", organism, "/", normaliseMethod))){
    dir.create(paste0("intermediateData/", organism, "/", normaliseMethod), recursive = T)
  }
  
  intermediateFolder <- paste0("intermediateData/", organism, "/", normaliseMethod)
  
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
  
  message("Done normalisation!")
  
  # Add kResults as an object we would like to keep
  
  #
  # Numbers - Kruskal Wallis on the number of taxa per taxonomic level and factor (i.e. "Management", "Region", "Season" and "Year_of_program").
  #
  
  ### Approach
  # Calculate the number of taxa by the following taxonomic levels: "OTUs", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus" and "Species".
  # Use Kruskal-Wallis to inspect potential interactions between the following factors: "Management", "Region", "Season" and "Year_of_program".
  # For the interactions found, inspect for additional interactions within each factor. 
  
  ### Setup
  # List all the taxonomic levels in biomData
  print(rank_names(biomDataNormalised))
  
  # The object taxonomicLevels contains all the taxonomic levels we will test.
  
  if( organism == "fungi"){
    # We would like to test all the taxonomic levels in biomData, with the exception of "Rank1".
    (taxonomicLevels <- c(rank_names(biomDataNormalised)[rank_names(biomDataNormalised) != "Kingdom"], "OTUs"))
  } 
  
  if( organism == "bacteria"){
    # We would like to test all the taxonomic levels in biomData, with the exception of "Rank1" and "Kingdom.
    (taxonomicLevels <- c(rank_names(biomDataNormalised)[rank_names(biomDataNormalised) != "Kingdom"], "OTUs"))
  }
  
  if( organism == "metazoa"){
    # We would like to test all the taxonomic levels in biomData, with the exception of "Rank1".
    (taxonomicLevels <- c("OTUs"))
  }
  
  if( organism == "all"){
    # We would like to test all the taxonomic levels in biomData, with the exception of "Rank1".
    (taxonomicLevels <- c("OTUs"))
  }
  
  taxonomicLevelsToTest <- taxonomicLevels
  
  
  
  # List the factors we would like to test for interaction. Factors will also be tested individually. 
  factorsToTest <- c("Management", "Region", "Season", "Year_of_program")
  factorLevels <- c("C", "F", "HB", "MB", "Budburst", "Veraison", "Harvest", "1", "2", "3")
  # Add kResults as an object we would like to keep
  toKeep <- c(toKeep, "biomDataNormalised", "factorsToTest", "factorLevels", "taxonomicLevels", "taxonomicLevelsToTest")
  
  
  ### Test for interaction at top level
  # Run Kruskal-Wallis test for the number of taxa per sample, 
  # per taxonomic level (i.e. "OTUs", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus" and "Species") and 
  # interaction between "Management", "Region", "Season" and "Year_of_program", individually, 2 x 2, 3 x 3 and 4 x 4.
  
  
  ## Number of taxa
  # The following script will: 
  # * filter biomData for each taxonomic level defined in taxonomicLevels
  # * test all potential interactions or combinations of the factors defined in factorsToTest
  source("supportingFunctions/analysis_test_interaction_kruskal_wallis_all_levels.R", local = environment())
  
  # Remove all objects we don't need anymore
  rm( list = ls()[!(ls() %in% toKeep)])
  
  
  
  ## Shannon
  # The following script will: 
  # * filter biomData for each taxonomic level defined in taxonomicLevels
  # * test all potential interactions or combinations of the factors defined in factorsToTest
  source("supportingFunctions/analysis_test_interaction_kruskal_wallis_all_levels_ShannonAndSimpson.R", local = environment())
  
  # Remove all objects we don't need anymore
  rm( list = ls()[!(ls() %in% toKeep)])
  
  
  ### Abundances and Types
  # Run Adonis test per taxonomic level (i.e. "OTUs", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus" and "Species") and
  # individually, 2 x 2, 3 x 3 and 4 x 4.
  
  # The following script will:
  # * filter biomData for each taxonomic level defined in taxonomicLevels
  # * run Adonis on the factors defined in factorsToTest
  source("supportingFunctions/analysis_test_interaction_adonis_all_levels.R", local = environment())
  
  
  # Remove all objects we don't need anymore
  rm( list = ls()[!(ls() %in% toKeep)])
  
}