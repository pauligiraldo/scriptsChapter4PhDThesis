#
##
###
####
##### Pathogenics
####
###
##
#

rm(list = ls())  

# The following script will: 
# * load the required libraries
# * load the input data "data/fungiASVTable.csv" | "data/metadata_fungi_coding2.csv" | "data/fungiTaxonomy.csv"
# * merge the input data into a phyloseq object

# Define the organisms to analyse: "bacteria", "fungi", "metazoa", "all"
targetOrganisms <- c("fungi")
# targetOrganisms <- c("all")

# Define years to use
yearsToUse <- 1:3


#
# Choose the normalisation method to use: "rarefy", "bc" [bias correction] or "notNormalised"
#
normaliseMethod <- "notNormalised"

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
           # taxonomicLevelsToTest <- c("OTUs")
           taxonomyDataFile <- "data/bacteriaTaxonomy.csv"
           readsValueForRarefying <- 2600
         },
         fungi = {
           dataFile <- "data/itsASVTable.csv"
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomicLevelsToTest <- c("Genus", "Species")
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
    dplyr::select( SampleID, Variety, Management, Region) %>%
    tidyr::gather( Variable, Value, -SampleID  ) %>%
    dplyr::group_by( Variable, Value ) %>%
    dplyr::summarise( Samples = n(),
                      organism = organism,
                      normaliseMethod = normaliseMethod,
                      time = "before")
  
  if(!file.exists(paste0("output/", organism, "/pathogenics/", normaliseMethod))){
    dir.create(paste0("output/", organism, "/pathogenics/", normaliseMethod), recursive = T)
  }
  
  outputFolder <- paste0("output/", organism, "/pathogenics/", normaliseMethod)
  
  if(!file.exists(paste0("intermediateData/", organism, "/pathogenics/", normaliseMethod))){
    dir.create(paste0("intermediateData/", organism, "/pathogenics/", normaliseMethod), recursive = T)
  }
  
  intermediateFolder <- paste0("intermediateData/", organism, "/pathogenics/", normaliseMethod)
  
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
  # Numbers - Kruskal Wallis on the number of taxa per taxonomic level and factor (i.e. "Management", "Region", "Variety" and "Year_of_program").
  #
  
  ### Approach
  # Calculate the number of taxa by the following taxonomic levels: "OTUs", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus" and "Species".
  # Use Kruskal-Wallis to inspect potential interactions between the following factors: "Management", "Region", "Variety" and "Year_of_program".
  # For the interactions found, inspect for additional interactions within each factor. 
  
  ### Setup
  # List all the taxonomic levels in biomData
  print(rank_names(biomDataNormalised))
  
  taxonomicLevels <- taxonomicLevelsToTest
  
  
  
  # List the factors we would like to test for interaction. Factors will also be tested individually. 
  factorsToTest <- c("Management", "Region", "Season", "Year_of_program")
  factorLevels <- c("C", "F", "HB", "MB", "Budburst", "Veraison", "Harvest", "1", "2", "3")
  # Add kResults as an object we would like to keep
  toKeep <- c(toKeep, "biomDataNormalised", "factorsToTest", "factorLevels", "taxonomicLevels")
  
  
  ### Analysis by OTU
  source("supportingFunctions/runKWallisByOTU.R", local = environment())
  
  # Remove all objects we don't need anymore
  rm( list = ls()[!(ls() %in% toKeep)])
  
  
}