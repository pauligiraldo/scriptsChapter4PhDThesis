#
##
###
####
##### Indicative species above family
####
###
##
#

rm(list = ls())  

library("MiscMetabar")

# The following script will: 
# * load the required libraries
# * load the input data "data/itsASVTable.csv" | "data/metadata_fungi_coding2.csv" | "data/itsTaxonomy.csv"
# * merge the input data into a phyloseq object

# Define the organisms to analyse: "bacteria", "fungi", "metazoa", "all"
targetOrganisms <- c("bacteria", "fungi", "metazoa", "all")
# targetOrganisms <- c("bacteria")

#
# Choose the normalisation method to use: "rarefy", "bc" [bias correction] or "notNormalised"
#
normaliseMethod <- "notNormalised"

# Should we use normalised or non-normalised data for this analysis?
useNormalisedData <- FALSE


# Define objects to save
toKeep <- ls()
toKeep <- ls()
toKeep <- c(toKeep)

for(organism in targetOrganisms){
  
  message(paste0("Analysing data for ", organism, "...."))
  
  outputFolderToLook <- paste0("output/", organism, "/", normaliseMethod, "/a_indicativeSpecies")
  
  message(paste0(" the output folder we will look at is: ", outputFolderToLook))
  
  # Load main data
  intermediateFolderToLook <- paste0("intermediateData/", organism, "/", normaliseMethod, "/a_indicativeSpecies")
  
  # Create folder for report
  if(!file.exists(paste0("output/", organism, "/", normaliseMethod, "/a_indicativeSpecies/a_indicSpeciesTaxReport" ))){
    dir.create(paste0("output/", organism, "/", normaliseMethod, "/a_indicativeSpecies/a_indicSpeciesTaxReport"), recursive = T)
  }
  
  destinationReport <- paste0("output/", organism, "/", normaliseMethod, "/a_indicativeSpecies/a_indicSpeciesTaxReport" )
  
  
  if(useNormalisedData){
    
    load(paste0(intermediateFolderToLook, "/biomDataNormalised.RData"))
    
    if(normaliseMethod == "rarefy"){
      
      dataToAnalyse <- biomDataNormalised
      
    }
    
    if(normaliseMethod == "bc"){
      
      otuData <- otu_table(bcData$feature_table, taxa_are_rows = T)
      
      load(paste0(intermediateFolderToLook, "/biomData.RData"))
      
      biomDataNormalised <- phyloseq(otuData, sample_data(biomData), tax_table(biomData))
      
      dataToAnalyse <- biomDataNormalised
      
    }
      

    
  } else {
    
    load(paste0(intermediateFolderToLook, "/biomData.RData"))
    
    dataToAnalyse <- biomData
    
  }
  
  
  # Convert data to binary
  dataToAnalyseBinary <- as_binary_otu_table(dataToAnalyse)
  
  
  # Get results files to analyse
  filesInOutputFolder <- dir(outputFolderToLook, full.names = T, recursive = F)
  
  filesToAnalyse <- data.frame( fileNames = filesInOutputFolder[grep("analysis_indicativeSpecies_OTUs", filesInOutputFolder)], stringsAsFactors = F ) 
  
  for(filesToInspect in filesToAnalyse$fileNames){

    rawData <- read.csv(filesToInspect) 
    
    if(nrow(rawData) > 0) {
      
      # Filter presence multipatt
      rawData <- rawData %>%
        dplyr::filter( presenceFromMultipatt == 1 )
    
      taxIdsToFilter <- unique(rawData$taxaID)
      
      # Prune phyloseq data to keep only selected taxa
      filteredDataToAnalyse <- prune_taxa(taxIdsToFilter, dataToAnalyseBinary)
      
      #tax_datatable(filteredDataToAnalyse)
      
      report <- tax_table(filteredDataToAnalyse) %>%
        data.frame() %>%
        dplyr::add_rownames( "taxaID" ) %>%
        dplyr::left_join( rawData %>% 
                            dplyr::select( taxaID, p.value, stat, subfactorFromMultipatt, factor, subfactor ),
                          by = "taxaID" )
      
      # Get path to save
      #write.csv(report, file = gsub(".csv", "_taxReport.csv", filesToInspect))
      
      if(file.exists(paste0(destinationReport, "/indicSpeciesTaxReport.csv"))){
        write.table(report, file = paste0(destinationReport, "/indicSpeciesTaxReport.csv"), sep = ",", col.names = !file.exists(paste0(destinationReport, "/indicSpeciesTaxReport.csv")), append = T, row.names = F, )
      } else {
        write.csv(report, file = paste0(destinationReport, "/indicSpeciesTaxReport.csv"), row.names = F )
      }
    
    }
    
  }
  
}
  