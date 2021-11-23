# Purpose: load data and convert to Biom file for Phyloseq; 
# Date: 07/03/2020
# Type: load data
# Author: Paulina Giraldo Perez


#
# Get OTU table
#

if(organism == "all"){
  
  mergedDataSets <- list()
  
  for(organismsToGetData in 1:length(dataFile)){
    
    dataSamples <- read.csv(dataFile[organismsToGetData], colClasses = c("#OTU ID"="character"), check.names = F)
    
    # Fix name of first column as it comes with a # symbol
    names(dataSamples)[1] <- "IdOTU"
    row.names(dataSamples) <- dataSamples$IdOTU
    dataSamples[1] <- NULL
    # See how it looks in Qiime
    otuData <- otu_table(dataSamples, taxa_are_rows = T)
    
    #
    # Get metadata
    #
    metaData <- read.csv(mappingFile) %>%
      dplyr::mutate( Year_of_program = as.factor(Year_of_program) )
    
    
    # Check if names in Metadata match names in OTU table
    mapName <- data.frame( OTU = names(dataSamples)[order(names(dataSamples))],
                           metaData = metaData$SampleID[order(metaData$SampleID)] ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate( isEqual = ifelse(OTU == metaData, "Yes", "No") )
    
    if( sum(mapName$isEqual == "No") > 0 ){
      print(mapName %>% 
              dplyr::filter(isEqual == "No"))
      stop("There are differences between samples names in OTU and metadata. See differences above.")
    }
    
    # See how it looks in Qiime
    mapData <- sample_data(metaData)
    row.names(mapData) <- mapData$SampleID
    
    #
    # Get Taxonomic data
    #
    taxonomyDataCsv <- read.csv(taxonomyDataFile[organismsToGetData]) %>%
      dplyr::mutate( Taxon = gsub("d__", "k__", Taxon, fixed = T) ) %>%
      dplyr::filter( Taxon != "Unassigned" )
    
    # Build Qiime taxonomicdata
    taxData <- list()
    
    for(i in 2:nrow(taxonomyDataCsv)){
      taxData[i-1] <- list(parse_taxonomy_qiime(as.character(taxonomyDataCsv$Taxon[i])))
      names(taxData)[i-1] <- as.character(taxonomyDataCsv$Feature.ID[i])
    }
    
    taxonomyData <- build_tax_table(taxData)
    
    if( !any(colnames(taxonomyData)[1] %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) ){
        
      colnames(taxonomyData) <- taxLevels
    
    }
    
    #
    # Merge data in a list
    #
    
    mergedDataSets[[organismsToGetData]] <- phyloseq(otuData, mapData, taxonomyData)
    
    #
    # Make sure all the organisms have the same taxonomic levels
    #
    
    # Merge data
    if(!exists("biomDataPre")){

      # Compile the phyloseq data
      biomDataPre <- phyloseq(otuData, mapData, taxonomyData)

    } else {

      dataToBind <- phyloseq(otuData, mapData, taxonomyData)

      biomDataPre <- biomDataPre %>%
        merge_phyloseq(dataToBind)

    }
    
  }
  
  # Find all ranks and keep only the ones present everywhere
  
  rankMergedDataSets <- c()
  
  for(i in 1:length(dataFile)){
  
    rankMergedDataSets <- c(rankMergedDataSets, rank_names(mergedDataSets[[i]]))
  
  }
  
  rankedNames <- table(rankMergedDataSets)
  
  selectedRanks <- rankedNames[rankedNames == max(rankedNames, na.rm = T)]
  
  for(i in 1:length(dataFile)){
    
    if(!exists("biomDataPre")){
      
      # Get the data
      dataToFilter <- mergedDataSets[[i]]
      
      # Get the taxa we need
      dataToMerge <- tax_glom(dataToFilter, taxrank = names(selectedRanks), NArm=F)
      
      # Compile the phyloseq data
      biomDataPre <- phyloseq(otuData, mapData, taxonomyData)
      
    } else {
      
      # Get the data
      dataToFilter <- mergedDataSets[[i]]
      
      # Get the taxa we need
      dataToMerge <- tax_glom(dataToFilter, taxrank = names(selectedRanks), NArm=F)
      
      biomDataPre <- biomDataPre %>%
        merge_phyloseq(dataToMerge)
      
    }
    
  }
  
} else {
  
  dataSamples <- read.csv(dataFile, colClasses = c("#OTU ID"="character"), check.names = F)
  
  # Fix name of first column as it comes with a # symbol
  names(dataSamples)[1] <- "IdOTU"
  row.names(dataSamples) <- dataSamples$IdOTU
  dataSamples[1] <- NULL
  # See how it looks in Qiime
  otuData <- otu_table(dataSamples, taxa_are_rows = T)
  
  #
  # Get metadata
  #
  metaData <- read.csv(mappingFile) %>%
    dplyr::mutate( Year_of_program = as.factor(Year_of_program) )
  
  
  # Check if names in Metadata match names in OTU table
  mapName <- data.frame( OTU = names(dataSamples)[order(names(dataSamples))],
                         metaData = metaData$SampleID[order(metaData$SampleID)] ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate( isEqual = ifelse(OTU == metaData, "Yes", "No") )
  
  if( sum(mapName$isEqual == "No") > 0 ){
    print(mapName %>% 
            dplyr::filter(isEqual == "No"))
    stop("There are differences between samples names in OTU and metadata. See differences above.")
  }
  
  # See how it looks in Qiime
  mapData <- sample_data(metaData)
  row.names(mapData) <- mapData$SampleID
  
  
  # taxonomyFile <- data.frame(Feature.ID = row.names(dataSamples),
  #            Taxon = "k__Eukaryota_2759;p__Ascomycota_4890;c__Sordariomycetes_147550;o__Hypocreales_5125;f__Peronosporaceae_4777;g__Phytophthora_4783;s__PhytophthoraTrem_4783",
  #            Confidence = 1)
  # 
  # write.csv(taxonomyFile, file = "data/metazoaTaxonomy.csv", row.names = F)
  # 
  #
  # Get Taxonomic data
  #
  taxonomyDataCsv <- read.csv(taxonomyDataFile) %>%
    dplyr::mutate( Taxon = gsub("d__", "k__", Taxon, fixed = T) ) %>%
    dplyr::filter( Taxon != "Unassigned" )
  
  # Build Qiime taxonomicdata
  taxData <- list()
  
  for(i in 2:nrow(taxonomyDataCsv)){
    taxData[i-1] <- list(parse_taxonomy_qiime(as.character(taxonomyDataCsv$Taxon[i])))
    names(taxData)[i-1] <- as.character(taxonomyDataCsv$Feature.ID[i])
  }
  
  taxonomyData <- build_tax_table(taxData)
  
  if( !is.na(taxLevels[1]) ){
    colnames(taxonomyData) <- taxLevels
  }
  
  # Compile the phyloseq data
  biomDataPre <- phyloseq(otuData, mapData, taxonomyData)
  
}

# Remove not required objects; keeping only biomData, the phyloseq object
rm( list = ls()[!(ls() %in% toKeep)] )

