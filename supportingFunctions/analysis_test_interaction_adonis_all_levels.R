# Purpose: run adoniss testing for all potential interactions
# Date: 07/09/2019
# Type: it is part of the script "analyses/fourYears/analysis_four_years_pipeline.R"
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
  
  # Save data per taxon to use in the future, e.g. to build the box plots later
  save(taxData, file = paste0(intermediateFolder, "/dataTypesAndAbundancesOfTaxaLevelTaxon_", taxonomicLevels[i], ".RData"))
  
  ## Create data profile showing number of samples per factor level
  dataProfile <- sample_data(taxData) %>%
    data.frame() %>%
    dplyr::select( SampleID, Region, Season, Management, Year_of_program) %>%
    tidyr::gather( Variable, Value, -SampleID  ) %>%
    dplyr::group_by( Variable, Value ) %>%
    dplyr::summarise( Samples = n() ) %>%
    dplyr::mutate( subfactorTested = NA )
  
  if(file.exists(paste0(outputFolder, "/sampleProfilePerComparison.csv"))){
    write.table(dataProfile, file = paste0(outputFolder, "/sampleProfilePerComparison.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/sampleProfilePerComparison.csv")), append = T, row.names = F)
  } else {
    write.csv(dataProfile, file = paste0(outputFolder, "/sampleProfilePerComparison.csv"), row.names = F )
  }
  
  if(any(dataProfile$Samples < 3)){
    
    flagSampleSize <- "!!!! less than 3 samples per factor level !!!!!"
    
  } else {
    
    flagSampleSize <- "Ok. More than 3 samples per factor level."
    
  }
  
  #
  # Run adonis at taxonomic level using the factors selected for testing 
  #
  
  # Generate distance matrix
  distanceDf <- distance(taxData, method="jaccard", binary = F)
  
  # Generate formula
  formulaToUse <- as.formula(paste("distanceDf", paste(factorsToTest, collapse=" * "), sep=" ~ "))
  
  # Generate distance matrix binary
  distanceDfBinary <- distance(taxData, method="jaccard", binary = T)
  
  # Generate formula binary
  formulaToUseBinary <- as.formula(paste("distanceDfBinary", paste(factorsToTest, collapse=" * "), sep=" ~ "))
  
  # Run adonis 
  set.seed(4321)
  message( "Start Adonis....")
  adonisGeneral <- data.frame((adonis2( formulaToUse, 
                                       data = data.frame(sample_data(taxData))))) %>%
    rownames_to_column( "factors" ) %>%
    dplyr::mutate( taxa = taxonomicLevelsToTest[i],
                   type = "jaccard-abundance",
                   subfactorTested = NA,
                   flag = flagSampleSize ) %>%
    dplyr::bind_rows( data.frame((adonis2( formulaToUseBinary, 
                                          data = data.frame(sample_data(taxData))))) %>%
                        rownames_to_column( "factors" ) %>%
                        dplyr::mutate( taxa = taxonomicLevelsToTest[i],
                                       type = "jaccard-types",
                                       subfactorTested = NA,
                                       flag = flagSampleSize ) )
  
  if(file.exists(paste0(outputFolder, "/analysis_abundance_and_type_adonis.csv"))){
    write.table(adonisGeneral, file = paste0(outputFolder, "/analysis_abundance_and_type_adonis.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_abundance_and_type_adonis.csv")), append = T, row.names = F)
  } else {
    write.csv(adonisGeneral, file = paste0(outputFolder, "/analysis_abundance_and_type_adonis.csv"), row.names = F )
  }
  
  # Generate NMDS graph
  source("supportingFunctions/getDataForNMDS.R")
  set.seed(4321)
  nmdsData <- getDataForNMDS( data = taxData,
                              distance = "jaccard",
                              minCounts = 1,
                              minSamplePercent = 0.05,
                              binary = F) %>% 
    filter( id.type == "Samples" ) %>%
    dplyr::select( NMDS1, NMDS2, one_of(factorsToTest) ) %>%
    tidyr::gather( factor, subfactor, -NMDS1, -NMDS2 ) %>%
    dplyr::mutate( subfactor = factor(subfactor, levels = factorLevels))
  
  save(nmdsData, file = paste0(intermediateFolder, "/NMDS_Taxon_", taxonomicLevelsToTest[i], "_Abundance.RData"))
  
  pNmds <- ggplot( nmdsData, 
                   aes(x = NMDS1, y = NMDS2, color = subfactor) ) +
    geom_point( size = 2 ) +
    scale_color_manual( values = c(alpha( brewer.pal(5, "Dark2"), 0.8 ), alpha( brewer.pal(6, "Spectral"), 0.8 ) ) ) +
    geom_hline( yintercept = 0, linetype = 2 ) +
    geom_vline( xintercept = 0, linetype = 2 ) +
    ggtitle( paste0(organism, " - ", taxonomicLevelsToTest[i], " - Jaccard | Abundance" ) ) +
    labs( x = "\nNMDS1", y = "NMDS2\n" ) +
    facet_wrap( ~factor, scales = "free" ) +
    theme( legend.position = "right",
           panel.background = element_blank(),
           legend.key=element_blank(),
           plot.title = element_text(hjust = 0.5) )
  
  pdf(file=paste0(outputFolder, "/NMDSGraph_", taxonomicLevelsToTest[i], "_abundance.pdf"), width=15, height=5)
  print(pNmds)
  dev.off()
  
  # Binary
  nmdsDataBinary <- getDataForNMDS( data = taxData,
                                    distance = "jaccard",
                                    minCounts = 1,
                                    minSamplePercent = 0.05,
                                    binary = T) %>% 
    filter( id.type == "Samples" ) %>%
    dplyr::select( NMDS1, NMDS2, one_of(factorsToTest) ) %>%
    tidyr::gather( factor, subfactor, -NMDS1, -NMDS2 ) %>%
    dplyr::mutate( subfactor = factor(subfactor, levels = factorLevels))
  
  save(nmdsData, file = paste0(intermediateFolder, "/NMDS_Taxon_", taxonomicLevelsToTest[i], "_Types.RData"))
  
  pNmdsBinary <- ggplot( nmdsDataBinary, 
                         aes(x = NMDS1, y = NMDS2, color = subfactor) ) +
    geom_point( size = 2 ) +
    scale_color_manual( values = c(alpha( brewer.pal(5, "Dark2"), 0.8 ), alpha( brewer.pal(6, "Spectral"), 0.8 ) ) ) +
    geom_hline( yintercept = 0, linetype = 2 ) +
    geom_vline( xintercept = 0, linetype = 2 ) +
    ggtitle( paste0(organism, " - ", taxonomicLevelsToTest[i], " - Jaccard | Types ") ) +
    labs( x = "\nNMDS1", y = "NMDS2\n" ) +
    facet_wrap( ~factor, scales = "free" ) +
    theme( legend.position = "right",
           panel.background = element_blank(),
           legend.key=element_blank(),
           plot.title = element_text(hjust = 0.5) )
  
  pdf(file=paste0(outputFolder, "/NMDSGraph_", taxonomicLevelsToTest[i], "_types.pdf"), width=15, height=5)
  print(pNmdsBinary)
  dev.off()
  
  
  # Build list for comparison to do at level 2
  factorsToTestLevel2 <- factorsToTest
  
  compList <- list()
  
  for(j in 1:length(factorsToTestLevel2)){
    
    compList[factorsToTestLevel2[j]] <- list(unique(as.character(data.frame(sample_data(taxData))[, factorsToTestLevel2[j]])))
    
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
  
  # Convert data to long format
  dataLong <- sample_data(taxData) %>%
    data.frame() %>%
    dplyr::select( SampleID, one_of(factorsToTest) ) %>%
    tidyr::gather( factor, subfactor, -SampleID )
  
  dataSelectedFactors <- sample_data(taxData) %>%
    data.frame() %>%
    dplyr::select( SampleID, one_of(factorsToTest) )
  
  # Filter data by subfactors and repeat test for interactions and pairwises
  source("supportingFunctions/filterDtaBySubfactors.R")
  
  for(j in 1:length(comparisonsToDo)){
    
    print(paste0("J Value: ", j))
    
    kValues <- 1:nrow(comparisonsToDo[[j]])
    
    for(k in kValues){
      
      print(paste0("J Value: ", j))
      print(paste0("K Value: ", k))
      
      if(ncol(comparisonsToDo[[j]]) == 1){
        
        subfactorsToFilter <- as.character(comparisonsToDo[[j]][k,])
        print(subfactorsToFilter)
        
      } else {
        
        subfactorsToFilter <- (comparisonsToDo[[j]][k,] %>% 
                                 gather( var, value ))$value
        
        print(subfactorsToFilter)
        
      }
      
      filteredData <- filterDataBySubfactors(dataSelectedFactors, dataLong, subfactorsToFilter )
      
      filteredBiomDataPre <- subset_samples(taxData, SampleID %in% filteredData$SampleID)
      
      filteredBiomData <- prune_taxa(taxa_sums(filteredBiomDataPre) > 0, filteredBiomDataPre)
      
      dataProfile <- sample_data(filteredBiomData) %>%
        data.frame() %>%
        dplyr::select( SampleID, Year_of_program, Season, Management, Region) %>%
        tidyr::gather( Variable, Value, -SampleID  ) %>%
        dplyr::group_by( Variable, Value ) %>%
        dplyr::summarise( Samples = n() ) %>%
        dplyr::mutate( subfactorTested = paste(subfactorsToFilter, collapse = "_") )
      
      if(file.exists(paste0(outputFolder, "/sampleProfilePerComparison.csv"))){
        write.table(dataProfile, file = paste0(outputFolder, "/sampleProfilePerComparison.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/sampleProfilePerComparison.csv")), append = T, row.names = F)
      } else {
        write.csv(dataProfile, file = paste0(outputFolder, "/sampleProfilePerComparison.csv"), row.names = F )
      }
      
      message("Going to Adonis based on the profile below:\n")
      
      print(dataProfile)
      
      if(any(dataProfile$Samples < 3)){
        
        flagSampleSize <- "!!!! less than 3 samples per factor level !!!!!"
        
      } else {
        
        flagSampleSize <- "Ok. More than 3 samples per factor level."
        
      }
      
      
      # Save data per taxon to use in the future, e.g. to build the box plots later
      # save(filteredBiomData, file = paste0(intermediateFolder, "/dataAbundanceAndTypeTaxaLevelTaxon_", taxonomicLevels[i], "_", paste(subfactorsToFilter, collapse = "_"), ".RData"))
      
      # Remove the factor under analysis
      factorsToBeRemoved <- dataLong %>%
        dplyr::filter( subfactor %in% subfactorsToFilter )
      
      newFactorList <- factorsToTest[!(factorsToTest %in% unique(factorsToBeRemoved$factor))]
      
      # Generate distance matrix
      distanceDf <- distance(filteredBiomData, method="jaccard", binary = F)
      
      # Generate formula
      formulaToUse <- as.formula(paste("distanceDf", paste(newFactorList, collapse=" * "), sep=" ~ "))
      
      # Generate distance matrix binary
      distanceDfBinary <- distance(filteredBiomData, method="jaccard", binary = T)
      
      # Generate formula binary
      formulaToUseBinary <- as.formula(paste("distanceDfBinary", paste(newFactorList, collapse=" * "), sep=" ~ "))
      
      # Run adonis
      set.seed(4321)
      message( "Start Adonis....")
      adonisGeneral <- data.frame((adonis2( formulaToUse,
                                           data = data.frame(sample_data(filteredBiomData))))) %>%
        rownames_to_column( "factors" ) %>%
        dplyr::mutate( taxa = taxonomicLevelsToTest[i],
                       type = "jaccard-abundance",
                       subfactorTested = paste(subfactorsToFilter, collapse = "_"),
                       flag = flagSampleSize ) %>%
        dplyr::bind_rows( data.frame((adonis2( formulaToUseBinary,
                                              data = data.frame(sample_data(filteredBiomData))))) %>%
                            rownames_to_column( "factors" ) %>%
                            dplyr::mutate( taxa = taxonomicLevelsToTest[i],
                                           type = "jaccard-types",
                                           subfactorTested = paste(subfactorsToFilter, collapse = "_"),
                                           flag = flagSampleSize ) )
      
      if(file.exists(paste0(outputFolder, "/analysis_abundance_and_type_adonis.csv"))){
        write.table(adonisGeneral, file = paste0(outputFolder, "/analysis_abundance_and_type_adonis.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/analysis_abundance_and_type_adonis.csv")), append = T, row.names = F)
      } else {
        write.csv(adonisGeneral, file = paste0(outputFolder, "/analysis_abundance_and_type_adonis.csv"), row.names = F )
      }
      
      
      # Generate NMDS graph
      source("supportingFunctions/getDataForNMDS.R")
      
      pdfSize <- switch( length(subfactorsToFilter),
                         c(12, 6),
                         c(10, 4),
                         c(4, 4) )
      
      set.seed(4321)
      message("NMDS non-binary...")
      message(paste("Testing taxonomic level:", taxonomicLevels[i]))
      nmdsData <- getDataForNMDS( data = filteredBiomData,
                                  distance = "jaccard",
                                  minCounts = 1,
                                  minSamplePercent = 0.05,
                                  binary = F) %>% 
        filter( id.type == "Samples" ) %>%
        dplyr::select( NMDS1, NMDS2, one_of(newFactorList) ) %>%
        tidyr::gather( factor, subfactor, -NMDS1, -NMDS2 ) %>%
        dplyr::mutate( subfactor = factor(subfactor, levels = factorLevels))
      
      
      pNmds <- ggplot( nmdsData, 
                       aes(x = NMDS1, y = NMDS2, color = subfactor) ) +
        geom_point( size = 2 ) +
        scale_color_manual( values = c(alpha( brewer.pal(5, "Dark2"), 0.8 ), alpha( brewer.pal(6, "Spectral"), 0.8 ) ) ) +
        geom_hline( yintercept = 0, linetype = 2 ) +
        geom_vline( xintercept = 0, linetype = 2 ) +
        ggtitle( paste0(targetOrganisms, " - ", taxonomicLevels[i], " - Jaccard | Abundance | Factor = ", paste(subfactorsToFilter, collapse = "_")) ) +
        labs( x = "\nNMDS1", y = "NMDS2\n" ) +
        facet_wrap( ~factor, scales = "free", nrow = 1) +
        theme( legend.position = "right",
               panel.background = element_blank(),
               legend.key=element_blank(),
               plot.title = element_text(hjust = 0.5) )
      
      pdf(file=paste0(outputFolder, "/NMDSGraph_", taxonomicLevels[i], "_", paste(subfactorsToFilter, collapse = "_"), "_abundance.pdf"), width=pdfSize[1], height=pdfSize[2])
      print(pNmds)
      dev.off()
      
      # Binary
      
      if( ntaxa(filteredBiomData) > 3 ){
        message("NMDS Binary...")
        message(paste("Testing taxonomic level:", taxonomicLevels[i]))
        
        testnmds <- try(getDataForNMDS( data = filteredBiomData,
                                 distance = "jaccard",
                                 minCounts = 1,
                                 minSamplePercent = 0.05,
                                 binary = T), silent = T)
        
        if(class(testnmds) != "try-error"){
          
          message(paste("Testing taxonomic level:", taxonomicLevels[i], "- ERROR!"))
          
          nmdsDataBinary <- getDataForNMDS( data = filteredBiomData,
                                            distance = "jaccard",
                                            minCounts = 1,
                                            minSamplePercent = 0.05,
                                            binary = T) %>% 
            filter( id.type == "Samples" ) %>%
            dplyr::select( NMDS1, NMDS2, one_of(newFactorList) ) %>%
            tidyr::gather( factor, subfactor, -NMDS1, -NMDS2 ) %>%
            dplyr::mutate( subfactor = factor(subfactor, levels = factorLevels))
          
          pNmdsBinary <- ggplot( nmdsDataBinary, 
                                 aes(x = NMDS1, y = NMDS2, color = subfactor) ) +
            geom_point( size = 2 ) +
            scale_color_manual( values = c(alpha( brewer.pal(5, "Dark2"), 0.8 ), alpha( brewer.pal(6, "Spectral"), 0.8 ) ) ) +
            geom_hline( yintercept = 0, linetype = 2 ) +
            geom_vline( xintercept = 0, linetype = 2 ) +
            ggtitle( paste0(targetOrganisms, " - ", taxonomicLevels[i], " - Jaccard | Types | Factor = ", paste(subfactorsToFilter, collapse = "_")) ) +
            labs( x = "\nNMDS1", y = "NMDS2\n" ) +
            facet_wrap( ~factor, scales = "free", nrow = 1) +
            theme( legend.position = "right",
                   panel.background = element_blank(),
                   legend.key=element_blank(),
                   plot.title = element_text(hjust = 0.5) )
          
          pdf(file=paste0(outputFolder, "/NMDSGraph_", taxonomicLevels[i], "_", paste(subfactorsToFilter, collapse = "_"), "_types.pdf"), width=pdfSize[1], height=pdfSize[2])
          print(pNmdsBinary)
          dev.off()
          
        } else {
          
          pdf(file=paste0(outputFolder, "/NMDSGraph_", taxonomicLevels[i], "_", paste(subfactorsToFilter, collapse = "_"), "_types_ERROR.pdf"), width=pdfSize[1], height=pdfSize[2])
          print("NMDS error for not enough taxa after filtering: Error in monoMDS(dist, y = cmdscale(dist, k = k), k = k, maxit = maxit. 'dist' cannot be all zero (all points are identical)")
          dev.off()
          
        }
      }
    }
  }
}







