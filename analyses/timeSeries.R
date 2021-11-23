
#
##
###
####
##### Time series
####
###
##
#

# 1 - All variables are run at the same time: Number of tax, Shannon, Simpson, Air temperature, etc...

rm(list = ls())  

#
# Choose the normalisation method to use: "rarefy", "bc" [bias correction] or "notNormalised".
#
normaliseMethod <- "rarefy"

# Organisms to be analysed:
organismsToBuildTimeSeries <- c("all", "bacteria", "fungi", "metazoa")

# For each organism 
for( organism in organismsToBuildTimeSeries ) {
  
  # Create output folder
  if(!file.exists(paste0("output/", organism, "/", normaliseMethod))){
    dir.create(paste0("output/", organism, "/", normaliseMethod), recursive = T)
  }
  
  outputFolder <- paste0("output/", organism, "/", normaliseMethod)
  
  # Create a special folder for time series
  if(!file.exists(paste0(outputFolder, "/a_timeSeries/"))){
    dir.create(paste0(outputFolder, "/a_timeSeries/"), recursive = T)
  }
  
  # Create intermediate folder
  if(!file.exists(paste0("intermediateData/", organism, "/", normaliseMethod))){
    dir.create(paste0("intermediateData/", organism, "/", normaliseMethod), recursive = T)
  }
  
  intermediateFolder <- paste0("intermediateData/", organism, "/", normaliseMethod)
  
  # Create a special folder for time series
  if(!file.exists(paste0(intermediateFolder, "/a_timeSeries/"))){
    dir.create(paste0(intermediateFolder, "/a_timeSeries/"), recursive = T)
  }
  
  
  switch(organism,
         bacteria = {
           taxonomicLevelsToBuildTimeSeries <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "OTUs")
         },
         fungi = {
           taxonomicLevelsToBuildTimeSeries <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "OTUs")
         },
         metazoa = {
           taxonomicLevelsToBuildTimeSeries <- c("OTUs")
         },
         all = {
           taxonomicLevelsToBuildTimeSeries <- c("OTUs")
         } 
  )
  
  # For each taxonomic level
  for( taxonomicLevel in taxonomicLevelsToBuildTimeSeries ) {
    
    # Get data number of taxa
    load(paste0(intermediateFolder, "/dataNumberOfTaxaLevelTaxon_", taxonomicLevel, ".RData"))
    
    # Get data 
    load(paste0(intermediateFolder, "/datashannonAndSimpsonLevelTaxon_", taxonomicLevel, ".RData"))
    
    
    
    #
    #
    #
    # Make sure the geom_smooth is ploting the same model as we are using for ANOVA
    #
    #
    #
    
    # # Get data only for HB and C
    # testData1 <- dataNumberOfTaxa %>%
    #   dplyr::filter( Region == "HB",
    #                  Management == "C" )
    # 
    # # Build model
    # model1 <- glm(numberOfTaxa ~ Timepoint, data = testData1, family = quasipoisson(link = "log"))
    # summary(model1)
    # 
    # # Get predicted values for the Timepoint we have, 1 to 9
    # modelToPlot <- data.frame( y = predict(model1, list(Timepoint = 1:9), type="response"), x = 1:9) %>%
    #   dplyr::mutate( Management = "Model" )
    # 
    # # Plot data with: data points, average and glm model from function geom_smooth
    # (p <- ggplot(dataNumberOfTaxa %>%
    #               dplyr::filter( Region == "HB"),
    #        aes(y=as.numeric(numberOfTaxa), x=Timepoint, color=Management) ) +
    #   geom_point() +
    #   geom_line( data = dataNumberOfTaxa %>%
    #                dplyr::filter( Region == "HB" ) %>%
    #                dplyr::group_by( Timepoint, Management ) %>%
    #                dplyr::summarise( numberOfTaxa = mean(numberOfTaxa, na.rm = T)), 
    #              aes(y=as.numeric(numberOfTaxa), x=Timepoint, color=Management), linetype=2, size = 2 ) +
    #   geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
    #   scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) )
    # 
    # # Now, we plot the data from the model we manually produced to run ANOVA - Note that the data we used here are for HB and C only. 
    # # This means that the line for the manual model should fit precisely on top of the glm model from geom_smooth for Management type C. 
    # (p2 <- p +
    #   geom_line(data = modelToPlot, aes(y = y, x=x, color=Management), linetype = 2 ) ) 
    
    
    
    #
    #
    #
    #
    # Run actual time series tests
    #
    #
    #
    #
    
    
    # Merge Number of Taxa, Simpson, Shannon
    
    dataToRun <- dataNumberOfTaxa %>%
      dplyr::left_join( shannonAndSimpson %>%
                          dplyr::select( SampleID, Shannon, Simpson, Chao1),
                        by = "SampleID")
    
    # Generate a statistical summary for each potential combimation between Management, Region and Year of program
    factorsToTest <- c("Management", "Region", "Season", "Year_of_program")
    variablesToProcess <- c("numberOfTaxa", "Shannon", "Simpson", "Chao1", "AnMiN", "OlsP", "pH", "PotN", "TotC", "TotN", "VolW",
                            "AvgTempAirMax", "AvgTempAirMin", "AvgTempAirRange", "AvgTempAir", "AvgRainfall")
    
    compList <- list()
    
    for(j in 1:length(factorsToTest)){
      
      compList[factorsToTest[j]] <- list(unique(as.character(dataToRun[, factorsToTest[j]])))
      
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
    dataLong <- dataToRun %>%
      dplyr::select( SampleID, one_of(variablesToProcess), one_of(factorsToTest) ) %>%
      tidyr::gather( factor, subfactor, -SampleID, -one_of(variablesToProcess) )
    
    dataSelected <- dataToRun %>%
      dplyr::select( SampleID, one_of(variablesToProcess), one_of(factorsToTest) )
    
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
        
        filteredData <- filterDataBySubfactors(dataSelected, dataLong, subfactorsToFilter)
        
        
        # Generate stats
        fullStats <- round(stat.desc( filteredData %>%
                                        dplyr::select( all_of(variablesToProcess) ) ), 2) %>%
          data.frame() %>%
          rownames_to_column("Statistic") %>%
          dplyr::mutate( taxon = taxonomicLevel,
                         organism = organism,
                         subfactorTested = paste(subfactorsToFilter, collapse = "_") ) 
        
        if(file.exists(paste0(outputFolder, "/fullStatsByFactor.csv"))){
          write.table(fullStats, file = paste0(outputFolder, "/fullStatsByFactor.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/fullStatsByFactor.csv")), append = T, row.names = F)
        } else {
          write.csv(fullStats, file = paste0(outputFolder, "/fullStatsByFactor.csv"), row.names = F )
        }
        
        
        
        
      }
      
    }
    
    
    
    
    
    
    
    
    
    
    
    #
    # Comparison 1: Per region, does management make any difference?
    #
    
    comparingHere <- c("HB", "MB")
    
    # Defined above
    # variablesToProcess <- c("numberOfTaxa", "Shannon", "Simpson", "AnMiN", "OlsP", "pH", "PotN", "TotC", "TotN", "VolW",
    #                         "AvgTempAirMax", "AvgTempAirMin", "AvgTempAirRange", "AvgTempAir", "AvgRainfall")
    
    dataToRunSelected <- dataToRun %>%
      dplyr::select( all_of(c(variablesToProcess, "Region", "Management", "Timepoint")) )
    
    # For each area, HB or MB
    for(areaToCollect in comparingHere){
      
      testData1 <- dataToRunSelected %>%
        dplyr::filter( Region == areaToCollect )
      
      # Produce plot for all variablesToProcess
      dataToPlot <- testData1 %>%
        # Select the variables we want to plot
        dplyr::select( Management, Timepoint, all_of(variablesToProcess) ) %>%
        # Make it long format
        tidyr::gather( Variable, Value, -Management, -Timepoint ) %>%
        dplyr::mutate( Variable = factor(Variable, levels = variablesToProcess) )
      
      pAllFactors <- ggplot( dataToPlot, aes(x = Timepoint, y = as.numeric(Value), color = Management) ) +
        geom_point(alpha = 0.5) +
        geom_line( data = dataToPlot %>%
                     dplyr::group_by( Timepoint, Management, Variable ) %>%
                     dplyr::summarise( Value = mean(Value, na.rm = T)),
                   aes(y=as.numeric(Value), x=Timepoint, color=Management), linetype=2, size = 1 ) +
        geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
        scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
        labs( y = "Variable in test\n", x = "Timepoint" ) +
        theme_minimal() +
        theme( legend.position = "top" ) + 
        facet_wrap( ~ Variable, scales = "free" ) 
      
      # Save dataset used to build the plot
      save(dataToPlot, file = paste0(intermediateFolder, "/a_timeSeries/dataToPlotAllFactors_Area_", areaToCollect, "_Taxon_", taxonomicLevel, ".RData"))
      
      # Save dataset used to build the plot
      save(pAllFactors, file = paste0(intermediateFolder, "/a_timeSeries/plotAllFactors_Area_", areaToCollect, "_Taxon_", taxonomicLevel, ".RData"))
      
      # Save plot to PDF
      pdf(paste0(outputFolder, "/a_timeSeries/plotAllFactors_Area_", areaToCollect, "_Taxon_", taxonomicLevel, ".pdf"), width = 15, height = 10)
      
      print(pAllFactors)
      
      dev.off()
      
      
      for(variableInTest in variablesToProcess){
        
        # Build model 1 - by time only
        formulaToUse1 <- formula(paste0(variableInTest, " ~ Timepoint"))
        model1 <- glm(formula = formulaToUse1, data = testData1, family = quasipoisson(link = "log"))
        
        # Build model 2 - by time and management
        formulaToUse2 <- formula(paste0(variableInTest, " ~ Timepoint + Management"))
        model2 <- glm(formula = formulaToUse2, data = testData1, family = quasipoisson(link = "log"))
        
        # Calculate ANOVA
        resultsAnova <- anova(model1, model2, test = "F")
        
        # Get numbers that we will save
        statsNumbersToSave <- data.frame(organism = organism,
                                         taxonomicLevel = taxonomicLevel,
                                         samplesFromSpecificArea = areaToCollect,
                                         variableInTest = variableInTest,
                                         pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                         fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                         deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                         df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                         residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                         residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
        
        # Save stats
        if(file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv"))){
          write.table(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv")), append = T, row.names = F, )
        } else {
          write.csv(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), row.names = F )
        }
        
        # Save the plot to use later
        dataToPlotVariable <- testData1 %>%
          dplyr::select( Management, Timepoint, "Value" = all_of(variableInTest) )
        
        pIndividual <- ggplot(dataToPlotVariable,
                              aes(y=as.numeric(Value), x=Timepoint, color=Management) ) +
          geom_point() +
          geom_line( data = dataToPlotVariable %>%
                       dplyr::group_by( Timepoint, Management ) %>%
                       dplyr::summarise( Value = mean(Value, na.rm = T)), 
                     aes(y=as.numeric(Value), x=Timepoint, color=Management), linetype=2, size = 1 ) +
          geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
          scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
          labs( y = as.character(variableInTest), x = "" ) +
          theme_minimal( ) +
          theme( legend.position = "top" )
        
        
        # Save plot to PDF
        pdf(paste0(outputFolder, "/a_timeSeries/plotIndividualFactor_", variableInTest, "_Area_", areaToCollect, "_Taxon_", taxonomicLevel, ".pdf"), width = 15, height = 10)
        
        print(pIndividual)
        
        dev.off()
        
        # Save dataset used to build the plot
        save(dataToPlotVariable, file = paste0(intermediateFolder, "/a_timeSeries/dataToplotIndividualFactor_", variableInTest, "_Area_", areaToCollect, "_Taxon_", taxonomicLevel, ".RData"))
        
        # Save dataset used to build the plot
        save(pIndividual, file = paste0(intermediateFolder, "/a_timeSeries/plotIndividualFactor_", variableInTest, "_Area_", areaToCollect, "_Taxon_", taxonomicLevel, ".RData"))
        
      }
      
    }
    
    
    
    
    
    
    
    #
    # Comparison 2: Does Region makes any difference?
    #
    # variablesToProcess <- c("numberOfTaxa", "Shannon", "Simpson", "AnMiN", "OlsP", "pH", "PotN", "TotC", "TotN", "VolW",
    #                         "AvgTempAirMax", "AvgTempAirMin", "AvgTempAirRange", "AvgTempAir", "AvgRainfall")
    
    testData1 <- dataToRunSelected
    
    # Produce plot for all variablesToProcess
    dataToPlot <- testData1 %>%
      # Select the variables we want to plot
      dplyr::select(Region, Timepoint, all_of(variablesToProcess) ) %>%
      # Make it long format
      tidyr::gather( Variable, Value, -Region, -Timepoint ) %>%
      dplyr::mutate( Variable = factor(Variable, levels = variablesToProcess) )
    
    pAllFactors <- ggplot( dataToPlot, aes(x = Timepoint, y = as.numeric(Value), color = Region) ) +
      geom_point(alpha = 0.5) +
      geom_line( data = dataToPlot %>%
                   dplyr::group_by( Timepoint, Region, Variable ) %>%
                   dplyr::summarise( Value = mean(Value, na.rm = T)),
                 aes(y=as.numeric(Value), x=Timepoint, color=Region), linetype=2, size = 1 ) +
      geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
      scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
      labs( y = "Variable in test\n", x = "Timepoint" ) +
      theme_minimal() +
      theme( legend.position = "top" ) + 
      facet_wrap( ~ Variable, scales = "free" ) 
    
    # Save dataset used to build the plot
    save(dataToPlot, file = paste0(intermediateFolder, "/a_timeSeries/dataToPlotAllFactors_CompareRegion_Taxon_", taxonomicLevel, ".RData"))
    
    # Save dataset used to build the plot
    save(pAllFactors, file = paste0(intermediateFolder, "/a_timeSeries/plotAllFactors_CompareRegion_Taxon_", taxonomicLevel, ".RData"))
    
    # Save plot to PDF
    pdf(paste0(outputFolder, "/a_timeSeries/plotAllFactors_CompareRegion_Taxon_", taxonomicLevel, ".pdf"), width = 15, height = 10)
    
    print(pAllFactors)
    
    dev.off()
    
    
    for(variableInTest in variablesToProcess){
      
      # Build model 1 - by time only
      formulaToUse1 <- formula(paste0(variableInTest, " ~ Timepoint"))
      model1 <- glm(formula = formulaToUse1, data = testData1, family = quasipoisson(link = "log"))
      
      # Build model 2 - by time and management
      formulaToUse2 <- formula(paste0(variableInTest, " ~ Timepoint + Region"))
      model2 <- glm(formula = formulaToUse2, data = testData1, family = quasipoisson(link = "log"))
      
      # Calculate ANOVA
      resultsAnova <- anova(model1, model2, test = "F")
      
      # Get numbers that we will save
      statsNumbersToSave <- data.frame(organism = organism,
                                       taxonomicLevel = taxonomicLevel,
                                       samplesFromSpecificArea = "Region",
                                       variableInTest = variableInTest,
                                       pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                       fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                       deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                       df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                       residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                       residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
      
      # Save stats
      if(file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv"))){
        write.table(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv")), append = T, row.names = F, )
      } else {
        write.csv(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), row.names = F )
      }
      
      # Save the plot to use later
      dataToPlotVariable <- testData1 %>%
        dplyr::select( Region, Timepoint, "Value" = all_of(variableInTest) )
      
      pIndividual <- ggplot(dataToPlotVariable,
                            aes(y=as.numeric(Value), x=Timepoint, color=Region) ) +
        geom_point() +
        geom_line( data = dataToPlotVariable %>%
                     dplyr::group_by( Timepoint, Region ) %>%
                     dplyr::summarise( Value = mean(Value, na.rm = T)), 
                   aes(y=as.numeric(Value), x=Timepoint, color=Region), linetype=2, size = 1 ) +
        geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
        scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
        labs( y = as.character(variableInTest), x = "" ) +
        theme_minimal( ) +
        theme( legend.position = "top" )
      
      
      # Save plot to PDF
      pdf(paste0(outputFolder, "/a_timeSeries/plotIndividualFactor_", variableInTest, "_CompareRegion_Taxon_", taxonomicLevel, ".pdf"), width = 15, height = 10)
      
      print(pIndividual)
      
      dev.off()
      
      # Save dataset used to build the plot
      save(dataToPlotVariable, file = paste0(intermediateFolder, "/a_timeSeries/dataToplotIndividualFactor_", variableInTest, "_CompareRegion_Taxon_", taxonomicLevel, ".RData"))
      
      # Save dataset used to build the plot
      save(pIndividual, file = paste0(intermediateFolder, "/a_timeSeries/plotIndividualFactor_", variableInTest, "_CompareRegion_Taxon_", taxonomicLevel, ".RData"))
      
    }
    
    
    
    
    
    
    #
    # Comparison 3: Does Management makes any difference?
    #
    # variablesToProcess <- c("numberOfTaxa", "AnMiN", "OlsP", "pH", "PotN", "TotC", "TotN", "VolW",
    #                         "AvgTempAirMax", "AvgTempAirMin", "AvgTempAirRange", "AvgTempAir", "AvgRainfall")
    
    testData1 <- dataToRunSelected
    
    # Produce plot for all variablesToProcess
    dataToPlot <- testData1 %>%
      # Select the variables we want to plot
      dplyr::select(Management, Timepoint, all_of(variablesToProcess) ) %>%
      # Make it long format
      tidyr::gather( Variable, Value, -Management, -Timepoint ) %>%
      dplyr::mutate( Variable = factor(Variable, levels = variablesToProcess) )
    
    pAllFactors <- ggplot( dataToPlot, aes(x = Timepoint, y = as.numeric(Value), color = Management) ) +
      geom_point(alpha = 0.5) +
      geom_line( data = dataToPlot %>%
                   dplyr::group_by( Timepoint, Management, Variable ) %>%
                   dplyr::summarise( Value = mean(Value, na.rm = T)),
                 aes(y=as.numeric(Value), x=Timepoint, color=Management), linetype=2, size = 1 ) +
      geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
      scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
      labs( y = "Variable in test\n", x = "Timepoint" ) +
      theme_minimal() +
      theme( legend.position = "top" ) + 
      facet_wrap( ~ Variable, scales = "free" ) 
    
    # Save dataset used to build the plot
    save(dataToPlot, file = paste0(intermediateFolder, "/a_timeSeries/dataToPlotAllFactors_CompareManagement_Taxon_", taxonomicLevel, ".RData"))
    
    # Save dataset used to build the plot
    save(pAllFactors, file = paste0(intermediateFolder, "/a_timeSeries/plotAllFactors_CompareManagement_Taxon_", taxonomicLevel, ".RData"))
    
    # Save plot to PDF
    pdf(paste0(outputFolder, "/a_timeSeries/plotAllFactors_CompareManagement_Taxon_", taxonomicLevel, ".pdf"), width = 15, height = 10)
    
    print(pAllFactors)
    
    dev.off()
    
    
    for(variableInTest in variablesToProcess){
      
      # Build model 1 - by time only
      formulaToUse1 <- formula(paste0(variableInTest, " ~ Timepoint"))
      model1 <- glm(formula = formulaToUse1, data = testData1, family = quasipoisson(link = "log"))
      
      # Build model 2 - by time and management
      formulaToUse2 <- formula(paste0(variableInTest, " ~ Timepoint + Management"))
      model2 <- glm(formula = formulaToUse2, data = testData1, family = quasipoisson(link = "log"))
      
      # Calculate ANOVA
      resultsAnova <- anova(model1, model2, test = "F")
      
      # Get numbers that we will save
      statsNumbersToSave <- data.frame(organism = organism,
                                       taxonomicLevel = taxonomicLevel,
                                       samplesFromSpecificArea = "Management",
                                       variableInTest = variableInTest,
                                       pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                       fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                       deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                       df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                       residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                       residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
      
      # Save stats
      if(file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv"))){
        write.table(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv")), append = T, row.names = F, )
      } else {
        write.csv(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), row.names = F )
      }
      
      # Save the plot to use later
      dataToPlotVariable <- testData1 %>%
        dplyr::select( Management, Timepoint, "Value" = all_of(variableInTest) )
      
      pIndividual <- ggplot(dataToPlotVariable,
                            aes(y=as.numeric(Value), x=Timepoint, color=Management) ) +
        geom_point() +
        geom_line( data = dataToPlotVariable %>%
                     dplyr::group_by( Timepoint, Management ) %>%
                     dplyr::summarise( Value = mean(Value, na.rm = T)), 
                   aes(y=as.numeric(Value), x=Timepoint, color=Management), linetype=2, size = 1 ) +
        geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
        scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
        labs( y = as.character(variableInTest), x = "" ) +
        theme_minimal( ) +
        theme( legend.position = "top" )
      
      
      # Save plot to PDF
      pdf(paste0(outputFolder, "/a_timeSeries/plotIndividualFactor_", variableInTest, "_CompareManagement_Taxon_", taxonomicLevel, ".pdf"), width = 15, height = 10)
      
      print(pIndividual)
      
      dev.off()
      
      # Save dataset used to build the plot
      save(dataToPlotVariable, file = paste0(intermediateFolder, "/a_timeSeries/dataToplotIndividualFactor_", variableInTest, "_CompareManagement_Taxon_", taxonomicLevel, ".RData"))
      
      # Save dataset used to build the plot
      save(pIndividual, file = paste0(intermediateFolder, "/a_timeSeries/plotIndividualFactor_", variableInTest, "_CompareManagement_Taxon_", taxonomicLevel, ".RData"))
      
    }
    
    
    #
    ## Build Temporal Graph from Jaccard 
    #
    
    # Load file
    load(file = paste0(intermediateFolder, "/dataTypesAndAbundancesOfTaxaLevelTaxon_", taxonomicLevel, ".RData"))
    
    # Get data for comparison 1 and 2: (1) Per region, does management make any difference?; (2) Is there difference between regions?
    comparingHere <- c("HB", "MB")
    
    timePointsToUse <- c(1:9)
    
    for(areaToCollect in comparingHere){
      
      filteredTaxDataArea <- subset_samples(taxData, Region == areaToCollect )
      
      for(selectedPoint in timePointsToUse){
        
        filteredTaxDataAreaTimePoint <- subset_samples(filteredTaxDataArea, Timepoint == selectedPoint )
        
        distanceDf <- as.data.frame(as.matrix(distance(filteredTaxDataAreaTimePoint, method="jaccard", binary = F)))
        
        dissimilarityData <- data.frame( organism = organism, taxonomicLevel = taxonomicLevel,
                                         samplesFromSpecificArea = areaToCollect,
                                         timePoint = selectedPoint, management = "all",
                                         dissIndex = distanceDf[-1,1], type = "noBinary")
        
        # Save index
        if(file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"))){
          write.table(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv")), append = T, row.names = F, )
        } else {
          write.csv(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), row.names = F )
        }
        
        # Binary
        distanceDf <- as.data.frame(as.matrix(distance(filteredTaxDataAreaTimePoint, method="jaccard", binary = T)))
        
        dissimilarityData <- data.frame( organism = organism, taxonomicLevel = taxonomicLevel,
                                         samplesFromSpecificArea = areaToCollect,
                                         timePoint = selectedPoint, management = "all",
                                         dissIndex = distanceDf[-1,1], type = "binary")
        
        # Save index
        if(file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"))){
          write.table(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv")), append = T, row.names = F, )
        } else {
          write.csv(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), row.names = F )
        }
        
        for(managementToUse in unique(sample_data(filteredTaxDataAreaTimePoint)$Management)){
          
          filteredTaxDataAreaTimePointManagement <- subset_samples(filteredTaxDataAreaTimePoint, Management == managementToUse )
          
          distanceDf <- as.data.frame(as.matrix(distance(filteredTaxDataAreaTimePointManagement, method="jaccard", binary = F)))
          
          dissimilarityData <- data.frame( organism = organism, taxonomicLevel = taxonomicLevel,
                                           samplesFromSpecificArea = areaToCollect,
                                           timePoint = selectedPoint, management = managementToUse,
                                           dissIndex = distanceDf[-1,1], type = "noBinary")
          
          # Save index
          if(file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"))){
            write.table(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv")), append = T, row.names = F, )
          } else {
            write.csv(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), row.names = F )
          }
          
          # Binary
          distanceDf <- as.data.frame(as.matrix(distance(filteredTaxDataAreaTimePointManagement, method="jaccard", binary = T)))
          
          dissimilarityData <- data.frame( organism = organism, taxonomicLevel = taxonomicLevel,
                                           samplesFromSpecificArea = areaToCollect,
                                           timePoint = selectedPoint, management = managementToUse,
                                           dissIndex = distanceDf[-1,1], type = "binary")
          
          # Save index
          if(file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"))){
            write.table(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv")), append = T, row.names = F, )
          } else {
            write.csv(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), row.names = F )
          }
          
          
        }
        
      }
      
    }
    
    
    # Get data for comparison 3: (3) Is there difference between managements?
    timePointsToUse <- c(1:9)
    
    for(managementToUse in unique(sample_data(taxData)$Management)){
      
      filteredTaxDataManagement <- subset_samples(taxData, Management == managementToUse )
      
      for(selectedPoint in timePointsToUse){
        
        filteredTaxDataManagementTimePoint <- subset_samples(filteredTaxDataManagement, Timepoint == selectedPoint )
        
        distanceDf <- as.data.frame(as.matrix(distance(filteredTaxDataManagementTimePoint, method="jaccard", binary = F)))
        
        dissimilarityData <- data.frame( organism = organism, taxonomicLevel = taxonomicLevel,
                                         samplesFromSpecificArea = "Management",
                                         timePoint = selectedPoint, management = managementToUse,
                                         dissIndex = distanceDf[-1,1], type = "noBinary")
        
        # Save index
        if(file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"))){
          write.table(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv")), append = T, row.names = F, )
        } else {
          write.csv(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), row.names = F )
        }
        
        # Binary
        distanceDf <- as.data.frame(as.matrix(distance(filteredTaxDataManagementTimePoint, method="jaccard", binary = T)))
        
        dissimilarityData <- data.frame( organism = organism, taxonomicLevel = taxonomicLevel,
                                         samplesFromSpecificArea = "Management",
                                         timePoint = selectedPoint, management = managementToUse,
                                         dissIndex = distanceDf[-1,1], type = "binary")
        
        # Save index
        if(file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"))){
          write.table(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv")), append = T, row.names = F, )
        } else {
          write.csv(dissimilarityData, file = paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), row.names = F )
        }
        
      }
      
    }
    
  }
  
}




# 3 - Jaccard Abundance
# The analysis is done using the file /a_timeSeries/dissimiralityIndexes.csv
rm(list = ls())

# Choose normalised method
normaliseMethod <- "rarefy"

# Organisms to be analysed:
organismsToBuildTimeSeries <- c("all", "bacteria", "fungi", "metazoa")
# organismsToBuildTimeSeries <- c("bacteria")

# Taxonomic levels to be analysed
#taxonomicLevelsToBuildTimeSeries <- c("OTUs", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#taxonomicLevelsToBuildTimeSeries <- c("OTUs")


# For each organism 
for( organism in organismsToBuildTimeSeries ) {
  
  # Create output folder
  if(!file.exists(paste0("output/", organism, "/", normaliseMethod))){
    dir.create(paste0("output/", organism, "/", normaliseMethod), recursive = T)
  }
  
  outputFolder <- paste0("output/", organism, "/", normaliseMethod)
  
  # Create a special folder for time series
  if(!file.exists(paste0(outputFolder, "/a_timeSeries/"))){
    dir.create(paste0(outputFolder, "/a_timeSeries/"), recursive = T)
  }
  
  # Create intermediate folder
  if(!file.exists(paste0("intermediateData/", organism, "/", normaliseMethod))){
    dir.create(paste0("intermediateData/", organism, "/", normaliseMethod), recursive = T)
  }
  
  intermediateFolder <- paste0("intermediateData/", organism, "/", normaliseMethod)
  
  # Create a special folder for time series
  if(!file.exists(paste0(intermediateFolder, "/a_timeSeries/"))){
    dir.create(paste0(intermediateFolder, "/a_timeSeries/"), recursive = T)
  }
  
  switch(organism,
         bacteria = {
           taxonomicLevelsToBuildTimeSeries <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "OTUs")
         },
         fungi = {
           taxonomicLevelsToBuildTimeSeries <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "OTUs")
         },
         metazoa = {
           taxonomicLevelsToBuildTimeSeries <- c("OTUs")
         },
         all = {
           taxonomicLevelsToBuildTimeSeries <- c("OTUs")
         } 
  )
  
  # For each taxonomic level
  for( taxonomicLevel in taxonomicLevelsToBuildTimeSeries ) {
    
    #
    # Read data
    #
    
    #
    ## Analysis - Compare C and F for HB and MB
    
    areasToCompare <- c("HB", "MB")
    
    for(areaToUse in areasToCompare){
      
      for(typeToUse in c("binary", "noBinary")){
        
        rawData <- read.csv(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), colClasses = "character") %>%
          dplyr::rename( "organismInAnalysis" = organism,
                         "taxLevelInAnalysis" = taxonomicLevel ) %>%
          dplyr::filter( organismInAnalysis == organism,
                         taxLevelInAnalysis == taxonomicLevel,
                         samplesFromSpecificArea == areaToUse,
                         management %in% c("C", "F"),
                         type == typeToUse ) %>%
          dplyr::mutate( timePoint = as.numeric(timePoint),
                         dissIndex = as.numeric(dissIndex) )
        
        # Plot data
        p <- ggplot( rawData, aes(x = timePoint, y = dissIndex, colour = management) ) +
          geom_point(alpha = 0.5) +
          geom_line( data = rawData %>%
                       dplyr::group_by(timePoint, management) %>%
                       dplyr::summarise( dissIndex = mean(dissIndex, na.rm = T) ), linetype=2, size = 1 ) +
          geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
          scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
          labs( y = "Dissimilarity - Jaccard\n", x = "" ) +
          theme_minimal( ) +
          theme( legend.position = "top" )
        
        # Save plot to PDF
        pdf(paste0(outputFolder, "/a_timeSeries/jaccardIndex_organism_", 
                   organism, "_taxa_", taxonomicLevel, "_Area_", areaToUse, "_type_", typeToUse, ".pdf"), width = 15, height = 10)
        
        print(p)
        
        dev.off()
        
        # Run GLM
        
        # Build model 1 - by time only
        formulaToUse1 <- formula("dissIndex ~ timePoint")
        model1 <- glm(formula = formulaToUse1, data = rawData, family = quasipoisson(link = "log"))
        
        # Build model 2 - by time and management
        formulaToUse2 <- formula("dissIndex ~ timePoint + management")
        model2 <- glm(formula = formulaToUse2, data = rawData, family = quasipoisson(link = "log"))
        
        # Calculate ANOVA
        resultsAnova <- anova(model1, model2, test = "F")
        
        # Get numbers that we will save
        statsNumbersToSave <- data.frame(organism = organism,
                                         taxonomicLevel = taxonomicLevel,
                                         samplesFromSpecificArea = areaToUse,
                                         variableInTest = paste0("JaccardIndex", typeToUse),
                                         pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                         fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                         deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                         df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                         residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                         residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
        
        # Save stats
        if(file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv"))){
          write.table(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv")), append = T, row.names = F, )
        } else {
          write.csv(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), row.names = F )
        }
        
      }
      
    }
    
    
    for(typeToUse in c("binary", "noBinary")){
      # Analysis - HB vs. MB
      rawData <- read.csv(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), colClasses = "character") %>%
        dplyr::rename( "organismInAnalysis" = organism,
                       "taxLevelInAnalysis" = taxonomicLevel ) %>%
        dplyr::filter( organismInAnalysis == organism,
                       taxLevelInAnalysis == taxonomicLevel,
                       samplesFromSpecificArea %in% c("HB", "MB"),
                       management == "all",
                       type == typeToUse ) %>%
        dplyr::mutate( timePoint = as.numeric(timePoint),
                       dissIndex = as.numeric(dissIndex) )
      
      # Plot data
      p <- ggplot( rawData, aes(x = timePoint, y = dissIndex, colour = samplesFromSpecificArea) ) +
        geom_point(alpha = 0.5) +
        geom_line( data = rawData %>%
                     dplyr::group_by(timePoint, samplesFromSpecificArea) %>%
                     dplyr::summarise( dissIndex = mean(dissIndex, na.rm = T) ), linetype=2, size = 1 ) +
        geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
        scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
        labs( y = "Dissimilarity - Jaccard\n", x = "" ) +
        theme_minimal( ) +
        theme( legend.position = "top" )
      
      # Save plot to PDF
      pdf(paste0(outputFolder, "/a_timeSeries/jaccardIndex_organism_", 
                 organism, "_taxa_", taxonomicLevel, "_compareRegion_Type_", typeToUse, ".pdf"), width = 15, height = 10)
      
      print(p)
      
      dev.off()
      
      # Run GLM
      
      # Build model 1 - by time only
      formulaToUse1 <- formula("dissIndex ~ timePoint")
      model1 <- glm(formula = formulaToUse1, data = rawData, family = quasipoisson(link = "log"))
      
      # Build model 2 - by time and management
      formulaToUse2 <- formula("dissIndex ~ timePoint + samplesFromSpecificArea")
      model2 <- glm(formula = formulaToUse2, data = rawData, family = quasipoisson(link = "log"))
      
      # Calculate ANOVA
      resultsAnova <- anova(model1, model2, test = "F")
      
      # Get numbers that we will save
      statsNumbersToSave <- data.frame(organism = organism,
                                       taxonomicLevel = taxonomicLevel,
                                       samplesFromSpecificArea = "Region",
                                       variableInTest = paste0("JaccardIndex", typeToUse),
                                       pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                       fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                       deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                       df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                       residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                       residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
      
      # Save stats
      if(file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv"))){
        write.table(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv")), append = T, row.names = F, )
      } else {
        write.csv(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), row.names = F )
      }
      
      
      # Analysis - C vs. F
      rawData <- read.csv(paste0(outputFolder, "/a_timeSeries/dissimiralityIndexes.csv"), colClasses = "character") %>%
        dplyr::rename( "organismInAnalysis" = organism,
                       "taxLevelInAnalysis" = taxonomicLevel ) %>%
        dplyr::filter( organismInAnalysis == organism,
                       taxLevelInAnalysis == taxonomicLevel,
                       samplesFromSpecificArea == "Management",
                       type == typeToUse ) %>%
        dplyr::mutate( timePoint = as.numeric(timePoint),
                       dissIndex = as.numeric(dissIndex) )
      
      # Plot data
      p <- ggplot( rawData, aes(x = timePoint, y = dissIndex, colour = management) ) +
        geom_point(alpha = 0.5) +
        geom_line( data = rawData %>%
                     dplyr::group_by(timePoint, management) %>%
                     dplyr::summarise( dissIndex = mean(dissIndex, na.rm = T) ), linetype=2, size = 1 ) +
        geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
        scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
        labs( y = "Dissimilarity - Jaccard\n", x = "" ) +
        theme_minimal( ) +
        theme( legend.position = "top" )
      
      # Save plot to PDF
      pdf(paste0(outputFolder, "/a_timeSeries/jaccardIndex_organism_", 
                 organism, "_taxa_", taxonomicLevel, "_compareManagement_Type_", typeToUse, ".pdf"), width = 15, height = 10)
      
      print(p)
      
      dev.off()
      
      # Run GLM
      
      # Build model 1 - by time only
      formulaToUse1 <- formula("dissIndex ~ timePoint")
      model1 <- glm(formula = formulaToUse1, data = rawData, family = quasipoisson(link = "log"))
      
      # Build model 2 - by time and management
      formulaToUse2 <- formula("dissIndex ~ timePoint + management")
      model2 <- glm(formula = formulaToUse2, data = rawData, family = quasipoisson(link = "log"))
      
      # Calculate ANOVA
      resultsAnova <- anova(model1, model2, test = "F")
      
      # Get numbers that we will save
      statsNumbersToSave <- data.frame(organism = organism,
                                       taxonomicLevel = taxonomicLevel,
                                       samplesFromSpecificArea = "Management",
                                       variableInTest = paste0("JaccardIndex", typeToUse),
                                       pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                       fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                       deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                       df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                       residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                       residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
      
      # Save stats
      if(file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv"))){
        write.table(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/a_timeSeries/anovaResults.csv")), append = T, row.names = F, )
      } else {
        write.csv(statsNumbersToSave, file = paste0(outputFolder, "/a_timeSeries/anovaResults.csv"), row.names = F )
      }
      
    }
    
    
  }
  
}










#
##
### All metrics
##
#

# In this chart we will combine types, numbers and abundances for bacteria, fungi and metazoa.

rm(list = ls())

# Choose normalised method - "bc", "rarefy", "notNormalised"
normaliseMethod <- "rarefy"

# Organisms to be included in the plot: "all", "bacteria", "fungi", "metazoa"
organismsToBuildTimeSeries <- c("bacteria", "fungi", "metazoa")

# Taxonomic levels to be analysed
taxonomicLevelsToBuildTimeSeries <- c("OTUs")

# Regions to include in the chart - "HB", "MB"
regionsToUse <- c("MB")

# Set output folder
outputFolder <- "output"

# Set intermediate folder
intermediateFolder <- "intermediateData"


#
# Collect and merge all the required data
#

################### Numbers

numbersData <- data.frame(stringsAsFactors = F)

# Get data for all organisms
for(organism in organismsToBuildTimeSeries){

  # Get data number of taxa
  load(paste0(intermediateFolder, "/", organism, "/", normaliseMethod, "/dataNumberOfTaxaLevelTaxon_", taxonomicLevelsToBuildTimeSeries, ".RData"))
  
  # Get data 
  # load(paste0(intermediateFolder, "/", organism, "/", normaliseMethod, "/datashannonAndSimpsonLevelTaxon_", taxonomicLevelsToBuildTimeSeries, ".RData"))

  numbersData <- numbersData %>%
    dplyr::bind_rows( dataNumberOfTaxa %>%
                        dplyr::mutate( Organism = organism ) %>%
                        dplyr::filter( Region %in% regionsToUse ) %>%
                        dplyr::select( Organism, value = numberOfTaxa, Management, Timepoint ) %>%
                        dplyr::mutate( measure = "numbers") )
  
}

################### Types and abundances

typesAndAbundanceData <- data.frame(stringsAsFactors = F)

# Get data for all organisms
for(organism in organismsToBuildTimeSeries){
  
  rawData <- read.csv(paste0(outputFolder, "/", organism, "/", normaliseMethod, "/a_timeSeries/dissimiralityIndexes.csv"), colClasses = "character") %>%
    dplyr::rename( "Organism" = organism ) %>%
    dplyr::filter( Organism == organism,
                   taxonomicLevel == "OTUs",
                   management %in% c("C", "F"),
                   samplesFromSpecificArea %in% regionsToUse ) %>%
    dplyr::mutate( timePoint = as.numeric(timePoint),
                   dissIndex = as.numeric(dissIndex),
                   measure = type ) %>%
    dplyr::select( Organism, value = dissIndex, Management = management, Timepoint = timePoint, measure ) %>%
    dplyr::mutate( measure = ifelse(measure == "noBinary", "abundances", "types") )
  
  
  typesAndAbundanceData <- typesAndAbundanceData %>%
    dplyr::bind_rows( rawData )
  
}

### Merge data

allData <- numbersData %>%
  dplyr::bind_rows(typesAndAbundanceData) %>%
  dplyr::mutate( measure = factor(measure, levels = c("numbers", "types", "abundances")),
                 Organism = factor(Organism, levels = c("bacteria", "fungi", "metazoa", "all")))

pAllData <- ggplot( allData, aes(x = Timepoint, y = as.numeric(value), color = Management) ) +
  geom_point(alpha = 0.5) +
  geom_line( data = allData %>%
               dplyr::group_by( Timepoint, Management, measure, Organism ) %>%
               dplyr::summarise( value = mean(value, na.rm = T)) %>%
               dplyr::ungroup(),
             aes(y=as.numeric(value), x=Timepoint, color=Management), linetype=2, size = 1 ) +
  geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
  scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
  labs( y = "Variable in test\n", x = "Timepoint" ) +
  theme_minimal() +
  theme( legend.position = "top" ) + 
  facet_wrap( measure ~ Organism, scales = "free" ) 

# Save plot to PDF
pdf(paste0(outputFolder, "/numbersTypesAndAbundancesForEveryone_", normaliseMethod, "_regionsIncluded_", paste0(regionsToUse, collapse = "_"), ".pdf"), width = 15, height = 10)

print(pAllData)

dev.off()
















