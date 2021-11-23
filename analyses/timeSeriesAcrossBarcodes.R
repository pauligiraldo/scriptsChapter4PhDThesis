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

library(RColorBrewer)


#
# Choose the normalisation method to use: "rarefy", "bc" [bias correction] or "notNormalised".
#
normaliseMethod <- "rarefy"

# Organisms to be analysed:
organismsToBuildTimeSeries <- c("bacteria", "fungi", "metazoa")

# Create output folder
if(!file.exists(paste0("output/", "timeSeriesBarCodes", "/", normaliseMethod))){
  dir.create(paste0("output/", "timeSeriesBarCodes", "/", normaliseMethod), recursive = T)
}
  
outputFolder <- paste0("output/", "timeSeriesBarCodes", "/", normaliseMethod)
  
# Create intermediate folder
if(!file.exists(paste0("intermediateData/", "timeSeriesBarCodes", "/", normaliseMethod))){
  dir.create(paste0("intermediateData/", "timeSeriesBarCodes", "/", normaliseMethod), recursive = T)
}
  
intermediateFolder <- paste0("intermediateData/", "timeSeriesBarCodes", "/", normaliseMethod)
  
taxonomicLevel <- c("OTUs")


# Get data for each organism and merge it
allDataTypesAndAbundances <- data.frame(stringsAsFactors = F)
allDataNumbers <- data.frame(stringsAsFactors = F)

for(organism in organismsToBuildTimeSeries){
  
    # Get data types and abundances
    intermediateDataFolder <- paste0("intermediateData/", organism, "/", normaliseMethod)
    
    load(file = paste0(intermediateDataFolder, "/dataTypesAndAbundancesOfTaxaLevelTaxon_", taxonomicLevel, ".RData"))
    
    # Time points to use
    timePointsToUse <- c(1:9)
    
    for(selectedPoint in timePointsToUse){
      
      filteredTaxDataAreaTimePoint <- subset_samples(taxData, Timepoint == selectedPoint )
      
      # Calculate abundances - Non binary
      distanceDf <- as.data.frame(as.matrix(distance(filteredTaxDataAreaTimePoint, method="jaccard", binary = F)))
      
      dissimilarityDataNonBinary <- data.frame( organism = organism, taxonomicLevel = taxonomicLevel,
                                       timePoint = selectedPoint,
                                       dissIndex = distanceDf[-1,1], type = "noBinary")
      
      
      # Binary
      distanceDf <- as.data.frame(as.matrix(distance(filteredTaxDataAreaTimePoint, method="jaccard", binary = T)))
      
      dissimilarityDataBinary <- data.frame( organism = organism, taxonomicLevel = taxonomicLevel,
                                       timePoint = selectedPoint,
                                       dissIndex = distanceDf[-1,1], type = "binary")
      
      allDataTypesAndAbundances <- allDataTypesAndAbundances %>%
        bind_rows( dissimilarityDataNonBinary %>%
                     dplyr::bind_rows( dissimilarityDataBinary ) )
      
    }
    
    
    
    # Get data numbers
    load(paste0(intermediateDataFolder, "/dataNumberOfTaxaLevelTaxon_", taxonomicLevel, ".RData"))
    # Get data diversity
    load(paste0(intermediateDataFolder, "/datashannonAndSimpsonLevelTaxon_", taxonomicLevel, ".RData"))
    
    # Merge Number of Taxa, Simpson, Shannon
    
    dataToRun <- dataNumberOfTaxa %>%
      dplyr::left_join( shannonAndSimpson %>%
                          dplyr::select( SampleID, Shannon, Simpson, Chao1),
                        by = "SampleID") %>%
      dplyr::mutate( Organism = organism )
    
    allDataNumbers <- allDataNumbers %>%
      dplyr::bind_rows( dataToRun )

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

}




#
##
### Generate ANOVA types and abundances
##
#
growthRates <- data.frame(stringsAsFactors = F)

for(typeToUse in c("binary", "noBinary")){
  
  rawData <- allDataTypesAndAbundances %>%
    dplyr::rename( "taxLevelInAnalysis" = taxonomicLevel ) %>%
    dplyr::filter( taxLevelInAnalysis == taxonomicLevel,
                   type == typeToUse ) %>%
    dplyr::mutate( timePoint = as.numeric(timePoint),
                   dissIndex = as.numeric(dissIndex) ) %>%
    dplyr::mutate( organism = ifelse(organism == "bacteria", "16S",
                                     ifelse(organism == "fungi", "ITS2",
                                            ifelse(organism == "metazoa", "COI", organism) ))) %>%
    dplyr::group_by(timePoint, organism) %>%
    dplyr::summarise( dissIndex = mean(dissIndex, na.rm = T) )
  
  
  # Get growth rates
  for(organismToRate in c("16S", "ITS2", "COI")){
    
    formulaToUse <- formula("dissIndex ~ timePoint")
    
    model <- glm(formula = formulaToUse, data = rawData %>%
                   dplyr::filter(organism == organismToRate), 
                 family = quasipoisson(link = "log"))
    
    summary(model)
    
    growthRatesPre <- data.frame(Organism = organismToRate, DataIndexed = "No", 
                                 Variable = ifelse(typeToUse == "binary", "Types", "Abundances"), 
                                 GrowthRate = as.numeric(model$coefficients[2]))
    
    growthRates <- growthRatesPre %>%
      dplyr::bind_rows(growthRates)
    
  }
  
  
  labelToUse <- ifelse(typeToUse == "binary", "Binary", "Non-Binary")
  
  # Plot data
  p <- ggplot( rawData, aes(x = timePoint, y = dissIndex, colour = organism) ) +
    geom_point(alpha = 0.5) +
    geom_line( data = rawData, linetype=2, size = 1 ) +
    geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
    scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
    labs( y = paste0("Dissimilarity [Jaccard - ", labelToUse, "]\n"), x = "Time points", colour = "" ) +
    scale_color_brewer(palette = "Set2") +
    theme_minimal() +
    theme( legend.position = "top" )
  
  # Save plot to PDF
  pdf(paste0(outputFolder, "/jaccardIndex_All_", 
             "taxa_", taxonomicLevel, "_type_", labelToUse, ".pdf"), width = 12, height = 8)
  
  print(p)
  
  dev.off()
  
  # Run GLM
  
  # Build model 1 - by time only
  formulaToUse1 <- formula("dissIndex ~ timePoint")
  model1 <- glm(formula = formulaToUse1, data = rawData, family = quasipoisson(link = "log"))
  
  # Build model 2 - by time and organism
  formulaToUse2 <- formula("dissIndex ~ timePoint + organism")
  model2 <- glm(formula = formulaToUse2, data = rawData, family = quasipoisson(link = "log"))
  
  # Calculate ANOVA
  resultsAnova <- anova(model1, model2, test = "F")
  
  # Get numbers that we will save
  statsNumbersToSave <- data.frame(taxonomicLevel = taxonomicLevel,
                                   barcode = "All",
                                   variableInTest = paste0("JaccardIndex", typeToUse),
                                   pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                   fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                   deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                   df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                   residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                   residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
  
  # Save stats
  if(file.exists(paste0(outputFolder, "/anovaResults.csv"))){
    write.table(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/anovaResults.csv")), append = T, row.names = F, )
  } else {
    write.csv(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), row.names = F )
  }
  
  
  #
  # Pairwise ANOVA
  #
  pairsToUse <- list( c("bacteria", "fungi"),
                      c("bacteria", "metazoa"),
                      c("fungi", "metazoa") )
  
  for(pairs in 1:length(pairsToUse)) {
  
    rawData <- allDataTypesAndAbundances %>%
      dplyr::rename( "taxLevelInAnalysis" = taxonomicLevel ) %>%
      dplyr::filter( taxLevelInAnalysis == taxonomicLevel,
                     type == typeToUse,
                     organism %in% pairsToUse[[pairs]]) %>%
      dplyr::mutate( timePoint = as.numeric(timePoint),
                     dissIndex = as.numeric(dissIndex) ) %>%
      dplyr::mutate( organism = ifelse(organism == "bacteria", "16S",
                                       ifelse(organism == "fungi", "ITS2",
                                              ifelse(organism == "metazoa", "COI", organism) ))) %>%
      dplyr::group_by(timePoint, organism) %>%
      dplyr::summarise( dissIndex = mean(dissIndex, na.rm = T) )
    
    # Plot data
    p <- ggplot( rawData, aes(x = timePoint, y = dissIndex, colour = organism) ) +
      geom_point(alpha = 0.5) +
      geom_line( data = rawData, linetype=2, size = 1 ) +
      geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
      scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
      labs( y = paste0("Dissimilarity [Jaccard - ", labelToUse, "]\n"), x = "Time points", colour = "" ) +
      scale_color_brewer(palette = "Set2") +
      theme_minimal() +
      theme( legend.position = "top" )
    
    # Save plot to PDF
    pdf(paste0(outputFolder, "/jaccardIndex_", paste0(pairsToUse[[pairs]], collapse = "-"),
               "_taxa_", taxonomicLevel, "_type_", labelToUse, ".pdf"), width = 12, height = 8)
    
    print(p)
    
    dev.off()
    
    # Run GLM
    
    # Build model 1 - by time only
    formulaToUse1 <- formula("dissIndex ~ timePoint")
    model1 <- glm(formula = formulaToUse1, data = rawData, family = quasipoisson(link = "log"))
    
    # Build model 2 - by time and organism
    formulaToUse2 <- formula("dissIndex ~ timePoint + organism")
    model2 <- glm(formula = formulaToUse2, data = rawData, family = quasipoisson(link = "log"))
    
    # Calculate ANOVA
    resultsAnova <- anova(model1, model2, test = "F")
    
    # Get numbers that we will save
    statsNumbersToSave <- data.frame(taxonomicLevel = taxonomicLevel,
                                     barcode = paste0(pairsToUse[[pairs]], collapse = "-"),
                                     variableInTest = paste0("JaccardIndex", typeToUse),
                                     pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                     fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                     deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                     df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                     residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                     residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
    
    # Save stats
    if(file.exists(paste0(outputFolder, "/anovaResults.csv"))){
      write.table(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/anovaResults.csv")), append = T, row.names = F, )
    } else {
      write.csv(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), row.names = F )
    }
  }
}


#
##
### Generate ANOVA numbers
##
#

# Convert data to long format
dataToAnova <- allDataNumbers %>%
  dplyr::select( SampleID, "numberOfTaxa", "Organism", "Timepoint" ) %>%
  dplyr::group_by(Timepoint, Organism) %>%
  dplyr::summarise( numberOfTaxa = mean(numberOfTaxa, na.rm = T) ) %>%
  dplyr::mutate( Organism = ifelse(Organism == "bacteria", "16S",
                                   ifelse(Organism == "fungi", "ITS2",
                                          ifelse(Organism == "metazoa", "COI", Organism) )))

# Get growth rates
for(organismToRate in c("16S", "ITS2", "COI")){
  
  formulaToUse <- formula("numberOfTaxa ~ Timepoint")
  
  model <- glm(formula = formulaToUse, data = dataToAnova %>%
                 dplyr::filter(Organism == organismToRate), 
               family = quasipoisson(link = "log"))
  
  summary(model)
  
  growthRatesPre <- data.frame(Organism = organismToRate, DataIndexed = "No", Variable = "Numbers", GrowthRate = as.numeric(model$coefficients[2]))
  
  growthRates <- growthRatesPre %>%
    dplyr::bind_rows(growthRates)
  
}


# Plot data
p <- ggplot( dataToAnova, aes(x = Timepoint, y = numberOfTaxa, colour = Organism) ) +
  geom_point(alpha = 0.5) +
  geom_line( data = dataToAnova, linetype=2, size = 1 ) +
  geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
  scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
  labs( y = paste0("Average number of taxa\n"), x = "Time points", colour = "" ) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  theme( legend.position = "top" )

# Save plot to PDF
pdf(paste0(outputFolder, "/numerOdfTaxa_All_", 
           "taxa_", taxonomicLevel, ".pdf"), width = 12, height = 8)

print(p)

dev.off()

# Run GLM

# Build model 1 - by time only
formulaToUse1 <- formula("numberOfTaxa ~ Timepoint")
model1 <- glm(formula = formulaToUse1, data = dataToAnova, family = quasipoisson(link = "log"))

# Build model 2 - by time and organism
formulaToUse2 <- formula("numberOfTaxa ~ Timepoint + Organism")
model2 <- glm(formula = formulaToUse2, data = dataToAnova, family = quasipoisson(link = "log"))

# Calculate ANOVA
resultsAnova <- anova(model1, model2, test = "F")

# Get numbers that we will save
statsNumbersToSave <- data.frame(taxonomicLevel = taxonomicLevel,
                                 barcode = "All",
                                 variableInTest = "numberOfTaxa",
                                 pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                 fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                 deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                 df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                 residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                 residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )

# Save stats
if(file.exists(paste0(outputFolder, "/anovaResults.csv"))){
  write.table(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/anovaResults.csv")), append = T, row.names = F, )
} else {
  write.csv(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), row.names = F )
}


#
# Pairwise ANOVA
#
pairsToUse <- list( c("bacteria", "fungi"),
                    c("bacteria", "metazoa"),
                    c("fungi", "metazoa") )

for(pairs in 1:length(pairsToUse)) {
  
  dataToAnova <- allDataNumbers %>%
    dplyr::filter( Organism %in% pairsToUse[[pairs]] ) %>%
    dplyr::select( SampleID, "numberOfTaxa", "Organism", "Timepoint" ) %>%
    dplyr::group_by(Timepoint, Organism) %>%
    dplyr::summarise( numberOfTaxa = mean(numberOfTaxa, na.rm = T) ) %>%
    dplyr::mutate( Organism = ifelse(Organism == "bacteria", "16S",
                                     ifelse(Organism == "fungi", "ITS2",
                                            ifelse(Organism == "metazoa", "COI", Organism) )))
  
  
  # Plot data
  p <- ggplot( dataToAnova, aes(x = Timepoint, y = numberOfTaxa, colour = Organism) ) +
    geom_point(alpha = 0.5) +
    geom_line( data = dataToAnova, linetype=2, size = 1 ) +
    geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
    scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
    labs( y = paste0("Average number of taxa\n"), x = "Time points", colour = "" ) +
    scale_color_brewer(palette = "Set2") +
    theme_minimal() +
    theme( legend.position = "top" )
  
  # Save plot to PDF
  pdf(paste0(outputFolder, "/numerOdfTaxa_", paste0(pairsToUse[[pairs]], collapse = "-"), 
             "_taxa_", taxonomicLevel, ".pdf"), width = 12, height = 8)
  
  print(p)
  
  dev.off()
  
  # Run GLM
  
  # Build model 1 - by time only
  formulaToUse1 <- formula("numberOfTaxa ~ Timepoint")
  model1 <- glm(formula = formulaToUse1, data = dataToAnova, family = quasipoisson(link = "log"))
  
  # Build model 2 - by time and organism
  formulaToUse2 <- formula("numberOfTaxa ~ Timepoint + Organism")
  model2 <- glm(formula = formulaToUse2, data = dataToAnova, family = quasipoisson(link = "log"))
  
  # Calculate ANOVA
  resultsAnova <- anova(model1, model2, test = "F")
  
  # Get numbers that we will save
  statsNumbersToSave <- data.frame(taxonomicLevel = taxonomicLevel,
                                   barcode = paste0(pairsToUse[[pairs]], collapse = "-"),
                                   variableInTest = "numberOfTaxa",
                                   pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                   fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                   deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                   df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                   residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                   residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
  
  # Save stats
  if(file.exists(paste0(outputFolder, "/anovaResults.csv"))){
    write.table(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/anovaResults.csv")), append = T, row.names = F, )
  } else {
    write.csv(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), row.names = F )
  }
}


################################################################################################################
############################################# Indexed values ###################################################
################################################################################################################

#
##
### Generate ANOVA types and abundances
##
#

for(typeToUse in c("binary", "noBinary")){
  
  rawDataPre <- allDataTypesAndAbundances %>%
    dplyr::rename( "taxLevelInAnalysis" = taxonomicLevel ) %>%
    dplyr::filter( taxLevelInAnalysis == taxonomicLevel,
                   type == typeToUse ) %>%
    dplyr::mutate( timePoint = as.numeric(timePoint),
                   dissIndex = as.numeric(dissIndex) ) %>%
    dplyr::mutate( organism = ifelse(organism == "bacteria", "16S",
                                     ifelse(organism == "fungi", "ITS2",
                                            ifelse(organism == "metazoa", "COI", organism) ))) %>%
    dplyr::group_by(timePoint, organism) %>%
    dplyr::summarise( dissIndex = mean(dissIndex, na.rm = T) ) %>%
    dplyr::ungroup()
  
  baseValue <- rawDataPre %>%
    dplyr::filter( timePoint == 1 )
  
  rawData <- rawDataPre %>%
    dplyr::mutate( dissIndex = ifelse( organism == "16S", 100 * (dissIndex/(baseValue %>% filter( organism == "16S"))$dissIndex),
                                       ifelse( organism == "ITS2", 100 * (dissIndex/(baseValue %>% filter( organism == "ITS2"))$dissIndex),
                                               ifelse( organism == "COI", 100 * (dissIndex/(baseValue %>% filter( organism == "COI"))$dissIndex), dissIndex)
                                               )))
  
  # Get growth rates
  for(organismToRate in c("16S", "ITS2", "COI")){
    
    formulaToUse <- formula("dissIndex ~ timePoint")
    
    model <- glm(formula = formulaToUse, data = rawData %>%
                   dplyr::filter(organism == organismToRate), 
                 family = quasipoisson(link = "log"))
    
    summary(model)
    
    growthRatesPre <- data.frame(Organism = organismToRate, DataIndexed = "Yes", 
                                 Variable = ifelse(typeToUse == "binary", "Types", "Abundances"), 
                                 GrowthRate = as.numeric(model$coefficients[2]))
    
    growthRates <- growthRatesPre %>%
      dplyr::bind_rows(growthRates)
    
  }
  
  labelToUse <- ifelse(typeToUse == "binary", "Binary", "Non-Binary")
  
  # Plot data
  p <- ggplot( rawData, aes(x = timePoint, y = dissIndex, colour = organism) ) +
    geom_point(alpha = 0.5) +
    geom_line( data = rawData, linetype=2, size = 1 ) +
    geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
    scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
    labs( y = paste0("Indexed dissimilarity (%) [Jaccard - ", labelToUse, "]\n"), x = "Time points", colour = "" ) +
    scale_color_brewer(palette = "Set2") +
    theme_minimal() +
    theme( legend.position = "top" )
  
  # Save plot to PDF
  pdf(paste0(outputFolder, "/jaccardIndex_All_", 
             "taxa_", taxonomicLevel, "_type_", labelToUse, "_IndexedValues.pdf"), width = 12, height = 8)
  
  print(p)
  
  dev.off()
  
  # Run GLM
  
  # Build model 1 - by time only
  formulaToUse1 <- formula("dissIndex ~ timePoint")
  model1 <- glm(formula = formulaToUse1, data = rawData, family = quasipoisson(link = "log"))
  
  # Build model 2 - by time and organism
  formulaToUse2 <- formula("dissIndex ~ timePoint + organism")
  model2 <- glm(formula = formulaToUse2, data = rawData, family = quasipoisson(link = "log"))
  
  # Calculate ANOVA
  resultsAnova <- anova(model1, model2, test = "F")
  
  # Get numbers that we will save
  statsNumbersToSave <- data.frame(taxonomicLevel = taxonomicLevel,
                                   barcode = "All",
                                   variableInTest = paste0("JaccardIndex", typeToUse, "_Indexed"),
                                   pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                   fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                   deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                   df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                   residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                   residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
  
  # Save stats
  if(file.exists(paste0(outputFolder, "/anovaResults.csv"))){
    write.table(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/anovaResults.csv")), append = T, row.names = F, )
  } else {
    write.csv(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), row.names = F )
  }
  
  
  #
  # Pairwise ANOVA
  #
  pairsToUse <- list( c("bacteria", "fungi"),
                      c("bacteria", "metazoa"),
                      c("fungi", "metazoa") )
  
  for(pairs in 1:length(pairsToUse)) {
    
    rawDataPre <- allDataTypesAndAbundances %>%
      dplyr::rename( "taxLevelInAnalysis" = taxonomicLevel ) %>%
      dplyr::filter( taxLevelInAnalysis == taxonomicLevel,
                     type == typeToUse,
                     organism %in% pairsToUse[[pairs]]) %>%
      dplyr::mutate( timePoint = as.numeric(timePoint),
                     dissIndex = as.numeric(dissIndex) ) %>%
      dplyr::mutate( organism = ifelse(organism == "bacteria", "16S",
                                       ifelse(organism == "fungi", "ITS2",
                                              ifelse(organism == "metazoa", "COI", organism) ))) %>%
      dplyr::group_by(timePoint, organism) %>%
      dplyr::summarise( dissIndex = mean(dissIndex, na.rm = T) )
    
    baseValue <- rawDataPre %>%
      dplyr::filter( timePoint == 1 )
    
    rawData <- rawDataPre %>%
      dplyr::mutate( dissIndex = ifelse( organism == "16S", 100 * (dissIndex/(baseValue %>% filter( organism == "16S"))$dissIndex),
                                         ifelse( organism == "ITS2", 100 * (dissIndex/(baseValue %>% filter( organism == "ITS2"))$dissIndex),
                                                 ifelse( organism == "COI", 100 * (dissIndex/(baseValue %>% filter( organism == "COI"))$dissIndex), dissIndex)
                                         )))
    
    # Plot data
    p <- ggplot( rawData, aes(x = timePoint, y = dissIndex, colour = organism) ) +
      geom_point(alpha = 0.5) +
      geom_line( data = rawData, linetype=2, size = 1 ) +
      geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
      scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
      labs( y = paste0("Indexed dissimilarity (%) [Jaccard - ", labelToUse, "]\n"), x = "Time points", colour = "" ) +
      scale_color_brewer(palette = "Set2") +
      theme_minimal() +
      theme( legend.position = "top" )
    
    # Save plot to PDF
    pdf(paste0(outputFolder, "/jaccardIndex_", paste0(pairsToUse[[pairs]], collapse = "-"),
               "_taxa_", taxonomicLevel, "_type_", labelToUse, "_IndexedValues.pdf"), width = 12, height = 8)
    
    print(p)
    
    dev.off()
    
    # Run GLM
    
    # Build model 1 - by time only
    formulaToUse1 <- formula("dissIndex ~ timePoint")
    model1 <- glm(formula = formulaToUse1, data = rawData, family = quasipoisson(link = "log"))
    
    # Build model 2 - by time and organism
    formulaToUse2 <- formula("dissIndex ~ timePoint + organism")
    model2 <- glm(formula = formulaToUse2, data = rawData, family = quasipoisson(link = "log"))
    
    # Calculate ANOVA
    resultsAnova <- anova(model1, model2, test = "F")
    
    # Get numbers that we will save
    statsNumbersToSave <- data.frame(taxonomicLevel = taxonomicLevel,
                                     barcode = paste0(pairsToUse[[pairs]], collapse = "-"),
                                     variableInTest = paste0("JaccardIndex", typeToUse, "_Indexed"),
                                     pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                     fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                     deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                     df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                     residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                     residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
    
    # Save stats
    if(file.exists(paste0(outputFolder, "/anovaResults.csv"))){
      write.table(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/anovaResults.csv")), append = T, row.names = F, )
    } else {
      write.csv(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), row.names = F )
    }
  }
}


#
##
### Generate ANOVA numbers
##
#

# Convert data to long format
dataToAnovaPre <- allDataNumbers %>%
  dplyr::select( SampleID, "numberOfTaxa", "Organism", "Timepoint" ) %>%
  dplyr::group_by(Timepoint, Organism) %>%
  dplyr::summarise( numberOfTaxa = mean(numberOfTaxa, na.rm = T) ) %>%
  dplyr::mutate( Organism = ifelse(Organism == "bacteria", "16S",
                                   ifelse(Organism == "fungi", "ITS2",
                                          ifelse(Organism == "metazoa", "COI", Organism) )))

baseValue <- dataToAnovaPre %>%
  dplyr::filter( Timepoint == 1 )

dataToAnova <- dataToAnovaPre %>%
  dplyr::mutate( numberOfTaxa = ifelse( Organism == "16S", 100 * (numberOfTaxa/(baseValue %>% filter( Organism == "16S"))$numberOfTaxa),
                                     ifelse( Organism == "ITS2", 100 * (numberOfTaxa/(baseValue %>% filter( Organism == "ITS2"))$numberOfTaxa),
                                             ifelse( Organism == "COI", 100 * (numberOfTaxa/(baseValue %>% filter( Organism == "COI"))$numberOfTaxa), numberOfTaxa)
                                     )))

# Get growth rates
for(organismToRate in c("16S", "ITS2", "COI")){
  
  formulaToUse <- formula("numberOfTaxa ~ Timepoint")
  
  model <- glm(formula = formulaToUse, data = dataToAnova %>%
                  dplyr::filter(Organism == organismToRate), 
                family = quasipoisson(link = "log"))
  
  summary(model)
  
  growthRatesPre <- data.frame(Organism = organismToRate, DataIndexed = "Yes", Variable = "Numbers", GrowthRate = as.numeric(model$coefficients[2]))
  
  growthRates <- growthRatesPre %>%
    dplyr::bind_rows(growthRates)
  
}


# Plot data
p <- ggplot( dataToAnova, aes(x = Timepoint, y = numberOfTaxa, colour = Organism) ) +
  geom_point(alpha = 0.5) +
  geom_line( data = dataToAnova, linetype=2, size = 1 ) +
  geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
  scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
  labs( y = paste0("Indexed average number of taxa (%)\n"), x = "Time points", colour = "" ) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  theme( legend.position = "top" )

# Save plot to PDF
pdf(paste0(outputFolder, "/numerOdfTaxa_All_", 
           "taxa_", taxonomicLevel, "_IndexedValues.pdf"), width = 12, height = 8)

print(p)

dev.off()

# Run GLM

# Build model 1 - by time only
formulaToUse1 <- formula("numberOfTaxa ~ Timepoint")
model1 <- glm(formula = formulaToUse1, data = dataToAnova, family = quasipoisson(link = "log"))

# Build model 2 - by time and organism
formulaToUse2 <- formula("numberOfTaxa ~ Timepoint + Organism")
model2 <- glm(formula = formulaToUse2, data = dataToAnova, family = quasipoisson(link = "log"))

# Calculate ANOVA
resultsAnova <- anova(model1, model2, test = "F")

# Get numbers that we will save
statsNumbersToSave <- data.frame(taxonomicLevel = taxonomicLevel,
                                 barcode = "All",
                                 variableInTest = "numberOfTaxa_Indexed",
                                 pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                 fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                 deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                 df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                 residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                 residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )

# Save stats
if(file.exists(paste0(outputFolder, "/anovaResults.csv"))){
  write.table(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/anovaResults.csv")), append = T, row.names = F, )
} else {
  write.csv(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), row.names = F )
}


#
# Pairwise ANOVA
#
pairsToUse <- list( c("bacteria", "fungi"),
                    c("bacteria", "metazoa"),
                    c("fungi", "metazoa") )

for(pairs in 1:length(pairsToUse)) {
  
  dataToAnovaPre <- allDataNumbers %>%
    dplyr::filter( Organism %in% pairsToUse[[pairs]] ) %>%
    dplyr::select( SampleID, "numberOfTaxa", "Organism", "Timepoint" ) %>%
    dplyr::group_by(Timepoint, Organism) %>%
    dplyr::summarise( numberOfTaxa = mean(numberOfTaxa, na.rm = T) ) %>%
    dplyr::mutate( Organism = ifelse(Organism == "bacteria", "16S",
                                     ifelse(Organism == "fungi", "ITS2",
                                            ifelse(Organism == "metazoa", "COI", Organism) )))
  
  baseValue <- dataToAnovaPre %>%
    dplyr::filter( Timepoint == 1 )
  
  dataToAnova <- dataToAnovaPre %>%
    dplyr::mutate( numberOfTaxa = ifelse( Organism == "16S", 100 * (numberOfTaxa/(baseValue %>% filter( Organism == "16S"))$numberOfTaxa),
                                          ifelse( Organism == "ITS2", 100 * (numberOfTaxa/(baseValue %>% filter( Organism == "ITS2"))$numberOfTaxa),
                                                  ifelse( Organism == "COI", 100 * (numberOfTaxa/(baseValue %>% filter( Organism == "COI"))$numberOfTaxa), numberOfTaxa)
                                          )))
  
  
  # Plot data
  p <- ggplot( dataToAnova, aes(x = Timepoint, y = numberOfTaxa, colour = Organism) ) +
    geom_point(alpha = 0.5) +
    geom_line( data = dataToAnova, linetype=2, size = 1 ) +
    geom_smooth( method = "glm", method.args = list(family = quasipoisson(link = "log")) ) +
    scale_x_continuous( "", breaks = 1:9, labels = c("B15", "V16", "H16", "B16", "V17", "H17", "B17", "V18", "H18")  ) +
    labs( y = paste0("Indexed average number of taxa (%)\n"), x = "Time points", colour = "" ) +
    scale_color_brewer(palette = "Set2") +
    theme_minimal() +
    theme( legend.position = "top" )
  
  # Save plot to PDF
  pdf(paste0(outputFolder, "/numerOdfTaxa_", paste0(pairsToUse[[pairs]], collapse = "-"), 
             "_taxa_", taxonomicLevel, "_IndexedValues.pdf"), width = 12, height = 8)
  
  print(p)
  
  dev.off()
  
  # Run GLM
  
  # Build model 1 - by time only
  formulaToUse1 <- formula("numberOfTaxa ~ Timepoint")
  model1 <- glm(formula = formulaToUse1, data = dataToAnova, family = quasipoisson(link = "log"))
  
  # Build model 2 - by time and organism
  formulaToUse2 <- formula("numberOfTaxa ~ Timepoint + Organism")
  model2 <- glm(formula = formulaToUse2, data = dataToAnova, family = quasipoisson(link = "log"))
  
  # Calculate ANOVA
  resultsAnova <- anova(model1, model2, test = "F")
  
  # Get numbers that we will save
  statsNumbersToSave <- data.frame(taxonomicLevel = taxonomicLevel,
                                   barcode = paste0(pairsToUse[[pairs]], collapse = "-"),
                                   variableInTest = "numberOfTaxa_Indexed",
                                   pValue = resultsAnova$`Pr(>F)`[!is.na(resultsAnova$`Pr(>F)`)],
                                   fValue = resultsAnova$F[!is.na(resultsAnova$F)],
                                   deviance = resultsAnova$`Deviance`[!is.na(resultsAnova$`Deviance`)],
                                   df = resultsAnova$`Df`[!is.na(resultsAnova$`Df`)],
                                   residualDeviation = paste(resultsAnova$`Resid. Dev`[!is.na(resultsAnova$`Resid. Dev`)], collapse = " | "),
                                   residualDf = paste(resultsAnova$`Resid. Df`[!is.na(resultsAnova$`Resid. Df`)], collapse = " | ") )
  
  # Save stats
  if(file.exists(paste0(outputFolder, "/anovaResults.csv"))){
    write.table(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/anovaResults.csv")), append = T, row.names = F, )
  } else {
    write.csv(statsNumbersToSave, file = paste0(outputFolder, "/anovaResults.csv"), row.names = F )
  }
}


# Save stats
write.csv(growthRates, file = paste0(outputFolder, "/growthRates.csv"), row.names = F )


    