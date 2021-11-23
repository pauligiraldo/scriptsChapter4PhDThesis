# Purpose: run Kruskal Wallis testing for all potential interactions on sub-levels
# Date: 07/09/2019
# Type: it is part of the script "analyses/fourYears/analysis_four_years_pipeline.R"
# Author: Paulina Giraldo Perez





#
# Kruskal-Wallis with interaction testing per taxonomic level on sub-levels of "Region", "Season" and "Year_of_program"
#

# The taxonomic levels that will be tested are defined through the object taxonomicLevels,
# which is defined in the script analyses/fourYears/analysisFourYearsPipeline.R

# The factors that will be tested are defined through the object factorsToTest,
# which is defined in the script analyses/fourYears/analysisFourYearsPipeline.R

forSubFactorToTest <- factorsToTest[factorsToTest != "Management"]

# kResults is the object where we will store the results of the loop through taxonomic levels and interactions
kResultsSubLevel <- data.frame()
# sampleStatsSubLevel is the object where we will store the stats (e.g. mean, median) of each comparison
sampleStatsSubLevel <- data.frame()
# dataForBoxPlot is the object where we will store the values to build a boxPlot
dataForBoxPlot <- data.frame()

# For each factor, isolate data by the sub-levels of each factor
for(i in 1:length(forSubFactorToTest)){
  
  message(paste("Testing factor:", forSubFactorToTest[i]))
  
  # Find sub-levels for this factor
  subLevels <- unique(sample_data(biomDataRarefied)[[forSubFactorToTest[i]]])
  
  # Test each sub-level
  for(j in 1:length(subLevels)){
    
    message(paste("Testing sublevel:", subLevels[j]))
    
    var_values <- sample_data(biomDataRarefied)[[forSubFactorToTest[i]]]
    
    subLevelData <- prune_samples( var_values == subLevels[j], biomDataRarefied)
  
    # Create all the possible interactions for all the sub-factors defined in factorsToTest
    interactionsToTest <- list()
    
    forInterectionTest <- factorsToTest[factorsToTest != forSubFactorToTest[i]]
    
    for(k in 1:length(forInterectionTest)){
      
      interactionsToTest <- append( interactionsToTest, combn( forInterectionTest, m = k, simplify = F ) )
      
    }
    
    # Start loop by taxonomic levels and interactions
    for(k in 1:length(taxonomicLevels)) {
      
      message(paste("Testing taxonomic level:", taxonomicLevels[k]))
      
      # Check if the analysis should be done by OTUs only or by a taxonomic level
      if(taxonomicLevels[k] == "OTUs"){
        # No need to filter data if analysis for OTUs only
        taxData <- subLevelData
      } else {
        # Filter the data for the taxonomic level we need
        taxData <- tax_glom(subLevelData, taxrank = taxonomicLevels[k], NArm=F)
      }
  
      #
      # Get the sample data and add to it a column counting the number of taxa by sample for the specific taxonomic level we are testing
      #

      dataNumberOfTaxa <- 
        # Get and groom OTU table
        otu_table(taxData) %>%
        data.frame(check.names = F) %>%
        rownames_to_column( var = "taxCode" ) %>%
        tidyr::gather( SampleID, isPresent, -taxCode ) %>%
        # Change data to binary or presence and absence
        dplyr::mutate( isPresent = ifelse(isPresent > 0, 1, isPresent) ) %>%
        # Group by sample and count number of taxa
        dplyr::group_by( SampleID ) %>%
        dplyr::summarise( numberOfTaxa = sum(isPresent, na.rm = T) ) %>%
        data.frame() %>%
        # Add the columns with factor values
        dplyr::left_join( sample_data(taxData) %>%
                          data.frame(check.names = F),
                        by = "SampleID" )
      
      # Get data for boxPlot
      dataForBoxPlot <- dataForBoxPlot %>%
        dplyr::bind_rows( dataNumberOfTaxa %>% 
                            dplyr::select(SampleID, numberOfTaxa, Year_of_program, Season, Management, Region) %>%
                            dplyr::mutate(taxLevel = taxonomicLevels[k],
                                          factorInTest = as.character(forSubFactorToTest[i]),
                                          subFactor = as.character(subLevels[j]) ) )
  
      # Start loop through possible interactions
      for(m in 1:length(interactionsToTest)) {
    
        message(paste("Testing interaction:", paste(interactionsToTest[[m]], collapse = " | ")))
    
        if(length(interactionsToTest[[m]]) == 1){
      
          kTest <- kruskal.test(dataNumberOfTaxa$numberOfTaxa ~ dataNumberOfTaxa[,interactionsToTest[[m]][1]])
          epSquared <- epsilonSquared(x = dataNumberOfTaxa$numberOfTaxa, g = dataNumberOfTaxa[,interactionsToTest[[m]][1]])
          dunnTestResults <- dunnTest(dataNumberOfTaxa$numberOfTaxa, dataNumberOfTaxa[,interactionsToTest[[m]][1]], method = "bonferroni")
      
          kResultsPre <- data.frame( taxLevel = taxonomicLevels[k],
                                     factor = forSubFactorToTest[i],
                                     subfactor = subLevels[j],
                                     interaction = paste0(interactionsToTest[[m]], collapse = " | "),
                                     comparison = dunnTestResults$res$Comparison,
                                     pvalue = kTest$p.value,
                                     pvalueAdjusted = dunnTestResults$res$P.adj,
                                     chiSquared = as.numeric(kTest$statistic),
                                     epsilonSquared = as.numeric(epSquared) )
          
          # Prepare stats for sub-levels of each factor
          statsSubLevels <- dataNumberOfTaxa %>%
            dplyr::group_by_at( vars(one_of(interactionsToTest[[m]])) ) %>%
            dplyr::summarise( minNumberOfTaxa = min(numberOfTaxa, na.rm = T),
                              meanNumberOfTaxa = mean(numberOfTaxa, na.rm = T),
                              medianNumberOfTaxa = median(numberOfTaxa, na.rm = T),
                              maxNumberOfTaxa = max(numberOfTaxa, na.rm = T),
                              sdNumberOfTaxa = stats::sd(numberOfTaxa, na.rm = T),
                              nSamples = n() ) %>%
            dplyr::rename( subFactor = interactionsToTest[[m]] ) %>%
            dplyr::mutate( factor = interactionsToTest[[m]],
                           taxLevel = taxonomicLevels[k] ) %>%
            dplyr::select( factor, subFactor, taxLevel, minNumberOfTaxa, meanNumberOfTaxa, 
                           medianNumberOfTaxa, maxNumberOfTaxa, sdNumberOfTaxa, nSamples)
          
          sampleStatsSubLevel <- sampleStatsSubLevel %>%
            dplyr::bind_rows( statsSubLevels )
      
        } 
    
        if(length(interactionsToTest[[m]]) > 1){
      
          interactionObject <- interaction(dataNumberOfTaxa[,interactionsToTest[[m]]])

          kTest <- kruskal.test(dataNumberOfTaxa$numberOfTaxa ~ interactionObject)
      
          kResultsPre <- data.frame( taxLevel = taxonomicLevels[k],
                                     factor = forSubFactorToTest[i],
                                     subfactor = subLevels[j],
                                     interaction = paste0(interactionsToTest[[m]], collapse = " | "),
                                     comparison = "-",
                                     pvalue = kTest$p.value,
                                     pvalueAdjusted = NA,
                                     chiSquared = as.numeric(kTest$statistic),
                                     epsilonSquared = NA )
      
        } 
    
        kResultsSubLevel <- kResultsSubLevel %>%
          dplyr::bind_rows( kResultsPre )
    
      }
  
    }
    
  }
  
}

# Save results to .RData
save(kResultsSubLevel, file = "intermediateData/fourYears/kResultsSubLevel.RData")
save(sampleStatsSubLevel, file = "intermediateData/fourYears/sampleStatsSubLevel.RData")
save(dataForBoxPlot, file = "intermediateData/fourYears/dataForBoxPlotSubLevel.RData")

# Save kResults as csv
write.csv(kResultsSubLevel, file = "output/fourYears/analysis_four_years_sub_level_kruskal_wallis.csv", row.names = F)
message("A CSV file with the results of Kruskal Wallis on sublevels was saved at \noutput/fourYears/analysis_four_years_sub_level_kruskal_wallis.csv")

# Prepare data for plotting the results of Kruskal-Wallis for the individual factors with no interaction
testIndividualFactors <- kResultsSubLevel %>%
  # The pvalueAdjusted is only collected for tests for each individual factor with no interaction 
  dplyr::filter( !is.na(pvalueAdjusted) ) %>%
  # Re-order taxonomic levels to "OTUs", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus" and "Species"
  dplyr::mutate( taxLevel = factor(taxLevel, levels = c("OTUs", rank_names(biomData)[rank_names(biomData) != "Rank1"]))) %>%
  # Mark the tests that returned adjusted pvalues lower than 0.05. It is basically the asterix we will place at the top of the bars. 
  dplyr::mutate( sigLabel = ifelse(pvalueAdjusted < 0.05, "*", "") )

# Plot results of Kruskal Wallis for individual factors  
pNoInteraction <- ggplot(testIndividualFactors,
                         # Set taxonomic levels on X axis; pvalues adjusted at Y axis; and fill bars by taxonomic levels.
                         aes(x = taxLevel, y = pvalueAdjusted, fill = taxLevel) ) +
  geom_bar( stat = "identity" ) +
  # Add asterix at the top of bars with adjusted pvalue lower than 0.05
  geom_text( aes(label = sigLabel), size = 10, color = "red", vjust = 0 ) +
  # Set the colors to fill bars
  scale_fill_manual( values = brewer.pal(8, "Dark2") ) +
  # Draw a horizontal line at 0.05
  geom_hline( yintercept = 0.05, linetype = 2 ) +
  # Plot facets for interaction and comparison within each level
  facet_wrap( ~factor + subfactor + interaction + comparison, ncol = 2, scales = "free" ) +
  # Define labels of X and Y axes
  labs( x = "", y = "Adjusted pvalue - Bonferroni\n" ) +
  # Remove background color and stripes
  theme( legend.position = "none",
         panel.background = element_blank(),
         legend.key=element_blank(),
         strip.background =element_rect(fill="white") )

# Save graph as PDF
while(length(dev.list()) != 0){ dev.off() } 
pdf(file = "output/fourYears/analysis_four_years_numbers_kruskal_sub_level_individual_factors.pdf", width = 10, height = 80)
print(pNoInteraction)
dev.off()
message("A PDF file with the plot of results from Kruskal Wallis on individual factors was saved at \noutput/fourYears/analysis_four_years_numbers_kruskal_sub_level_individual_factors.pdf")


# Prepare data for plotting the results of Kruskal-Wallis on potential interactions
testInteractions <- kResultsSubLevel %>%
  # The pvalueAdjusted is only collected for tests for each individual factor with no interaction 
  dplyr::filter( is.na(pvalueAdjusted) ) %>%
  # Re-order taxonomic levels to "OTUs", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus" and "Species"
  dplyr::mutate( taxLevel = factor(taxLevel, levels = c("OTUs", rank_names(biomData)[rank_names(biomData) != "Rank1"]))) %>%
  # Mark the tests that returned adjusted pvalues lower than 0.05. It is basically the asterix we will place at the top of the bars. 
  dplyr::mutate( sigLabel = ifelse(pvalue < 0.05, "*", "") )

pInteraction <- ggplot( testInteractions,
                        # Set taxonomic levels on X axis; pvalues at Y axis; and fill bars by taxonomic levels.
                        aes(x = taxLevel, y = pvalue, fill = taxLevel) ) +
  geom_bar( stat = "identity" ) +
  # Add asterix at the top of bars with adjusted pvalue lower than 0.05
  geom_text( aes(label = sigLabel), size = 10, color = "red", vjust = 0 ) +
  # Set the colors to fill bars
  scale_fill_manual( values = brewer.pal(8, "Dark2") ) +
  # Draw a horizontal line at 0.05
  geom_hline( yintercept = 0.05, linetype = 2 ) +
  # Plot facets for interactions between factors tested
  facet_wrap( ~factor + subfactor + interaction, ncol = 2, scales = "free"  ) +
  # Define labels of X and Y axes
  labs( x = "", y = "pvalue\n" ) +
  # Remove background color and stripes
  theme( legend.position = "none",
         panel.background = element_blank(),
         legend.key=element_blank(),
         strip.background =element_rect(fill="white") )

# Save graph as PDF
while(length(dev.list()) != 0){ dev.off() } 
pdf(file = "output/fourYears/analysis_four_years_numbers_kruskal_sub_level_interactions.pdf", width = 10, height = 35)
print(pInteraction)
dev.off()
message("A PDF file with the plot of results from Kruskal Wallis on interactions was saved at \noutput/fourYears/analysis_four_years_numbers_kruskal_sub_level_interactions.pdf")

# Add standard error to sampleStats
sampleStatsSubLevel <- sampleStatsSubLevel %>%
  dplyr::mutate( se = sdNumberOfTaxa / sqrt(nSamples) )

# Save sampleStats as csv
write.csv(sampleStatsSubLevel, file = "output/fourYears/analysis_four_years_sub_level_stats.csv", row.names = F)
message("\n\n\nA CSV file with the stats (e.g. mean, median and min) of each factor was saved at \noutput/fourYears/analysis_four_years_sub_level_stats.csv")

# Build graph for stats (e.g. mean, median, sd, etc...)
pBoxPlot <- ggplot( dataForBoxPlot %>%
                      tidyr::gather( factor, value, -SampleID, -numberOfTaxa, -taxLevel, -subFactor, -factorInTest ) %>%
                      dplyr::mutate(taxLevel = factor(taxLevel, levels = taxonomicLevels),
                                    value = factor(value, 
                                                       levels = rev(c("C", "F", "HB", "MB", "Budburst", "Veraison", "Harvest", "1", "2", "3", "4"))) ), 
                    aes(x = value, y = numberOfTaxa, color = factor) ) +
  geom_boxplot() +
  coord_flip() +
  #scale_color_manual( values = brewer.pal(11, "Dark2") ) +
  facet_wrap( ~subFactor + taxLevel, scales = "free", ncol = 2 ) +
  # Define labels of X and Y axes
  labs( x = "", y = "\nNumber of taxa" ) +
  # Remove background color and stripes
  theme( legend.position = "none",
         panel.background = element_blank(),
         legend.key=element_blank(),
         strip.background =element_rect(fill="white") )

# Save graph as PDF
while(length(dev.list()) != 0){ dev.off() } 
pdf(file = "output/fourYears/analysis_four_years_numbers_sub_level_stats.pdf", width = 10, height = 35)
print(pBoxPlot)
dev.off()
message("\n\n\nA PDF file with the plot with standard stats by taxonomic level was saved at \noutput/fourYears/analysis_four_years_numbers_sub_level_stats.pdf")

# Add kResults as an object we would like to keep
toKeep <- c(toKeep, "kResultsSubLevel", "sampleStatsSubLevel")

# Remove all objects we don't need anymore
rm( list = ls()[!(ls() %in% toKeep)])






