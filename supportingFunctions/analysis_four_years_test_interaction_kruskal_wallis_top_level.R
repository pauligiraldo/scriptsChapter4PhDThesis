# Purpose: run Kruskal Wallis testing for all potential interactions
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

# Create all the possible interactions for all the factors defined in factorsToTest
interactionsToTest <- list()

for(i in 1:length(factorsToTest)){
  
  interactionsToTest <- append( interactionsToTest, combn( factorsToTest, m = i, simplify = F ) )
  
}

# kResults is the object where we will store the results of the loop through taxonomic levels and interactions
kResults <- data.frame()
# sampleStats is the object where we will store the stats (e.g. mean, median) of each comparison
sampleStats <- data.frame()
# dataForBoxPlot is the object where we will store the values to build a boxPlot
dataForBoxPlot <- data.frame()

# Start loop by taxonomic levels and interactions
for(i in 1:length(taxonomicLevels)) {
  
  message(paste("Testing taxonomic level:", taxonomicLevels[i]))
  
  # Check if the analysis should be done by OTUs only or by a taxonomic level
  if(taxonomicLevels[i] == "OTUs"){
    # No need to filter data if analysis for OTUs only
    taxData <- biomDataRarefied
  } else {
    # Filter the data for the taxonomic level we need
    taxData <- tax_glom(biomDataRarefied, taxrank = taxonomicLevels[i], NArm=F)
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
                        dplyr::mutate(taxLevel = taxonomicLevels[i]) )
  
  # Start loop through possible interactions
  for(j in 1:length(interactionsToTest)) {
    
    message(paste("Testing interaction:", paste(interactionsToTest[[j]], collapse = " | ")))
    
    if(length(interactionsToTest[[j]]) == 1){
      
      kTest <- kruskal.test(dataNumberOfTaxa$numberOfTaxa ~ dataNumberOfTaxa[,interactionsToTest[[j]][1]])
      epSquared <- epsilonSquared(x = dataNumberOfTaxa$numberOfTaxa, g = dataNumberOfTaxa[,interactionsToTest[[j]][1]])
      dunnTestResults <- dunnTest(dataNumberOfTaxa$numberOfTaxa, dataNumberOfTaxa[,interactionsToTest[[j]][1]], method = "bonferroni")
      
      kResultsPre <- data.frame( taxLevel = taxonomicLevels[i], 
                                 interaction = paste0(interactionsToTest[[j]], collapse = " | "),
                                 comparison = dunnTestResults$res$Comparison,
                                 pvalue = kTest$p.value,
                                 pvalueAdjusted = dunnTestResults$res$P.adj,
                                 chiSquared = as.numeric(kTest$statistic),
                                 epsilonSquared = as.numeric(epSquared) ) 
      
      # Prepare stats for sub-levels of each factor
      statsSubLevels <- dataNumberOfTaxa %>%
        dplyr::group_by_at( vars(one_of(interactionsToTest[[j]])) ) %>%
        dplyr::summarise( minNumberOfTaxa = min(numberOfTaxa, na.rm = T),
                          meanNumberOfTaxa = mean(numberOfTaxa, na.rm = T),
                          medianNumberOfTaxa = median(numberOfTaxa, na.rm = T),
                          maxNumberOfTaxa = max(numberOfTaxa, na.rm = T),
                          sdNumberOfTaxa = stats::sd(numberOfTaxa, na.rm = T),
                          nSamples = n() ) %>%
        dplyr::rename( subFactor = interactionsToTest[[j]] ) %>%
        dplyr::mutate( factor = interactionsToTest[[j]],
                       taxLevel = taxonomicLevels[i] ) %>%
        dplyr::select( factor, subFactor, taxLevel, minNumberOfTaxa, meanNumberOfTaxa, 
                       medianNumberOfTaxa, maxNumberOfTaxa, sdNumberOfTaxa, nSamples)
      
      sampleStats <- sampleStats %>%
        dplyr::bind_rows( statsSubLevels )
      
    } 
    
    if(length(interactionsToTest[[j]]) > 1){
      
      interactionObject <- interaction(dataNumberOfTaxa[,interactionsToTest[[j]]])

      kTest <- kruskal.test(dataNumberOfTaxa$numberOfTaxa ~ interactionObject)
      
      kResultsPre <- data.frame( taxLevel = taxonomicLevels[i], 
                                 interaction = paste0(interactionsToTest[[j]], collapse = " | "),
                                 comparison = "-",
                                 pvalue = kTest$p.value,
                                 pvalueAdjusted = NA,
                                 chiSquared = as.numeric(kTest$statistic),
                                 epsilonSquared = NA )
      
    } 
    
    kResults <- kResults %>%
      dplyr::bind_rows( kResultsPre )
    
  }
  
}

# Save results to .RData
save(kResults, file = "intermediateData/fourYears/kResults.RData")
save(sampleStats, file = "intermediateData/fourYears/sampleStats.RData")
save(dataForBoxPlot, file = "intermediateData/fourYears/dataForBoxPlot.RData")

# Save kResults as csv
write.csv(kResults, file = "output/fourYears/analysis_four_years_top_level_kruskal_wallis.csv", row.names = F)
message("\n\n\nA CSV file with the results of Kruskal Wallis was saved at \noutput/fourYears/analysis_four_years_top_level_kruskal_wallis.csv")

# Prepare data for plotting the results of Kruskal-Wallis for the individual factors with no interaction
testIndividualFactors <- kResults %>%
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
  facet_wrap( ~interaction + comparison, ncol = 2, scales = "free" ) +
  # Define labels of X and Y axes
  labs( x = "", y = "Adjusted pvalue - Bonferroni\n" ) +
  # Remove background color and stripes
  theme( legend.position = "none",
         panel.background = element_blank(),
         legend.key=element_blank(),
         strip.background =element_rect(fill="white") )

# Save graph as PDF
while(length(dev.list()) != 0){ dev.off() } 
pdf(file = "output/fourYears/analysis_four_years_numbers_top_level_kruskal_individual_factors.pdf", width = 10, height = 25)
print(pNoInteraction)
dev.off()
message("\n\n\nA PDF file with the plot of results from Kruskal Wallis on individual factors was saved at \noutput/fourYears/analysis_four_years_numbers_top_level_kruskal_individual_factors.pdf")


# Prepare data for plotting the results of Kruskal-Wallis on potential interactions
testInteractions <- kResults %>%
  # The pvalueAdjusted is only collected for tests for each individual factor with no interaction 
  dplyr::filter( is.na(pvalueAdjusted) ) %>%
  # Re-order taxonomic levels to "OTUs", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus" and "Species"
  dplyr::mutate( taxLevel = factor(taxLevel, levels = c("OTUs", rank_names(biomData)[rank_names(biomData) != "Rank1"]))) %>%
  # Mark the tests that returned adjusted pvalues lower than 0.05. It is basically the asterix we will place at the top of the bars. 
  dplyr::mutate( sigLabel = ifelse(pvalueAdjusted < 0.05, "*", "") )

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
  facet_wrap( ~interaction, ncol = 2, scales = "free"  ) +
  # Define labels of X and Y axes
  labs( x = "", y = "pvalue\n" ) +
  # Remove background color and stripes
  theme( legend.position = "none",
         panel.background = element_blank(),
         legend.key=element_blank(),
         strip.background =element_rect(fill="white") )

# Save graph as PDF
while(length(dev.list()) != 0){ dev.off() }
pdf(file = "output/fourYears/analysis_four_years_numbers_top_level_kruskal_interactions.pdf", width = 10, height = 25)
print(pInteraction)
dev.off()
message("\n\n\nA PDF file with the plot of results from Kruskal Wallis on interactions was saved at \noutput/fourYears/analysis_four_years_numbers_top_level_kruskal_interactions.pdf")

# Add standard error to sampleStats
sampleStats <- sampleStats %>%
  dplyr::mutate( se = sdNumberOfTaxa / sqrt(nSamples) )

# Save sampleStats as csv
write.csv(sampleStats, file = "output/fourYears/analysis_four_years_numbers_top_level_stats.csv", row.names = F)
message("\n\n\nA CSV file with the stats (e.g. mean, median and min) of each factor was saved at \noutput/fourYears/analysis_four_years_numbers_top_level_stats.csv")

# Build graph for stats (e.g. mean, median, sd, etc...)
pBoxPlot <- ggplot( dataForBoxPlot %>%
          tidyr::gather( factor, subFactor, -SampleID, -numberOfTaxa, -taxLevel ) %>%
          dplyr::mutate(taxLevel = factor(taxLevel, levels = taxonomicLevels),
                        subFactor = factor(subFactor, 
                                           levels = rev(c("C", "F", "HB", "MB", "Budburst", "Veraison", "Harvest", "1", "2", "3", "4"))) ), 
        aes(x = subFactor, y = numberOfTaxa, color = subFactor) ) +
  geom_boxplot() +
  coord_flip() +
  #scale_color_manual( values = brewer.pal(11, "Dark2") ) +
  facet_wrap( ~taxLevel, scales = "free", ncol = 2 ) +
  # Define labels of X and Y axes
  labs( x = "", y = "\nNumber of taxa" ) +
  # Remove background color and stripes
  theme( legend.position = "none",
         panel.background = element_blank(),
         legend.key=element_blank(),
         strip.background =element_rect(fill="white") )

# Save graph as PDF
while(length(dev.list()) != 0){ dev.off() }
pdf(file = "output/fourYears/analysis_four_years_numbers_top_level_stats.pdf", width = 10, height = 15)
print(pBoxPlot)
dev.off()
message("\n\n\nA PDF file with the plot with standard stats by taxonomic level was saved at \noutput/fourYears/analysis_four_years_numbers_top_level_stats.pdf")

# Add kResults as an object we would like to keep
toKeep <- c(toKeep, "kResults", "sampleStats")

# Remove all objects we don't need anymore
rm( list = ls()[!(ls() %in% toKeep)])






