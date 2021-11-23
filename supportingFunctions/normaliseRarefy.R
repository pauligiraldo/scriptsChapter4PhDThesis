message("Normalising data using rarefy")

#
# Find best sample size a.k.a depth for rarefying
#

# The rarefying process is performed by normalising or indexing the number of reads per sample.
# During rarefying, we can define the minimun number of reads per sample that is required to maintain
# a sample in the dataset. Samples showing a number of reads lower than the number we define 
# will be removed from the analysis.

if(shouldGeneratePlotNReads) {
  # The folowing code allows you to visualise and list sequencing depths of each sample.
  rarecurve(t(otu_table(biomData)), step=50, cex = 0.5) 
  (depthsList <- sample_sums(biomData)[order(as.numeric(sample_sums(biomData)))])
  
  pdf(file = paste0(outputFolder, "/analysis_rarefaction_curve.pdf"), width = 15, height = 10)
  rarecurve(t(otu_table(biomData)), step=50, cex = 0.5)
  dev.off()
}


# The following script will: 
# * test the number of reads that offers the most stable proportion of samples accross classes after rarefying.
# * the best value of number of reads to use in rarefying is defined as:
# * The number of reads that, when used for rarefying, will show a the same proportion of samples across classes.
# * For example, suppose your experiment is testing for Management type and Region. These are two classes.
# * The class Management has 2 values, "Conservation" and "Future".
# * The class Rgeion has 3 values, "Hawk's Bay", "Marlborough" and "Nelson".
# * The ideal resulting dataset after rarefying would contain 50% of the sample from "Conservation" and 50% from "Future"; 
# * and 33.3% of samples from "Hawk's Bay", 33% from "Marlborough" and 33% from "Nelson".
# * We don't want to end up with classes containing too few samples. To make sure it doesn't happen, we can set minimun proportion
# * of sample a class should have.

if(is.na(readsValueForRarefying)){
  
  minimumProp <- 0.25
  
  # * It will show the depth value selected and the resulting proportion of samples per class. 
  
  source("supportingFunctions/analysis_four_years_test_sample_size_values.R", local = environment())
  
  # Outputs in intermediateData folder
  
  print(numberOfReadsToUse)
  
  # Copy rarefying log to output
  file.copy(file.path(intermediateFolderRarefying, paste0("LogForRarefyingValue_", numberOfReadsToUse, ".csv") ), outputFolder)
  
} else {
  
  numberOfReadsToUse <- readsValueForRarefying
  
}
#
# Apply rarefying using the sampleSizeValue defined above
#

# Set the seed for ramdom number
set.seed(4002)

# Run rarefying function !! Using the numberOfReadsToUse object we defined and tested above !!
biomDataNormalised <- rarefy_even_depth(biomData, rngseed = F, sample.size = numberOfReadsToUse,
                                        replace = F, trimOTUs = T)

# Save biomDataRarefied as a .RData
save(biomData, file = paste0(intermediateFolder, "/biomData.RData"))
save(biomDataNormalised, file = paste0(intermediateFolder, "/biomDataNormalised.RData"))

# Generate data profile after rarefying
#biomDataRarefied

dataProfileAfter <- sample_data(biomDataNormalised) %>%
  data.frame() %>%
  dplyr::select( SampleID, Year_of_program, Season, Management, Region) %>%
  tidyr::gather( Variable, Value, -SampleID  ) %>%
  dplyr::group_by( Variable, Value ) %>%
  dplyr::summarise( Samples = n(),
                    organism = organism,
                    normaliseMethod = normaliseMethod,
                    time = "after")

dataSizeLog <- dataSizeLog %>%
  dplyr::bind_rows( dataProfile %>%
                      # Remove rarefying column
                      dplyr::bind_rows( dataProfileAfter ) %>%
                      tidyr::spread( time, Samples ) )

# Save data profile
if(file.exists("output/dataProfileBeforAndAfterRarefying.csv")){
  write.table(dataSizeLog, file = paste0(outputFolder, "/dataProfileBeforAndAfterNormalising.csv"), sep = ",", col.names = !file.exists(paste0(outputFolder, "/dataProfileBeforAndAfterNormalising.csv")), append = T, row.names = F, )
} else {
  write.csv(dataSizeLog, file = paste0(outputFolder, "/dataProfileBeforAndAfterNormalising.csv"), row.names = F )
}