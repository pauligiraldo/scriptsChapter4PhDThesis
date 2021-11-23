# Purpose: test the best sample.size value to be used in rarefying
# Date: 04/09/2019
# Type: it is part of the script "analyses/fourYears/analysis_four_years_pipeline.R"
# Author: Paulina Giraldo Perez

#
# Find best sample size for rarefying
#

if(!file.exists(file.path(intermediateFolder, "rarefyingLogs"))){
  dir.create(file.path(intermediateFolder, "rarefyingLogs"))
}

intermediateFolderRarefying <- file.path(intermediateFolder, "rarefyingLogs")

toKeep <- c(toKeep, "intermediateFolderRarefying")

depthRange <- seq(from = depthRangeFrom, to = depthRangeTo, by = depthRangeBy)

# The value of depthRange is defined at the top of the script integrate.R

disimilarityIndexResults <- data.frame(stringsAsFactors = F)

for( depthRangeValue in depthRange ) {
  
  message(paste0("Testing rarefying when using the number of reads of: ", depthRangeValue))
  
  # List samoples with their sum of reads for each sample
  reads <- sample_sums(biomData)[order(as.numeric(sample_sums(biomData)))] %>%
    data.frame() %>%
    dplyr::rename( "numberOfReads" = "." ) %>%
    # Get sample names
    tibble::rownames_to_column( var = "SampleID" )
  
  # Filter data according to the depthRangeValue
  samplesThatWillBeRemoved <- reads %>%
    # Join to metadata so we know the Region, Season etc.. associated with each sample
    dplyr::left_join( data.frame(sample_data(biomData)), by = "SampleID" ) %>%
    # Select the columns we need to inspect
    dplyr::select( numberOfReads, SampleID, Management, Region, Season, Year_of_program ) %>%
    # Define the number of reads that would be used as cutoff
    dplyr::filter( numberOfReads < depthRangeValue )
  
  # Get proportion of remaining samples per class
  propSamplesThatWillBeRemoved <- samplesThatWillBeRemoved %>%
    tidyr::gather( Variable, Value, -numberOfReads, -SampleID ) %>%
    dplyr::group_by( Variable, Value ) %>%
    dplyr::summarise( SamplesRemoved = n() ) %>%
    dplyr::mutate( FreqOfRemovedSamples = SamplesRemoved / sum(SamplesRemoved),
                   TotalNumberOfSamplesRemoved = sum(SamplesRemoved) )
  
  # Filter data according to the sampleSizeValue defined above
  samplesThatWillRemain <- reads %>%
    # Join to metadata so we know the Region, Season etc.. associated with each sample
    dplyr::left_join( data.frame(sample_data(biomData)), by = "SampleID" ) %>%
    # Select the columns we need to inspect
    dplyr::select( numberOfReads, SampleID, Management, Region, Season, Year_of_program ) %>%
    # Define the number of reads that would be used as cutoff
    dplyr::filter( numberOfReads >= depthRangeValue )
  
  # Get proportion of remaining samples per class
  propSamplesThatWillRemain <- samplesThatWillRemain %>%
    tidyr::gather( Variable, Value, -numberOfReads, -SampleID ) %>%
    dplyr::group_by( Variable, Value ) %>%
    dplyr::summarise( SamplesRemained = n() ) %>%
    dplyr::mutate( FreqOfRemainingSamples = SamplesRemained / sum(SamplesRemained),
                   TotalNumberOfSamplesRemaining = sum(SamplesRemained) ) %>%
    dplyr::full_join( propSamplesThatWillBeRemoved, 
                      by = c("Variable", "Value") ) %>%
    dplyr::mutate( TotalSamples = sum(c(SamplesRemained, SamplesRemoved), na.rm = T) ) %>%
    dplyr::ungroup( )
  
  
  # Calculate differences from perfect proportion
  fullStats <- propSamplesThatWillRemain %>%
    dplyr::left_join( propSamplesThatWillRemain %>%
                          dplyr::group_by( Variable ) %>%
                          dplyr::summarise( subVariables = n() ),
                      by = "Variable" ) %>%
    dplyr::mutate( perfectRatio = (1/subVariables) )
  
  write.csv(fullStats, 
            file = file.path(intermediateFolderRarefying, paste0("LogForRarefyingValue_", depthRangeValue, ".csv") ))
    
  # Define the similarity index
  disimilarityIndexResults <- disimilarityIndexResults %>%
    dplyr::bind_rows( data.frame( NumberOfReads = depthRangeValue,
                                  disimilarityIndex = sum(abs(fullStats$FreqOfRemainingSamples - fullStats$perfectRatio), na.rm = F),
                                  minProportionCheck = any( fullStats$FreqOfRemainingSamples < minimumProp ),
                                  stringsAsFactors = F) )
  
}

print(ggplotly(ggplot(disimilarityIndexResults,
       aes(x = NumberOfReads, y = disimilarityIndex, colour = minProportionCheck) ) +
  geom_line( ) +
  geom_point() +
  scale_color_manual( "It didn't pass min proportion check", values = c("black", "red") ) +
  theme( legend.position = "bottom" )))

message(paste0("You are analysing organism: ", organism))
message(paste0("The full logs of the rarefying options are in folder: \n", intermediateFolderRarefying))

numberOfReadsToUse <- select.list(depthRange, title = "Enter the number of reads to be used for rarefying: \n")

toKeep <- c(toKeep, "numberOfReadsToUse")

rm( list = ls()[!(ls() %in% toKeep)])
