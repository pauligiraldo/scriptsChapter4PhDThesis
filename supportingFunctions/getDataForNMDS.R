# Purpose: a function to generate NMDS data so we can produce graphs manually
# Date: 25/08/2019
# Author: Paulina Giraldo Perez

getDataForNMDS <- function(data,
                           distance = "jaccard",
                           minCounts = 1,
                           minSamplePercent = 0.05,
                           binary = F) {
  
  # Test removing OTUs that do not appear more than minCounts times in more than minSamplePercent of the samples
  howToFilter <- genefilter_sample(data, filterfun_sample(function(x) x > minCounts), A=minSamplePercent*nsamples(data))
  filteredDataForNMDS <- prune_taxa(howToFilter, data)
  
  ordinatedObj <- ordinate(physeq = filteredDataForNMDS, method = "NMDS", distance = distance, binary = binary)
  
  plotData <- plot_ordination(filteredDataForNMDS, ordinatedObj, type="split", color=NULL, title="", justDF = T)
  
  return(plotData)
  
}

# distanceMatrixBinary <- distance(taxData, method="jaccard", binary = T)
# 
# ordinatedObjBinary <- ordinate(taxData, method = "NMDS", distance = "jaccard", binary = T)
# 
# distanceMatrix <- distance(taxData, method="jaccard", binary = F)
# 
# ordinatedObj <- ordinate(taxData, method = "NMDS", distance = "jaccard")
# 
# plotData <- plot_ordination(data, ordinatedObj, type="split", color=NULL, title="", justDF = T)
# 
