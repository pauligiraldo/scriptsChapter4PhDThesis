#
##
###
####
##### Report OTUs by organism and proportion of kindom and phylum
####
###
##
#

# Clean environment and load libraries
source("supportingFunctions/initialSetup.R")

library(RColorBrewer)

rm(list = ls())  

# Define the organisms to analyse: "bacteria", "fungi", "metazoa", "all"
targetOrganisms <- c("bacteria", "fungi", "metazoa")
# targetOrganisms <- c("all")

# Define years to use
yearsToUse <- 1:3

# The object dataSizeLog is used to store the profile of samples per factor level for each organism before and after rarefying.
dataSizeLog <- data.frame()

# Define objects to save
toKeep <- ls()
toKeep <- ls()
toKeep <- c(toKeep, "biomData", "biomDataPre", "organism", "shouldTestValuesForRarefying", "readsValueForRarefying",
            "outputFolder", "intermediateFolder", "shouldGeneratePlotNReads", "dataSizeLog",
            "mappingFile", "taxonomyDataFile", "normaliseMethod", "taxLevels")

for(organism in targetOrganisms){
  
  message(paste0("Analysing data for ", organism, "...."))
  
  switch(organism,
         bacteria = {
           dataFile <- "data/bacteriaASVTable.csv"
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomyDataFile <- "data/bacteriaTaxonomy.csv"
           readsValueForRarefying <- 2600
         },
         fungi = {
           dataFile <- "data/itsASVTable.csv"
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomyDataFile <- "data/itsTaxonomy.csv"
           readsValueForRarefying <- 2000
         },
         metazoa = {
           dataFile <- "data/metazoaASVTable.csv"
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomyDataFile <- "data/metazoaTaxonomy.csv"
           readsValueForRarefying <- 2200
         },
         all = {
           dataFile <- c("data/bacteriaASVTable.csv", "data/itsASVTable.csv", "data/metazoaASVTable.csv") 
           mappingFile <- "data/metadata_4Years.csv"
           taxLevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
           taxonomyDataFile <- c("data/bacteriaTaxonomy.csv", "data/itsTaxonomy.csv", "data/metazoaTaxonomy.csv")
           readsValueForRarefying <- 2000
         } 
  )
  
  source("supportingFunctions/loadSamples.R")
  
  # Loaded samples
  message("Loaded samples.....")
  
  # Select samples from the desired years
  biomData <- subset_samples(biomDataPre, Year_of_program %in% yearsToUse)
  taxData <- prune_taxa(taxa_sums(biomData) > 0, biomData)
  
  allTaxData <- data.frame(tax_table(taxData), stringsAsFactors = F)
  
  # Number of OTUs
  # Percentage of each item in Kingdom, Phylum, Class, etc...
  message(organism)
  message(paste0("Total number of OTUs: ", nrow(allTaxData)))
  
  for(i in 1:length(taxLevels)){
    
    message(paste0("Percentage of taxa by ", taxLevels[i], ":"))
    
    frequencies <- as.data.frame(table(allTaxData[,taxLevels[i]]))
    
    totalOTUs <- sum(frequencies$Freq)
    
    frequenciesAsTable <- frequencies %>%
      dplyr::rename( "TaxName" = Var1,
                     "Percent" = Freq ) %>%
      dplyr::mutate( Percent = (Percent/totalOTUs)*100,
                     TaxLevel = taxLevels[i],
                     Organism = organism,
                     TotalNumberOfOTUs = totalOTUs ) %>%
      dplyr::arrange( desc(Percent) )
    
    print(paste0("Sum of frequencies: ", sum(frequenciesAsTable$Percent)))
    
    if(file.exists("output/summaryOTUsByOrganismAndTaxLevel.csv")){
      write.table(frequenciesAsTable, file = "output/summaryOTUsByOrganismAndTaxLevel.csv", sep = ",", col.names = !file.exists("output/summaryOTUsByOrganismAndTaxLevel.csv"), append = T, row.names = F)
    } else {
      write.csv(frequenciesAsTable, file = "output/summaryOTUsByOrganismAndTaxLevel.csv", row.names = F )
    }
    
  }
  
}


#
## Produce a bar plot with the results
#

rawData <- read.csv("output/summaryOTUsByOrganismAndTaxLevel.csv", stringsAsFactors = F)

# Bacteria Phylum
groomedDataBacteriaPhylum <- rawData %>%
  dplyr::filter( TaxLevel %in% c("Phylum"),
                 Organism %in% c("bacteria") ) %>%
  dplyr::arrange( (Percent) ) %>%
  dplyr::mutate( TaxName = gsub("[", "", TaxName, fixed = T),
                 TaxName = gsub("]", "", TaxName, fixed = T),
                 TaxName = factor(TaxName, levels = TaxName) )
  

p1 <- ggplot( groomedDataBacteriaPhylum, aes(x = Percent, y = TaxName, fill = TaxName, label = percent(Percent/100, accuracy = 0.001) ) ) +
  geom_bar(stat="identity", color = "black") +
  labs(title="16S - Phyla",
       x ="Percent (%)", y = "") +
  geom_text( nudge_x = 2, size = 5 ) +
  #ggrepel::geom_text_repel( position = position_stack(vjust = 0.5), max.overlaps = 40 ) +
  theme_bw() +
  theme( legend.position = "none",
         axis.text.y = element_text(size = 15),
         axis.text.x = element_text(size = 15),
         axis.title.x = element_text(size = 15),
         plot.title = element_text(size = 20) ) 

pdf(file = "output/summaryBacteriaPhyla.pdf", width = 11, height = 10)
print(p1)
dev.off()

# Bacteria Kingdom
groomedDataBacteriaKingdom <- rawData %>%
  dplyr::filter( TaxLevel %in% c("Kingdom"),
                 Organism %in% c("bacteria") ) %>%
  dplyr::arrange( (Percent) ) %>%
  dplyr::mutate( TaxName = factor(TaxName, levels = TaxName) )


p2 <- ggplot( groomedDataBacteriaKingdom, aes(x = Percent, y = TaxName, fill = TaxName, label = percent(Percent/100, accuracy = 0.001) ) ) +
  geom_bar(stat="identity", color = "black") +
  labs(title="16S - Kingdom",
       x ="Percent (%)", y = "") +
  geom_text( nudge_x = 5, size = 6 ) +
  #ggrepel::geom_text_repel( position = position_stack(vjust = 0.5), max.overlaps = 40 ) +
  theme_bw() +
  theme( legend.position = "none",
         axis.text.y = element_text(size = 25),
         axis.text.x = element_text(size = 25),
         axis.title.x = element_text(size = 20),
         plot.title = element_text(size = 25) ) 

pdf(file = "output/summaryBacteriaKindom.pdf", width = 15, height = 10)
print(p2)
dev.off()



# Bacteria Phylum
groomedDataFungiPhylum <- rawData %>%
  dplyr::filter( TaxLevel %in% c("Phylum"),
                 Organism %in% c("fungi") ) %>%
  dplyr::arrange( (Percent) ) %>%
  dplyr::mutate( TaxName = gsub("unidentified", "Unassigned", TaxName, fixed = T),
                 TaxName = factor(TaxName, levels = TaxName) )


p3 <- ggplot( groomedDataFungiPhylum, aes(x = Percent, y = TaxName, fill = TaxName, label = percent(Percent/100, accuracy = 0.001) ) ) +
  geom_bar(stat="identity", color = "black") +
  labs(title="ITS2 - Phyla",
       x ="Percent (%)", y = "") +
  geom_text( nudge_x = 2.5, size = 5 ) +
  #ggrepel::geom_text_repel( position = position_stack(vjust = 0.5), max.overlaps = 40 ) +
  theme_bw() +
  theme( legend.position = "none",
         axis.text.y = element_text(size = 20),
         axis.text.x = element_text(size = 20),
         axis.title.x = element_text(size = 15),
         plot.title = element_text(size = 20) ) 

pdf(file = "output/summaryFungiPhylum.pdf", width = 12, height = 10)
print(p3)
dev.off()

# Fungi Kingdom
groomedDataFungiKingdom <- rawData %>%
  dplyr::filter( TaxLevel %in% c("Kingdom"),
                 Organism %in% c("fungi") ) %>%
  dplyr::arrange( (Percent) ) %>%
  dplyr::mutate( TaxName = gsub("unidentified", "Unassigned", TaxName, fixed = T),
                 TaxName = factor(TaxName, levels = TaxName) )


p4 <- ggplot( groomedDataFungiKingdom, aes(x = Percent, y = TaxName, fill = TaxName, label = percent(Percent/100, accuracy = 0.001) ) ) +
  geom_bar(stat="identity", color = "black") +
  labs(title="ITS2 - Kindom",
       x ="Percent (%)", y = "") +
  geom_text( nudge_x = 4, size = 5 ) +
  #ggrepel::geom_text_repel( position = position_stack(vjust = 0.5), max.overlaps = 40 ) +
  theme_bw() +
  theme( legend.position = "none",
         axis.text.y = element_text(size = 20),
         axis.text.x = element_text(size = 20),
         axis.title.x = element_text(size = 15),
         plot.title = element_text(size = 20) ) 

pdf(file = "output/summaryFungiKingdom.pdf", width = 12, height = 8)
print(p4)
dev.off()




    