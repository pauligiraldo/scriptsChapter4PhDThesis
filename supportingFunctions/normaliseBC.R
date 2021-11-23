message("Normalising data using BC or Bias Correction")

#
# Apply Bias correction
#

# Set the seed for ramdom number
set.seed(4002)

# Run ancomb
bcData <- ancombc(phyloseq = biomData, formula = "Season * Management * Region", 
                  p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
                  group = "Management", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                  max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

# Save ancombc results to a .RData file
save(biomData, file = paste0(intermediateFolder, "/biomData.RData"))
save(bcData, file = paste0(intermediateFolder, "/biomDataNormalised.RData"))

# re-Build phyloseq object from BC's feature table
otuData <- otu_table(bcData$feature_table, taxa_are_rows = T)

biomDataNormalised <- phyloseq(otuData, sample_data(biomData), tax_table(biomData))

# Generate data profile after rarefying
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
  print("test fin bc")
  write.csv(dataSizeLog, file = paste0(outputFolder, "/dataProfileBeforAndAfterNormalising.csv"), row.names = F )
}