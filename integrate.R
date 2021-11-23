# Purpose: report on differences between Management, Region and Season 
# Date: 07/03/2020
# Type: this is the integrate. It connects all analyses.
# Author: Paulina Giraldo-Perez

# How it works:
# This script connects all the analyses performed in chapter 4.
# Run the initial setup script found in line 18, then scroll down to find the script use in each 
# specific analysis.
# At the top of each one of these scripts you will find parameters or variables you can modify to
# fit your analyses.



#
# Initial setup | Sample loading | Produce phyloseq object
#

# Clean environment and load libraries
source("supportingFunctions/initialSetup.R")



##
### Analysis Types and Abundances
##
#

source("analyses/typesAndAbundances.R")

# !! NOTE !!
# The following parameters can be modified directly accessing the script "analyses/typesAndAbundances.R"

# # Define the organisms to analyse. Options: "bacteria", "fungi", "metazoa" and/or "all".

# targetOrganisms <- c("bacteria", "fungi", "metazoa", "all")
 
# # Define years to use. The dataset used in this analyses contained 4 years of data, however, only years 1 to 3 were used.

# yearsToUse <- 1:3

# # Choose the normalisation method to use. Options are: "rarefy", "bc" [bias correction] or "notNormalised".
# 
# normaliseMethod <- "rarefy"
# 
# # If the normalisation method chosen is "rarefy", use the options below for further 
# # configure how rarefying will be applied:
# 
# # Test sample size for rarefication:
# # The variable or objects below depthRangeFrom, 
# # depthRangeTo and depthRangeBy can be used to define the range of reads or sample depths that should be tested.
# # The default values of depthRangeFrom = 500, depthRangeTo = 10000 and depthRangeBy = 500. These values would result
# # in testing from 500 reads to 10000 reads by 500 reads, or the ranges below:
# 
# # 500  1000  1500  2000  2500  3000  3500  4000  4500  5000  5500  6000  6500  7000  7500  8000  8500  9000  9500 10000
# 
# # If you don't want to rarefy, just sets depthRangeFrom, depthRangeTo and depthRangeBy to 0, or zero.
# 
# # Your turn:
# depthRangeFrom = 500
# depthRangeTo = 10000
# depthRangeBy = 100
# 
# #### Finishes choosing normalisation method ##############
# 
# # Define minimum proportion os samples per data class
# # For example, after rarefying, my final dataset cannot contain a data class, like Management:Future, represented by
# # less than the minimumPropToUse:
# minimumPropToUse <- 0.25 
# 
# # There is an option to create a plot showing the impat of rarefying
# shouldGeneratePlotNReads <- F









##
### Analysis Time Series
##

source("analyses/timeSeries.R" )

source("analyses/timeSeriesAcrossBarcodes.R" )

# !! NOTE !!
# The following parameters can be modified directly accessing the scripts "analyses/timeSeries.R" and "analyses/timeSeriesAcrossBarcodes.R" 

## !!!! Time series analyses depend on the following intermediate objects created by the script "analyses/typesAndAbundances.R":
# paste0(intermediateFolder, "/dataNumberOfTaxaLevelTaxon_", taxonomicLevel, ".RData")
# paste0(intermediateFolder, "/datashannonAndSimpsonLevelTaxon_", taxonomicLevel, ".RData")
# paste0(intermediateFolder, "/dataTypesAndAbundancesOfTaxaLevelTaxon_", taxonomicLevel, ".RData")
# The intermediate folder is:
# paste0("intermediateData/", organism, "/", normaliseMethod)

# normaliseMethod <- "rarefy"
# 
# # Organisms to be analysed:
# organismsToBuildTimeSeries <- c("all", "bacteria", "fungi", "metazoa")



##
### Analysis Indicative Species
##

source("analyses/indicativeSpecies.R" )

# Analysis indicative species above species
source("analyses/indicativeSpeciesFamilyAndAbove.R" )

# !! NOTE !!
# The following parameters can be modified directly accessing the scripts "analyses/indicativeSpecies.R" and "analyses/indicativeSpeciesFamilyAndAbove.R" 

# # Define the organisms to analyse: "bacteria", "fungi", "metazoa", "all"
# targetOrganisms <- c("bacteria", "fungi", "metazoa", "all")
# 
# # Define years to use
# yearsToUse <- 1:3
# 
# # Permutation value to use in the multipatti function
# permutationValue <- 100
# 
# #
# # Choose the normalisation method to use: "rarefy", "bc" [bias correction] or "notNormalised"
# #
# normaliseMethod <- "rarefy"
# 
# # If the normalisation method chosen is "rarefy", use the options below for further 
# # configure how rarefying will be applied:
# 
# # Test sample size for rarefication:
# # The variable or objects below depthRangeFrom, 
# # depthRangeTo and depthRangeBy can be used to define the range of reads or sample depths that should be tested.
# # The default values of depthRangeFrom = 500, depthRangeTo = 10000 and depthRangeBy = 500. These values would result
# # in testing from 500 reads to 10000 reads by 500 reads, or the ranges below:
# 
# # 500  1000  1500  2000  2500  3000  3500  4000  4500  5000  5500  6000  6500  7000  7500  8000  8500  9000  9500 10000
# 
# # If you don't want to rarefy, just sets depthRangeFrom, depthRangeTo and depthRangeBy to 0, or zero.
# 
# # Your turn:
# depthRangeFrom = 500
# depthRangeTo = 10000
# depthRangeBy = 100
# 
# #### Finishes choosing normalisation method ##############
# 
# # Define minimum proportion os samples per data class
# # For example, after rarefying, my final dataset cannot contain a data class, like Management:Future, represented by
# # less than the minimumPropToUse:
# minimumPropToUse <- 0.25 
# 
# shouldGeneratePlotNReads <- F







##
### Analysis Pathogenics
##
# 

source("analyses/pathogenics.R")

  


