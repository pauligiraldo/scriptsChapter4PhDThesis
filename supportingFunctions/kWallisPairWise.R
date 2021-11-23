# Test pairwisely
kWallisPairWise <- function(data, nameColumnWithNumberOfTaxa, factor) {
  
  # data = otuData
  # nameColumnWithNumberOfTaxa = "abundance"
  # factor = factor
  
  resultsData <- data.frame()
  
  for(factorIndex in 1:length(factor)){
    
    # message("KWallis pairwise....")
    # message(paste0("Factor to test: ", factor[factorIndex]))
    # message("Data:")
    # print(data[,factor[factorIndex]])
    
    kTest <- try(kruskal.test(data[,nameColumnWithNumberOfTaxa] ~ data[,factor[factorIndex]]), T)
    epSquared <- try(epsilonSquared(x = data[,nameColumnWithNumberOfTaxa], g = data[,factor[factorIndex]]), T)
    dunnTestResults <- try(dunnTest(data[,nameColumnWithNumberOfTaxa], data[,factor[factorIndex]], method = "bonferroni"), T)
    # Shapiro test
    shapiroTest <- try(shapiro.test( data[,nameColumnWithNumberOfTaxa] ), T)
    
    if( class(dunnTestResults) == "try-error" ) {
      
      kResultsPre <- data.frame( comparison = "Error",
                                 pvalueKTest = ifelse( class(kTest) == "try-error",  "Error", as.character(kTest$p.value)),
                                 pvalueAdjustedDunnTest = "Error",
                                 chiSquaredKTest = ifelse( class(kTest) == "try-error",  "Error", as.character(as.numeric(kTest$statistic))),
                                 epsilonSquared = ifelse( class(epSquared) == "try-error",  "Error", as.character(as.numeric(epSquared))),
                                 factor = factor[factorIndex], 
                                 shapiroPValue = ifelse( class(shapiroTest) == "try-error",  "Error", as.character(as.numeric(shapiroTest$p.value))), 
                                 shapiroW = ifelse( class(shapiroTest) == "try-error",  "Error", as.character(as.numeric(shapiroTest$statistic))) )
      
    } else {
      
      kResultsPre <- data.frame( comparison = ifelse( rep( class(dunnTestResults) == "try-error", length(dunnTestResults$res$Comparison) ),  "Error", as.character(dunnTestResults$res$Comparison)),
                                 pvalueKTest = ifelse( class(kTest) == "try-error",  "Error", as.character(kTest$p.value)),
                                 pvalueAdjustedDunnTest = ifelse( rep( class(dunnTestResults) == "try-error", length(dunnTestResults$res$P.adj) ),  "Error", as.character(dunnTestResults$res$P.adj)),
                                 chiSquaredKTest = ifelse( class(kTest) == "try-error",  "Error", as.character(as.numeric(kTest$statistic))),
                                 epsilonSquared = ifelse( class(epSquared) == "try-error",  "Error", as.character(as.numeric(epSquared))),
                                 factor = factor[factorIndex], 
                                 shapiroPValue = ifelse( class(shapiroTest) == "try-error",  "Error", as.character(as.numeric(shapiroTest$p.value))), 
                                 shapiroW = ifelse( class(shapiroTest) == "try-error",  "Error", as.character(as.numeric(shapiroTest$statistic))) )
    }
    
    resultsData <- resultsData %>%
      dplyr::bind_rows( kResultsPre )
  }
  
  return( resultsData %>%
           dplyr::mutate( measure = nameColumnWithNumberOfTaxa ) )
  
}