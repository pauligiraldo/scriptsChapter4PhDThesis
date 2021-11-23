# Test pairwisely
kWallisInteractions <- function(data, nameColumnWithNumberOfTaxa, factor) {
  
  # data = otuData
  # nameColumnWithNumberOfTaxa = "abundance"
  # factor = factor
  
  resultsData <- data.frame()
  
  interactionsToTest <- list()
  
  for(i in 2:length(factor)){
    
    interactionsToTest <- append( interactionsToTest, combn( factor, m = i, simplify = F ) )
    
  }
  
  for(interactionIndex in 1:length(interactionsToTest)){
  
    interactionObject <- interaction(data[,interactionsToTest[[interactionIndex]]])
  
    kTest <- try(kruskal.test(data[, nameColumnWithNumberOfTaxa] ~ interactionObject), T)
  
    kResultsPre <- data.frame( interaction = paste0(interactionsToTest[[interactionIndex]], collapse = " | "),
                               pvalue = ifelse( class(kTest) == "try-error",  "Error", as.character(kTest$p.value)),
                               chiSquared = ifelse( class(kTest) == "try-error",  "Error", as.character(kTest$statistic))
    )
  
    resultsData <- resultsData %>%
      dplyr::bind_rows( kResultsPre )
  }
  
  return( resultsData %>%
           dplyr::mutate( measure = nameColumnWithNumberOfTaxa ) )
  
}