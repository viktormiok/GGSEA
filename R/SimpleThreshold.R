#' Sub-Function for Simple Threshold Enrichment
#' @description Internal use. Function  to performed Fisher exact test from predefined list
#' @keywords internal
#' @importFrom broom glance
#' @importFrom dplyr rename
#' @return Fisher exact tests results for the given pathways and FeatureID as a data frame
EnrichFisherPathway <- function(FAThres,FBThres,PathwayList,FeatureID){
  #Construct table for Chisquare test
  #Features Before Threshold part of the Pathway of interest
  FeaturesBT_InP <- length(FBThres[FBThres %in% PathwayList])
  #Features Before Threshold not part of the Pathway of interest
  FeaturesBT_NotInP <- length(FBThres) - FeaturesBT_InP
  #Features After Threshold part of the Pathway of interest
  FeaturesAT_InP <- length(FAThres[FAThres %in% PathwayList])
  #Features After Threshold not part of the Pathway of interest
  FeaturesAT_NotInP <- length(FAThres) - FeaturesAT_InP
  #Construct table
  TMatrix <- matrix(c(FeaturesBT_InP,FeaturesBT_NotInP,
                      FeaturesAT_InP,FeaturesAT_NotInP),
                    ncol = 2, nrow = 2)
  rownames(TMatrix) <- c("InP","NotInP")
  colnames(TMatrix) <- c("InL","NotInL")
  
  #Fisher's exact test
  FResultsPathway <- fisher.test(TMatrix, alternative = "greater")
  FResultsPathwayT <- glance(FResultsPathway)
  FResultsPathwayT <- FResultsPathwayT %>% dplyr::rename(oddsRatio = estimate)
  FResultsPathwayT$FeatureID <- FeatureID
  return(FResultsPathwayT)
}

#' Simple Threshold
#' @description This function tests for every gene of a distance matrix, if the group of closest genes (defined by threshold) 
#' are enriched for the given pathways
#' @param DistM Distance Matrix
#' @param PathwayL List of pathways
#' @param desc logical. If TRUE sorts distances row descending
#' @param thres_quant Quantile threshold to define the list to compare
#' @return Enrichment results for all genes and pathways as a data frame
SimpleThreshold <- function(DistM,PathwayL,desc,thres_quant, ...){
  ResultsRows <- list() #store results later
  for(rw in 1: ncol(DistM)){
    # print progress
    progress(rw, ncol(DistM))
    
    #Feature ID
    FID <- colnames(DistM)[rw]
    
    #Get Row names
    rwIterColnames <- colnames(DistM[rw,-rw,drop=FALSE])
    #Get values row
    rwIter <- DistM[rw,-rw] #Remove distance to itself
    
    #Name the vector
    names(rwIter) <- rwIterColnames
    #Sort descending
    rwIter <- sort(rwIter,decreasing = desc)
    
    #Threshold
    threshold <- unname(quantile(rwIter,thres_quant))
    
    #Before threshold
    FeaturesBThres <- names(rwIter[rwIter < threshold])
    #After threshold
    FeaturesAThres <- names(rwIter[rwIter >= threshold])
    
    #Loop pathways and apply enrichment
    PathwaysGene <- lapply(PathwayL, function(pathway){
      EnrichFisherPathway(FAThres = FeaturesAThres,
                          FBThres = FeaturesBThres, 
                          PathwayList = pathway,
                          FeatureID = FID)
    })
    
    #Join results
    PathwaysGene <- do.call(rbind,PathwaysGene)
    #Move rownames to column
    PathwaysGene <- PathwaysGene %>% rownames_to_column(var = "Pathway")
    #Store results
    ResultsRows[[rw]] <- PathwaysGene
    
  }
  ResultEnrichment <- do.call(rbind,ResultsRows)
  ResultEnrichment$p.adj <- p.adjust(ResultEnrichment$p.value, ...)
  ResultEnrichment <- ResultEnrichment %>% relocate(FeatureID, oddsRatio, p.value, p.adj)
  return(ResultEnrichment)
}