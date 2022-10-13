#'Simple enrichment
#'@description This function takes as input a distance matrix and a list of specific pathway IDs.
#'The function performs an enrichment test of the closest elements of every feature.  
#'The closest elements are selected based on a given quantile 
#'@param matrix Distance matrix
#'@param PathwayL List of all IDs of a given pathway. IDs between the list and the matrix should be compatible
#'@param desc logical. When TRUE the distances are ordered in a descending manner NOTE_DEP: probably not necessary as parameter... probably remove
#'@param thres_quant Quantile used as threshold to select the closest elements
SEnrich <- function(DistM, PathwayL, TestType, desc = TRUE, thres_quant=0.05, ...){
  #####
  #Idea TO DO Function Improvement:
  #Input as list... alternatively we could provide a df with multiple ids... 
  #Create some function that automatically detects the type of ID... 
  #loops first entry of all columns and finds the ID that matches the format of the input matrix, 
  #if not found throw error to inform the user
  
  #Remove NAs if there are any in the pathways
  PathwayL <- sapply(PathwayL,function(pathway){
    pathway[!is.na(pathway)]
  })
  
  REnrichment <- switch(TestType, 
                        SimpleThreshold = {
                          SimpleThreshold(DistM = DistM,PathwayL = PathwayL,
                                          desc = desc,thres_quant = thres_quant)
                        },
                        DistributionBased={
                          ###############################
                          # qickndirty implementation of a distribution based test
                          # needs to go into own loop with switch like function
                          
                          # intersection if matrix genes with actual PW list
                          ixPWAct <- intersect(colnames(DistM), PathwayL)
                          
                          # get only relevant distances
                          DistMPW <- DistM[,ixPWAct]
                          rwMax <- rowMax(DistM)
                          dist2maxact <- 1:nrow(DistMPW)
                          # go for each row/gene in dist matrix
                          for(rw in 1:nrow(DistMPW)){
                            # moving window size k
                            nmax= rwMax[rw]
                            k <- nmax/100
                            scorek <- 1:nmax
                            ddistact <- DistMPW[rw,]
                            # go along distance calculate local score
                            # this must be dome easier <- use some local density gaussian stuff here
                            for(i in 1:nmax){
                              scorek[i] <- sum(ddistact < i+k & ddistact > i-k)
                            }
                            # get the distance to the maximum k-density
                            dist2maxact[rw] <- which.max(scorek)
                          }
                          
                          # create null model <- takes ages, need to speed up
                          nRand <- 100
                          dist2maxsample = 1:nRand
                          for(j in 1:nRand){
                            # print progress
                            progress(j, nRand)
                            
                            ixPWsmpl = sample.int(ncol(DistM), length(ixPWAct))
                            DistMsmpl <- DistM[,ixPWsmpl]
                            rwMax <- rowMax(DistM)
                            dist2maxsmpl <- 1:nrow(DistMsmpl)
                            
                            # !! ATTENTION need to outsource this together with upstream test REDUNDANT
                            for(rw in 1:nrow(DistMsmpl)){
                              # moving window size k
                              nmax= rwMax[rw]
                              k <- nmax/100
                              scoresmpl <- 1:nmax
                              ddistsmpl <- DistMsmpl[rw,]
                              # go along distance calculate local score
                              # this must be dome easier <- use some local density gaussian stuff here
                              for(i in 1:nmax){
                                scoresmpl[i] <- sum(ddistsmpl < i+k & ddistsmpl > i-k)
                              }
                              dist2maxsmpl[rw] <- which.max(scoresmpl)
                            }
                            if (j == 1) {
                              NullDistMax =  dist2maxsmpl
                            } else {
                              NullDistMax = cbind(NullDistMax, dist2maxsmpl)
                            }
                          }
                          
                          # Quantile testing
                          # not working well yet some genes show generally higher distances than overall mean...
                          for(rw in 1:nrow(NullDistMax)){
                            q <- unname(quantile(NullDistMax[rw,], 0.05))
                            # do calculate manual p value?
                            if(dist2maxact[rw] < q){
                              # jep, we have a hit!
                              hist(NullDistMax[rw,])
                            }
                          }
                          
                          
                        },
                        {
                          print('Option not recognized')
                        }
  )
  
  
  
  #TO DO:
  #When creating density plots it should match row order of the heatmap
  
  return(REnrichment)
}
