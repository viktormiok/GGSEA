## Load libraries 
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)

annotation_file <- data.table::fread("./Rdata/data/MouseGeneMapping.csv") %>% tibble::column_to_rownames(var = "Gene stable ID")

## get the ENTRZID for each ENSEMB and remove duplicates in both 
annotation_file <- clusterProfiler::bitr(geneID = as.character(row.names(annotation_file)),
                                         fromType="ENSEMBL",
                                         toType = "ENTREZID",
                                         OrgDb = org.Mm.eg.db,
                                         drop = TRUE) %>% 
  dplyr::distinct(ENSEMBL,.keep_all = TRUE) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames(var = "ENSEMBL") %>% dplyr::distinct(ENTREZID)

## count number of duplicates in the ENTERZID if present 
sum(duplicated(annotation_file$ENTREZID)) 


## read the expression data
expr_matrix <- data.table::fread("./Rdata/data/GSE112582_raw_count.csv") %>% tibble::column_to_rownames(var = "ID")

## Match the rownames of the annotation_file with the expr_matrix DF
all(row.names(expr_matrix) %in% row.names(annotation_file))
length(which(row.names(expr_matrix) %in% row.names(annotation_file))) # 29226 rownames are matched

##take the expression values of matched IDs
final_data_expre <- expr_matrix[rownames(annotation_file),] %>% tidyr::drop_na()  

## change the rowname of expression data from ENSEMBL TO ENTEREZID
expression_rownames <- clusterProfiler::bitr(geneID = as.character(row.names(final_data_expre)),
                                             fromType="ENSEMBL",
                                             toType = "ENTREZID",
                                             OrgDb = org.Mm.eg.db,
                                             drop = TRUE) %>%
  dplyr::distinct(ENSEMBL,.keep_all = TRUE) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames(var = "ENSEMBL") %>% dplyr::distinct(ENTREZID)

##make sure that the rownames of both are identical 
identical(row.names(final_data_expre), row.names(expression_rownames))
row.names(final_data_expre) <- as.character(expression_rownames$ENTREZID)

## meta data(sample information)
sample_data <- data.frame(replicate= as.factor(names(final_data_expre)),
                          group = as.factor(rep(c("BAT", "BAT_cold", "WAT", "WAT_cold"), c(4,4,4,4))))

## save the file 
write.csv(sample_data, "sample_data.csv")

## Pre-filtering 
dds <- DESeqDataSetFromMatrix(countData = final_data_expre, 
                              colData = sample_data, 
                              design =  ~ group)


## include genes higher than 0 read on each sample 
keep <- rowSums(counts(dds)) > 0   
dds <- dds[keep,]
nrow(dds) ## 31784

## Further I include genes that have a total of 4000 reads in at least 4 samples 
keep <- rowSums(counts(dds) >= 9000) >= 8
dds <- dds[keep,]
nrow(dds)# 4342

## Extract the count matrix from dds object and pass the count matrix to biofig  plot_ExploratoryPlots function 
## for visualization 
filter_data<-assay(dds)
# library(BioFig)
# plot_ExploratoryPlots(expression=filter_data,
#                       savePlots = F,
#                       group=sample_data$group,
#                       samplenames=colnames(filter_data),
#                       plottype="all",
#                       point.size=3,
#                       size_gglab = 4)


## estimating Dispersion
dds <- DESeq(dds)

## transform the count
vsd_BAT_WAT_Data <- vst(dds, blind = FALSE) 
#write.csv(as.data.frame(assay(vsd_BAT_WAT_Data)), file = "./RData/data/GSE112582_processed_count.csv")
 
#Get expression matrix
expr_matrix <- data.table::fread("./RData/data/GSE112582_processed_count.csv") %>% tibble::column_to_rownames(var = "V1")

# plot expression matrix
heatmapObj <- pheatmap(expr_matrix, cluster_cols = FALSE)

##################
#Distance Matrix
distanceMatrixBAT_WAT <- as.matrix(dist(expr_matrix))

## Function to read the .gmt file and save the result as list of list
readGeneSet_gmt <- function(file){
  if(!file.exists(file)){
    stop('File',file,' not available\n')
  }
  
  if (!grepl("*.gmt", file)[1]){
    stop("Pathway information must be a .gmt file")
  }
  geneSetData = strsplit(readLines(file), "\t")
  names(geneSetData) = sapply(geneSetData, "[", 1)
  geneSetData = lapply(geneSetData, "[", -1:-2)
  geneSetData = lapply(geneSetData, function(x) {
    x[which(x != "")]
  })
  return(geneSetData)
  
}

KEGG_GeneSets <- readGeneSet_gmt("c2.cp.kegg.v7.4.symbols.gmt")

#Read new Pathway Info
library(tidyverse)
library(pathwayPCA)
KEGG_GeneSets <- read_gmt("./RData/data/GeneSets/c2.cp.kegg.v7.4.symbols.gmt")
KEGG_l <- KEGG_GeneSets$pathways
names(KEGG_l) <- KEGG_GeneSets$TERMS
KEGG_l[[1]]

#Read Data with the pathway information
Alzheimer <- read.csv("./RData/data/ListPathways/Alzheimer_disease.csv",row.names = 1)
DiabetesCardio <- read.csv("./RData/data/ListPathways/Diabetic_cardiomyopathy.csv", row.names = 1)
PathwaysTest <- list(Alzheimer = Alzheimer$ENSEMBL,DiabetesCardio = DiabetesCardio$ENSEMBL)


################
# do the DistMEnrichment

source("./R/SEnrich.R")
source("./R/SimpleThreshold.R")
debug(SEnrich)
ERes <- SEnrich(DistM = distanceMatrixBAT_WAT,TestType = "SimpleThreshold", PathwayL = PathwaysTest)

#plot(ERes$p.value[ord],nrow(ERes):1, type="s", xlim=c(0,1), col="green")

# plot the heatmap of the samples
pheatmap(expr_matrix, cluster_cols = FALSE, show_rownames=F)
# plot the heatmap of the ERes for pathways
pheatmap(cbind(-log10(ERes_dc$p.value[ord]),-log10(ERes_a$p.value[ord])), cluster_cols = FALSE, cluster_rows = FALSE, show_rownames=F)




##############################
#Transition probability Matrix
tprobabilityMatrixBAT_WAT<-apply(distanceMatrixBAT_WAT,1,function(r){
  totSumRow <- sum(r)
  1 - r/totSumRow
})
max(tprobabilityMatrixBAT_WAT)


###########
#Option 1 Simple Enrichment/ Guilt by association 1

## get the list of mouse KEEG pathways


## retrieve gene sets from KEGG database for a mouse
genesets_kegg_pathways_mm <- EnrichmentBrowser::getGenesets(org = "mmu", db = "kegg", cache = TRUE, return.type="list")

## retrieve gene sets involved in Alzheimer disease  
Alzheimer_disease <- data.frame(genesets_kegg_pathways_mm$mmu05010_Alzheimer_disease) 
Alzheimer_disease <- clusterProfiler::bitr(geneID = as.character(Alzheimer_disease$genesets_kegg_pathways_mm.mmu05010_Alzheimer_disease),
                                       fromType="ENTREZID",
                                       toType = c("SYMBOL","ENSEMBL"),
                                       OrgDb = org.Mm.eg.db,
                                       drop = TRUE)
write.csv(Alzheimer_disease, "./RData/data/ListPathways/Alzheimer_disease.csv")

##gene sets involved in the diabetic cardiomyopathy
DiabetesCardio <- data.frame(genesets_kegg_pathways_mm$mmu05415_Diabetic_cardiomyopathy)
DiabetesCardio <- clusterProfiler::bitr(geneID = as.character(DiabetesCardio$genesets_kegg_pathways_mm.mmu05415_Diabetic_cardiomyopathy),
                                           fromType="ENTREZID",
                                           toType = c("SYMBOL","ENSEMBL"),
                                           OrgDb = org.Mm.eg.db,
                                           drop = TRUE)

write.csv(DiabetesCardio, "./RData/data/ListPathways/Diabetic_cardiomyopathy.csv")

Alzheimer <- read.csv("./RData/data/ListPathways/Alzheimer_disease.csv",row.names = 1)
DiabetesCardio <- read.csv("./RData/data/ListPathways/Diabetic_cardiomyopathy.csv", row.names = 1)

#TO DO:
#Get ENTREZ ID for the whole expression matrix

#TO DO: 
#Extract the Heatmap function of BioFig... and get row order

#TO DO:
#Matrix should actually match the heatmap row - order


#TO DO:
#add description parameters roxygen style
# SEnrich <- function(matrix, PathwayL, desc = TRUE, thres_quant=0.05, ...){
#   #####
#   #Idea TO DO Function Improvement:
#   #Input as list... alternatively we could provide a df with multiple ids... 
#   #Create some function that automatically detects the type of ID... 
#   #loops first entry of all columns and finds the ID that matches the format of the input matrix, 
#   #if not found throw error to inform the user
#   
#   #Remove NAs if there are any in the list
#   PathwayL <- PathwayL[!is.na(PathwayL)]
#   
#   ResultsRows <- list() #store results later
#   for(rw in 1: ncol(matrix)){
#     #Get Row names
#     rwIterColnames <- colnames(matrix[rw,-rw,drop=FALSE])
#     #Get values row
#     rwIter <- matrix[rw,-rw] #Remove distance to itself
#     #Name the vector
#     names(rwIter) <- rwIterColnames
#     #Sort descending
#     rwIter <- sort(rwIter,decreasing = desc)
#     
#     #Threshold
#     threshold <- unname(quantile(rwIter,thres_quant))
#     
#     #Before threshold
#     FeaturesBThres <- names(rwIter[rwIter < threshold])
#     #After threshold
#     FeaturesAThres <- names(rwIter[rwIter >= threshold])
#     
#     #Construct table for Chisquare test
#     #Features Before Threshold part of the Pathway of interest
#     FeaturesBT_InP <- length(FeaturesBThres[FeaturesBThres %in% PathwayL])
#     #Features Before Threshold not part of the Pathway of interest
#     FeaturesBT_NotInP <- length(FeaturesBThres) - FeaturesBT_InP
#     #Features After Threshold part of the Pathway of interest
#     FeaturesAT_InP <- length(FeaturesAThres[FeaturesAThres %in% PathwayL])
#     #Features After Threshold not part of the Pathway of interest
#     FeaturesAT_NotInP <- length(FeaturesAThres) - FeaturesAT_InP
#     #Construct table
#     TMatrix <- matrix(c(FeaturesBT_InP,FeaturesBT_NotInP,
#                          FeaturesAT_InP,FeaturesAT_NotInP),
#                        ncol = 2, nrow = 2)
#     rownames(TMatrix) <- c("InP","NotInP")
#     colnames(TMatrix) <- c("InL","NotInL")
#     
#     #Fisher's exact test
#     FResultsRow <- fisher.test(TMatrix, alternative = "greater")
#     FResultsRowT <- glance(FResultsRow)
#     FResultsRowT <- FResultsRowT %>% dplyr::rename(oddsRatio = estimate)
#     FResultsRowT$FeatureID <- colnames(matrix)[rw]
#     ResultsRows[[rw]] <- FResultsRowT
#     # phyper(FeaturesBT_InP,length(PathwayL),(length(rwIter)-length(PathwayL)),length(FeaturesBThres),lower.tail = FALSE)
#   }
#   ResultEnrichment <- do.call(rbind,ResultsRows)
#   ResultEnrichment$p.adj <- p.adjust(ResultEnrichment$p.value, ...)
#   ResultEnrichment <- ResultEnrichment %>% relocate(FeatureID, oddsRatio, p.value, p.adj)
#   
#   #TO DO:
#   #When creating density plots it should match row order of the heatmap
#   
#   return(ResultEnrichment)
# }



#To solve... What Ids to use as default?
#For now as the data we are using has Ensembl we will use those. Probably entrez ids are better... 

undebug(SEnrich)
Enrich_Alzheimer <- SEnrich(distanceMatrixBAT_WAT, thres_quant = 0.01, Alzheimer$ENSEMBL,desc = F)


#Work in progress
# Enrich_Alzheimer <- Enrich_Alzheimer %>% mutate(OneMpadj = 1 - p.adj)
# Enrich_Alzheimer <- Enrich_Alzheimer %>% arrange(p.adj)
# 
# resl<-list()
# increase<-10
# previous_padj<-0
# a<-0
# for(i in 1:nrow(Enrich_Alzheimer)){
#   if(i>1){
#     if(previous_padj<0.05){
#       increase <- increase * 2
#     }else{
#       increase <- 1
#     }
#   }
#   
#   if(Enrich_Alzheimer$p.adj[i]<0.05){
#     increase <- 10
#     resl[[i]] <- a + increase
#     a <- a + increase
#     previous_padj <- Enrich_Alzheimer$p.adj[i]
#   }else{
#     resl[[i]] <- 0
#     a <- 0
#     previous_padj <- Enrich_Alzheimer$p.adj[i]
#   }
# }
# resl<-unlist(resl)
# 
# Enrich_Alzheimer$N <- resl
# Enrich_Alzheimer$N2 <- 1:nrow(Enrich_Alzheimer)
# 
# ggplot(data = Enrich_Alzheimer,aes(x=N2, y = N))+
#   stat_smooth(aes(x = N2, y = N), method = "lm",
#               formula = y ~ poly(x, 21), se = FALSE)
# 
# 
# ggplot(d) + 
#   geom_point(aes(x = hour, y = impressions, colour = cvr), size = 3) +
#   stat_smooth(aes(x = hour, y = impressions), method = "lm",
#               formula = y ~ poly(x, 21), se = FALSE) +
#   coord_cartesian(ylim = c(0, 1.5e7))


###########
#Option 2 GSEA

###########
#Option 3

###########
#Explore GSEA function
library(clusterProfiler)

clusterProfiler::GSEA
clusterProfiler:::GSEA_internal
clusterProfiler:::.GSEA
file.path(find.package("GSEA"))
find.package("clusterProfiler")
lazyLoad(filebase = "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/clusterProfiler/R/clusterProfiler", envir = parent.frame())

###############################
## Distance matrix as in WGCNA
library(dplyr)
library(WGCNA)
library(igraph)
## Calculate distance for the gene expression and convert it to weighted adjacency matrix using
## make_adjmatrix_graph() function or  adjacency() function in igraph or WGCNA package
gene_dist_expr_matrix_Wadjacency <- WGCNA::adjacency(t(expr_matrix),type = "distance") %>% as.matrix()
#gene_dist_expr_matrix_Wadjacency <- igraph::make_adjmatrix_graph(t(expr_matrix)) 
class(gene_dist_expr_matrix_Wadjacency)

adjacencyM<-gene_dist_expr_matrix_Wadjacency



#########
#Create a network

library('igraph')
#Function to create network from adjacency matrix
export_network_to_graphml <- function (adj_mat, NetworkFileName=NULL, weighted=TRUE,
                                       threshold=0, upperthreshold = NULL,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE) {
  
  
  # Determine filename to use
  if (is.null(NetworkFileName)) {
    NetworkFileName<-'network.graphml'
  }
  
  # Remove edges with weights lower than the cutoff
  adj_mat[abs(adj_mat) < threshold] <- 0
  
  #I added this part to the function to create a network based on a range
  #Create network from a specific range e.g. 0.5 - 0.8
  if (!is.null(upperthreshold)) {
    adj_mat[abs(adj_mat) > upperthreshold] <- 0
  }
  
  
  if (verbose) {
    message(sprintf("Outputting matrix with %d nodes and %d edges", 
                    nrow(adj_mat), sum(adj_mat > 0)-(ncol(adj_mat))))
    #message(sprintf("The final threshold was %f",threshold))
  }
  
  # Create a new graph and add vertices
  # Weighted graph
  if (weighted) {
    g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
  } else {
    adj_mat[adj_mat != 0] <- 1
    g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
  }
  
  # Add single node annotation from vector
  if (!is.null(nodeAttr)) {
    g <- set.vertex.attribute(g, "attr", value=nodeAttr)
  }
  
  # Add node one or more node annotations from a data frame
  if (!is.null(nodeAttrDataFrame)) {
    for (colname in colnames(nodeAttrDataFrame)) {
      g <- set.vertex.attribute(g, colname, value=as.character(nodeAttrDataFrame[,colname]))
    }
  }

  # edge_correlation_negative <- c()
  # 
  # # neg_correlations[edge_list]
  # edge_list <- get.edgelist(g)
  # 
  # for (i in 1:nrow(edge_list)) {
  #   from <- edge_list[i, 1]    
  #   to   <- edge_list[i, 2]    
  # }
  
  # Save graph to a file
  write.graph(g, NetworkFileName, format='graphml')
  
  # return igraph
  return(g)
}

g <- export_network_to_graphml(adjacencyM, NetworkFileName ='./RData/NetworkFiles/Network_Test.graphml',
                               threshold=0.5, verbose=TRUE )
library(dnet)

Test<-dRWR(g = g,normalise = "none",restart = 0)












































