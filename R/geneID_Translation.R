## Function for gene ID translation using biomart DB

geneID_TR <- function(df, keys, columns, keytype){
  require(biomaRt)
  require(tidyverse)
  require(EnsDb.Mmusculus.v79)
  require(EnsDb.Hsapiens.v86)
  require(ensembldb)
  require(org.Mm.eg.db)
  require(org.Hs.eg.db)
  require(AnnotationDbi)
  
  for(i in  1:nrow(df)) {  
    if(grepl("ENSMUSG", row.names(df)[i])){
      # get the moues ensembl gene symbols from the biomart
      mart_mmg <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl") 
      expr_data <- biomaRt::getBM(attributes = c("ensembl_gene_id","mgi_symbol"), 
                                  filters = "ensembl_gene_id", values = rownames(df),mart = mart_mmg) %>%
        dplyr::distinct(mgi_symbol, .keep_all = T) %>% tidyr::drop_na()  %>%
        left_join(tibble::rownames_to_column(df), by = c("ensembl_gene_id" = "rowname")) %>%
        dplyr::select(-ensembl_gene_id) %>%
        tibble::column_to_rownames(var = "mgi_symbol")
      return(expr_data)
      
    } 
    if(grepl("ENSG", row.names(df)[i])){ 
      # to convert from ensembl.gene to gene.symbol for hsg 
      mart_hsg <- biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
      expr_data <- biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), 
                                  filters = "ensembl_gene_id", values = rownames(df),mart = mart_hsg) %>%
        dplyr::distinct(mgi_symbol, .keep_all = T) %>% tidyr::drop_na()  %>%
        left_join(tibble::rownames_to_column(df), by = c("ensembl_gene_id" = "rowname")) %>%
        dplyr::select(-ensembl_gene_id) %>%
        tibble::column_to_rownames(var = "hgnc_symbol")
      return(expr_data)
    }
   
  if(!grepl("ENSMUSG", row.names(df)[i]) & is.character(row.names(df))){
    # to convert from moues gene symbol to ensembl gene ID 
    mart_mmg <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl") 
    expr_data <- biomaRt::getBM(attributes = c("ensembl_gene_id","mgi_symbol"), 
                                filters = "mgi_symbol", values = rownames(df),mart = mart_mmg) %>%
      left_join(tibble::rownames_to_column(df), by = c("mgi_symbol" = "rowname")) %>%
      dplyr::select(-mgi_symbol) %>%
      dplyr::distinct(ensembl_gene_id, .keep_all = T) %>% tidyr::drop_na()  %>%
      tibble::column_to_rownames(var = "ensembl_gene_id")
    return(expr_data)
}

  if(!grepl("ENSG", row.names(df)[i]) & is.character(row.names(df))){ 
    # to convert human gene symbol to  ensembl gene ID 
    mart_hsg <- biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    expr_data <- biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), 
                                filters = "hgnc_symbol", values = rownames(df),mart = mart_hsg) %>%
      left_join(tibble::rownames_to_column(df), by = c("hgnc_symbol" = "rowname")) %>%
      dplyr::select(-hgnc_symbol) %>%
      dplyr::distinct(ensembl_gene_id, .keep_all = T) %>% tidyr::drop_na()  %>%
      tibble::column_to_rownames(var = "ensembl_gene_id")
    return(expr_data)
  }
}
}


