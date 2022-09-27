#' function to read the .gmt file 
#' @return save the pathway genesets as list of list
readGSet_gmt <- function(file){
  if(!file.exists(file)){
    stop('File',file,' not available\n')
  }
  
  if (!grepl("*.gmt", file)[1]){
    stop("Pathway information must be a .gmt file")
  }
  geneSet = strsplit(readLines(file), "\t")
  names(geneSet) = sapply(geneSet, "[", 1)
  geneSet = lapply(geneSet, "[", -c(1:2))
  geneSet = lapply(geneSet, function(data) {
    data[which(! data == "")]
  })
  return(geneSet)
  
}
