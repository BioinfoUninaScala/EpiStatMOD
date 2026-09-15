#' Calculating the number of epialleles with just one observation
#'
#' @param matrix Epiallele binary matrix
#' @return An integer indicating the number of epiallele with a unique observation in a given interval
#' @export

singleton <- function(matrix){
  matrix$epi= apply(matrix, 1, function(x) paste(x, collapse = ""))
  matrix= as.data.frame(table(matrix$epi))
  singleton= length(matrix$Freq[matrix$Freq==1])
  return(singleton)
}


#' Calculating the number of epialleles with just one observation
#'
#' @param list_epi_cpg List containing the epiallele dataframe and the covered CpG sites in a aspecific genomic range.
#' @return An integer indicating the number of epiallele with a unique observation in a given interval
#' @export

singleton_mhap <- function(list_epi_cpg){
  epi_mat = list_epi_cpg$epi
  singleton= length(epi_mat$nreads [epi_mat$nreads==1])
  return(singleton)
}
