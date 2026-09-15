#' Calculating the number of different epialleles observed for one genomic region
#'
#' @param matrix Epiallele binary matrix.
#' @return A numeric indicating the number of observed epialleles in one interval.
#' @export

epi <- function(matrix){
  matrix$epi= apply(matrix, 1, function(x) paste(x, collapse = ""))
  matrix= as.data.frame(table(matrix$epi))
  epi=nrow(matrix)
  return(epi)
}


#' Calculating the number of different epialleles observed for one genomic region
#'
#' @param list_epi_cpg List containing the epiallele dataframe and the covered CpG sites in a aspecific genomic range.
#' @return A numeric indicating the number of observed epialleles in one interval.
#' @export

epi_mhap <- function(list_epi_cpg){
  epi_mat = list_epi_cpg$epi
  epi=length(unique(epi_mat$epi))
  return(epi)
}
