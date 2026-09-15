#' Counting the number of CpG sites contained in one interval
#'
#' @param matrix Epialleles binary matrix.
#' @return An integer indicating the number of CpGs displayed for one genomic interval.
#' @export

numCG <- function(matrix){
  num_cg=ncol(matrix)
  return(num_cg)
}


#' Counting the number of CpG sites contained in one interval
#'
#' @param list_epi_cpg List containing the epiallele dataframe and the covered CpG sites in a aspecific genomic range.
#' @return An integer indicating the number of CpGs displayed for one genomic interval.
#' @export

numCG_mhap <- function(list_epi_cpg){
  epi_mat = list_epi_cpg$epi
  num_cg = unique(epi_mat$n_cpg)[1]
  return(num_cg)
}