#' Counting the number of reads mapping on one genomic interval.
#'
#' @param matrix Epialleles binary matrix.
#' @return An integer indicating the coverage observed in the genomic region.
#' @export

num_reads= function(matrix){
  num_reads=nrow(matrix)
  return(num_reads)
}


#' Counting the number of reads mapping on one genomic interval.
#'
#' @param list_epi_cpg List containing the epiallele dataframe and the covered CpG sites in a aspecific genomic range.
#' @return An integer indicating the coverage observed in the genomic region.
#' @export

num_reads_mhap= function(list_epi_cpg){
  epi_mat = list_epi_cpg$epi
  num_reads=sum(epi_mat$nreads, na.rm = TRUE)
  return(num_reads)
}
