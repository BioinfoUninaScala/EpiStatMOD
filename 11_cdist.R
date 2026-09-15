#' Calculating the mean distance among the CpG sites contained in the analysed region
#'
#' @param matrix Epiallele binary matrix for one genomic region
#' @return Integer indicating the mean distance between CpG in the interval
#' @export

cdist <- function(matrix){
  Cpos=as.numeric(names(matrix))
  dist=round(mean(Cpos[-1]-Cpos[-length(Cpos)]),2)
  return(dist)
}


#' Calculating the mean distance among the CpG sites contained in the analysed region from colnames of sparse matrix.
#'
#' @param matrix Epiallele matrix for one genomic region
#' @return Integer indicating the mean distance between CpG in the interval
#' @export

cdist_mod <- function(matrix){
  Cpos=as.numeric(colnames(matrix))
  dist=round(mean(Cpos[-1]-Cpos[-length(Cpos)]),2)
  return(dist)
}


#' Calculating the mean distance among the CpG sites contained in the analysed region from a mhap file.
#'
#' @param list_epi_cpg List containing the epiallele dataframe and the covered CpG sites in a aspecific genomic range.
#' @return Integer indicating the mean distance between CpG in the interval
#' @export

cdist_mhap <- function(list_epi_cpg){
  Cpos=as.numeric(list_epi_cpg$cpg_pos)
  dist=round(mean(Cpos[-1]-Cpos[-length(Cpos)]),2)
  return(dist)
}

