#' Calculating the average DNA methylation for one interval
#'
#' @param matrix Epialleles binary matrix.
#' @return A numeric indicating the average DNA methylation of one genomic region
#' @export

meanMeth <- function(matrix){
  mean_met=round(sum(as.matrix(matrix),na.rm=T)/(dim(matrix)[1]*dim(matrix)[2]),2)
  return(mean_met)
}


#' Calculating the average DNA methylation for a not binary epialleles matrix.
#'
#' @param matrix Epialleles matrix.
#' @param param A list of two numeric values indicating the codes for methylated and unmethylated status in the epialleles matrix.
#' @return A numeric indicating the average DNA methylation of one genomic region
#' @export

meanMeth_mod <- function(matrix, param = list(meth_code = 1L, unmeth_code = 0L)){
  meth_code = param[[1]]
  unmeth_code = param[[2]]
  new_mat <-as.matrix(matrix)
  methylated_sites <- sum(new_mat == meth_code,na.rm=T)
  total_sites <- sum(new_mat %in% c(meth_code, unmeth_code),na.rm=T)
  mean_met=round(methylated_sites/total_sites,2)
  return(mean_met)
}



#' Calculating the average DNA methylation for a not binary epialleles matrix.
#'
#' @param list_epi_cpg List containing the epiallele dataframe and the covered CpG sites in a aspecific genomic range.
#' @param param A list of two numeric values indicating the codes for methylated and unmethylated status in the epialleles matrix.
#' @return A numeric indicating the average DNA methylation of one genomic region
#' @export
 

meanMeth_mhap <- function(
    list_epi_cpg,
    param = list(meth_code = 1L, unmeth_code = 0L)
) {
  
  meth_code <- param[[1]]
  unmeth_code <- param[[2]]
  
  epi_dt <- data.table::as.data.table(list_epi_cpg$epi)
  
  if (!"epi" %in% colnames(epi_dt)) {
    stop("`list_epi_cpg$epi` must contain an `epi` column.")
  }
  
  count_col <- NULL
  
  if ("Freq" %in% colnames(epi_dt)) {
    count_col <- "Freq"
  } else if ("nreads" %in% colnames(epi_dt)) {
    count_col <- "nreads"
  } else {
    stop("`list_epi_cpg$epi` must contain either `Freq` or `nreads`.")
  }
  
  weights <- as.numeric(epi_dt[[count_col]])
  patterns <- as.character(epi_dt$epi)
  
  mat_list <- strsplit(patterns, "", fixed = TRUE)
  
  n_cpg <- unique(lengths(mat_list))
  
  if (length(n_cpg) != 1L) {
    stop("All epialleles must have the same number of CpGs.")
  }
  
  mat <- matrix(
    as.integer(unlist(mat_list, use.names = FALSE)),
    nrow = length(mat_list),
    ncol = as.integer(n_cpg),
    byrow = TRUE
  )
  
  meth_mat <- mat == meth_code
  valid_mat <- (mat == meth_code) | (mat == unmeth_code)
  
  methylated_sites <- sum(weights * rowSums(meth_mat), na.rm = TRUE)
  total_sites <- sum(weights * rowSums(valid_mat), na.rm = TRUE)
  
  if (total_sites == 0) {
    return(NA_real_)
  }
  
  mean_met <- round(methylated_sites / total_sites, 2)
  
  return(mean_met)
}