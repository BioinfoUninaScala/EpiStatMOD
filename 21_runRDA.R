#' Running principal component analysis for one genomic interval
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr filter select full_join
#' @importFrom vegan rda
#' @param samples A list object containing the epiallele composition matrices from all the samples of the dataset.
#' @param region A character indicating the regionID wanted to perform the analysis.
#' @param printData Logical indicating whether the epiallele matrix should be printed in the standard output or not (Default = FALSE).
#' @param rmUnmeth Logical indicating if 0 methylated epialleles should be removed or not from the analysis.
#' @return Principal Component analysis results and plot.
#' @export

runRDA <- function(samples, region, printData = FALSE, rmUnmeth = FALSE){
  data <- getEpimatrix(samples, region)
  if(rmUnmeth == TRUE){
    data <- data[grep("1", colnames(data))]
    data = data %>% dplyr::filter(!rowSums(.) == 0)
  }
  if(length(colnames(data)) <= 1){
    print("Plotting is not possible with just one epiallele specie")
  } else {
    myrda <- vegan::rda(data)
    p <- stats::biplot(myrda, scaling = "symmetric")
    if(printData == TRUE){
      print(data)
    }
    return(p)
  }
}
