#' Finding the epiallele species which are contributing most to the dissimilarity among groups
#'
#' @description
#' The function returns the information regarding which epialleles are the ones contributing for the dissimilarity of samples among groups.
#' @importFrom purrr map reduce
#' @importFrom dplyr filter select full_join
#' @importFrom vegan simper
#' @param samples A list object containing the epiallele composition matrices from all the samples of the dataset.
#' @param region A string containing the regionID wanted to perform the analysis.
#' @param metadata A dataframe object containing samples metadata. Dataframe should contain dedicated columns for samples IDs and the one indicating the group they belong to.
#' @return A list containing the dissimilarity results.
#' @export



epidiv <- function(samples, region, metadata){
  data <- getEpimatrix(samples, region)
  metadata = metadata %>%
    dplyr::filter(Samples %in% rownames(data))
  sim <- with(metadata, suppressMessages(vegan::simper(data, Group, permutations = 999)))
  return(sim)
}




