#' Visualizing the ordination plot for a selected region
#'
#' @description
#' Function used to display the ordination plot of samples relative to specific genomic regions.
#' @importFrom purrr map reduce
#' @importFrom dplyr filter select full_join
#' @importFrom vegan cca ordiellipse
#' @param samples A list object containing the epiallele composition matrices from all the samples of the dataset.
#' @param region A string containing the regionID wanted to perform the analysis.
#' @param metadata A dataframe object containing samples metadata. Dataframe should contain dedicated columns for samples IDs and the one indicating the group they belong to.
#' @param printData Logical indicating whether the epiallele matrix should be printed in the standard output or not.
#' @param rmUnmeth Logical indicating whether unmethylated epialleles should be discarded from the analysis.
#' @return A list object storing the ordination results and an ordination plot.
#' @export

runCCA <- function(samples, region, metadata, printData = FALSE, rmUnmeth = FALSE){
  
  data <- getEpimatrix(samples, region)
  
  if(rmUnmeth == TRUE){
    data <- data[, grep("1", colnames(data)), drop = FALSE]
    data <- data %>% dplyr::filter(!rowSums(.) == 0)
  }
  
  if(ncol(data) <= 1){
    
    print("Plotting is not possible with just one epiallele specie")
    return(NULL)
    
  } else {
    
    metadata <- metadata %>%
      dplyr::filter(Samples %in% rownames(data)) %>%
      dplyr::mutate(Group = factor(Group))
    
    metadata <- metadata[match(rownames(data), metadata$Samples), ]
    
    if(any(is.na(metadata$Group))){
      warning("Some samples have NA Group after metadata matching")
    }
    
    mod <- vegan::cca(data ~ Group, data = metadata)
    
    plot(
      mod,
      type = "n",
      display = "sites",
      main = region
    )
    
    site_scores <- vegan::scores(mod, display = "sites", choices = 1:2)
    
    graphics::text(
      x = site_scores[, 1],
      y = site_scores[, 2],
      labels = as.character(metadata$Group),
      col = as.numeric(metadata$Group)
    )
    
    pl <- vegan::ordiellipse(
      mod,
      metadata$Group,
      kind = "se",
      conf = 0.95,
      lwd = 2,
      draw = "polygon",
      col = seq_along(levels(metadata$Group)),
      border = seq_along(levels(metadata$Group)),
      alpha = 63
    )
    
    if(printData == TRUE){
      print(data)
    }
    
    return(pl)
  }
}