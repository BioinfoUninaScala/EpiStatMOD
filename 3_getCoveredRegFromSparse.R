#' Selecting regions with a minimum user-defined coverage
#'
#' @description
#' The functions allows the selection of genomic regions that satisfy a minimum coverage threshold defined by the user as a parameter.
#' @import dplyr
#' @param sparseMat Sparse matrix generated with loadInputFromModFile().
#' @param threshold Integer indicating the minimum coverage required to keep genomic regions.
#' @param minsize Integer indicating the minimum size of the genomic range. 
#'  This value is also depending from the method of choice that will be used by 
#'  the dedicated functions makeBins() and makeWindows(). 
#'  If you use makeBins() the minimum size will be equal to the length of n(CpG sites to be included in the target region)x2.
#'  Instead, if you use makeWindows, the minimum size of the genomic regions 
#'  used to build target regions will be equal to the size of the sliding window that will be used next.
#' @return filtered GRanges object.
#' @export


getCoveredRegFromSparse <- function(sparseMat,
                          threshold = 10,
                          minsize = 50
                          )
{
  cat("Filtering positions with minimum coverage", threshold, "\n")
  
  # 1) Classifying sites based on coverage threshold
  cov_counts <- Matrix::colSums(sparseMat != 0)
  covered_vec <- cov_counts >= threshold
  
  # 2) Extract column names
  col_keys <- colnames(sparseMat)
  if (is.null(col_keys)) {
    stop("colnames(sparseMat) are NULL - cannot derive coordinates.")
  }
  
  parts <- data.table::tstrsplit(col_keys, "_", fixed = TRUE)
  if (length(parts) != 2L) {
    stop("Column names must be of the form 'chr_pos', e.g. 'chr21_50808277'.")
  }
  
  chr_vec <- parts[[1]]
  pos_vec <- as.integer(parts[[2]])
  
  if (any(is.na(pos_vec))) {
    stop("Some positions could not be converted to integer. Check column names.")
  }
  
  dt <- data.table::data.table(
    idx     = seq_along(col_keys),
    chr     = chr_vec,
    pos     = pos_vec,
    covered = covered_vec
  )
  
  out_list <- list()
  
  # 3) Merging sites per chromosome
  for (chr_name in unique(dt$chr)) {
    dt_chr <- dt[chr == chr_name]
    
    data.table::setorder(dt_chr, pos)
    
    v      <- dt_chr$covered
    coords <- dt_chr$pos
    
    if (length(v) == 0L) next
    
    r <- rle(v)
    run_end   <- cumsum(r$lengths)
    run_start <- run_end - r$lengths + 1
    
    true_runs <- which(r$values)
    if (length(true_runs) == 0L) next
    
    true_start_indexes <- run_start[true_runs]
    true_end_indexes <- run_end[true_runs]
    
    chr_df <- data.frame(
      chr   = chr_name,
      start = coords[true_start_indexes],
      end   = coords[true_end_indexes],
      stringsAsFactors = FALSE
    )
    
    out_list[[chr_name]] <- chr_df
  }
  
  if (length(out_list) == 0L) {
    cat("Any site has a coverage >=", threshold, "\n")
    return(
      data.frame(chr = character(),
                 start = integer(),
                 end = integer(),
                 stringsAsFactors = FALSE)
    )
  }
  
  bed <- do.call(rbind, out_list)
  
  bed <- bed[order(bed$chr, bed$start), ]
  rownames(bed) <- NULL
  
  cat("Number of covered regions:", nrow(bed), "\n")
  
  filt_bed <- bed %>% filter(end - start >= minsize)
    
  cat("Number of covered regions after filtering for minsize =", minsize, ':', nrow(filt_bed), "\n")
  
  gr_bed <- GenomicRanges::makeGRangesFromDataFrame(filt_bed, keep.extra.columns = TRUE)
  cat("done\n")
  
  return(gr_bed)
}

