#' Designing fixed-width windows from a sparse CpG matrix
#'
#' @description
#' This function builds flexible-width genomic bins directly from the target base
#' positions encoded in the columns of a sparse matrix (as returned by
#' \code{loadInputFromMod() or loadInputFromMhap()}). Each column name is expected to be of the form
#' \code{"chr_pos"} (e.g. \code{"chr21_5011700"}). Candidate bins are created by
#'  generating not-redundant genomic ranges that contain n target base positions.
#' All bins are then filtered according to the minimum and maximum bin size
#' (between \code{min.binsize} and \code{max.binsize}).
#'
#' Optionally, the number of reads overlapping each bin can be added
#' using a \code{GRanges} object with read alignments.
#'
#' @importFrom data.table tstrsplit data.table setorder
#' @importFrom GenomicRanges GRanges countOverlaps findOverlaps
#' @importFrom IRanges IRanges
#'
#' @param sparseMat A sparse matrix (typically a \code{Matrix::dgCMatrix})
#'   with rows corresponding to reads and columns corresponding to CpG
#'   positions. Column names must be of the form \code{"chr_pos"},
#'   for example \code{"chr21_5011700"}.
#' @param cov_regs A \code{GRanges} object specifying the genomic regions
#'   where windows are allowed. Only windows entirely contained within
#'   these regions are kept.
#' @param read_coord_gr Optional \code{GRanges} object describing read
#'   genomic coordinates (e.g. from a BED/BAM conversion). When provided,
#'   the number of reads overlapping each window is computed and stored
#'   in the \code{num_reads} metadata column.
#' @param n An integer indicating the number of CpGs sites to be contained in the target region.
#' @param min.binsize An integer indicating the minimum length size of each bin
#' @param max.binsize An integer indicating the maximum length size of each bin
#' 
#' @return A \code{GRanges} object containing the selected windows. Each
#'   range has at least the metadata column \code{num_cpg}, indicating
#'   the number of CpGs within the window. If \code{read_coord_gr} is
#'   provided, an additional column \code{num_reads} is added with the
#'   number of overlapping reads.
#' @export
#' @examples
#' \dontrun{
#'   # sparseMat: dgCMatrix with colnames like "chr21_5011700"
#'   # cov_regs: GRanges with high-coverage regions
#'   # reads_gr: GRanges with read alignments
#'
#'   bin_gr <- makeBinsFromSparse(
#'     sparseMat    = sparseMat,
#'     cov_regs     = cov_regs,
#'     read_coord_gr = reads_gr,
#'     n = 4,
#'     min.binsize  = 50L,
#'     max.binsize  = 150L
#'   )
#' }
#'


makeBinsFromSparse <- function(sparseMat,
                                  cov_regs,
                                  read_coord_gr = NULL,
                                  n = 4,
                                  min.binsize  = 50L,
                                  max.binsize  = 150L) {
  stopifnot(!is.null(colnames(sparseMat)))
  if (!methods::is(cov_regs, "GRanges")) {
    stop("`cov_regs` must be a GRanges object.")
  }
  
  # ----- 1) Estrai chr e pos dalle colonne -----
  parts <- data.table::tstrsplit(colnames(sparseMat), "_", fixed = TRUE)
  if (length(parts) != 2L) {
    stop("Column names must be of the form 'chr_pos', e.g. 'chr21_5011700'.")
  }
  chr_vec <- parts[[1]]
  pos_vec <- as.integer(parts[[2]])
  
  if (any(is.na(pos_vec))) {
    stop("Some positions could not be converted to integer. Check column names.")
  }
  
  dt <- data.table::data.table(chr = chr_vec, pos = pos_vec)
  data.table::setorder(dt, chr, pos)
  
  # ----- 2) Costruisci finestre fisse per cromosoma -----
  bins_2cpg <- makeCpgCountBins(
    cpg_df = dt,
    n = n,
    min.binsize  = min.binsize,
    max.binsize  = max.binsize
  )
  
  if (!is.null(read_coord_gr)) {
    read_per_bin <- GenomicRanges::countOverlaps(bins_2cpg, read_coord_gr)
    S4Vectors::mcols(bins_2cpg)$num_reads <- read_per_bin
  }
  
  # ----- 4) Tieni solo finestre dentro cov_regs -----
  hits <- GenomicRanges::findOverlaps(bins_2cpg, cov_regs, type = "within")
  if (length(hits) == 0L) {
    message("No windows fall completely within `cov_regs`.")
    return(GenomicRanges::GRanges())
  }
  contained_idx <- unique(S4Vectors::queryHits(hits))
  final_bins <- bins_2cpg[contained_idx]
  
  return(final_bins)
}


makeCpgCountBins <- function(cpg_df, n = 2L, step = 1L, min.binsize = 50L, max.binsize = 150L) {
  
  data.table::setDT(cpg_df)
  
  if (!all(c("chr", "pos") %in% colnames(cpg_df))) {
    stop("`cpg_df` must contain columns `chr` and `pos`.")
  }
  
  if (n <= 0L) {
    stop("`n` must be > 0.")
  }
  
  if (step <= 0L) {
    stop("`step` must be > 0.")
  }
  
  cpg_df <- cpg_df[
    !is.na(chr) &
      !is.na(pos)
  ]
  
  cpg_df <- unique(cpg_df)
  
  data.table::setorder(cpg_df, chr, pos)
  
  out_list <- cpg_df[
    ,
    {
      if (.N < n) {
        NULL
      } else {
        
        starts_idx <- seq(
          from = 1L,
          to = .N - n + 1L,
          by = step
        )
        
        starts <- pos[starts_idx]
        ends <- pos[starts_idx + n - 1L]
        
        data.table::data.table(
          start = starts,
          end = ends,
          num_cpg = n,
          len = ends - starts
        )
      }
    },
    by = chr
  ]
  
  out_list <- out_list[
    len >= min.binsize &
      len <= max.binsize
  ]
  
  if (nrow(out_list) == 0L) {
    return(GenomicRanges::GRanges())
  }
  
  GenomicRanges::makeGRangesFromDataFrame(
    as.data.frame(out_list),
    keep.extra.columns = TRUE
  )
}
