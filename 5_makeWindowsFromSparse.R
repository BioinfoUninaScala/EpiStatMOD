#' Designing fixed-width windows from a sparse CpG matrix
#'
#' @description
#' This function builds fixed-width genomic windows directly from the target base
#' positions encoded in the columns of a sparse matrix (as returned by
#' \code{loadInputFromMod() or loadInputFromMhap()}). Each column name is expected to be of the form
#' \code{"chr_pos"} (e.g. \code{"chr21_5011700"}). For each CpG position,
#' two candidate windows are created:
#' \itemize{
#'   \item a left-anchored window: \code{[pos - window + 1, pos]}
#'   \item a right-anchored window: \code{[pos, pos + window - 1]}
#' }
#' All windows are then filtered according to the number of CpGs they contain
#' (between \code{min.C} and \code{max.C}), and only windows fully contained
#' in the user-provided \code{cov_regs} regions are retained.
#'
#' Optionally, the number of reads overlapping each window can be added
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
#' @param window Integer scalar indicating the window size (in base pairs)
#'   used to construct left- and right-anchored windows around each CpG
#'   position. Default is \code{50}.
#' @param min.C Integer scalar indicating the minimum number of CpGs
#'   required for a window to be retained. Default is \code{2}.
#' @param max.C Integer scalar indicating the maximum number of CpGs
#'   allowed in a window. Windows with more than \code{max.C} CpGs are
#'   discarded. Default is \code{50}.
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
#'   win_gr <- makeWindowsFromSparse(
#'     sparseMat    = sparseMat,
#'     cov_regs     = cov_regs,
#'     read_coord_gr = reads_gr,
#'     window       = 50L,
#'     min.C        = 2L,
#'     max.C        = 50L
#'   )
#' }
#'


makeWindowsFromSparse <- function(sparseMat,
                                  cov_regs,
                                  read_coord_gr = NULL,
                                  window = 50L,
                                  min.C  = 2L,
                                  max.C  = 50L) {
  stopifnot(!is.null(colnames(sparseMat)))
  if (!methods::is(cov_regs, "GRanges")) {
    stop("`cov_regs` must be a GRanges object.")
  }
  
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
  
  cpg_gr <- GenomicRanges::GRanges(
    seqnames = dt$chr,
    ranges   = IRanges::IRanges(start = dt$pos, width = 1L)
  )
  
  chr_levels <- unique(dt$chr)
  win_list   <- vector("list", length(chr_levels))
  names(win_list) <- chr_levels
  
  for (chr_name in chr_levels) {
    dt_chr <- dt[chr == chr_name]
    pos    <- dt_chr$pos
    if (length(pos) == 0L) next
    
    left_start  <- pmax(pos - window + 1L, 1L)
    left_end    <- pos
    right_start <- pos
    right_end   <- pos + window - 1L
    
    starts <- c(left_start, right_start)
    ends   <- c(left_end,   right_end)
    
    seq_levels <- unique(dt$chr)
    gr_chr <- GenomicRanges::GRanges(
      seqnames = chr_name,
      ranges   = IRanges::IRanges(start = starts, end = ends),
      seqinfo  = GenomeInfoDb::Seqinfo(seqnames = seq_levels)
    )
    
    gr_chr <- BiocGenerics::unique(gr_chr)
    
    win_list[[chr_name]] <- gr_chr
  }
  
  win_list <- win_list[!vapply(win_list, is.null, logical(1))]
  if (length(win_list) == 0L) {
    message("No windows generated from sparse matrix.")
    return(GenomicRanges::GRanges())
  }
  
  all_win <- do.call(c, base::unname(win_list))
  
  cpg_per_win <- GenomicRanges::countOverlaps(all_win, cpg_gr)
  keep <- (cpg_per_win >= min.C) & (cpg_per_win <= max.C)
  
  if (!any(keep)) {
    message("No windows with CpG count in [", min.C, ", ", max.C, "].")
    return(GenomicRanges::GRanges())
  }
  
  all_win_filt <- all_win[keep]
  S4Vectors::mcols(all_win_filt)$num_cpg <- cpg_per_win[keep]
  
  if (!is.null(read_coord_gr)) {
    read_per_win <- GenomicRanges::countOverlaps(all_win_filt, read_coord_gr)
    S4Vectors::mcols(all_win_filt)$num_reads <- read_per_win
  }
  
  hits <- GenomicRanges::findOverlaps(all_win_filt, cov_regs, type = "within")
  if (length(hits) == 0L) {
    message("No windows fall completely within `cov_regs`.")
    return(GenomicRanges::GRanges())
  }
  contained_idx <- unique(S4Vectors::queryHits(hits))
  final_windows <- all_win_filt[contained_idx]
  
  final_windows
}
