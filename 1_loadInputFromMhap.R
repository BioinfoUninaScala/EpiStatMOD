#' Load mHap input and optionally build a read-by-CpG sparse matrix
#'
#' `loadInputFromMhap()` loads methylation haplotypes in `.mhap` format and
#' prepares them for epiallele analysis. The function can either return a
#' lightweight mHap-based object, build a sparse read-by-CpG matrix compatible
#' with the sparse backend, or return both objects.
#'
#' The function uses two CpG coordinate sets internally:
#'
#' \itemize{
#'   \item `cpg_df_all`: all CpG positions used to correctly map each character
#'   of the mHap epiallele string to its genomic CpG coordinate.
#'   \item `cpg_df_target`: the CpG positions retained as columns of the sparse
#'   matrix. If `region_file = NULL`, this corresponds to all CpGs. If
#'   `region_file` is provided, only CpGs falling inside the selected regions
#'   are retained as matrix columns.
#' }
#'
#' If `region_file` is provided, mHap records are first filtered to retain only
#' records overlapping at least one selected genomic region. CpG columns are then
#' restricted to CpGs falling inside the same regions. Importantly, the complete
#' mHap pattern is still interpreted using the full CpG coordinate set, so that
#' partial regional selection preserves the correct correspondence between
#' epiallele characters and CpG positions.
#'
#' @param mhap_file Character string. Path to the `.mhap` or `.mhap.gz` file.
#' The file is expected to contain at least six columns:
#' chromosome, start, end, epiallele pattern, number of reads, and strand.
#'
#' @param cpg_file Optional character string. Path to a CpG coordinate file.
#' The file should contain at least two columns: chromosome and CpG position.
#' If `NULL`, CpG coordinates are extracted from the genome using
#' `genome_fasta` or `genome_assembly`.
#'
#' @param genome_fasta Optional character string. Path to a genome FASTA file.
#' Used only when `cpg_file = NULL`.
#'
#' @param genome_assembly Character string or genome object. Genome assembly
#' used to retrieve CpG coordinates when `cpg_file = NULL`. This can be a
#' supported assembly name, such as `"hg38"`, or an object accepted by
#' `load_genome()`.
#'
#' @param mhap_file_1based Logical. Whether coordinates in `mhap_file` are
#' already 1-based inclusive. If `FALSE`, the mHap start coordinate is converted
#' from BED-like 0-based coordinates by adding 1. Default is `TRUE`.
#'
#' @param cpg_file_1based Logical. Whether CpG positions in `cpg_file` are
#' already 1-based. If `FALSE`, 1 is added to the CpG position column.
#' Default is `TRUE`.
#'
#' @param region_file Optional character string. Path to a BED-like file
#' containing genomic regions to restrict the analysis. The file must contain
#' at least three columns: chromosome, start, and end. Additional columns are
#' ignored. If `NULL`, all CpGs and mHap records on the selected chromosomes are
#' retained.
#'
#' @param region_file_1based Logical. Whether the region file already uses
#' 1-based inclusive coordinates. If `FALSE`, regions are interpreted as
#' BED-like 0-based half-open intervals and the start coordinate is converted by
#' adding 1. Default is `FALSE`.
#'
#' @param sel_chr Optional character vector specifying chromosomes to retain,
#' for example `"chr22"` or `c("chr21", "chr22")`. If `NULL`, all chromosomes
#' present in the input are processed.
#'
#' @param merge_regions Logical. Whether overlapping or adjacent regions in
#' `region_file` should be merged before filtering mHap records and CpG
#' positions. Default is `TRUE`.
#'
#' @param bases Character vector specifying the motif used to define target
#' bases when CpG coordinates are extracted from the genome. Default is `"CG"`.
#' Used only when `cpg_file = NULL`.
#'
#' @param code_map Named integer vector mapping mHap pattern characters to
#' integer values stored in the sparse matrix. The default maps `"0"` to `2`
#' and `"1"` to `3`, with `default = 1L` for unexpected characters.
#'
#' @param rm_zero_col Logical. Whether to remove columns with no non-zero
#' entries from the sparse matrix. Default is `FALSE`.
#'
#' @param output_type Character string specifying the output type. One of:
#' \itemize{
#'   \item `"mhap"`: return a lightweight mHap object only.
#'   \item `"sparse"`: build and return a read-by-CpG sparse matrix only.
#'   \item `"both"`: return both the mHap object and the sparse matrix.
#' }
#'
#' @param min_cpg_per_read Integer. Minimum number of target CpGs that an mHap
#' record must overlap to be retained for sparse matrix construction. Default
#' is `1`.
#'
#' @param output_file Optional character string. Path to an `.rds` file where
#' the returned object should be saved. If `NULL`, the object is returned only
#' in memory.
#'
#' @param nThread Integer. Number of threads used by `data.table::fread()`.
#' Default is `8`.
#'
#' @param install_if_missing Logical. Passed to `load_genome()`. If `TRUE`,
#' missing genome resources may be installed automatically when supported.
#' Default is `FALSE`.
#'
#' @param parallel Logical. Whether to build the sparse matrix using genomic
#' CpG chunks and parallel workers. Default is `TRUE`.
#'
#' @param n_workers Integer. Number of parallel workers used when
#' `parallel = TRUE`. Default is `4`.
#'
#' @param chunks_per_core Integer. Number of genomic CpG chunks generated per
#' worker. Larger values may improve load balancing but increase scheduling
#' overhead. Default is `10`.
#'
#' @param strand_sign Logical. If `TRUE`, sparse matrix values corresponding to
#' mHap records on the negative strand are multiplied by `-1`. Default is
#' `FALSE`.
#'
#' @param strict Logical. If `TRUE`, the function stops when the length of an
#' mHap epiallele pattern does not match the number of CpGs inferred between
#' the mHap start and end coordinates. If `FALSE`, inconsistent records are
#' skipped with a warning. Default is `TRUE`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{sparse_matrix}{A sparse read-by-CpG matrix of class
#'   `Matrix::dgCMatrix`, or `NULL` when `output_type = "mhap"`. Rows are
#'   artificial read identifiers obtained by expanding the mHap `nreads` field.
#'   Columns are target CpG positions.}
#'
#'   \item{read_coord_gr}{A `GenomicRanges::GRanges` object containing the
#'   genomic coordinates of the artificial reads used as matrix rows, or `NULL`
#'   when no sparse matrix is built. Coordinates correspond to the original mHap
#'   record span.}
#'
#'   \item{mhap_obj}{A lightweight mHap object containing the filtered mHap
#'   table, target CpG coordinates, full CpG coordinates, selected regions, and
#'   backend name. Returned when `output_type` is `"mhap"` or `"both"`.}
#'
#'   \item{mhap_df_used_for_sparse}{The mHap records retained for sparse matrix
#'   construction after target-CpG filtering, or `NULL` when no sparse matrix is
#'   built.}
#' }
#'
#' @details
#' The `.mhap` format stores methylation haplotypes in compressed form by
#' aggregating identical epiallele patterns and recording their read count.
#' When `output_type = "sparse"` or `"both"`, this function expands each mHap
#' record according to its `nreads` value to generate artificial read rows.
#' This mode is intended for compatibility with sparse-matrix workflows such as
#' `makeWindowsFromSparse()` and `epiAnalysisFromSparse()`.
#'
#' For large datasets, expanding mHap counts into read-level rows may still be
#' memory intensive. In those cases, using `output_type = "mhap"` together with
#' a direct mHap backend, such as `makeWindowsFromMhap()` and
#' `epiAnalysisFromMhap()`, is recommended.
#'
#' @examples
#' \dontrun{
#' res <- loadInputFromMhap(
#'   mhap_file = "sample.mhap.gz",
#'   cpg_file = "hg38_CpG.gz",
#'   sel_chr = "chr22",
#'   region_file = "regions.bed",
#'   output_type = "both",
#'   min_cpg_per_read = 2,
#'   parallel = TRUE,
#'   n_workers = 4
#' )
#'
#' sparseMat <- res$sparse_matrix
#' read_coord_gr <- res$read_coord_gr
#' mhap_obj <- res$mhap_obj
#' }
#'
#' @import Matrix
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @import foreach
#' @import doParallel
#' @importFrom data.table fread setDT setnames setorder getDTthreads :=
#' @export
#' 

loadInputFromMhap <- function(
    mhap_file,
    cpg_file = NULL,
    genome_fasta = NULL,
    genome_assembly = "hg38",
    mhap_file_1based = TRUE,
    cpg_file_1based = TRUE,
    region_file = NULL,
    region_file_1based = FALSE,
    sel_chr = NULL,
    merge_regions = TRUE,
    bases = "CG",
    code_map = c(default = 1L, `0` = 2L, `1` = 3L),
    rm_zero_col = FALSE,
    output_type = c("mhap", "sparse", "both"),
    min_cpg_per_read = 1L,
    output_file = NULL,
    nThread = 8,
    install_if_missing = FALSE,
    parallel = TRUE,
    n_workers = 4,
    chunks_per_core = 10,
    strand_sign = FALSE,
    strict = TRUE
) {
  
  output_type <- match.arg(output_type)
  
  message("Loading mHap input...")
  
  mhap_df <- load_mhap(
    mhap_file = mhap_file,
    mhap_file_1based = mhap_file_1based,
    sel_chr = sel_chr,
    nThread = nThread
  )
  
  regions_gr <- NULL
  
  if (!is.null(region_file) && !region_file %in% c("", "NA", "NULL", "null", "None", "none")) {
    
    message("Loading regions: ", region_file)
    
    regions_gr <- load_regions_bed(
      region_file = region_file,
      sel_chr = sel_chr,
      merge_regions = merge_regions,
      region_file_1based = region_file_1based,
      nThread = nThread
    )
    
    message("Filtering mHap records by regions...")
    
    mhap_df <- filter_mhap_by_regions(
      mhap_df = mhap_df,
      regions_gr = regions_gr
    )
    
    if (nrow(mhap_df) == 0L) {
      stop("No mHap records remain after region filtering.")
    }
  }
  
  message("Loading CpG coordinates...")
  
  if (!is.null(cpg_file)) {
    
    cpg_df_all <- load_cpg_file(
      cpg_file = cpg_file,
      cpg_file_1based = cpg_file_1based,
      sel_chr = sel_chr,
      nThread = nThread
    )
    
  } else {
    
    genome <- load_genome(
      genome_assembly = genome_assembly,
      sel_chr = sel_chr,
      genome_fasta = genome_fasta,
      install_if_missing = install_if_missing
    )
    
    cpg_df_all <- get_coordinates(
      genome = genome,
      bases = bases,
      sel_chr = sel_chr,
      regions = NULL,
      nThread = nThread
    )
    
    cpg_df_all <- normalize_cpg_df(cpg_df_all)
  }
  
  if (nrow(cpg_df_all) == 0L) {
    stop("No CpG coordinates available.")
  }
  
  if (!is.null(regions_gr)) {
    
    message("Filtering CpG coordinates by regions...")
    
    cpg_df_target <- filter_cpg_by_regions(
      cpg_df = cpg_df_all,
      regions_gr = regions_gr
    )
    
  } else {
    
    cpg_df_target <- cpg_df_all
  }
  
  if (nrow(cpg_df_target) == 0L) {
    stop("No target CpG coordinates remain after region filtering.")
  }
  
  sparse_matrix <- NULL
  read_coord_gr <- NULL
  mhap_df_used_for_sparse <- NULL
  
  if (output_type %in% c("sparse", "both")) {
    
    message("Building sparse matrix from mHap...")
    
    sparse_res <- build_sparse_matrix_from_mhap_parallel_chunks(
      mhap_df = mhap_df,
      cpg_df_all = cpg_df_all,
      cpg_df_target = cpg_df_target,
      code_map = code_map,
      min_cpg_per_read = min_cpg_per_read,
      rm_zero_col = rm_zero_col,
      strand_sign = strand_sign,
      strict = strict,
      cores = if (parallel) n_workers else 1L,
      chunks_per_core = chunks_per_core
    )
    
    sparse_matrix <- sparse_res$sparse_matrix
    read_coord_gr <- sparse_res$read_coord_gr
    mhap_df_used_for_sparse <- sparse_res$mhap_df_used
  }
  
  mhap_obj <- NULL
  
  if (output_type %in% c("mhap", "both")) {
    
    mhap_obj <- list(
      mhap = mhap_df,
      cpg_df = cpg_df_target,
      cpg_df_all = cpg_df_all,
      regions = regions_gr,
      backend = "mhap"
    )
  }
  
  res <- list(
    sparse_matrix = sparse_matrix,
    read_coord_gr = read_coord_gr,
    mhap_obj = mhap_obj #,
    # mhap_df_used_for_sparse = mhap_df_used_for_sparse
  )
  
  if (!is.null(output_file)) {
    message("Saving output object to: ", output_file)
    saveRDS(res, output_file)
  }
  
  message("Completed!")
  
  return(res)
}


.fread_maybe_gz <- function(
    file,
    ...,
    nThread = data.table::getDTthreads()
) {
  
  if (grepl("\\.gz$", file)) {
    data.table::fread(
      cmd = paste("zcat", shQuote(file)),
      ...,
      nThread = nThread
    )
  } else {
    data.table::fread(
      file,
      ...,
      nThread = nThread
    )
  }
}


load_mhap <- function(
    mhap_file,
    mhap_file_1based = TRUE,
    sel_chr = NULL,
    nThread = data.table::getDTthreads()
) {
  
  if (is.null(mhap_file) || !file.exists(mhap_file)) {
    stop("mhap_file does not exist: ", mhap_file)
  }
  
  message("Reading mHap file: ", mhap_file)
  
  mhap_df <- .fread_maybe_gz(
    mhap_file,
    header = FALSE,
    nThread = nThread,
    showProgress = TRUE
  )
  
  if (ncol(mhap_df) < 6L) {
    stop("mHap file must contain at least 6 columns: chr, start, end, epi, nreads, strand.")
  }
  
  mhap_df <- mhap_df[, 1:6]
  
  data.table::setnames(
    mhap_df,
    c("chr", "start", "end", "epi", "nreads", "strand")
  )
  
  mhap_df[, chr := as.character(chr)]
  mhap_df[, start := as.integer(start)]
  mhap_df[, end := as.integer(end)]
  mhap_df[, epi := as.character(epi)]
  mhap_df[, nreads := as.integer(nreads)]
  mhap_df[, strand := as.character(strand)]
  
  if (!mhap_file_1based) {
    mhap_df[, start := start + 1L]
  }
  
  mhap_df <- mhap_df[
    !is.na(chr) &
      !is.na(start) &
      !is.na(end) &
      start <= end &
      !is.na(nreads) &
      nreads > 0L &
      !is.na(epi) &
      nchar(epi) > 0L
  ]
  
  if (nrow(mhap_df) == 0L) {
    stop("mhap_file has 0 valid rows - nothing to process.")
  }
  
  if (!is.null(sel_chr)) {
    
    mhap_df <- mhap_df[chr %in% sel_chr]
    
    if (nrow(mhap_df) == 0L) {
      stop(
        "No mHap records on selected chromosomes: ",
        paste(sel_chr, collapse = ", ")
      )
    }
  }
  
  data.table::setorder(mhap_df, chr, start, end)
  
  mhap_df
}


load_cpg_file <- function(
    cpg_file,
    cpg_file_1based = TRUE,
    sel_chr = NULL,
    nThread = data.table::getDTthreads()
) {
  
  if (is.null(cpg_file) || !file.exists(cpg_file)) {
    stop("cpg_file does not exist: ", cpg_file)
  }
  
  message("Reading CpG file: ", cpg_file)
  
  cpg_raw <- .fread_maybe_gz(
    cpg_file,
    header = FALSE,
    nThread = nThread,
    showProgress = TRUE
  )
  
  if (ncol(cpg_raw) < 2L) {
    stop("cpg_file must contain at least two columns: chr and position.")
  }
  
  cpg_df <- data.table::data.table(
    chr = as.character(cpg_raw[[1]]),
    pos = as.integer(cpg_raw[[2]])
  )
  
  if (!cpg_file_1based) {
    cpg_df[, pos := pos + 1L]
  }
  
  cpg_df <- cpg_df[!is.na(chr) & !is.na(pos)]
  
  if (!is.null(sel_chr)) {
    cpg_df <- cpg_df[chr %in% sel_chr]
  }
  
  cpg_df <- unique(cpg_df)
  
  data.table::setorder(cpg_df, chr, pos)
  
  cpg_df
}


normalize_cpg_df <- function(cpg_df) {
  
  data.table::setDT(cpg_df)
  
  if (all(c("chr", "pos") %in% names(cpg_df))) {
    
    out <- cpg_df[, .(
      chr = as.character(chr),
      pos = as.integer(pos)
    )]
    
  } else if (all(c("chr", "ref_pos") %in% names(cpg_df))) {
    
    out <- cpg_df[, .(
      chr = as.character(chr),
      pos = as.integer(ref_pos)
    )]
    
  } else if (all(c("seqnames", "pos") %in% names(cpg_df))) {
    
    out <- cpg_df[, .(
      chr = as.character(seqnames),
      pos = as.integer(pos)
    )]
    
  } else if (all(c("chr", "start") %in% names(cpg_df))) {
    
    out <- cpg_df[, .(
      chr = as.character(chr),
      pos = as.integer(start)
    )]
    
  } else {
    
    stop("Cannot infer CpG coordinate columns. Expected chr/pos, chr/ref_pos, seqnames/pos, or chr/start.")
  }
  
  out <- out[!is.na(chr) & !is.na(pos)]
  out <- unique(out)
  
  data.table::setorder(out, chr, pos)
  
  out
}


load_regions_bed <- function(
    region_file,
    sel_chr = NULL,
    merge_regions = TRUE,
    region_file_1based = FALSE,
    nThread = data.table::getDTthreads()
) {
  
  if (is.null(region_file) || !file.exists(region_file)) {
    stop("region_file does not exist: ", region_file)
  }
  
  bed <- .fread_maybe_gz(
    region_file,
    header = FALSE,
    nThread = nThread,
    showProgress = TRUE
  )
  
  if (ncol(bed) < 3L) {
    stop("region_file must contain at least 3 columns: chr, start, end.")
  }
  
  regions <- data.table::data.table(
    chr = as.character(bed[[1]]),
    start = as.integer(bed[[2]]),
    end = as.integer(bed[[3]])
  )
  
  if (!region_file_1based) {
    regions[, start := start + 1L]
  }
  
  regions <- regions[
    !is.na(chr) &
      !is.na(start) &
      !is.na(end) &
      start <= end
  ]
  
  if (!is.null(sel_chr)) {
    regions <- regions[chr %in% sel_chr]
  }
  
  if (nrow(regions) == 0L) {
    stop("No regions remain after chromosome filtering.")
  }
  
  regions_gr <- GenomicRanges::GRanges(
    seqnames = regions$chr,
    ranges = IRanges::IRanges(
      start = regions$start,
      end = regions$end
    )
  )
  
  if (merge_regions) {
    regions_gr <- GenomicRanges::reduce(regions_gr)
  }
  
  regions_gr
}


filter_mhap_by_regions <- function(
    mhap_df,
    regions_gr
) {
  
  if (is.null(regions_gr)) {
    return(mhap_df)
  }
  
  data.table::setDT(mhap_df)
  
  mhap_gr <- GenomicRanges::GRanges(
    seqnames = mhap_df$chr,
    ranges = IRanges::IRanges(
      start = mhap_df$start,
      end = mhap_df$end
    )
  )
  
  ov <- GenomicRanges::findOverlaps(
    mhap_gr,
    regions_gr,
    ignore.strand = TRUE
  )
  
  if (length(ov) == 0L) {
    return(mhap_df[0])
  }
  
  mhap_df[unique(S4Vectors::queryHits(ov))]
}


filter_cpg_by_regions <- function(
    cpg_df,
    regions_gr
) {
  
  if (is.null(regions_gr)) {
    return(cpg_df)
  }
  
  data.table::setDT(cpg_df)
  
  cpg_gr <- GenomicRanges::GRanges(
    seqnames = cpg_df$chr,
    ranges = IRanges::IRanges(
      start = cpg_df$pos,
      end = cpg_df$pos
    )
  )
  
  ov <- GenomicRanges::findOverlaps(
    cpg_gr,
    regions_gr,
    ignore.strand = TRUE
  )
  
  if (length(ov) == 0L) {
    return(cpg_df[0])
  }
  
  cpg_df[unique(S4Vectors::queryHits(ov))]
}


filter_mhap_by_target_cpg <- function(
    mhap_df,
    cpg_df_target,
    min_cpg_per_read = 1L
) {
  
  data.table::setDT(mhap_df)
  data.table::setDT(cpg_df_target)
  
  mhap_gr <- GenomicRanges::GRanges(
    seqnames = mhap_df$chr,
    ranges = IRanges::IRanges(
      start = mhap_df$start,
      end = mhap_df$end
    )
  )
  
  cpg_gr <- GenomicRanges::GRanges(
    seqnames = cpg_df_target$chr,
    ranges = IRanges::IRanges(
      start = cpg_df_target$pos,
      end = cpg_df_target$pos
    )
  )
  
  ov <- GenomicRanges::findOverlaps(
    mhap_gr,
    cpg_gr,
    ignore.strand = TRUE
  )
  
  if (length(ov) == 0L) {
    return(mhap_df[0])
  }
  
  n_target_cpg <- tabulate(
    S4Vectors::queryHits(ov),
    nbins = nrow(mhap_df)
  )
  
  mhap_df[, n_target_cpg := as.integer(n_target_cpg)]
  
  mhap_df[n_target_cpg >= min_cpg_per_read]
}


build_sparse_matrix_from_mhap_parallel_chunks <- function(
    mhap_df,
    cpg_df_all,
    cpg_df_target = cpg_df_all,
    code_map = c(default = 1L, `0` = 2L, `1` = 3L),
    min_cpg_per_read = 1L,
    rm_zero_col = FALSE,
    strand_sign = FALSE,
    strict = TRUE,
    cores = 4,
    chunks_per_core = 10
) {
  
  message("Building sparse matrix from mHap by target CpG chunks...")
  
  data.table::setDT(mhap_df)
  data.table::setDT(cpg_df_all)
  data.table::setDT(cpg_df_target)
  
  cpg_df_all <- normalize_cpg_df(cpg_df_all)
  cpg_df_target <- normalize_cpg_df(cpg_df_target)
  
  cpg_df_all[, cpg_name := paste(chr, pos, sep = "_")]
  cpg_df_target[, cpg_name := paste(chr, pos, sep = "_")]
  
  cpg_df_all <- unique(cpg_df_all)
  cpg_df_target <- unique(cpg_df_target)
  
  chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM", "MT")
  
  cpg_df_all[, chr_order := match(chr, chr_levels)]
  cpg_df_all[is.na(chr_order), chr_order := length(chr_levels) + 1L]
  data.table::setorder(cpg_df_all, chr_order, chr, pos)
  cpg_df_all[, chr_order := NULL]
  
  cpg_df_target[, chr_order := match(chr, chr_levels)]
  cpg_df_target[is.na(chr_order), chr_order := length(chr_levels) + 1L]
  data.table::setorder(cpg_df_target, chr_order, chr, pos)
  cpg_df_target[, chr_order := NULL]
  
  if (nrow(cpg_df_target) == 0L) {
    stop("cpg_df_target has 0 rows. No target CpG columns available.")
  }
  
  mhap_df[, mhap_row_id := seq_len(.N)]
  mhap_df[, nreads := as.integer(nreads)]
  mhap_df <- mhap_df[!is.na(nreads) & nreads > 0L]
  
  if (nrow(mhap_df) == 0L) {
    stop("No valid mHap records with nreads > 0.")
  }
  
  message("Filtering mHap records with at least one target CpG...")
  
  mhap_df <- filter_mhap_by_target_cpg(
    mhap_df = mhap_df,
    cpg_df_target = cpg_df_target,
    min_cpg_per_read = min_cpg_per_read
  )
  
  if (nrow(mhap_df) == 0L) {
    stop("No mHap records overlap target CpGs after filtering.")
  }
  
  mhap_df[, row_end := cumsum(nreads)]
  mhap_df[, row_start := row_end - nreads + 1L]
  
  total_reads <- max(mhap_df$row_end)
  
  message("Total artificial reads: ", total_reads)
  message("Total target CpG columns: ", nrow(cpg_df_target))
  
  row_names <- make_mhap_read_names(
    mhap_id = mhap_df$mhap_row_id,
    nreads = mhap_df$nreads
  )
  
  target_base_key <- cpg_df_target$cpg_name
  
  cpg_pos_by_chr <- split(cpg_df_all$pos, cpg_df_all$chr)
  
  n_chunks <- max(1L, cores * chunks_per_core)
  
  chunk_list <- prepare_mhap_sparse_chunks(
    cpg_df_target = cpg_df_target,
    mhap_df = mhap_df,
    n_chunks = n_chunks
  )
  
  message("Number of chunks: ", length(chunk_list))
  
  if (cores > 1L) {
    
    cl <- parallel::makeCluster(cores, type = "PSOCK")
    doParallel::registerDoParallel(cl)
    
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
    
    mat_list <- tryCatch(
      {
        foreach::foreach(
          chunk = chunk_list,
          .packages = c("Matrix", "data.table"),
          .export = c(
            "build_mhap_sparse_one_chunk",
            "get_cpg_between_mhap"
          )
        ) %dopar% {
          build_mhap_sparse_one_chunk(
            chunk = chunk,
            cpg_pos_by_chr = cpg_pos_by_chr,
            total_reads = total_reads,
            code_map = code_map,
            strand_sign = strand_sign,
            strict = strict
          )
        }
      },
      finally = {
        try(parallel::stopCluster(cl), silent = TRUE)
        foreach::registerDoSEQ()
      }
    )
    
  } else {
    
    mat_list <- lapply(
      chunk_list,
      build_mhap_sparse_one_chunk,
      cpg_pos_by_chr = cpg_pos_by_chr,
      total_reads = total_reads,
      code_map = code_map,
      strand_sign = strand_sign,
      strict = strict
    )
  }
  
  message("Combining sparse mHap chunks...")
  
  mat_list <- mat_list[!vapply(mat_list, is.null, logical(1))]
  
  if (length(mat_list) == 0L) {
    stop("No mHap chunks produced a sparse matrix.")
  }
  
  mat_base <- Reduce(Matrix::cbind2, mat_list)
  mat_base <- mat_base[, target_base_key, drop = FALSE]
  
  rownames(mat_base) <- row_names
  
  if (rm_zero_col) {
    keep_cols <- Matrix::colSums(mat_base != 0) > 0
    mat_base <- mat_base[, keep_cols, drop = FALSE]
  }
  
  read_coord_gr <- build_mhap_read_granges(
    mhap_df = mhap_df,
    row_names = row_names
  )
  
  list(
    sparse_matrix = mat_base,
    read_coord_gr = read_coord_gr,
    mhap_df_used = mhap_df,
    cpg_df_target = cpg_df_target
  )
}


prepare_mhap_sparse_chunks <- function(
    cpg_df_target,
    mhap_df,
    n_chunks
) {
  
  data.table::setDT(cpg_df_target)
  data.table::setDT(mhap_df)
  
  cpg_by_chr <- split(cpg_df_target, cpg_df_target$chr)
  
  chunk_list <- list()
  k <- 1L
  
  for (chr_i in names(cpg_by_chr)) {
    
    cpg_chr <- cpg_by_chr[[chr_i]]
    
    n_chr_chunks <- max(
      1L,
      round(nrow(cpg_chr) / nrow(cpg_df_target) * n_chunks)
    )
    
    chr_split <- split(
      cpg_chr,
      cut(
        seq_len(nrow(cpg_chr)),
        breaks = n_chr_chunks,
        labels = FALSE
      )
    )
    
    for (block_cpg in chr_split) {
      
      block_start <- min(block_cpg$pos)
      block_end <- max(block_cpg$pos)
      
      block_mhap <- mhap_df[
        chr == chr_i &
          end >= block_start &
          start <= block_end
      ]
      
      chunk_list[[k]] <- list(
        chunk_id = k,
        chr = chr_i,
        start = block_start,
        end = block_end,
        block_cpg = block_cpg,
        block_mhap = block_mhap
      )
      
      k <- k + 1L
    }
  }
  
  chunk_list
}


build_mhap_sparse_one_chunk <- function(
    chunk,
    cpg_pos_by_chr,
    total_reads,
    code_map = c(default = 1L, `0` = 2L, `1` = 3L),
    strand_sign = FALSE,
    strict = TRUE
) {
  
  block_cpg <- chunk$block_cpg
  block_mhap <- chunk$block_mhap
  
  empty_mat <- function() {
    Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = c(total_reads, nrow(block_cpg)),
      dimnames = list(NULL, block_cpg$cpg_name)
    )
  }
  
  if (nrow(block_cpg) == 0L) {
    return(empty_mat())
  }
  
  if (nrow(block_mhap) == 0L) {
    return(empty_mat())
  }
  
  block_col_index <- stats::setNames(
    seq_len(nrow(block_cpg)),
    block_cpg$cpg_name
  )
  
  default_code <- base::unname(code_map["default"])
  
  if (is.na(default_code)) {
    stop("code_map must contain a 'default' entry.")
  }
  
  entries_i <- list()
  entries_j <- list()
  entries_x <- list()
  
  out_idx <- 0L
  chr_i <- chunk$chr
  
  if (!chr_i %in% names(cpg_pos_by_chr)) {
    return(empty_mat())
  }
  
  pos_vec_chr <- cpg_pos_by_chr[[chr_i]]
  
  for (k in seq_len(nrow(block_mhap))) {
    
    chr_k <- block_mhap$chr[k]
    start_k <- block_mhap$start[k]
    end_k <- block_mhap$end[k]
    epi_k <- block_mhap$epi[k]
    nreads_k <- as.integer(block_mhap$nreads[k])
    strand_k <- block_mhap$strand[k]
    row_start_k <- block_mhap$row_start[k]
    row_end_k <- block_mhap$row_end[k]
    mhap_row_id_k <- block_mhap$mhap_row_id[k]
    
    if (is.na(nreads_k) || nreads_k <= 0L) {
      next
    }
    
    pattern_vec <- strsplit(epi_k, "", fixed = TRUE)[[1]]
    
    full_cpg_pos <- get_cpg_between_mhap(
      pos_vec = pos_vec_chr,
      start = start_k,
      end = end_k
    )
    
    if (length(full_cpg_pos) != length(pattern_vec)) {
      
      msg <- paste0(
        "mHap pattern length does not match number of CpGs for mHap row ",
        mhap_row_id_k,
        " (", chr_k, ":", start_k, "-", end_k, "). ",
        "nchar(epi) = ", length(pattern_vec),
        "; n CpGs = ", length(full_cpg_pos), "."
      )
      
      if (strict) {
        stop(msg)
      } else {
        warning(msg)
        next
      }
    }
    
    full_cpg_names <- paste(chr_k, full_cpg_pos, sep = "_")
    
    local_col_idx <- base::unname(block_col_index[full_cpg_names])
    keep <- !is.na(local_col_idx)
    
    if (!any(keep)) {
      next
    }
    
    local_col_idx <- local_col_idx[keep]
    local_pattern <- pattern_vec[keep]
    
    local_values <- base::unname(code_map[local_pattern])
    local_values[is.na(local_values)] <- default_code
    local_values <- as.integer(local_values)
    
    if (strand_sign && strand_k == "-") {
      local_values <- -local_values
    }
    
    new_rows <- row_start_k:row_end_k
    
    out_idx <- out_idx + 1L
    
    entries_i[[out_idx]] <- rep(new_rows, each = length(local_col_idx))
    entries_j[[out_idx]] <- rep(local_col_idx, times = nreads_k)
    entries_x[[out_idx]] <- rep(local_values, times = nreads_k)
  }
  
  if (length(entries_i) == 0L) {
    return(empty_mat())
  }
  
  i <- unlist(entries_i, use.names = FALSE)
  j <- unlist(entries_j, use.names = FALSE)
  x <- unlist(entries_x, use.names = FALSE)
  
  Matrix::sparseMatrix(
    i = i,
    j = j,
    x = x,
    dims = c(total_reads, nrow(block_cpg)),
    dimnames = list(NULL, block_cpg$cpg_name),
    giveCsparse = TRUE
  )
}


get_cpg_between_mhap <- function(
    pos_vec,
    start,
    end
) {
  
  if (length(pos_vec) == 0L) {
    return(integer(0))
  }
  
  lo <- findInterval(start - 1L, pos_vec) + 1L
  hi <- findInterval(end, pos_vec)
  
  if (lo > hi) {
    return(integer(0))
  }
  
  pos_vec[lo:hi]
}


make_mhap_read_names <- function(
    mhap_id,
    nreads
) {
  
  unlist(
    Map(
      function(id, n) {
        paste0("mhap_", id, "_", seq_len(n))
      },
      mhap_id,
      nreads
    ),
    use.names = FALSE
  )
}


build_mhap_read_granges <- function(
    mhap_df,
    row_names
) {
  
  strand_vec <- mhap_df$strand
  strand_vec[!strand_vec %in% c("+", "-")] <- "*"
  
  read_coord_gr <- GenomicRanges::GRanges(
    seqnames = rep(mhap_df$chr, mhap_df$nreads),
    ranges = IRanges::IRanges(
      start = rep(mhap_df$start, mhap_df$nreads),
      end = rep(mhap_df$end, mhap_df$nreads)
    ),
    strand = rep(strand_vec, mhap_df$nreads)
  )
  
  names(read_coord_gr) <- row_names
  
  read_coord_gr
}