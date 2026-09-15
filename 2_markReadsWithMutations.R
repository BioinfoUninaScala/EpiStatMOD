#' Mark reads carrying SNP, insertion, or deletion events
#'
#' @description
#' This function annotates the reads represented in a sparse CpG matrix according
#' to whether they carry at least one mutation event. Mutation events are read
#' from an external mutation table and summarized per read as three logical
#' flags: \code{SNPevent}, \code{INS}, and \code{DEL}.
#'
#' The sparse matrix is used to define the set of reads to be annotated, using
#' \code{rownames(sparseMat)} as read identifiers. The mutation file is expected
#' to contain one row per mutation event and to include at least the columns
#' \code{read_id}, \code{chrom}, \code{ref_pos}, \code{event}, and
#' \code{event_len}.
#'
#' If a BED file with regions of interest is provided, only mutation events
#' overlapping those regions are considered. In this case, a read is marked as
#' \code{TRUE} for a given event type only if it carries at least one event of
#' that type inside the regions of interest.
#'
#' @param sparseMat A sparse matrix, typically a \code{Matrix::dgCMatrix}, with
#'   rows corresponding to reads. The row names of \code{sparseMat} must contain
#'   the read identifiers to be annotated.
#' @param mutation_file Character string or \code{data.frame}/\code{data.table}.
#'   Path to the mutation-event file, or an already loaded mutation table. The
#'   expected file format contains at least the following columns:
#'   \code{read_id}, \code{chrom}, \code{ref_pos}, \code{event}, and
#'   \code{event_len}. Optional columns such as \code{base_qual}, \code{mapq},
#'   and \code{strand} can be used for filtering or retained in the input.
#' @param events_to_be_excluded vector of characters containing the mutation events 
#'  that must be excluded from the read attributes. These events must be provided 
#'  in the following format: event-on-the-reference_event-on-the-read, e.g. C_T.
#'  This exclusion vector is essential when working with BS and Illumina 5-base 
#'  sequencing data, where C->T events indicate methylation status of cytosines
#'  rather than the presence of SNPs.
#' @param roi_bed Optional character string. Path to a BED file containing
#'   regions of interest. If provided, only mutation events overlapping these
#'   regions are considered when setting the output logical flags. If
#'   \code{NULL}, mutation events are considered genome-wide.
#' @param mut_file_1based Logical. Whether the genomic coordinates in
#'   \code{mutation_file}, specifically \code{ref_pos}, are already 1-based.
#'   If \code{FALSE}, positions are converted from 0-based to 1-based by adding
#'   1. Default is \code{TRUE}.
#' @param roi_bed_1based Logical. Whether the coordinates in \code{roi_bed} are
#'   already 1-based inclusive. Standard BED files are 0-based, half-open, so the
#'   default is \code{FALSE}. When \code{FALSE}, the BED start coordinate is
#'   converted to 1-based by adding 1.
#' @param min_base_qual Optional numeric. Minimum base quality required for a
#'   mutation event to be considered. If \code{NULL}, no filtering by base
#'   quality is applied. This filter requires a \code{base_qual} column in the
#'   mutation table.
#' @param min_mapq Optional numeric. Minimum mapping quality required for a
#'   mutation event to be considered. If \code{NULL}, no filtering by mapping
#'   quality is applied. This filter requires a \code{mapq} column in the
#'   mutation table.
#' @param nThread Integer. Number of threads to use when reading large mutation
#'   files with \code{data.table::fread}. Default is \code{4}.
#'
#' @return A \code{data.frame} with one row per read in \code{sparseMat} and
#'   four columns:
#'   \describe{
#'     \item{\code{readID}}{Read identifier, taken from \code{rownames(sparseMat)}.}
#'     \item{\code{SNPevent}}{Logical flag indicating whether the read carries
#'       at least one SNP event.}
#'     \item{\code{INS}}{Logical flag indicating whether the read carries at
#'       least one insertion event.}
#'     \item{\code{DEL}}{Logical flag indicating whether the read carries at
#'       least one deletion event.}
#'   }
#'
#' @details
#' Mutation events are interpreted from the \code{event} column of the mutation
#' table. The expected values are \code{"SNP"}, \code{"INS"}, and
#' \code{"DEL"}.
#'
#' For SNP and insertion events, \code{ref_pos} is treated as the event anchor
#' position. For deletion events, the genomic interval is defined as
#' \code{ref_pos} to \code{ref_pos + event_len - 1}. This allows deletions to
#' be correctly overlapped with regions of interest.
#'
#' When \code{roi_bed} is provided, mutation events are first converted to
#' genomic ranges and overlapped with the regions in the BED file using
#' \code{GenomicRanges::findOverlaps}. Only overlapping events contribute to
#' the final per-read flags.
#'
#' Reads present in \code{sparseMat} but absent from the mutation table are
#' retained in the output and assigned \code{FALSE} for all three mutation
#' classes.
#'
#' @importFrom data.table fread data.table as.data.table is.data.table setDTthreads
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   mut_flags <- markReadsWithMutations(
#'     sparseMat = sparseMat,
#'     mutation_file = "sample.mutations.tsv.gz"
#'   )
#'
#'   mut_flags_roi <- markReadsWithMutations(
#'     sparseMat = sparseMat,
#'     mutation_file = "sample.mutations.tsv.gz",
#'     events_to_be_excluded = c('C_T'),
#'     roi_bed = "regions_of_interest.bed",
#'     roi_bed_1based = FALSE,
#'     min_base_qual = 20,
#'     min_mapq = 30,
#'     nThread = 8
#'   )
#' }


markReadsWithMutations <- function(
    sparseMat,
    mutation_file,
    events_to_be_excluded = c('C_T'),
    roi_bed = NULL,
    mut_file_1based = TRUE,
    roi_bed_1based = FALSE,
    min_base_qual = NULL,
    min_mapq = NULL,
    nThread = 4
) {
  
  requireNamespace("data.table")
  requireNamespace("GenomicRanges")
  requireNamespace("IRanges")
  requireNamespace("S4Vectors")
  

  # 1. Read IDs dalla sparse matrix

  
  read_ids <- rownames(sparseMat)
  
  if (is.null(read_ids)) {
    stop("`sparseMat` must have rownames corresponding to read IDs.")
  }
  
  read_ids <- as.character(read_ids)
  
  out <- data.table::data.table(
    readID = read_ids,
    SNP = FALSE,
    INS = FALSE,
    DEL = FALSE
  )
  

  # 2. Caricamento file mutazioni

  
  if (data.table::is.data.table(mutation_file) | is.data.frame(mutation_file)) {
    mut_dt <- data.table::as.data.table(mutation_file)
  } else {
    if (!file.exists(mutation_file)) {
      stop("`mutation_file` does not exist: ", mutation_file)
    }
    
    # mut_dt <- freadMutationsForSparseReads(
    #   mutation_file = mutation_file,
    #   sparseMat = sparseMat
    # )
    
    mut_dt <- data.table::fread(mutation_file, nThread = nThread)
  }
  
  required_cols <- c(
    "read_id",
    "chrom",
    "ref_pos",
    "event",
    "event_len"
  )
  
  missing_cols <- base::setdiff(required_cols, colnames(mut_dt))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in `mutation_file`: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (nrow(mut_dt) == 0L) {
    return(as.data.frame(out))
  }
  
  mut_dt <- mut_dt[read_id %in% read_ids]
  
  mut_dt <- mut_dt[
    ,
    .(
      readID = as.character(read_id),
      chr = as.character(chrom),
      ref_pos = as.integer(ref_pos),
      event = toupper(as.character(event)),
      event_len = as.integer(event_len),
      event_type = paste(toupper(ref_base), toupper(read_base ), sep = '_'),
      base_qual = if ("base_qual" %in% colnames(mut_dt)) as.numeric(base_qual) else NA_real_,
      mapq = if ("mapq" %in% colnames(mut_dt)) as.numeric(mapq) else NA_real_
    )
  ]
  
  mut_dt <- mut_dt[
    !is.na(readID) &
      !is.na(chr) &
      !is.na(ref_pos) &
      !is.na(event) &
      !event_type %in% events_to_be_excluded
  ]
  

  
  if (nrow(mut_dt) == 0L) {
    return(as.data.frame(out))
  }
  

  # 3. Filtri opzionali su qualità
  if (!is.null(min_base_qual) && "base_qual" %in% colnames(mut_dt)) {
    mut_dt <- mut_dt[
      !is.na(base_qual) &
        base_qual >= min_base_qual
    ]
  }
  
  if (!is.null(min_mapq) && "mapq" %in% colnames(mut_dt)) {
    mut_dt <- mut_dt[
      !is.na(mapq) &
        mapq >= min_mapq
    ]
  }
  
  if (nrow(mut_dt) == 0L) {
    return(as.data.frame(out))
  }
  

  # 4. Coordinate mutazioni
  # SNP: posizione singola
  # INS: uso ref_pos come anchor
  # DEL: uso ref_pos : ref_pos + event_len - 1
  
  if (!mut_file_1based) {
    mut_dt[, ref_pos := ref_pos + 1L]
  }
  
  mut_dt[
    is.na(event_len) | event_len < 1L,
    event_len := 1L
  ]
  
  mut_dt[
    event == "DEL",
    mut_end := ref_pos + event_len - 1L
  ]
  
  mut_dt[
    event != "DEL",
    mut_end := ref_pos
  ]
  

  # 5. Eventuale filtro sulle ROI
  if (!is.null(roi_bed)) {
    
    if (!file.exists(roi_bed)) {
      stop("`roi_bed` does not exist: ", roi_bed)
    }
    
    roi_dt <- data.table::fread(
      roi_bed,
      header = FALSE
    )
    
    if (ncol(roi_dt) < 3L) {
      stop("`roi_bed` must contain at least 3 columns: chr, start, end.")
    }
    
    roi_dt <- roi_dt[, 1:3]
    colnames(roi_dt) <- c("chr", "start", "end")
    
    roi_dt[
      ,
      `:=`(
        chr = as.character(chr),
        start = as.integer(start),
        end = as.integer(end)
      )
    ]
    
    # BED standard: start 0-based, end 1-based esclusivo.
    # Lo converto in coordinate genomiche 1-based inclusive.
    if (!roi_bed_1based) {
      roi_dt[, start := start + 1L]
    }
    
    roi_dt <- roi_dt[
      !is.na(chr) &
        !is.na(start) &
        !is.na(end) &
        end >= start
    ]
    
    if (nrow(roi_dt) == 0L) {
      stop("No valid intervals in `roi_bed`.")
    }
    
    mut_gr <- GenomicRanges::GRanges(
      seqnames = mut_dt$chr,
      ranges = IRanges::IRanges(
        start = mut_dt$ref_pos,
        end = mut_dt$mut_end
      )
    )
    
    roi_gr <- GenomicRanges::GRanges(
      seqnames = roi_dt$chr,
      ranges = IRanges::IRanges(
        start = roi_dt$start,
        end = roi_dt$end
      )
    )
    
    hits <- GenomicRanges::findOverlaps(
      mut_gr,
      roi_gr,
      ignore.strand = TRUE
    )
    
    if (length(hits) == 0L) {
      return(as.data.frame(out))
    }
    
    mut_dt <- mut_dt[
      unique(S4Vectors::queryHits(hits))
    ]
  }
  
  if (nrow(mut_dt) == 0L) {
    return(as.data.frame(out))
  }
  

  # 6. Booleani per read
  mut_by_read <- mut_dt[
    ,
    .(
      SNP = any(event == "SNP", na.rm = TRUE),
      INS = any(event == "INS", na.rm = TRUE),
      DEL = any(event == "DEL", na.rm = TRUE)
    ),
    by = readID
  ]
  
  out <- merge(
    out[, .(readID)],
    mut_by_read,
    by = "readID",
    all.x = TRUE
  )
  
  out[is.na(SNP), SNP := FALSE]
  out[is.na(INS), INS := FALSE]
  out[is.na(DEL), DEL := FALSE]
  
  out[
    ,
    .(
      readID = as.character(readID),
      SNP = as.logical(SNP),
      INS = as.logical(INS),
      DEL = as.logical(DEL)
    )
  ]
}


freadMutationsForSparseReads <- function(
    mutation_file,
    sparseMat,
    read_id_file = tempfile(fileext = ".read_ids.txt"),
    sep = "\t"
) {
  
  requireNamespace("data.table")
  
  read_ids <- rownames(sparseMat)
  
  if (is.null(read_ids)) {
    stop("`sparseMat` must have rownames corresponding to read IDs.")
  }
  
  read_ids <- unique(as.character(read_ids))
  
  data.table::fwrite(
    data.table::data.table(readID = read_ids),
    file = read_id_file,
    sep = "\t",
    col.names = FALSE
  )
  
  if (!file.exists(mutation_file)) {
    stop("`mutation_file` does not exist: ", mutation_file)
  }
  
  reader_cmd <- if (grepl("\\.gz$", mutation_file)) {
    paste("gzip -cd --", shQuote(mutation_file))
  } else {
    paste("cat --", shQuote(mutation_file))
  }
  
  awk_script <- paste0(
    "awk -v idsfile=", shQuote(read_id_file), " ",
    shQuote(
      paste0(
        "BEGIN { FS = OFS = \"\t\"; ",
        "while ((getline line < idsfile) > 0) ids[line] = 1 } ",
        "NR == 1 { print; next } ",
        "($1 in ids) { print }"
      )
    )
  )
  
  cmd <- paste(reader_cmd, "|", awk_script)
  
  data.table::fread(
    cmd = cmd,
    sep = sep,
    header = TRUE
  )
}