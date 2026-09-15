# ---- Modkit wrapper ---------------------------------------------------------------

#' Prepare Modkit input files
#'
#' Run the `prep_modkit.sh` script inside Conda environment.
#'
#' @param sample Sample identifier used in output filenames.
#' @param chromosome Chromosome or contig to process, for example `"chr1"`.
#' @param reference_fasta Path to the reference genome FASTA file.
#' @param bam_file Path to an indexed BAM file containing modified-base calls.
#' @param outdir Base output directory.
#' @param samtools_threads Number of threads used by samtools and bgzip.
#' @param readMod_multicore Value passed to the Bismark `--multicore` option.
#' @param env_name Name of the Conda environment containing Modkit, samtools and bedtools.
#' @param conda Optional path to the Conda executable.
#'
#' @return Log about the completed execution.
#'
#' @export

prepare_modkit_input <- function(
    sample,
    chromosome,
    reference_fasta,
    bam_file,
    outdir,
    samtools_threads = 4L,
    readMod_multicore = 1L,
    env_name = "modkit_env",
    conda = 'conda'
) {
  
  if (length(samtools_threads) != 1L | is.na(samtools_threads) | samtools_threads < 1L) {
    stop("`samtools_threads` must be a positive integer.", call. = FALSE)
  }
  
  if (length(readMod_multicore) != 1L | is.na(readMod_multicore) | readMod_multicore < 1L) {
    stop("`readMod_multicore` must be a positive integer.", call. = FALSE)
  }
  
  env_variables <- c(
    sprintf(
      "SAMTOOLS_THREADS=%d", samtools_threads
    ),
    sprintf(
      "MODKIT_MULTICORE=%d", readMod_multicore
    )
  )
  
  run_input_preparation(
    method = "modkit",
    sample = sample,
    chromosome = chromosome,
    reference_fasta = reference_fasta,
    alignment_file = bam_file,
    outdir = outdir,
    env_name = env_name,
    conda = conda,
    env_variables = env_variables
  )
}


# ---- Bismark wrapper --------------------------------------------------------------

#' Prepare Bismark input files
#'
#' Run the `prep_bismark.sh` script inside a Conda environment.
#'
#' @param sample Sample identifier used in output filenames.
#' @param chromosome Chromosome or contig to process, for example `"chr1"`.
#' @param reference_fasta Path to the reference genome FASTA file.
#' @param alignment_file Path to an indexed BAM or CRAM file containing Bismark tags.
#' @param outdir Base output directory.
#' @param samtools_threads Number of threads used by samtools and bgzip.
#' @param readMod_multicore Value passed to the Bismark `--multicore` option.
#' @param keep_intermediates Logical. Keep intermediate filtered and name-sorted CRAM files. (Default: FALSE)
#' @param env_name Name of the conda environment containing Bismark.
#' @param conda Optional path to the Conda executable.
#'
#' @return Log about the completed execution.
#'
#' @export

prepare_bismark_input <- function(
    sample,
    chromosome,
    reference_fasta,
    alignment_file,
    outdir,
    samtools_threads = 4L,
    readMod_multicore = 1L,
    keep_intermediates = FALSE,
    env_name = "bismark_env",
    conda = 'conda'
) {
  
  if (length(samtools_threads) != 1L | is.na(samtools_threads) | samtools_threads < 1L) {
    stop("`samtools_threads` must be a positive integer.", call. = FALSE)
  }
  
  if (length(readMod_multicore) != 1L | is.na(readMod_multicore) | readMod_multicore < 1L) {
    stop("`readMod_multicore` must be a positive integer.", call. = FALSE)
  }
  
  if (!is.logical(keep_intermediates) | length(keep_intermediates) != 1L | is.na(keep_intermediates)) {
    stop("`keep_intermediates` must be TRUE or FALSE.", call. = FALSE)
  }
  
  env_variables <- c(
    sprintf(
      "SAMTOOLS_THREADS=%d", samtools_threads
    ),
    sprintf(
      "BISMARK_MULTICORE=%d", readMod_multicore
    ),
    sprintf(
      "KEEP_INTERMEDIATES=%d", as.integer(keep_intermediates)
    )
  )
  
  run_input_preparation(
    method = "bismark",
    sample = sample,
    chromosome = chromosome,
    reference_fasta = reference_fasta,
    alignment_file = alignment_file,
    outdir = outdir,
    env_name = env_name,
    conda = conda,
    env_variables = env_variables
  )
}


# ---- Utilities -----------------------------------------------------------

epistats_file <- function(...) {
  path <- system.file(..., package = "epistats")
  
  if (!nzchar(path)) {
    stop(
      "The requested file was not found in the installed epistats package.",
      call. = FALSE
    )
  }
  
  return(path)
}


check_string <- function(x, name) {
  if (!is.character(x) | length(x) != 1L | is.na(x) | !nzchar(x)) {
    stop(
      sprintf("`%s` must be a single non-empty character string.", name),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


check_conda_environment <- function(
    conda_command,
    env_name
) {
  result <- suppressWarnings(
    system2(
      command = conda_command,
      args = c(
        "run",
        "--name",
        shQuote(env_name),
        "bash",
        "-c",
        shQuote("exit 0")
      ),
      stdout = TRUE,
      stderr = TRUE
    )
  )
  
  status <- attr(result, "status")
  
  if (is.null(status)) {
    status <- 0L
  }
  
  identical(as.integer(status), 0L)
}


run_input_preparation <- function(
    method,
    sample,
    chromosome,
    reference_fasta,
    alignment_file,
    outdir,
    env_name,
    conda = NULL,
    env_variables
) {
  method <- match.arg(method, c("modkit", "bismark"))
  
  check_string(sample, "sample")
  check_string(chromosome, "chromosome")
  check_string(reference_fasta, "reference_fasta")
  check_string(alignment_file, "alignment_file")
  check_string(outdir, "outdir")
  check_string(env_name, "env_name")
  
  if (grepl("[/\\\\]", sample)) {
    stop("`sample` cannot contain directory separators.", call. = FALSE)
  }
  
  if (!file.exists(reference_fasta)) {
    stop(sprintf("Reference FASTA not found: %s", reference_fasta), call. = FALSE)
  }
  
  if (!file.exists(alignment_file)) {
    stop(sprintf("Alignment file not found: %s", alignment_file), call. = FALSE)
  }
  
  extension <- tolower(tools::file_ext(alignment_file))
  
  if (method == "modkit" && extension != "bam") {
    stop("`prep_modkit.sh` currently requires an input BAM file.", call. = FALSE)
  }
  
  if (method == "bismark" && !extension %in% c("bam", "cram")) {
    stop("`prep_bismark.sh` requires a BAM or CRAM input file.", call. = FALSE)
  }
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  reference_fasta <- normalizePath(reference_fasta, mustWork = TRUE)
  alignment_file <- normalizePath(alignment_file, mustWork = TRUE)
  outdir <- normalizePath(outdir, mustWork = TRUE)
  
  script_name <- switch(
    method,
    modkit = "prep_modkit.sh",
    bismark = "prep_bismark.sh"
  )
  
  script <- epistats_file("exec", script_name)
  
  conda_command <- if (file.exists(conda)) {
    normalizePath(conda, mustWork = TRUE)
  } else {
    Sys.which(conda)
  }
  
  if (!nzchar(conda_command)) {
    stop(
      paste0(
        "Conda executable not found. ",
        "Ensure that Conda is available in PATH or provide its full path ",
        "through the `conda` argument."
      ),
      call. = FALSE
    )
  }
  
  
  if (!check_conda_environment(conda_command = conda_command, env_name = env_name)) {
    stop(
      sprintf(
        paste0(
          "Conda environment '%s' was not found or cannot be started.\n",
          "Create it first with:\n\n",
          "  create_input_environment(\n",
          "    method = \"%s\",\n",
          "    env_name = \"%s\",\n",
          "    conda = \"%s\"\n",
          "  )"
        ),
        env_name,
        method,
        env_name,
        conda_command
      ),
      call. = FALSE
    )
  }
  
  args <- c(
    "run",
    "--no-capture-output",
    "--name",
    shQuote(env_name),
    "bash",
    shQuote(script),
    shQuote(sample),
    shQuote(chromosome),
    shQuote(reference_fasta),
    shQuote(alignment_file),
    shQuote(outdir)
  )
  
  message(
    sprintf(
      "Running %s input preparation for sample '%s', chromosome '%s'.",
      method,
      sample,
      chromosome
    )
  )
  
  status <- system2(
    command = conda_command,
    args = args,
    stdout = "",
    stderr = "",
    env = env_variables
  )
  
  if (!identical(status, 0L)) {
    stop(
      sprintf(
        "%s input preparation failed with exit status %s.",
        method,
        status
      ),
      call. = FALSE
    )
  }
  
  message("Input preparation completed successfully.")
  
  invisible(
    list(
      method = method,
      sample = sample,
      chromosome = chromosome,
      environment = env_name,
      script = script,
      output_directory = outdir,
      status = status
    )
  )
}


# ---- Conda environment management -------------------------------------------------

#' Create a Conda environment for input preparation
#'
#' Create or update the Conda environment required by the Modkit or Bismark
#' input-preparation script distributed with the package.
#'
#' @param method Input-preparation method. One of `"modkit"` or `"bismark"`.
#' @param env_name Name assigned to the Conda environment. By default,
#'   `"modkit_env"` or `"bismark_env"` is used according to `method`.
#' @param update Logical. If `FALSE`, create a new environment. If `TRUE`,
#'   update an existing environment and remove dependencies that are no longer
#'   present in the YAML file.
#' @param conda Optional path to the Conda executable. By default, Conda is
#'   located using `Sys.which("conda")`.
#'
#' @return Invisibly returns the Conda exit status.
#'
#' @export

create_input_environment <- function(
    method = c("modkit", "bismark"),
    env_name = NULL,
    update = FALSE,
    conda = 'conda'
) {
  method <- match.arg(method)
  
  if (!is.logical(update) | length(update) != 1L | is.na(update)) {
    stop(
      "`update` must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  
  default_env_name <- switch(
    method,
    modkit = "modkit_env",
    bismark = "bismark_env"
  )
  
  if (is.null(env_name)) {
    env_name <- default_env_name
  }
  
  check_string(env_name, "env_name")
  
  yaml_name <- switch(
    method,
    modkit = "environment_modkit.yml",
    bismark = "environment_bismark.yml"
  )
  
  yaml_file <- epistats_file("conda", yaml_name)
  
  conda_command <- if (file.exists(conda)) {
    normalizePath(conda, mustWork = TRUE)
  } else {
    Sys.which(conda)
  }
  
  if (!nzchar(conda_command)) {
    stop(
      paste0(
        "Conda executable not found. ",
        "Ensure that Conda is available in PATH or provide its full path ",
        "through the `conda` argument."
      ),
      call. = FALSE
    )
  }
  
  env_exists <- check_conda_environment(
    conda_command = conda_command,
    env_name = env_name
  )
  
  if (env_exists && !update) {
    message(
      sprintf(
        "Conda environment '%s' already exists. Nothing to create.",
        env_name
      )
    )
    
    return(invisible(0L))
  }
  
  if (!env_exists && update) {
    stop(
      sprintf(
        paste0(
          "Conda environment '%s' does not exist and therefore cannot be updated. ",
          "Run `create_input_environment(..., update = FALSE)` to create it."
        ),
        env_name
      ),
      call. = FALSE
    )
  }
  
  if (update) {
    args <- c(
      "env",
      "update",
      "--name",
      shQuote(env_name),
      "--file",
      shQuote(yaml_file),
      "--prune"
    )
    
    action <- "Updating"
  } else {
    args <- c(
      "env",
      "create",
      "--name",
      shQuote(env_name),
      "--file",
      shQuote(yaml_file)
    )
    
    action <- "Creating"
  }
  
  message(
    sprintf(
      "%s Conda environment '%s' using %s.",
      action,
      env_name,
      basename(yaml_file)
    )
  )
  
  status <- system2(
    command = conda_command,
    args = args,
    stdout = "",
    stderr = ""
  )
  
  if (!identical(status, 0L)) {
    stop(
      sprintf("Conda environment operation failed with exit status %s.", status),
      call. = FALSE
    )
  }
  
  message(
    sprintf("Conda environment '%s' is ready.", env_name)
  )
  
  invisible(status)
}