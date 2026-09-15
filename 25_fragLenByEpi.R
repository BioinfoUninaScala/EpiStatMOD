#' Extracting the length of fragments for each epiallele observed in a genomic region.
#'
#' @import magrittr 
#' @param matrix Epiallele matrix.
#' @return A string reporting for each epiallele the length of fragments carring such an epiallele.
#' @export

fragLenByEpi <- function(matrix) {
  # fragLenByEpi
  frag_len <- attr(matrix, "frag_len")
  
  if (is.null(frag_len)){
    return(NA)
  } else {
    df <- gen_frag_len_mat(matrix)
    
    agg <- df %>% 
      dplyr::group_by(epi) %>% 
      dplyr::summarise(
        aggr_frag_len = paste(frag_len, collapse =','),
        .groups = "drop"
      ) %>% 
      dplyr::arrange(epi)
    
    res_str <- paste0(agg$epi, ":", agg$aggr_frag_len, collapse = ";")
    
    return(res_str)
  }
  
}


#' Calculating the average fragment length for each epiallele observed in a genomic region.
#'
#' @import magrittr 
#' @param matrix Epiallele matrix.
#' @return A string reporting for each epiallele the average length of fragments carrying such an epiallele.
#' @export

meanfragLenByEpi <- function(matrix) {
  frag_len <- attr(matrix, "frag_len")
  
  if (is.null(frag_len)){
    return(NA)
  } else {
    df <- gen_frag_len_mat(matrix)
    
    agg <- df %>% 
      dplyr::group_by(epi) %>% 
      dplyr::summarise(
        mean_frag_len = round(mean(frag_len, na.rm = TRUE), 1),
        .groups = "drop"
      ) %>% 
      dplyr::arrange(epi)
    
    res_str <- paste0(agg$epi, ":", agg$mean_frag_len, collapse = ",")
    
    return(res_str)
  }
}


gen_frag_len_mat <- function(matrix) {
  frag_len <- attr(matrix, "frag_len")
  if (is.null(frag_len)) {
    stop("The 'frag_len' attribute is not present in the epialleles matrix.")
  }
  
  if (is.null(names(frag_len))) {
    warning("L'attributo 'frag_len' non ha nomi: si assume che l'ordine corrisponda alle righe della matrice.")
    frag_len_row <- as.numeric(frag_len)
  } else {
    frag_len_row <- frag_len[rownames(matrix)]
  }
  
  if (length(frag_len_row) != nrow(matrix)) {
    stop("La lunghezza del vettore 'frag_len' non corrisponde al numero di righe della matrice.")
  }
  
  mat <- as.matrix(matrix)
  epi_str <- apply(mat, 1, function(x) paste0(x, collapse = ""))
  
  df <- data.frame(
    epi = epi_str,
    frag_len = as.numeric(frag_len_row),
    stringsAsFactors = FALSE
  )
  return(df)
}
