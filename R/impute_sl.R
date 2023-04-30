#' Impute missing values in shotgun lipidomics data matrix
#' 
#'
#' This function performs a sample minimum imputation of a shotgun lipidomics data matrix where missing values are represented by NA
#' @param m a matrix of shotgun lipidomics data where rownames are lipid compositions, structures, or any other unique identifier, columns are samples or replicates, and all entries of matrix are coercible to numeric.
#' @keywords normalization
#' @export
#' @examples
#' m <- df %>% dplyr::select(c("rel.intensity","sample","comp")) %>% 
#'             tidyr::pivot_wider(names_from = sample,
#'                                values_from = rel.intensity,
#'                                values_fn = mean))
#' 
#' m.norm <- normalize_sl(m)
#' 
#' m.final <- impute_sl(m.norm)

impute_sl <- function(m){
  if(!is.matrix(m)){
    m <- as.matrix(m)
  }
  if(!is.numeric(m)){
    m <- convert_to_numeric(m)
  }
  apply(m, 2, function(x) replace(x, is.na(x), (min(x, na.rm = TRUE)) ) )
}

convert_to_numeric <- function(m){
  m2 <- apply(m,2,as.numeric)
  colnames(m2) <- colnames(m)
  rownames(m2) <- rownames(m)
  return(m2)
}