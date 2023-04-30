#' Normalize shotgun lipidomics data matrix
#' 
#'
#' This function performs a quantile normalization (using preprocessCore::normalize.quantiles) of a shotgun lipidomics data matrix of 100 rows or more
#' @param m a matrix of shotgun lipidomics data where rownames are lipid compositions, structures, or any other unique identifier, columns are samples or replicates, and all entries of matrix are coercible to numeric.
#' @param log.normalize log2 transform data before normalization? Default = TRUE.
#' @export
#' @examples
#' m <- df %>% dplyr::select(c("rel.intensity","sample","comp")) %>% 
#'             tidyr::pivot_wider(names_from = sample,
#'                                values_from = rel.intensity,
#'                                values_fn = mean))
#' 
#' m.norm <- normalize_sl(m)

normalize_sl <- function(m, log.normalize = T){
  
  if(!is.matrix(m)){
    m <- as.matrix(m)
  }
  if(nrow(m)<100){
    abort("too few rows (= lipids) to normalize matrix (minimum = 100)")
  }
  if(!is.numeric(m)){
    m <- convert_to_numeric(m)
  }
  
  if(log.normalize){
    m2 <- preprocessCore::normalize.quantiles(log2(m))
    rownames(m2)<-rownames(m)
    colnames(m2)<-colnames(m)
    return(m2)
  }else{
    m2 <- preprocessCore::normalize.quantiles(m)
    rownames(m2)<-rownames(m)
    colnames(m2)<-colnames(m)
    return(m2)
  }
  
}

convert_to_numeric <- function(m){
  m2 <- apply(m,2,as.numeric)
  colnames(m2) <- colnames(m)
  rownames(m2) <- rownames(m)
  return(m2)
}