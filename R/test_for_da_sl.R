#' Test all lipid species for differential abundance
#' 
#'
#' This function performs tests all species for differential abundance across two biological conditions
#' @param m a matrix of shotgun lipidomics data where rownames are lipid compositions, structures, or any other unique identifier, columns correspond to samples or replicates, and all entries of matrix are coercible to numeric (matrix should already be normalized and have had missing values imputted)
#' @param a column names corresponding to control condition
#' @param b column names corresponding to experimental condition
#' @keywords testing
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
#' 
#' results <- test_for_da_sl(m.final, 
#'                           c("wt1","wt2","wt3"), 
#'                           c("ko1","ko2","ko3"))

test_for_da_sl <- function(m,a,b){
  pvals<-apply(m,
               1,
               function(x) t.test(as.numeric(x[a]),as.numeric(x[b]))$p.value )
  compwise.df<-data.frame(pvals,row.names = rownames(m))
  compwise.df$p.adj<-p.adjust(compwise.df$pvals,method = "BH")
  compwise.df$ratio<-apply(m,
                           1,
                           function(x) mean(as.numeric(x[a]))/mean(as.numeric(x[b])) )
  compwise.df$intensity<-apply(m,
                               1,
                               function(x) mean(c(as.numeric(x[a]),
                                                  as.numeric(x[b]))
                               ) )
  compwise.df[order(compwise.df$p.adj),]
}

  
