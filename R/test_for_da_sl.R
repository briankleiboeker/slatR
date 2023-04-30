#' Test shotgun lipidomic data for differential abundance of all lipid species
#' 
#'
#' This function tests all lipid species for differential abundance across two biological conditions in a shotgun lipidomics data matrix
#' @param m a matrix of shotgun lipidomics data where rownames are lipid compositions, structures, or any other unique identifier, columns are samples or replicates, and all entries of matrix are coercible to numeric.
#' @param a column names corresponding to control condition.
#' @param b column names corresponding to experimental condition.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param paired a logical indicating whether you want a paired t-test. Defaults to FALSE.
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used. Defaults to FALSE.
#' @param p.adj.method p-value correction method. One of c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"). Call ?p.adjust.methods for more information. Defaults to "BH".
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

test_for_da_sl <- function(m,
                           a,
                           b,
                           alternative = "two.sided",
                           paired = FALSE,
                           var.equal = FALSE,
                           p.adj.method = "BH"){
  
  drop <- unname(unlist(apply(m,1,function(x) ifelse(length(unique(x))==1,F,T))))
  
  if(length(which(!drop)) >= 1){
    warning(paste("One or more rows was dropped due to data being constant. Dropped row(s): ",which(!drop)))
  }
  
  m <- m[drop,]
  pvals<-apply(m,
               1,
               function(x) t.test(as.numeric(x[a]),
                                  as.numeric(x[b]),
                                  alternative = alternative,
                                  paired = paired,
                                  var.equal = var.equal)$p.value )
  compwise.df<-data.frame(pvals,row.names = rownames(m))
  compwise.df$p.adj<-p.adjust(compwise.df$pvals,method = p.adj.method)
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

  
