#' Guess exact structures from general structures
#' 
#'
#' This function guesses an exact lipid structure  (e.g. 16:0/18:1-PE) from general lipid structures (e.g. PE(34:1)) from the assign_structures() function 
#' @param structures list (or columns) containing the general structures (from assign_structures() function) for which exact structures are to be assigned.
#' @keywords structural assignment
#' @export
#' @examples
#' df$c <- extract_num_elements("C",df$composition)
#' df$h <- extract_num_elements("H",df$composition)
#' df$o <- extract_num_elements("O",df$composition)
#' df$n <- extract_num_elements("N",df$composition)
#' df$p <- extract_num_elements("P",df$composition)
#' df$na <- extract_num_elements("Na",df$composition)
#' df$cl <- extract_num_elements("Cl",df$composition)
#'
#' df$gen.structures <- assign_structures(df$c,df$h,df$o,df$n,df$p,df$na,df$cl,
#'                                        "neg.ion",
#'                                        "m.minus.h",
#'                                        "euk",
#'                                        c("pe","pa","pg","ffa"),
#'                                        6)
#'                                        
#' exact_structures_df <- guess_exact_structures(df$gen.structures)
#' 
#' df <- cbind(df,exact_structures_df)

guess_exact_structures <- function(structures){
  df <- data.frame("structures" = structures)
  
  if(any(structures %in% pl_database$gen.structure)){
    
    df$exact.structure.1 <- unlist(lapply(structures,function(x) ifelse(x %in% pl_database$gen.structure,unname(unlist(pl_database[match(x,pl_database$gen.structure),2:ncol(pl_database)]))[1],"" )) )
    df$exact.structure.2 <- unlist(lapply(structures,function(x) ifelse(x %in% pl_database$gen.structure,unname(unlist(pl_database[match(x,pl_database$gen.structure),2:ncol(pl_database)]))[2],"" )) )
    df$exact.structure.3 <- unlist(lapply(structures,function(x) ifelse(x %in% pl_database$gen.structure,unname(unlist(pl_database[match(x,pl_database$gen.structure),2:ncol(pl_database)]))[3],"" )) )
    df$exact.structure.4 <- unlist(lapply(structures,function(x) ifelse(x %in% pl_database$gen.structure,unname(unlist(pl_database[match(x,pl_database$gen.structure),2:ncol(pl_database)]))[4],"" )) )
    df$exact.structure.5 <- unlist(lapply(structures,function(x) ifelse(x %in% pl_database$gen.structure,unname(unlist(pl_database[match(x,pl_database$gen.structure),2:ncol(pl_database)]))[5],"" )) )
    df$exact.structure.6 <- unlist(lapply(structures,function(x) ifelse(x %in% pl_database$gen.structure,unname(unlist(pl_database[match(x,pl_database$gen.structure),2:ncol(pl_database)]))[6],"" )) )
    df$exact.structure.7 <- unlist(lapply(structures,function(x) ifelse(x %in% pl_database$gen.structure,unname(unlist(pl_database[match(x,pl_database$gen.structure),2:ncol(pl_database)]))[7],"" )) )
    
    if( all(df$exact.structure.7 == "") ){
      df<-df[, !(colnames(df) %in% c("exact.structure.7"))]
    }
    if( all(df$exact.structure.6 == "") ){
      df<-df[, !(colnames(df) %in% c("exact.structure.6"))]
    }
    if( all(df$exact.structure.5 == "") ){
      df<-df[, !(colnames(df) %in% c("exact.structure.5"))]
    }
    if( all(df$exact.structure.4 == "") ){
      df<-df[, !(colnames(df) %in% c("exact.structure.4"))]
    }
    if( all(df$exact.structure.3 == "") ){
      df<-df[, !(colnames(df) %in% c("exact.structure.3"))]
    }
    if( all(df$exact.structure.2 == "") ){
      df<-df[, !(colnames(df) %in% c("exact.structure.2"))]
    }
    if( all(df$exact.structure.1 == "") ){
      df<-df[, !(colnames(df) %in% c("exact.structure.1"))]
    }
    
  }
  
  df
}