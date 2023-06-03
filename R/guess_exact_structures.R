#' Guess exact structures from general structures
#' 
#'
#' This function guesses an exact lipid structure (e.g. 16:0/18:1-PE) from general lipid structures (e.g. PE(34:1)) from the assign_structures() function 
#' @param structures list (or dataframe column) containing the general structures (from assign_structures() function) for which exact structures are to be assigned.
#' @param strain which bacterial strain data originates from, if applicable. Strain options are "e.coli", "listeria", "streptococcus". Defaults to n/a (meaning data is either not from bacteria or is from an unsupported bacterial strain).
#' @keywords structural assignment
#' @export
#' @examples
#' df$gen.structures <- assign_structures(comp = df$composition,
#'                                        ion.mode = "neg.ion",
#'                                        domain = "euk",
#'                                        max.dbl.bnds = 8)
#'                                        
#' exact_structures_df <- guess_exact_structures(df$gen.structures)
#' 
#' df <- cbind(df,exact_structures_df)

guess_exact_structures <- function(structures,
                                   strain = "n/a"){
  
  if(!(strain %in% c("n/a", "e.coli","listeria","streptococcus"))){
    abort("invalid input for strain")
  }
  
  df <- data.frame("structure" = structures)
  
  if(any(structures %in% pl_database$gen.structure)){
    
    if(strain != "n/a"  && any(structures %in% pl_database[pl_database$strain == strain , colnames(pl_database) == "gen.structure"] ) ){
      pls2 <- pl_database[pl_database$strain == strain,]
      
      df$exact.structure.1 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[1],"" )) )
      df$exact.structure.2 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[2],"" )) )
      df$exact.structure.3 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[3],"" )) )
      df$exact.structure.4 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[4],"" )) )
      df$exact.structure.5 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[5],"" )) )
      df$exact.structure.6 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[6],"" )) )
      df$exact.structure.7 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[7],"" )) )
      
      df$strain.specific.assignment <- ifelse(df$exact.structure.1 != "","Y","N")
      
      df$exact.structure.1 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pl_database$gen.structure,
                                                                 unname(unlist(pl_database[match(y,pl_database$gen.structure),2:ncol(pl_database)]))[1],
                                                                 x)},
                                            df$exact.structure.1, # this is x
                                            df$structure # this is y
      ))
      df$exact.structure.2 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pl_database$gen.structure,
                                                                 unname(unlist(pl_database[match(y,pl_database$gen.structure),2:ncol(pl_database)]))[2],
                                                                 x)},
                                            df$exact.structure.2, # this is x
                                            df$structure # this is y
      ))
      df$exact.structure.3 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pl_database$gen.structure,
                                                                 unname(unlist(pl_database[match(y,pl_database$gen.structure),2:ncol(pl_database)]))[3],
                                                                 x)},
                                            df$exact.structure.3, # this is x
                                            df$structure # this is y
      ))
      df$exact.structure.4 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pl_database$gen.structure,
                                                                 unname(unlist(pl_database[match(y,pl_database$gen.structure),2:ncol(pl_database)]))[4],
                                                                 x)},
                                            df$exact.structure.4, # this is x
                                            df$structure # this is y
      ))
      df$exact.structure.5 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pl_database$gen.structure,
                                                                 unname(unlist(pl_database[match(y,pl_database$gen.structure),2:ncol(pl_database)]))[5],
                                                                 x)},
                                            df$exact.structure.5, # this is x
                                            df$structure # this is y
      ))
      df$exact.structure.6 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pl_database$gen.structure,
                                                                 unname(unlist(pl_database[match(y,pl_database$gen.structure),2:ncol(pl_database)]))[6],
                                                                 x)},
                                            df$exact.structure.6, # this is x
                                            df$structure # this is y
      ))
      df$exact.structure.7 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pl_database$gen.structure,
                                                                 unname(unlist(pl_database[match(y,pl_database$gen.structure),2:ncol(pl_database)]))[7],
                                                                 x)},
                                            df$exact.structure.7, # this is x
                                            df$structure # this is y
      ))
    }else{
      pls2 <- pl_database[pl_database$strain == "",]
      df$exact.structure.1 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[1],"" )) )
      df$exact.structure.2 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[2],"" )) )
      df$exact.structure.3 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[3],"" )) )
      df$exact.structure.4 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[4],"" )) )
      df$exact.structure.5 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[5],"" )) )
      df$exact.structure.6 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[6],"" )) )
      df$exact.structure.7 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[7],"" )) )
    }
    
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