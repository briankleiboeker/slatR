#' Assign elemental composition to m/z value
#' 
#'
#' This function assigns elemental composition to user-supplied m/z values in accordance with user-supplied elemental constraints. 
#' Atomic weights are referenced from: Coursey, J.S., Schwab, D.J., Tsai, J.J., and Dragoset, R.A. (2015), Atomic Weights and Isotopic Compositions (version 4.1). [Online] Available: http://physics.nist.gov/Comp [2023, June, 4]. National Institute of Standards and Technology, Gaithersburg, MD.
#' @param mz list of m/z values to which composition is to be assigned.
#' @param elements list of elements which may be included in composition (Exact element letter only, case sensitive).
#' @param isotope list of atomic number/isotope number of desired element (NOT atomic mass--atomic mass will be assigned from NIST table cited above).
#' @param mins.list list of minimum number of elements which are allowed in composition (must be in the same order as elements and isotope arguments). 
#' @param maxs.list list of maximum number of elements which are allowed in composition (must be in the same order as elements and isotope arguments). 
#' @param ionmode ion mode the samples were run in. Must be one of "neg.ion", "pos.ion", or "neutral". Used to account for mass from gain/loss of electron in ionization process.
#' @param error.ppm the maximum error (delta, units ppm) allowed for assigned compositions. Defaults to 5 ppm.
#' @param rdbrange the range of allowed RDBs (relative double bond equivalents) of assigned compositions. Defaults to -1 to 100 (denoted: c(-1,100)).
#' @param rdb.rule force assignment of elemental compositions to only those that have a certain RDB characteristic (e.g integer values or non-integer values)? Options are "none", "nonint" (non-integer, e.g. 0.5, 1.5, etc.), or "int" (integer, e.g. 0, 1, 2, etc.). Defaults to "none" meaning composition will be assigned as the closest matched m/z value regardless of RDB. 
#' @keywords elemental composition assignment
#' @export
#' @examples
#' assign_comp_from_mz(df$mz,
#'                     elements = c("C","H","O","N","P","Na","Cl"),
#'                     isotope = c(12,1,16,14,31,23,35),
#'                     mins.list = c(10,20,2,0,0,0,0),
#'                     maxs.list = c(60,160,16,2,2,0,1),
#'                     ionmode = "neg.ion",
#'                     error.ppm = 3,
#'                     rdb.rule = "nonint")

assign_comp_from_mz <- function(mz,
                                elements,
                                isotope,
                                mins.list,
                                maxs.list,
                                ionmode,
                                error.ppm = 5,
                                rdbrange = c(-1,100),
                                rdb.rule = "none"){
  
  mz <- as.numeric(mz)
  mins.list <- as.numeric(mins.list)
  maxs.list <- as.numeric(maxs.list)
  if(length(elements) != length(isotope) | length(isotope) != length(mins.list) | length(mins.list) != length(maxs.list)){
    rlang::abort("elements, isotope, mins.list, and maxs.list arguments are not all of the same length")
  }
  if(!is.numeric(error.ppm)){
    rlang::abort("error.ppm argument must be numeric")
  }
  if(length(rdbrange) != 2){
    rlang::abort("rdbrange argument must be a list of length 2, with first entry corresponding to lower limit for RDB and second entry corresponding to upper limit for RDB")
  }
  if(!(rdb.rule %in% c("none","int","nonint"))){
    rlang::abort("invalid rdb.rule argument. rdb.rule must be one of: none, int, or nonint")
  }
  if(!all(paste0(elements,isotope) %in% paste0(all.elements$element,all.elements$isotope))){
    rlang::abort(paste0("one or more invalid element/isotope combinations. Most likely invalid element(s) are: ",
                 paste(paste0(elements,isotope)[!(paste0(elements,isotope) %in% paste0(all.elements$element,all.elements$isotope))],collapse = ", ")
                 )
          )
  }
  if(!( ionmode %in% c("neutral","neg.ion","pos.ion") )){
    rlang::abort("invalid ionmode argument. ionmode must be one of: neutral, neg.ion, or pos.ion")
  }
  
  all.elements$full.element = paste0(all.elements$element,"[",all.elements$isotope,"]")
  atomic.mass <- all.elements$atomic.mass[match(paste0(elements,isotope),paste0(all.elements$element,all.elements$isotope))]
  
  if(rdb.rule == "none"){
    v <- c(atomic.mass,1/1837)
    
    ele.list <- c(mapply(paste0,elements,isotope),"e")
    
    if(ionmode == "neutral"){
      e.bound = 0
    }
    if(ionmode == "neg.ion"){
      e.bound = 1
    }
    if(ionmode == "pos.ion"){
      e.bound = -1
    }
    
    bounds <- list(lower = list(ind = c(1:length(v)), val = c(mins.list,e.bound)),
                   upper = list(ind = c(1:length(v)), val = c(maxs.list,e.bound)))
    
    rdb.row <- rep(0,length(v))
    rdb.row[ele.list %in% c("H1","D2","T3","Cl35","Cl37")] <- (-0.5)
    rdb.row[ele.list %in% c("C12","C13","C14")] <- 1
    rdb.row[ele.list %in% c("N14","N15","P31")] <- 0.5
    if(error.ppm == "" | is.null(error.ppm) | is.na(error.ppm)){
      error.ppm <- 5
    }
    
    clean_elements <- ifelse(paste0(elements,"[",isotope,"]") %in% all.elements[all.elements$isotopic.composition > 0.9,colnames(all.elements) == "full.element"],elements,paste0(elements,"[",isotope,"]"))
    order <- c("C","C[13]","C[14]","H","D[2]","T[3]","O","O[17]","O[18}","N","N[15]","P","Na","Cl[35]","Cl[37]")
    
    results <- lapply(mz,
             function(x) {
               a<-Rglpk::Rglpk_solve_LP(obj = v, mat = rbind(t(v),rdb.row,rdb.row), dir = c("<=","<=",">="), rhs = c(x,rdbrange[2]-1,rdbrange[1]-1), bounds = bounds,max = TRUE, types = "I")
               b<-Rglpk::Rglpk_solve_LP(obj = v, mat = rbind(t(v),rdb.row,rdb.row), dir = c(">=","<=",">="), rhs = c(x,rdbrange[2]-1,rdbrange[1]-1), bounds = bounds,max = FALSE, types = "I")
               if(((abs(a$optimum-x)/a$optimum)*10^6) <= error.ppm | ((abs(b$optimum-x)/b$optimum)*10^6) <= error.ppm){
                 if(abs(a$optimum-x) < abs(b$optimum-x)){
                   tempdf<-a$solution[!(ele.list %in% c("e","rdbvar"))]
                   clean_elements <- clean_elements[tempdf != 0]
                   tempdf <- tempdf[tempdf != 0]
                   tempdf[tempdf == 1] <- ""
                   return(
                     c(paste(mapply(paste0,
                                    clean_elements[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ],
                                    tempdf[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ]),
                             collapse = " "),
                       a$optimum,
                       ((x-a$optimum)/x)*10^6))
                 }else{
                   tempdf<-b$solution[!(ele.list %in% c("e","rdbvar"))]
                   clean_elements <- clean_elements[tempdf != 0]
                   tempdf <- tempdf[tempdf != 0]
                   tempdf[tempdf == 1] <- ""
                   return(
                     c(paste(mapply(paste0,
                                    clean_elements[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ],
                                    tempdf[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ]),
                             collapse = " "),
                       b$optimum,
                       ((b$optimum-x)/x)*10^6))
                 }
               }else{
                 return(c(rep("",3)))
               }
             }
      )
    return(data.frame("mz" = mz, 
                      "composition" = unlist(lapply(results,`[[`, 1)),
                      "theo.mass" = unlist(lapply(results,`[[`, 2)),
                      "delta" = unlist(lapply(results,`[[`, 3))
                      )
           )
    
    
    
  }else{
    if(rdb.rule == "nonint"){
      rhs.var = 0.5
    }else{
      #rdb rule = int
      rhs.var = 0
    }
    
    v <- c(atomic.mass,1/1837,0)
    ele.list <- c(mapply(paste0,elements,isotope),"e","rdbvar")
    
    if(ionmode == "neutral"){
      e.bound = 0
    }
    if(ionmode == "neg.ion"){
      e.bound = 1
    }
    if(ionmode == "pos.ion"){
      e.bound = -1
    }
    
    min.r = floor((rhs.var-rdbrange[2])/2)-1
    max.r = ceiling((rhs.var-rdbrange[1])/2)+1
    
    bounds <- list(lower = list(ind = c(1:length(v)), val = c(mins.list,e.bound,min.r)),
                   upper = list(ind = c(1:length(v)), val = c(maxs.list,e.bound,max.r)))
    
    rdb.row <- rep(0,length(v))
    rdb.row[ele.list %in% c("H1","D2","T3","Cl35","Cl37")] <- (-0.5)
    rdb.row[ele.list %in% c("C12","C13","C14")] <- 1
    rdb.row[ele.list %in% c("N14","N15","P31")] <- 0.5
    
    rdb.row2 <- rdb.row
    rdb.row2[length(rdb.row2)] <- 1
    
    if(error.ppm == "" | is.null(error.ppm) | is.na(error.ppm)){
      error.ppm <- 5
    }
    
    clean_elements <- ifelse(paste0(elements,"[",isotope,"]") %in% all.elements[all.elements$isotopic.composition > 0.9,colnames(all.elements) == "full.element"],elements,paste0(elements,"[",isotope,"]"))
    order <- c("C","C[13]","C[14]","H","D[2]","T[3]","O","O[17]","O[18}","N","N[15]","P","Na","Cl[35]","Cl[37]")
    
    
    results <- lapply(mz,function(x){
      a<-Rglpk::Rglpk_solve_LP(obj = v, mat = rbind(t(v),rdb.row,rdb.row,rdb.row2), dir = c("<=","<=",">=","=="), rhs = c(x,rdbrange[2]-1,rdbrange[1]-1,rhs.var), bounds,max = TRUE, types = "I")
      b<-Rglpk::Rglpk_solve_LP(obj = v, mat = rbind(t(v),rdb.row,rdb.row,rdb.row2), dir = c(">=","<=",">=","=="), rhs = c(x,rdbrange[2]-1,rdbrange[1]-1,rhs.var), bounds,max = FALSE, types = "I")
      
      if(((abs(a$optimum-x)/a$optimum)*10^6) <= error.ppm | ((abs(b$optimum-x)/b$optimum)*10^6) <= error.ppm){
        if(abs(a$optimum-x) < abs(b$optimum-x)){
          tempdf<-a$solution[!(ele.list %in% c("e","rdbvar"))]
          clean_elements <- clean_elements[tempdf != 0]
          tempdf <- tempdf[tempdf != 0]
          tempdf[tempdf == 1] <- ""
          return(c(paste(mapply(paste0,
                                clean_elements[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ],
                                tempdf[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ]),
                         collapse = " "),
                   a$optimum,
                   ((x-a$optimum)/x)*10^6))
        }else{
          tempdf<-b$solution[!(ele.list %in% c("e","rdbvar"))]
          clean_elements <- clean_elements[tempdf != 0]
          tempdf <- tempdf[tempdf != 0]
          tempdf[tempdf == 1] <- ""
          return(c(paste(mapply(paste0,
                                clean_elements[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ],
                                tempdf[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ]),
                         collapse = " "),
                   b$optimum,
                   ((b$optimum-x)/x)*10^6))
        }
      }else{
        return(c(rep("",3)))
      }})
    
    return(data.frame("mz" = mz, 
                      "composition" = unlist(lapply(results,`[[`, 1)),
                      "theo.mass" = unlist(lapply(results,`[[`, 2)),
                      "delta" = unlist(lapply(results,`[[`, 3))
                      )
           )
    
    
  }
  
}