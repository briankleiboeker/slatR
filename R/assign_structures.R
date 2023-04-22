#' Assign structures to a list of compositions
#' 
#'
#' This function assigns general lipid structure (e.g. PC(32:1)) to a user-defined elemental composition
#' @param c,h,o,n,p lists (or columns) containing the numbers of each respective elemennt in each respective composition (these should come from from extract_num_elements() function)
#' @param na,cl lists (or columns) containing the numbers of each respective elemennt in each respective composition. Default to zero.
#' @param ion.mode which ion mode the samples were run in. One of "neg.ion","pos.ion",or "neutral" (neutral searches for structures as exact compositions / neutral species)
#' @param adducts which adducts to search for compositions as. Defaults to c("m.plus.h","m.plus.ammonia","m.plus.sodium") in positive ion mode, c("m.minus.h","m.plus.chloride","m.minus.2h") in negative ion mode, and the inclusive union of these two lists in neutral mode
#' @param domain what domain of life the data originates from. "euk" for eukaryotic, "bact" for bacterial.
#' @param lois which lipids to look for. Defaults to c("pc","sm","tag","dag") in positive ion mode, c("pe","cer","cl","pi","pg","pa","ps","ffa") in negative ion mode, and the inclusive union of these two lists in neutral mode.
#' @param max.dbl.bnds maximum acyl+alkyl chain double bonds allowed in returned structures. Defaults to 14.
#' @keywords structureal assignment
#' @examples
#' assign_structures(df$c,df$h,df$o,df$n,df$p,df$na,df$cl,
#'                   "neg.ion",
#'                   "m.minus.h",
#'                   "euk",
#'                    c("pe","pa","pg","ffa"),
#'                    6)


assign_structures <- function(c,h,o,n,p,
                              na = NULL,
                              cl = NULL,
                              ion.mode,
                              adducts = NULL,
                              domain = "euk",
                              lois = NULL,
                              max.dbl.bnds = 14){
  if(!is.numeric(c) | !is.numeric(h) | !is.numeric(o) | !is.numeric(n) | !is.numeric(p) | (!is.null(na) & !is.numeric(na)) | (!is.null(cl) & !is.numeric(cl)) ){
    abort("one of elemental inputs is not numeric")
  }
  if(is.null(na)){
    na<-rep(0,length(c))
  }
  if(is.null(cl)){
    cl<-rep(0,length(c))
  }
  if(!(ion.mode %in% c("neg.ion","pos.ion","neutral")) | length(ion.mode) > 1){
    abort("problem in ion.mode argument (ion.mode must be one of neg.ion, pos.ion, or neutral)")
  }
  if(!is.null(adducts) & any(!(adducts %in% c("m.minus.h","m.plus.chloride","m.minus.2h","m.plus.h","m.plus.ammonia","m.plus.sodium")))){
    abort("one or more of inputted adducts are invalid (only valid adducts are m.minus.h,m.plus.chloride,m.minus.2h,m.plus.h,m.plus.ammonia,m.plus.sodium)")
  }
  if(!(domain %in% c("euk","bact")) | length(domain) > 1){
    abort("Invalid domain argument (domain must be one of euk or bact)")
  }
  if(!is.null(lois) & any(!(lois %in% c("pe","cer","cl","pi","pg","pa","ps","ffa","pc","sm","tag","dag")))){
    abort("one or more user-supplied input to lois argument is invalid (only valid lois are pe,cer,cl,pi,pg,pa,ps,ffa,pc,sm,tag,dag)")
  }
  if(!is.numeric(max.dbl.bnds) | length(max.dbl.bnds) > 1){
    abort("max.dbl.bnds argument not numeric")
  }
  # if(!(length(c) == length(h) == length(o) == length(n) == length(p))){
  #   abort("elemental columns/lists are not of equal length")
  # }
  
  #define adducts
  if(is.null(adducts)){
    if(ion.mode == "neg.ion"){
      adducts <- c("m.minus.h","m.plus.chloride","m.minus.2h")
    }
    if(ion.mode == "pos.ion"){
      adducts <- c("m.plus.h","m.plus.ammonia","m.plus.sodium")
    }
    if(ion.mode == "neutral"){
      adducts <- c("m.minus.h","m.plus.chloride","m.minus.2h","m.plus.h","m.plus.ammonia","m.plus.sodium")
    }
  }
  
  #define lois
  if(is.null(lois)){
    if(ion.mode == "neg.ion"){
      lois <- c("pe","cer","cl","pi","pg","pa","ps","ffa")
    }
    if(ion.mode == "pos.ion"){
      lois <- c("pc","sm","tag","dag")
    }
    if(ion.mode == "neutral"){
      lois <- c("pe","cer","cl","pi","pg","pa","ps","ffa","pc","sm","tag","dag")
    }
  }

  mapply(assign_species,c,h,o,n,p,na,cl,ion.mode,rep(list(adducts),length(c)),( c - (h/2) + ((n+p)/2) +1 ),domain,rep(list(lois),length(c)),max.dbl.bnds,SIMPLIFY = T)

}
#' @export

assign_species <- function(c,h,o,n,p,na,cl,ion.mode,adducts,rdb,domain,lois,max.dbl.bnds){
  if(is.na(rdb) | is.null(adducts)){
    return(NA)
  }
  if(ion.mode == "neutral"){
    return(structure_from_exact_comp(c,h,o,n,p,c("nonacylglycerols","acylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
  }
  if(ion.mode == "pos.ion"){
    if("m.plus.sodium" %in% adducts){
      if(na > 0){
        return(structure_from_exact_comp(c,h,o,n,p,c("nonacylglycerols","acylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
        # look for sodium adduct species -- ALL
        # assign exact composition by changing nothing
        # regardless if find one or not, return the result here
      }
      # else, keep going
    }
    if("m.plus.ammonia" %in% adducts){
      if(n > 0){
        result <- structure_from_exact_comp(c,h-4,o,n-1,p,c("acylglycerols"),domain,lois,as.numeric(max.dbl.bnds))
        if(!is.na(result)){
          return(result)
        }
        # look for ammonia adduct species -- only DAG/TAG 
        # do this by subtracting 4 H's and 1 N and seeing if the species could exist
        # if find a match, return it, else, go look for m.plus.h ion if ("m.plus.h" %in% adducts)
        # but if we do find one here, add it to a list of things to return
      }
      # else, keep going
    }
    if("m.plus.h" %in% adducts){
      return(structure_from_exact_comp(c,h-1,o,n,p,c("nonacylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
      # look for lipid matches
      # do this by subtracting 1H and seeing if the species could exist
      # if it does exist -> return m.plus.h
      # if it doesn't exist -> return some clever paste0 showing what we searched for, beginning with 'Species not detected (traceback:'
    }
    return(NA)
    # if we make it here, return some clever paste0 showing what we searched for, beginning with 'Species not detected (traceback:'
    
  }
  if(ion.mode == "neg.ion"){
    if("m.plus.chloride" %in% adducts){
      if(cl > 0){
        return(structure_from_exact_comp(c,h,o,n,p,c("nonacylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
        # look for chloride adduct species
        # assign exact composition by changing nothing
        # regardless if find one or not, return the result here
      }
      # else, keep going
    }
    if("m.minus.h" %in% adducts){
      result <- structure_from_exact_comp(c,h+1,o,n,p,c("nonacylglycerols"),domain,lois,as.numeric(max.dbl.bnds))
      if(!is.na(result)){
        return(result)
      }
      # look for lipid matches
      # do this by adding 1H and seeing if the species could exist
      # if it does exist return m.minus.h
      # if it doesn't exist -> return some clever paste0 showing what we searched for, beginning with 'Species not detected (traceback:'
      
    }
    if("m.minus.2h" %in% adducts){
      
      result <- structure_from_exact_comp(c,h+2,o,n,p,c("nonacylglycerols"),domain,lois,as.numeric(max.dbl.bnds))
      if(!is.na(result) && grepl("CL",result)){
        return(result)
      }
      # look for lipid matches
      # do this by adding 1H and seeing if the species could exist
      # if it does exist return m.minus.h
      # if it doesn't exist -> return some clever paste0 showing what we searched for, beginning with 'Species not detected (traceback:'
      
    }
    return(NA)
  }
}
assign.el.vs.plasmalogen.vs.lyso <- function(c,o,rdb,max.o,general.class,base.carbons,lyso.carbon.cutoff,base.rdb){
  if(o == max.o){
    #not an EL or lyso- species
    return(paste0(general.class,"("))
  }else{
    #either an EL or lyso-species (but not both, as that would be o-2)
    if((c - base.carbons) <= lyso.carbon.cutoff){
      #lyso species
      return(paste0("L",general.class,"("))
    }else{
      #an EL
      if(rdb == (base.rdb- 1)){
        #a non-plasmalogen EL
        #this is true no matter what
        return(paste0(general.class,"(O-"))
      }else{
        # there are one or more double bonds -- these can be adjacent to ether linkage (plasmalogen) or anywhere else
        if(general.class %in% c("PC","PI")){
          return(paste0(general.class,"(O-"))
        }else{
          return(paste0(general.class,"(P-"))
        }
        
      }
      
    }
  }
}
structure_from_exact_comp <- function(c,h,o,n,p,which_to_look_for,domain,lois,max.dbl.bnds){
  if(domain == "bact"){
    lyso.carbon.cutoff.value <- 19
  }else{
    lyso.carbon.cutoff.value <- 27
  }
  rdb.equiv = ( c - (h/2) + ((n+p)/2) +1 )
  if(rdb.equiv %% 1 != 0 | (((2*(c+n))-h)/2) > 14 ){
    return(NA)
  }
  
  if( "acylglycerols" %in% which_to_look_for){
    if(n == 0 && p == 0 && o == 5 && c >= 35){
      #DAG
      if(rdb.equiv < 2 | !("dag" %in% lois)){
        return(NA)
      }
      
      n.double.bonds <-  (rdb.equiv - 2)
      n.fa.carbons <- (c - 3)
      
      if(n.double.bonds > max.dbl.bnds){
        return(NA)
      }
      
      return(paste0("DAG(",n.fa.carbons,":",n.double.bonds,")"))
      #DAG has base 2 RDB , base 3 carbons
    }
    if(n == 0 && p == 0 && o == 6 && c >= 49){
      #TAG
      if(rdb.equiv < 3 | !("tag" %in% lois)){
        return(NA)
      }
      
      n.double.bonds <-  (rdb.equiv - 3)
      n.fa.carbons <- (c - 3)
      
      if(n.double.bonds > max.dbl.bnds){
        return(NA)
      }
      
      return(paste0("TAG(",n.fa.carbons,":",n.double.bonds,")"))
      #TG has base 3 RDB , base 3 carbons
    }
    if(length(which_to_look_for)==1){
      return(NA)
    }
  }
  if("nonacylglycerols" %in% which_to_look_for){ #look for everything but TAG and DAG
    if(n == 0){
      if(p == 1){
        #pa, pg, or pi
        if(o >=7 && o<= 8){
          #PA has base 2 rdb , base 3 carbons
          if((o == 8 && rdb.equiv < 2) | (o == 7 && rdb.equiv < 1) | !("pa" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 2)
          n.fa.carbons <- (c - 3)
          species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,8,"PA",3,lyso.carbon.cutoff.value,2)
          if(grepl("O-",species) | grepl("^L",species)){
            n.double.bonds <- n.double.bonds+1
          }
          
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          
          return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
          
        }
        if(o >=9 && o<= 10 ){
          #PG has base 2 rdb , base 6 carbons
          if((o == 10 && rdb.equiv < 2) | (o == 9 && rdb.equiv < 1) | !("pg" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 2)
          n.fa.carbons <- (c - 6)
          
          species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,10,"PG",6,lyso.carbon.cutoff.value,2)
          if(grepl("O-",species) | grepl("^L",species)){
            n.double.bonds <- n.double.bonds+1
          }
          
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          
          return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
          
        }
        if(o >=12 && o<= 13 ){
          #PI has base 3 rdb , base 9 carbons
          if((o == 13 && rdb.equiv < 3) | (o == 12 && rdb.equiv < 2) | !("pi" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 3)
          n.fa.carbons <- (c - 9)
          
          species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,13,"PI",9,lyso.carbon.cutoff.value,3)
          if(grepl("O-",species) | grepl("^L",species)){
            n.double.bonds <- n.double.bonds+1
          }
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          
          return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
          
        }
        return(NA)
      }
      if(c > 35 && p == 2){
        #some type of CL
        if(!("cl" %in% lois)){
          return(NA)
        }else{
          
          if(o == 15){
            #dilysoCL has 2 rdb, base 9 carbons
            if(rdb.equiv < 2){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 2)
            n.fa.carbons <- (c - 9)
            
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            
            return(paste0("dilysoCL(",n.fa.carbons,":",n.double.bonds,")"))
            
          }
          
          if(o == 16 && c > 50){
            #lysoCL has 3 rdb , base 9 carbons
            if(rdb.equiv < 3){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 3)
            n.fa.carbons <- (c - 9)
            
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            
            return(paste0("monolysoCL(",n.fa.carbons,":",n.double.bonds,")"))
            
          }
          
          if(o == 17 && c > 60){
            #CL has base 4 rdb , base 9 carbons
            if(rdb.equiv < 4){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 4)
            n.fa.carbons <- (c - 9)
            
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            
            return(paste0("CL(",n.fa.carbons,":",n.double.bonds,")"))
          }
          return(NA)
        }
      }
      if(p == 0 && o >= 2 && o <= 4 && c > 13){
        #some FFA probably
        if(o == 2 && c > 13 && c < 33){
          #plain old FFA
          if(rdb.equiv < 1 | !("ffa" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 1)
          n.fa.carbons <- (c - 0)
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          return(paste0(n.fa.carbons,":",n.double.bonds,"-FA"))
        }
        
        if(o == 3 && c > 13 && c < 33){
          #hFFA
          if(rdb.equiv < 1 | !("ffa" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 1)
          n.fa.carbons <- (c - 0)
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          return(paste0("h",n.fa.carbons,":",n.double.bonds,"-FA"))
        }
        
        if(o == 4 && c > 25 && c < 65){
          #FAHFA
          if(rdb.equiv < 2 | !("ffa" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 2)
          n.fa.carbons <- (c - 0)
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          return(paste0(n.fa.carbons,":",n.double.bonds,"-FAHFA"))
        }
        return(NA)
      }
      return(NA)
    }
    if(n == 2){
      if(p == 1 && o >= 5 && o<= 7){
        #SM has 1 rdb, base 5 carbons
        if(rdb.equiv < 1 | !("sm" %in% lois)){
          return(NA)
        }
        n.double.bonds <-  (rdb.equiv - 1)
        n.fa.carbons <- (c - 5)
        if(n.double.bonds > max.dbl.bnds){
          return(NA)
        }
        return(paste0("SM(",n.fa.carbons,":",n.double.bonds,")"))
        
      }
      return(NA)
    }
    if(n == 1){
      if(p == 0 && o >= 3 && o <= 5){
        #if o == 4 it could be MAG+NH4 adduct
        #ceramides have 1 rdb, base 0 carbons
        if(rdb.equiv < 1 | !("cer" %in% lois)){
          return(NA)
        }
        n.double.bonds <-  (rdb.equiv - 1)
        n.fa.carbons <- (c)
        if(n.double.bonds > max.dbl.bnds){
          return(NA)
        }
        
        if(o == 3){
          if(n.double.bonds == 0){
            return(paste0("dihydroCer(",n.fa.carbons,":",n.double.bonds,")"))
          }else{
            return(paste0("Ceramide(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }else{
          return(paste0("phytoCer(",n.fa.carbons,":",n.double.bonds,")"))
        }
      }
      if(p == 1){
        if(o >= 9 && o <= 10 ){
          #"PS has 3 rdb, base 6 carbons
          if((o == 10 && rdb.equiv < 3) | (o == 9 && rdb.equiv < 2) | !("ps" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 3)
          n.fa.carbons <- (c - 6)
          species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,10,"PS",6,lyso.carbon.cutoff.value,3)
          if(grepl("O-",species) | grepl("^L",species)){
            n.double.bonds <- n.double.bonds+1
          }
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
          
        }
        if(o>=7 && o<=8 ){
          if((c %% 2) == 0) {
            #PC has 2 rdb , base 8 carbons
            if((o == 8 && rdb.equiv < 2) | (o == 7 && rdb.equiv < 1)| !("pc" %in% lois)){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 2)
            n.fa.carbons <- (c - 8)
            species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,8,"PC",8,lyso.carbon.cutoff.value,2)
            if(grepl("O-",species) | grepl("^L",species)){
              n.double.bonds <- n.double.bonds+1
            }
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
            
          }else{
            #PE has 2 rdb , base 5 carbons
            if((o == 8 && rdb.equiv < 2) | (o == 7 && rdb.equiv < 1) | !("pe" %in% lois)){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 2)
            n.fa.carbons <- (c - 5)
            species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,8,"PE",5,lyso.carbon.cutoff.value,2)
            if(grepl("O-",species) | grepl("^L",species)){
              n.double.bonds <- n.double.bonds+1
            }
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
            
          }
        }
      }
      if(length(which_to_look_for)==1){
        return(NA)
      }
    }
  }
  return(NA)
}
