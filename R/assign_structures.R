#' Assign structures to a list of elemental compositions
#' 
#'
#' This function assigns general lipid structure (e.g. PC(32:1)) to a user-defined elemental composition (e.g. C40 H79 N O8 P), or a list/column of elemental compositions
#' @param comp list (or dataframe column) of elemental compositions for which structures are to be assigned.
#' @param ion.mode which ion mode the samples were run in. One of "neg.ion","pos.ion", or "neutral" ("neutral" searches for structures as exact compositions / neutral species).
#' @param adducts which ions/adducts to search for compositions as. Defaults to c("m.plus.h","m.plus.ammonia","m.plus.sodium") in positive ion mode, c("m.minus.h","m.plus.chloride","m.minus.2h") in negative ion mode, and the inclusive union of these two lists in neutral mode.
#' @param domain what domain of life the data originates from. "euk" for eukaryotic (default) or "bact" for bacterial.
#' @param lois which lipid classes to look for. Defaults to c("pc","sm","tag","dag","gdg") in positive ion mode, c("pe","cer","cl","pi","pg","pa","ps","ffa","gdg") in negative ion mode, and the inclusive union of these two lists in neutral mode (Diacylglycerol = "dag", Triacylglycerol = "tag", Phosphatidic acid = "pa", Phosphatidylglycerol = "pg", Phosphatidylinositol = "pi", Cardiolipin = "cl", Sphingomyelin = "sm", Ceramide = "cer", Phosphatidylcholine = "pc", Phosphatidylethanolamine = "pe", Phosphatidylserine = "ps", Free fatty acids = "ffa", Glycosyldiacylglycerols = "gdg")
#' @param max.dbl.bnds maximum number of acyl+alkyl chain double bonds allowed in returned structures. Defaults to 14.
#' @param carbon_isotope_symbol string which is used to denote carbon isotope used, if applicable. Defaults to "C\[13\]".
#' @param hydrogen_isotope_symbol string which is used to denote hydrogen isotope used, if applicable. Defaults to "D\[2\]".
#' @keywords structural assignment
#' @export
#' @examples
#' df$gen.structures <- assign_structures(comp = df$composition,
#'                                        ion.mode = "neg.ion",
#'                                        domain = "euk",
#'                                        max.dbl.bnds = 8)

assign_structures <- function(comp,
                              ion.mode,
                              adducts = NULL,
                              domain = "euk",
                              lois = NULL,
                              max.dbl.bnds = 14,
                              carbon_isotope_symbol = "C[13]",
                              hydrogen_isotope_symbol = "D[2]"){
  if(!(ion.mode %in% c("neg.ion","pos.ion","neutral")) | length(ion.mode) > 1){
    rlang::abort("problem in ion.mode argument (ion.mode must be one of neg.ion, pos.ion, or neutral)")
  }
  if(!is.null(adducts) & any(!(adducts %in% c("m.minus.h","m.plus.chloride","m.minus.2h","m.plus.h","m.plus.ammonia","m.plus.sodium")))){
    rlang::abort("one or more of inputted adducts are invalid (only valid adducts are m.minus.h,m.plus.chloride,m.minus.2h,m.plus.h,m.plus.ammonia,m.plus.sodium)")
  }
  if(!(domain %in% c("euk","bact")) | length(domain) > 1){
    rlang::abort("Invalid domain argument (domain must be one of euk or bact)")
  }
  if(!is.null(lois) & any(!(lois %in% c("pe","cer","cl","pi","pg","pa","ps","ffa","pc","sm","tag","dag")))){
    rlang::abort("one or more user-supplied input to lois argument is invalid (only valid lois are pe,cer,cl,pi,pg,pa,ps,ffa,pc,sm,tag,dag)")
  }
  if(!is.numeric(max.dbl.bnds) | length(max.dbl.bnds) > 1){
    rlang::abort("max.dbl.bnds argument not numeric")
  }
  
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
      lois <- c("pe","cer","cl","pi","pg","pa","ps","ffa","gdg")
    }
    if(ion.mode == "pos.ion"){
      lois <- c("pc","sm","tag","dag","gdg")
    }
    if(ion.mode == "neutral"){
      lois <- c("pe","cer","cl","pi","pg","pa","ps","ffa","pc","sm","tag","dag","gdg")
    }
  }

  c_tot <- (extract_num_elements_internal("C",comp) + extract_num_elements_internal(carbon_isotope_symbol,comp))
  h_tot <- (extract_num_elements_internal("H",comp) + extract_num_elements_internal(hydrogen_isotope_symbol,comp))
  n <- extract_num_elements_internal("N",comp)
  p <- extract_num_elements_internal("P",comp)
  
  result <- mapply(assign_species,
         c_tot,
         h_tot,
         extract_num_elements_internal("O",comp),
         n,
         p,
         extract_num_elements_internal("Na",comp),
         extract_num_elements_internal("Cl",comp),
         ion.mode,
         rep(list(adducts),length(c_tot)),
         ( c_tot - (h_tot/2) + ((n+p)/2) +1 ),
         domain,
         rep(list(lois),length(c_tot)),
         max.dbl.bnds,
         SIMPLIFY = T)
  if(all(is.na(result))){
    warning("No structures were assigned. Ensure that numbers of elements were extracted correctly from compositions: \n",print_and_capture(head(as.matrix(data.frame("c" = extract_num_elements_internal("C",comp),
                                                                                                                                                          "c_isotope" = extract_num_elements_internal(carbon_isotope_symbol,comp),
                                                                                                                                                          "h" = extract_num_elements_internal("H",comp),
                                                                                                                                                          "h_isotope" = extract_num_elements_internal(hydrogen_isotope_symbol,comp),
                                                                                                                                                          "o" = extract_num_elements_internal("O",comp),
                                                                                                                                                          "n" = extract_num_elements_internal("N",comp),
                                                                                                                                                          "p" = extract_num_elements_internal("P",comp),
                                                                                                                                                          "na" = extract_num_elements_internal("Na",comp),
                                                                                                                                                          "cl" = extract_num_elements_internal("Cl",comp)
                                                                                                                                                          )))))
    result
  }else{
    result
  }

}

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
        return(structure_from_exact_comp(c,h,o,n,p,c("nonacylglycerols","acylglycerols","galactosyldiacylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
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
      # regardless if find one or not, return the result here
    }
    return(NA)
    # if we make it here, return some clever paste0 showing what we searched for, beginning with 'Species not detected (traceback:'
    
  }
  if(ion.mode == "neg.ion"){
    if("m.plus.chloride" %in% adducts){
      if(cl > 0){
        return(structure_from_exact_comp(c,h,o,n,p,c("nonacylglycerols","galactosyldiacylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
        # look for chloride adduct species
        # assign exact composition by changing nothing
        # regardless if find one or not, return the result here
      }
      # else, keep going
    }
    if("m.minus.h" %in% adducts){
      result <- structure_from_exact_comp(c,h+1,o,n,p,c("nonacylglycerols","galactosyldiacylglycerols"),domain,lois,as.numeric(max.dbl.bnds))
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
  
  if(rdb.equiv %% 1 != 0 ){
    return(NA)
  }
  
  
  if( "acylglycerols" %in% which_to_look_for){
    if(n == 0 && p == 0 && o == 5 && c >= 35){
      #DAG
      #DAG has base 2 RDB , base 3 carbons
      if(rdb.equiv < 2 | !("dag" %in% lois)){
        return(NA)
      }
      n.double.bonds <-  (rdb.equiv - 2)
      n.fa.carbons <- (c - 3)
      if(n.double.bonds > max.dbl.bnds){
        return(NA)
      }
      return(paste0("DAG(",n.fa.carbons,":",n.double.bonds,")"))
    }
    if(n == 0 && p == 0 && o == 6 && c >= 49){
      #TAG
      #TAG has base 3 RDB , base 3 carbons
      if(rdb.equiv < 3 | !("tag" %in% lois)){
        return(NA)
      }
      n.double.bonds <-  (rdb.equiv - 3)
      n.fa.carbons <- (c - 3)
      if(n.double.bonds > max.dbl.bnds){
        return(NA)
      }
      return(paste0("TAG(",n.fa.carbons,":",n.double.bonds,")"))
    }
    if(length(which_to_look_for)==1){
      return(NA)
    }
  }
  if("galactosyldiacylglycerols" %in% which_to_look_for){
    if(n == 0 && p == 0 && c > 29 && (o == 10 | o == 15)){
      if(o == 15){
        if(rdb.equiv < 4 | !("gdg" %in% lois)){
          return(NA)
        }
        n.double.bonds <-  (rdb.equiv - 4)
        n.fa.carbons <- (c - 15)
        if(n.double.bonds > max.dbl.bnds){
          return(NA)
        }
        return(paste0("DGDG(",n.fa.carbons,":",n.double.bonds,")"))
      }
      if(o == 10){
        if(rdb.equiv < 3 | !("gdg" %in% lois)){
          return(NA)
        }
        n.double.bonds <-  (rdb.equiv - 3)
        n.fa.carbons <- (c - 9)
        if(n.double.bonds > max.dbl.bnds){
          return(NA)
        }
        return(paste0("MGDG(",n.fa.carbons,":",n.double.bonds,")"))
      }
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
        #some FFA 
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
        }
        if(o == 4){
          return(paste0("phytoCer(",n.fa.carbons,":",n.double.bonds,")"))
        }
        if(o == 5){
          return(paste0("phytoCer(h2,",n.fa.carbons,":",n.double.bonds,")"))
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

extract_num_elements_internal <- function(element_letter,column_of_df){
  sapply(column_of_df, function(x){
    if(element_letter=="N"|element_letter=="n"){
      element_letter<-"(?<!I|M|Z|S|R)N(?!a|s|d|e|p|i|b|o)"
    }
    if(element_letter=="C"|element_letter=="c"){
      element_letter<-"(?<!A|T|S)C(?!a|d|f|e|s|l|r|o|u|m)"
    }
    if(element_letter=="H"|element_letter=="h"){
      element_letter<-"(?<!T|R)H(?!g|f|s|e|o)"
    }
    if(element_letter=="O"|element_letter=="o"){
      element_letter<-"(?<!H|C|P|N)O(?!s)"
    }
    if(element_letter=="P"|element_letter=="p"){
      element_letter<-"(?<!N)P(?!b|d|t|u|o|r|m|a)"
    }
    
    #if there is a square bracket (i.e. if it's an isotope), substitute in the double brackets so it works with regex below
    if(grepl("\\[",element_letter)){
      element_letter<-gsub("\\[","\\\\[",element_letter)
      element_letter<-gsub("\\]","\\\\]",element_letter)
    }
    
    num <- ifelse(!grepl(element_letter,x,perl = T,ignore.case = T),
                  #paste0(element_letter,"(?=\\s)|",element_letter,"(?=[0-9]+)|",element_letter,"$")
                  0,  #return zero
                  ifelse(!grepl(paste0("(?<=",element_letter,")[0-9]+|(?<=",element_letter,"\\s)[0-9]+"),x,perl = T,ignore.case = T), #else, 
                         1,
                         as.numeric(regmatches(x,regexpr(paste0("(?<=",element_letter,")[0-9]+|(?<=",element_letter,"\\s)[0-9]+"),x,perl = T,ignore.case = T)))))
    return(num)
  }
  )
}
print_and_capture <- function(x){
  paste(capture.output(print(x)), collapse = "\n")
}