#' Extract the number of a given element from elemental compositions
#' 
#'
#' This function extracts the number of a single element (e.g. carbon, C) from each composition (e.g. C43 H67 O8 N P) in a column of compositions
#' @param element_letter The element letter you would like to extract the number of. Case does not matter. Supported options are c("C","H","O","N","P","Na","Cl","\[13C\]","\[2H\]"), but the regex should be pretty flexible--just double check your results if using something other than these 9 options explicitly.
#' @param column_of_df A list (e.g. a column of a dataframe) of compositions which you'd like to extract the number of elements from.
#' @keywords composition
#' @export
#' @examples
#' df$c <- extract_num_elements("C",df$composition)
#' df$c_isotope <-  extract_num_elements("[13C]",df$composition)
#' df$h <- extract_num_elements("H",df$composition)
#' df$h_isotope <- extract_num_elements("[2H]",df$composition)
#' df$o <- extract_num_elements("O",df$composition)
#' df$n <- extract_num_elements("N",df$composition)
#' df$p <- extract_num_elements("P",df$composition)
#' df$na <- extract_num_elements("Na",df$composition)
#' df$cl <- extract_num_elements("Cl",df$composition)

extract_num_elements <- function(element_letter,column_of_df){
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
