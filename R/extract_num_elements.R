#' Extract number of a given element
#' 
#'
#' This function extracts the number of a single element from each composition in a column of compositions
#' @param element_letter The element letter you would like to return the number of. Case does not matter. Allowed options are c("C","H","O","N","P","Na","Cl")
#' @param column_of_df A list (ideally, a column of a dataframe) which you'd like to extract the number of elements from
#' @keywords composition
#' @export
#' @examples
#' cat_function("C",df$composition)

extract_num_elements <- function(element_letter,column_of_df){
  sapply(column_of_df, function(x){
    if(element_letter=="N"|element_letter=="n"){
      element_letter<-"N(?!a)"
    }
    if(element_letter=="C"|element_letter=="c"){
      element_letter<-"C(?!l)"
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
