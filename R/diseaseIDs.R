#' Identify disease cases within the UK Biobank.
#'
#' @description
#' In the UK Biobank (UKBB) diseases are reported based on either in-hospital records (ICD-10 classification system; field id 41270),
#' or based on self-reported diseases (field id 20001 or 20002).
#'
#' The reported diseases are listed in large matrices where many of the cells are filled with NA's.
#' This utility function enables easy extraction of disease cases based on a vector of ICD-10 codes and/or a vector of self-reported diseases.
#'
#' This function requires two matrices 'self' and 'icd10', which are matrices with self-reported and ICD-10 diagnoses,
#' where each row correspond to one individual and the columns correspond to the array number (see UKBB showcase for details).
#'
#' @param ICD a vector of ICD-10 codes
#' @param SR a vector of self-reported disease codes.
#'
#' @author Palle Duun Rohde
#' @export
#'
getID <- function(ICD=NULL, SR=NULL){
    if(!is.null(SR)){
        srID <- apply(self, 2, function(x){
            case <- names(na.omit(x[x==SR]))
            return(case)
           })
        srID <- unique(unname(unlist(srID)))
    }
    if(!is.null(ICD)){
        icdID <- NULL
        for(i in 1:length(ICD)){
            tmp <- apply(icd10, 2, function(x){
            case <- names(na.omit(x[grep(ICD[i],x)]))
            return(case)
            })
            icdID <- c(icdID, unique(unname(unlist(tmp))))
        }
        icdID <- unique(icdID)
        }
    if(!is.null(SR)&!is.null(ICD)){ids <- unique(c(srID, icdID))}
    if(!is.null(SR)&is.null(ICD)){ids <- unique(c(srID))}
    if(is.null(SR)){ids <- unique(c(icdID))}
    return(ids)
    }
