#' Convert observed heritability to the liability scale
#'
#' @description
#' Whether the estiamted heritability is based on GREML or on summary statistics, if the trait is a disease trait,
#' the estimate may be biased due to assertainment bias. To account for this, one can convert the observed heritability
#' to the heritability on the liability scale by correcting for population disease prevalence and the proportion of cases in the sample
#'
#' @param h2 observed SNP heritability
#' @param K Population disease prevalence
#' @param P Proportion of cases in the sample
#'
#' @author Palle Duun Rohde
#'
#' @export
#'
obs2lia <- function(h2=NULL,K=NULL, P=NULL){
  X <- qnorm(K,lower.tail=FALSE)
  z <- (1/sqrt(2*pi))*(exp(-(X**2)/2))
  h2lia <- h2*(K*(1-K)*K*(1-K))/(P*(1-P)*(z**2))
  return(h2lia)
}
