obs2lia <- function(h2=NULL,K=NULL, P=NULL){
  # h2: observed SNP heritability
  # K: Population disease prevalence
  # P: Proportion of cases in the sample
  X <- qnorm(K,lower.tail=FALSE)
  z <- (1/sqrt(2*pi))*(exp(-(X**2)/2))
  h2lia <- h2*(K*(1-K)*K*(1-K))/(P*(1-P)*(z**2))
  return(h2lia)
}
