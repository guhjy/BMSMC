#' Probability Generating Function
#'
#' This function generates an approximation to the posterior probability corresponding to a particular partition
#'
#' @param y vector of data
#' @param grpvec pooling vector of interest as genereated by grp.asn()
#' @export
#' @examples
#' prob.compute()
prob.compute <- function(y, grpvec, ...){
  lun <- length(unique(grpvec))
  if (lun == 1) {
    mod <- glm(y ~ 1, ...) 
    bic <- BIC(mod) 
  } else {
    grpvec <- factor(grpvec)
    mod <- glm(y ~ grpvec, ...) 
    bic <- BIC(mod) 
  }
  return(bic)
}