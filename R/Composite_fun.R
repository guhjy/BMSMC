#' Composite Bayesian Multiple Comparisons Function
#'
#' This function performs the Baysian Model Selection approach to Multiple Comparisons Problems 
#'
#' @param y vector of data
#' @param npg vector describing the number of observations in each group
#' @param pairs k by 2 matrix of k pairs of interest
#' @param cutoff the cutoff value for display of model posterior probabilities 
#' @param p.probs Do you want a table of posterior probabilities for all generated models? Default is FALSE 
#' @export
#' @examples
#' BayMC.composite()
BayMC.composite <- function(y, npg, pairs = NA, cutoff = NA, p.probs = FALSE, ...) {
  totobs  <- length(y)
  maxgrps <- length(npg) 
  parts <- setparts(maxgrps)
  
  
  allpt.vec <- apply(parts, 2, grp.asn, npg = npg)
  allpt.pool <- parts

  modbics <- apply(allpt.vec, 2, prob.compute, y = y, ...)
  
  delta_modbics <- modbics - min(modbics)
  preprobs <- exp(-0.5*delta_modbics)
  probs <- preprobs/sum(preprobs)
  mod.bic <- modbics
  probs <- round(probs, 4)
  
  modcoefs <- matrix(-111, nrow = maxgrps, ncol = dim(parts)[2])
  bmest.pre <- matrix(-111, nrow = maxgrps, ncol = dim(parts)[2])
  for(i in 1:dim(parts)[2]) {
    modcoefs[,i] <- coef.compute(y = y, grpvec = allpt.vec[,i], pooling = allpt.pool[,i], ...)
    bmest.pre[,i] <- modcoefs[,i]*probs[i]
  }
  bmacoefs <- rowSums(bmest.pre)
  
  
  if(!all(is.na(pairs))){
    pw.bmaprob <- NULL
    plab <- NULL
    for (j in 1:dim(pairs)[1]){
      pairvec <- NULL
      for (i in 1:dim(parts)[2]){
        pairvec[i] <- (allpt.pool[pairs[j,1], i] == allpt.pool[pairs[j,2], i]) 
      }
      pw.bmaprob[j] <- sum(probs[pairvec])
      plab[j] <- paste(pairs[j,], collapse = "")
    }
    
    bmadf <- data.frame(Pairing = plab,
                        Prob    = pw.bmaprob)
    bmadf <- bmadf[order(bmadf[,2], decreasing = TRUE),]
  }
  
  if(all(is.na(pairs))){
    pairs <- t(combn(1:maxgrps,2))
    
    pw.bmaprob <- NULL
    plab <- NULL
    for (j in 1:dim(pairs)[1]){
      pairvec <- NULL
      for (i in 1:dim(parts)[2]){
        pairvec[i] <- (allpt.pool[pairs[j,1], i] == allpt.pool[pairs[j,2], i]) 
      }
      pw.bmaprob[j] <- sum(probs[pairvec])
      plab[j] <- paste(pairs[j,], collapse = "")
    }
    
    bmadf <- data.frame(Pairing = plab,
                        Prob    = pw.bmaprob)
    bmadf <- bmadf[order(bmadf[,2], decreasing = TRUE),]
  }
  
  
  char.mods <- apply(allpt.pool,2,paste,collapse = "")
  
  postprobs <- data.frame(model = char.mods, probs = probs, bic = mod.bic)
  postprobs <- postprobs[order(postprobs[,2], decreasing = TRUE),]
  rownames(postprobs) <- paste("M", 1:dim(parts)[2], sep = "")
  
  if (!is.na(cutoff)){
    postprobs <- postprobs[postprobs$probs > cutoff, ]
  }
  
  results <- list(p.pairprobs = bmadf,
                  BMA.coefs = bmacoefs)
  
  if (p.probs == TRUE){
    results <- list(p.probs = postprobs,
                    p.pairprobs = bmadf,
                    BMA.coefs = bmacoefs)
    return(results)
  } else{
    return(results)
  }
}