#' Composite Bayesian Multiple Comparisons Function
#'
#' This function performs analyses as described in "Bayesian Multiple Comparisons and Model Selection" by Andrew A. Neath,
#' Javier E. Flore, and Joseph E. Cavanaugh. Specifically, this function computes approximate posterior probabilities for all
#' models arising from partitions of a set of I elements. Additionally, posterior pairwise probabilities are computed for
#' pairwise comparisons of interest.
#'
#' @param y vector of data
#' @param npg vector describing the number of observations in each group
#' @param pairs k by 2 matrix of k pairs of interest
#' @param cutoff the cutoff value for display of model posterior probabilities
#' @param p.probs Do you want a table of posterior probabilities for all generated models? Default is FALSE
#' @param ... Specification of the family and link function as in the glm procedure
#' @export
#' @examples
#' Consider the following data:
#' Group 1: 12,15,19,11
#' Group 2: 4,7,9,3
#' Group 3: 18,11,17
#' Group 4: 22,24,23,29,27
#'
#' We first define our observation vector using these data.
#' ex.y <- c(12,13,19,11, 4,7,9,3, 18,11,17, 22,24,23,29,27)
#'
#' Next, we define the vector specifying the number per group.
#' ex.npg <- c(4,4,3,5)
#'
#' We are interested in obtaining the posterior pairwise equality probability for the means of groups 1 and 3.
#' The posterior pairwise equality probability for the means of groups 3 and 4 is also of interest. We next define
#' the matrix which specifies this. (Note: Should interest be in all pairwise equalities, do not provide an object to pairs)
#' ex.pairs <- matrix(c(1,3,3,4), nrow = 2, ncol = 2, byrow = TRUE)
#'
#' Feeding these objects into the function, we obtain the desired results. We specify "p.probs = TRUE" in order to obtain
#' a table of posterior probabilities for all models defined by each partition. A restricted table displaying models with
#' posterior probabilites larger than some cutoff, c, may be obtained by specifying "cutoff = c". The models fit are ANOVA type
#' models, so we specify the link function accordingly.
#'
#' BayMC.composite(y = ex.y, npg = ex.npg, pairs = ex.pairs, p.probs = TRUE, family=gaussian(link = "identity")))
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
