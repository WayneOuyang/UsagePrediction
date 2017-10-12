Forecast_comb <- function(obs, fhat, fhat_new= NULL, decay_factor=0.9, interval=26,
                          Averaging_scheme=c("simple", "ols", "robust", "cls", "variance based", "best")) {
  pckg = c("quantreg", "quadprog")
  temp <- unlist(lapply(pckg, require, character.only=T))
  if (!all(temp==1) ) {
    stop("This function relies on packages \"quadprog\" and \"quantreg\".
         Use ?install.packages if they are not yet installed. \n")
  }
  
  if(1 < decay_factor || decay_factor <= 0){
    stop("Warning: decay_factor is the statistical decay factor, and is between 0 and 1")
  }
  
  mat_err <- apply(fhat, 2, function(x) obs - x)
  TT <- NROW(fhat)
  TT_new <- NROW(fhat_new)
  p <- NCOL(fhat)
  pred <- NULL
  sq_er <- function(obs, pred) { mean( (obs - pred)^2 )  }
  
  if(length(Averaging_scheme) != 1) {
    stop("Pick only one of the following:
         c(\"simple\", \"ols\", \"robust\", \"variance based\", \"cls\", \"best\")")
  }
  
  ## Different forecast averaging schemes
  if(Averaging_scheme== "simple") {
    pred <- apply(fhat, 1, mean)
    weights <- matrix( 1/p, nrow = 1, ncol = p)
    if (!is.null(fhat_new)) { pred_new <- apply(fhat_new, 1, mean) }
  } else if (Averaging_scheme== "ols") {
    decay_list <- decay_factor^((TT-1):0)
    weights <- lm(obs ~ fhat, weights = decay_list)$coef
    pred <- t(weights %*% t(cbind(rep(1, TT), fhat)))
    if (!is.null(fhat_new)) { pred_new <- t(weights %*% t(cbind(rep(1, TT_new), fhat_new))) }
  } else if (Averaging_scheme== "cls") {
    cls1 = function(y, predictions){
      Rinv <- solve(chol(t(predictions) %*% predictions))
      C <- cbind(rep(1, NCOL(predictions)), diag(NCOL(predictions)))
      b = c(1, rep(0, NCOL(predictions)))
      d = t(y) %*% predictions
      qp1 = solve.QP(Dmat= Rinv, factorized= TRUE, dvec= d, Amat= C, bvec = b, meq = 1)
      weights = qp1$sol
      yhatcls = t(weights %*% t(predictions))
      list(yhat= yhatcls, weights= weights)
    }
    weights <- cls1(obs, fhat)$weights
    pred <- t(weights %*% t(fhat))
    if (!is.null(fhat_new)) { pred_new <- t(weights %*% t(fhat_new)) }
  } else if (Averaging_scheme== "robust") {
    decay_list <- decay_factor^((TT-1):0)
    weights <- rq(obs ~ fhat, weights = decay_list)$coef
    pred <- t(weights %*% t(cbind(rep(1,TT), fhat)))
    if (!is.null(fhat_new)) { pred_new <- t(weights %*% t(cbind(rep(1, TT_new), fhat_new))) }
  } else if (Averaging_scheme=="variance based") {
    decay_list <- decay_factor^((TT-1):0)
    temp = apply(mat_err^2, 2, weighted.mean, w=decay_list)/sum(apply(mat_err^2, 2, weighted.mean, w=decay_list))
    weights <- (1/temp)/sum(1/temp)
    pred <- t(weights %*% t(fhat))
    if ( !is.null(fhat_new) ) { pred_new <- t(weights %*% t(fhat_new)) }
  } else if (Averaging_scheme== "best") {
    temp <- apply(fhat, 2, sq_er, obs= obs)
    weights <- rep(0, p)
    weights[which.min(temp)] <- 1
    pred <- t(weights %*% t(fhat))
    if (!is.null(fhat_new)) { pred_new <- t(weights %*% t(fhat_new)) }
  }
  if (is.null(fhat_new)) { pred_new <- NULL }
  return( list(fitted= pred, pred = pred_new, weights = weights ) )
  }