### Author: Matthew Tyler
### MDY Date: 08-11-20

calc_logliks_DS <- function(pars, dat) {
  lalpha <- log(pars$alpha)
  lbeta <- log(pars$beta)

  post_Z <- array(NA, dim = c(dat$I, dat$K))
  logliks <- array(NA, dim = dat$I)

  for (i in 1:dat$I) {
    post_Z[i, ] <- lalpha
  }
  for (n in 1:dat$N) {
    post_Z[dat$ii[n], ] <- post_Z[dat$ii[n], ] + lbeta[dat$jj[n], , dat$yy[n]]
  }
  for (i in 1:dat$I) {
    logliks[i] <- logSumExp(post_Z[i, ])
    post_Z[i, ] <- softmax(post_Z[i, ])
  }

  res <- list(logliks = logliks, post_Z = post_Z)
  return(res)
}

pars_map_DS <- function(pars, dat, prior = 1.0) {
  logliks <- calc_logliks_DS(pars, dat)
  post_Z <- logliks$post_Z
  logliks <- logliks$logliks

  alpha <- rep(prior, dat$K)
  for (i in 1:dat$I)
    alpha <- alpha + post_Z[i,]
  alpha <- alpha / sum(alpha)

  beta <- array(prior, c(dat$J, dat$K, dat$K))
  for (n in 1:dat$N)
    for (k in 1:dat$K)
      beta[dat$jj[n], k, dat$yy[n]] <- beta[dat$jj[n], k, dat$y[n]] +
        post_Z[dat$ii[n], k]
  
  for (j in 1:dat$J)
    for (k in 1:dat$K)
      beta[j, k, ] <- beta[j, k, ] / sum(beta[j, k, ]);

  mean_LL <- mean(logliks)

  pars <- list(alpha = alpha, beta = beta, mean_LL = mean_LL)
  return(pars)
}