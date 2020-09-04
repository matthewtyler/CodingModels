### Author: Matthew Tyler
### MDY Date: 08-10-20

calc_logliks_BACE <- function(pars, dat) {
    lalpha <- log(pars$alpha)
    lbeta <- log(pars$beta)
    lgamma <- log(pars$gamma)

    ll_dict <- array(NA, dim = c(dat$J, dat$K, dat$K)) # loglik that coder j outputs k given that the Z == k
    for (j in 1:dat$J) { # coder j
      for (k in 1:dat$K) { # Z == k
        a <- (1 - pars$beta[j]) * pars$gamma[j, ]; # clear == 0
        b <- pars$beta[j]
        ll_dict[j, k, k] <- log(a[k] + b)
        ll_dict[j, k, -k] <- log(a[-k])
      }
    }

    lpy <- array(NA, dim = c(dat$I, dat$K))
    for (k in 1:dat$K) lpy[, k] <- lalpha[k]

    for (n in 1:dat$N) {
      i <- dat$ii[n]
      j <- dat$jj[n]
      val <- dat$yy[n]
      lpy[i, ] <- lpy[i, ] + ll_dict[j, , val]
    }

    post_Z <- array(NA, dim = c(dat$I, dat$K))
    logliks <- array(NA, dim = dat$I)
    for (i in 1:dat$I) {
      logliks[i] <- logSumExp(lpy[i, ])
      post_Z[i, ] <- softmax(lpy[i, ])
    }

    res <- list(logliks = logliks, post_Z = post_Z, ll_dict = ll_dict)
    return(res)
  }

pars_map_BACE <- function(pars, dat, prior = 1.0) {

  logliks <- calc_logliks_BACE(pars, dat)
  post_Z <- logliks$post_Z
  p_dict <- exp(logliks$ll_dict)
  logliks <- logliks$logliks

  alpha <- rep(prior, dat$K)
  for (i in 1:dat$I)
    alpha <- alpha + post_Z[i,]
  alpha <- alpha / sum(alpha)

  post_C_1 <- array(NA, dim = c(dat$I, dat$J))
  for (n in 1:dat$N) {
    i <- dat$ii[n]
    j <- dat$jj[n]
    val <- dat$yy[n]
    post_C_1[i, j] <- pars$beta[j] * post_Z[i, val] / p_dict[j, val, val]
    # because the Z[ii[n]] != yy[n] components contribute 0 mass
    # note that
    # post_Z[i, val] / p_dict[j, val, val]
    # = p(z = k | x_j = k) / p(x_j = k | z = k)
    # = p(z = k) / p(x = k) [multiply num. and denom. by p(z)*p(x)]
    # = p(z = k | c = 1) / p(x = k)
    # = p(x = k | z = k, c = 1) p(z = k | c = 1) / p(x = k)
    # = p(x = z = k | c = 1) / p(x = k)
    # = p(x = k | c = 1) / p(x = k)
  }
  post_C_0 <- 1 - post_C_1

  beta_1 <- beta_0 <- array(prior, dim = dat$J)  
  for (n in 1:dat$N) {
    i <- dat$ii[n]
    j <- dat$jj[n]
    beta_1[j] <- beta_1[j] + post_C_1[i, j]
    beta_0[j] <- beta_0[j] + post_C_0[i, j]
  }
  beta <- beta_1 / (beta_1 + beta_0)

  gamma <- array(prior, dim = c(dat$J, dat$K))
  for (n in 1:dat$N) {
    i <- dat$ii[n]
    j <- dat$jj[n]
    val <- dat$yy[n]
    gamma[j, val] <- gamma[j, val] + post_C_0[i, j]
  }
  for (j in 1:dat$J) gamma[j, ] <- gamma[j, ] / sum(gamma[j, ])

  mean_LL <- mean(logliks)

  pars <- list(alpha = alpha, beta = beta, gamma = gamma,
    mean_LL = mean_LL)
  return(pars)
}