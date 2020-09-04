### Author: Matthew Tyler
### Date (mdy): 08-10-20

softmax <- function(x) {
  M <- max(x)
  x <- x - M
  ex <- exp(x)
  return(ex / sum(ex))
}

to_one_hot <- function(Y, K) sapply(1:K, function(u) 1 * (Y == u))


extract_vectors <- function(long_labels) {
  ii <- as.integer(long_labels$ii)
  jj <- as.integer(long_labels$jj)
  yy <- as.integer(long_labels$yy)

  N <- length(ii)  
  stopifnot(identical(N, length(jj)))
  stopifnot(identical(N, length(yy)))

  stopifnot(sum(is.na(ii)) == 0)
  stopifnot(sum(is.na(jj)) == 0)

  ii[is.na(yy)] <- NA
  jj[is.na(yy)] <- NA
  ii <- na.omit(ii)
  jj <- na.omit(jj)
  yy <- na.omit(yy)

  return(list(ii = ii, jj = jj, yy = yy))
}

#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @import matrixStats
#' @export
#' @examples
#' cat_function()
fit_em <- function(long_labels, model = "BACE",
  max_iter = 1e3, tol = 1e-8,
  verbose = FALSE) {
  # long_labels should be a data frame with columns:
  # ii: object index
  # jj: coder index
  # yy: coder label
  # because it's in long format, there shouldn't be any missing data

  vectors <- extract_vectors(long_labels)
  ii <- vectors$ii
  jj <- vectors$jj
  yy <- vectors$yy

  N <- length(ii)
  I <- max(ii); stopifnot(identical(sort(unique(ii)), 1:I))
  J <- max(jj); stopifnot(identical(sort(unique(jj)), 1:J))
  K <- max(yy); stopifnot(identical(sort(unique(yy)), 1:K))

  dat <- list(N = N, I = I, J = J, K = K,
    ii = ii, jj = jj, yy = yy)

  if (model == "BACE") {
    a <- rep(1/K, K)
    b <- array(0.8, dim = J)
    g <- array(NA, dim = c(J, K))
    for (j in 1:J) {
      g[j, ] <- rep(1 / K, K)
    }
    pars <- list(alpha = a, beta = b, gamma = g)

    num_params <- (K - 1) + J + J * (K - 1)
    pars_map <- pars_map_BACE
    calc_logliks <- calc_logliks_BACE

  } else if (model == "MACE") {
    a <- rep(1/K, K)
    b <- array(0.8, dim = J)
    pars <- list(alpha = a, beta = b)

    num_params <- (K - 1) + J
    pars_map <- pars_map_MACE
    calc_logliks <- calc_logliks_MACE   

  } else if (model == "DS") {

    a <- rep(1/K, K)
    b <- array(NA, dim = c(J, K, K))
    for (j in 1:J) {
      for (k in 1:K) {
        b[j, k, k] <- 0.9 # optimistic starting value to suggest good coders
        b[j, k, -k] <- (1 - b[j, k, k]) / (K - 1)
      }
    }
    pars <- list(alpha = a, beta = b)

    num_params <- (K - 1) + J * K * (K - 1)
    pars_map <- pars_map_DS
    calc_logliks <- calc_logliks_DS   

  } else {
    stopifnot(model %in% c("BACE", "DS", "MACE"))
  }
 
  old_mean_LL <- -Inf
  mean_LL <- array(NA, dim = max_iter)
  
  for (iter in 1:max_iter) {
    pars <- pars_map(pars, dat)
    
    new_LL <- pars$mean_LL
    mean_LL[iter] <- new_LL
    change_LL <- new_LL - old_mean_LL
    old_mean_LL <- new_LL
    
    if (abs(change_LL) < tol) {
      message(paste0("CONGERGED @ ", iter, " iterations!"))
      break
    } else {
      if (verbose) print(paste0("Iter: ", iter, ", LL: ", round(mean_LL[iter], -log10(tol))))
    }
  }

  if (iter == max_iter) {
    message("Maximum iterations reached before convergence :(")
  }

  mean_LL_path <- mean_LL[1:iter]
  mean_LL <- mean_LL_path[iter]
  
  logliks <- calc_logliks(pars, dat)
  post_Z <- logliks$post_Z

  sum_LL <- sum(logliks$logliks)
  AIC <- 2 * (num_params - sum_LL)
  AIC_norm <- AIC / (-2 * I) # alt., (sum_ll - num_params) / I

  res <- list(pars = pars,
    post_Z = post_Z, mean_LL_path = mean_LL_path,
    AIC = AIC, AIC_norm = AIC_norm,
    num_params = num_params, sum_LL = sum_LL, mean_LL = mean_LL)

  return(res)
}