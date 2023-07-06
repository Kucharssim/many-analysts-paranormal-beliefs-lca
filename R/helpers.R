repeatOptim <- function(n = 10, parallel = TRUE, ...) {
  require(rstan)
  if (parallel == FALSE || .Platform$OS.type %in% c("Windows", "windows")) {
    repeater <- lapply
  } else {
    require(parallel)
    require(pbmcapply)
    repeater <- pbmcapply::pbmclapply
  }

  cat("Fitting the model...\n")
  fits <- repeater(seq_len(n), function(seed) {
    dots <- list(...)
    dots$as_vector <- FALSE
    dots$draws <- 0
    dots$seed <- seed
    dots$verbose <- FALSE
    f <- do.call(rstan::optimizing, dots)
    cat("Fit", seed, "done!\n")
    return(f)
  })
  return_code <- sapply(fits, \(f) f$return_code)
  if(all(return_code != 0L))
    warning("None of the fits converged normally!")

  logLik <- sapply(fits, \(f) f$par$logLik)

  fake_logLik <- logLik
  fake_logLik[return_code != 0] <- -Inf
  best <- which.max(fake_logLik)

  fit <- fits[[best]]
  init <- fit$par

  cat("Refitting, sampling...\n")
  fit <- rstan::optimizing(..., seed = best, init = init)
  fit$refits <- list(
    logLik = logLik,
    return_code = return_code
  )

  return(fit)
}



inv_logit <- function (x) {
  1/(1 + exp(-x))
}
