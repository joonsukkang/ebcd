#' Backfit an ebcd object
#'
#' @param ebcd_obj An ebcd object.
#' @param tol The convergence tolerance parameter.
#' @param maxiter The maximum number of iterations for the backfit procedure.
#'
#' @return An ebcd object.
#' @export
#'
ebcd_backfit <- function(ebcd_obj,
                         tol = 1e-6,
                         maxiter = 5000) {
  Kmax <- ncol(ebcd_obj$Z)
  ebcd_obj$KL <- rep(0, length = Kmax)
  ebcd_obj$obj.old <- -Inf
  ebcd_obj$vec.obj <- c()

  for (iter in 1:maxiter) {
    # Shrinkage step
    ebcd_obj$V <- matrix(0, nrow = nrow(ebcd_obj$EL), ncol = ncol(ebcd_obj$EL))

    for (k in 1:Kmax) {
      ebnm_fn <- ebcd_obj$ebnm_fn
      x <- c(crossprod(ebcd_obj$A, ebcd_obj$Z[, k])) / ebcd_obj$nrowA
      s <- rep(sqrt(1 / (ebcd_obj$N * ebcd_obj$tau)), times=length(x))
      e <- ebnm_fn(x = x, s = s)

      ebcd_obj$EL[, k] <- e$posterior$mean
      ebcd_obj$V[, k] <- e$posterior$sd^2
      ebcd_obj$KL[k] <- e$log_likelihood +
        - normal_means_loglik(x, s, ebcd_obj$EL[, k], ebcd_obj$EL[, k]^2 + ebcd_obj$V[, k])
    }

    # Rotation step
    ebcd_obj$Z <- sqrt(ebcd_obj$nrowA) * PolarU(ebcd_obj$A %*% ebcd_obj$EL)

    # Precision step
    ebcd_obj$tau <- prod(dim(ebcd_obj$A)) / (sum((ebcd_obj$A - ebcd_obj$Z %*% t(ebcd_obj$EL))^2) + ebcd_obj$nrowA * sum(ebcd_obj$V))

    # check convergence
    ebcd_obj$obj <- -ebcd_obj$N * ncol(ebcd_obj$A) / 2 * log(2 * pi / ebcd_obj$tau) +
      -(ebcd_obj$N * ebcd_obj$tau / 2) * (
        sum(ebcd_obj$A^2) / ebcd_obj$nrowA - 2 * sum(diag(t(ebcd_obj$A) %*% ebcd_obj$Z %*% t(ebcd_obj$EL))) / ebcd_obj$nrowA + sum(ebcd_obj$EL^2) + sum(ebcd_obj$V)
      ) +
      +sum(ebcd_obj$KL)

    ebcd_obj$vec.obj <- c(ebcd_obj$vec.obj, ebcd_obj$obj)
    if (iter >= 10 & ebcd_obj$obj - ebcd_obj$obj.old < tol) break
    ebcd_obj$obj.old <- ebcd_obj$obj
  }

  return(ebcd_obj)
}
