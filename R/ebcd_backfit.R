#' Title
#'
#' @param ebcd An ebcd object.
#' @param tol The convergence tolerance parameter.
#' @param maxiter The maximum number of iterations for the backfit procedure.
#'
#' @return An ebcd object.
#' @export
#'
ebcd_backfit <- function(ebcd,
                         tol = 1e-6,
                         maxiter = 5000) {
  Kmax <- ncol(ebcd$Z)
  ebcd$KL <- rep(0, length = Kmax)
  ebcd$obj.old <- -Inf
  ebcd$vec.obj <- c()

  for (iter in 1:maxiter) {
    # Shrinkage step
    ebcd$V <- matrix(0, nrow = nrow(ebcd$EL), ncol = ncol(ebcd$EL))

    for (k in 1:Kmax) {
      ebnm_fn <- ebcd$ebnm_fn
      x <- c(crossprod(ebcd$A, ebcd$Z[, k])) / ebcd$nrowA
      s <- sqrt(1 / (ebcd$N * ebcd$tau))
      e <- ebnm_fn(x = x, s = s)

      ebcd$EL[, k] <- e$posterior$mean
      ebcd$V[, k] <- e$posterior$sd^2
      ebcd$KL[k] <- e$log_likelihood +
        -flashier:::normal.means.loglik(x, s, ebcd$EL[, k], ebcd$EL[, k]^2 + ebcd$V[, k])
    }

    # Rotation step
    ebcd$Z <- sqrt(ebcd$nrowA) * polar(ebcd$A %*% ebcd$EL)

    # Precision step
    ebcd$tau <- prod(dim(ebcd$A)) / (sum((ebcd$A - ebcd$Z %*% t(ebcd$EL))^2) + ebcd$nrowA * sum(ebcd$V))

    # check convergence
    ebcd$obj <- -ebcd$N * ncol(ebcd$A) / 2 * log(2 * pi / ebcd$tau) +
      -(ebcd$N * ebcd$tau / 2) * (
        sum(ebcd$A^2) / ebcd$nrowA - 2 * sum(diag(t(ebcd$A) %*% ebcd$Z %*% t(ebcd$EL))) / ebcd$nrowA + sum(ebcd$EL^2) + sum(ebcd$V)
      ) +
      +sum(ebcd$KL)

    ebcd$vec.obj <- c(ebcd$vec.obj, ebcd$obj)
    if (iter >= 10 & ebcd$obj - ebcd$obj.old < ebcd$tol) break
    ebcd$obj.old <- ebcd$obj
  }

  return(ebcd)
}
