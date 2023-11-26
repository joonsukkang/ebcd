#' Greedily add factors to an ebcd object
#'
#' @param ebcd An ebcd object.
#' @param Kmax The maximum number of factors to be added.
#' @param ebnm_fn The ebnm function to be used.
#' @param tol  The convergence tolerance parameter.
#' @param maxiter The maximum number of iterations when optimizing the greedily added factor.
#'
#' @return An ebcd object.
#' @export
#'
ebcd_greedy <- function(ebcd,
                        Kmax = 1,
                        ebnm_fn = ebnm::ebnm_point_laplace,
                        tol = 1e-6,
                        maxiter = 500) {
  for (K in 1:Kmax) {
    R <- ebcd$A - tcrossprod(ebcd$Z, ebcd$EL)
    svd1 <- irlba::irlba(R, nv = 1, nu = 0)
    dv <- svd1$d * svd1$v

    l <- dv / sqrt(ebcd$nrowA)
    Rl <- R %*% l
    z <- sqrt(ebcd$nrowA) * Rl / sqrt(sum(Rl^2))

    ef1 <- list(
      A = R,
      Z = z,
      EL = l,
      maxiter = maxiter,
      ebnm_fn = ebnm_fn,
      N = ebcd$N,
      nrowA = ebcd$nrowA,
      tau = ebcd$tau,
      tol = ebcd$tol
    )
    ef1 <- ebcd_backfit(ef1)

    ebcd$EL <- cbind(ebcd$EL, ef1$EL)
    ebcd$Z <- sqrt(ebcd$nrowA) * PolarU(ebcd$A %*% ebcd$EL)
    ebcd$tau <- ef1$tau
  }

  ebcd$ebnm_fn <- ebnm_fn

  return(ebcd)
}
