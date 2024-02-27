#' Greedily add factors to an ebcd object
#'
#' @param ebcd_obj An ebcd object.
#' @param Kmax The maximum number of factors to be added.
#' @param ebnm_fn The ebnm function to be used.
#' @param tol  The convergence tolerance parameter.
#' @param maxiter The maximum number of iterations when optimizing the greedily added factor.
#'
#' @return An ebcd object.
#' @export
#'
ebcd_greedy <- function(ebcd_obj,
                        Kmax = 1,
                        ebnm_fn = ebnm::ebnm_point_laplace,
                        tol = 1e-6,
                        maxiter = 500) {

  if(length(ebnm_fn)==1){
    ebnm_fn <- rep(list(ebnm_fn), Kmax)
  }

  for (K in 1:Kmax) {
    R <- ebcd_obj$A - tcrossprod(ebcd_obj$Z, ebcd_obj$EL)
    svd1 <- irlba::irlba(R, nv = 1, nu = 0)
    dv <- svd1$d * svd1$v

    l <- dv / sqrt(ebcd_obj$nrowA)
    Rl <- R %*% l
    z <- sqrt(ebcd_obj$nrowA) * Rl / sqrt(sum(Rl^2))

    ef1 <- list(
      A = R,
      Z = z,
      EL = l,
      maxiter = maxiter,
      ebnm_fn = ebnm_fn[[K]],
      N = ebcd_obj$N,
      nrowA = ebcd_obj$nrowA,
      tau = ebcd_obj$tau,
      tol = ebcd_obj$tol
    )
    ef1 <- ebcd_backfit(ef1)

    ebcd_obj$EL <- cbind(ebcd_obj$EL, ef1$EL)
    ebcd_obj$Z <- sqrt(ebcd_obj$nrowA) * PolarU(ebcd_obj$A %*% ebcd_obj$EL)
    ebcd_obj$tau <- ef1$tau
  }

  ebcd_obj$ebnm_fn <- ebnm_fn

  return(ebcd_obj)
}
