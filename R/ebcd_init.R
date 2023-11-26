#' Initialize ebcd object
#'
#' Initialize an ebcd object with no factors, taking as an input one of
#' data matrix \eqn{X}, covariance matrix \eqn{S}, or compact target matrix \eqn{C}.
#' The algorithm runs on either \eqn{X} or \eqn{C}, which is internally saved as \eqn{A}.
#'
#' @param X N-by-P data matrix with N observations and P variables
#' @param S P-by-P covariance/Gram matrix \eqn{=X'X/N}
#' @param C P-by-P compact target matrix such that \eqn{C'C=X'X}
#' @param N number of observations
#' @param tol tolerance
#'
#' @return An initialized ebcd object
#' @export
#'
ebcd_init <- function(X = NULL,
                      S = NULL,
                      C = NULL,
                      N = NULL,
                      tol = 1e-6) {
  if (!is.null(C)) {
    A <- C
  } else if (!is.null(S)) {
    eigS <- eigen(S)
    C <- sqrt(nrow(S)) * eigS$vectors %*% diag(sqrt(pmax(eigS$values, 0))) %*% t(eigS$vectors)
    A <- C
  } else if (!is.null(X)) {
    N <- nrow(X)

    if (nrow(X) > ncol(X)) {
      S <- crossprod(X) / N
      eigS <- eigen(S)
      C <- sqrt(nrow(S)) * eigS$vectors %*% diag(sqrt(pmax(eigS$values, 0))) %*% t(eigS$vectors)
      A <- C
    } else {
      A <- X
    }
  }

  if (is.null(N)) {
    N <- nrow(A)
  }

  nrowA <- nrow(A)

  tau <- prod(dim(A)) / sum(A^2)
  Z <- matrix(0, nrow = nrowA, ncol = 0)
  EL <- matrix(0, nrow = ncol(A), ncol = 0)


  ebcd <- list(
    A = A, N = N, nrowA = nrowA, tol = tol,
    tau = tau, EL = EL, Z = Z
  )


  return(ebcd)
}
