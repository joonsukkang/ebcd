#' Initialize ebcd object
#'
#' Initialize an ebcd object with no factors, taking as an input one of
#' a data matrix \eqn{X}, a covariance matrix \eqn{S}, or a compact target matrix \eqn{C}.
#' The algorithm runs on either \eqn{X} or \eqn{C}, which is internally saved as \eqn{A}.
#'
#' @param X A N-by-P data matrix with N observations and P variables
#' @param S A P-by-P covariance/Gram matrix \eqn{=X'X/N}
#' @param C A P-by-P compact target matrix such that \eqn{C'C=X'X}
#' @param N The number of observations
#'
#' @return An initialized ebcd object
#' @export
#'
ebcd_init <- function(X = NULL,
                      S = NULL,
                      C = NULL,
                      N = NULL) {
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


  ebcd_obj <- list(
    A = A, N = N, nrowA = nrowA,
    tau = tau, Z = Z, EL = EL
  )


  return(ebcd_obj)
}
