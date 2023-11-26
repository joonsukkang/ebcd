#' Empirical Bayes Covariance Decomposition
#'
#' @param X A N-by-P data matrix with N observations and P variables
#' @param S A P-by-P covariance/Gram matrix \eqn{=X'X/N}
#' @param C A P-by-P compact target matrix such that \eqn{C'C=X'X}
#' @param N The number of observations
#' @param Kmax The maximum number of factors to be added.
#' @param ebnm_fn The ebnm function to be used.
#' @param tol_greedy  The convergence tolerance parameter when optimizing the greedily added factor.
#' @param maxiter_greedy The maximum number of iterations when optimizing the greedily added factor.
#' @param tol_backfit The convergence tolerance parameter for the backfit procedure.
#' @param maxiter_backfit The maximum number of iterations for the backfit procedure.
#'
#' @return An ebcd object.
#' @export
#'
ebcd <- function(X = NULL,
                 S = NULL,
                 C = NULL,
                 N = NULL,
                 Kmax = 5,
                 ebnm_fn = ebnm::ebnm_point_laplace,
                 tol_greedy = 1e-6,
                 maxiter_greedy = 500,
                 tol_backfit = 1e-6,
                 maxiter_backfit = 5000) {
  ebcd_obj <- ebcd_init(X = X, S = S, C = C, N = N) |>
    ebcd_greedy(
      Kmax = Kmax,
      ebnm_fn = ebnm_fn,
      tol = tol_greedy,
      maxiter = maxiter_greedy
    ) |>
    ebcd_backfit(
      tol = tol_backfit,
      maxiter = maxiter_backfit
    )

  return(ebcd_obj)
}
