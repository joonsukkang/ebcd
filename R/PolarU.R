PolarU <- function(A) {
  svdA <- svd(A)
  out <- svdA$u %*% t(svdA$v)
  return(out)
}
