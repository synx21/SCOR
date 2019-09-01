#' Smooth Approximations Of Empirical Hyper Volume Under Manifolds
#'
#' `SHUM` is a smooth approximation of EHUM
#'
#' @param beta The parameter we measure SHUM based on.
#' @param labels The labels of the Columns of the data matrix.
#' @param x_mat The Data Matrix
#' @param p p decides whether to use \eqn{s_n(x)} or \eqn{\phi_n(x)}. p = 1 stands for \eqn{\phi_n(x)} and p = 0 stands for \eqn{s_n(x)}
#'
#' @return Smooth approximation of the empirical Hyper-volume Under Manifolds Estimate
#'
#' @examples
#' estimate_SHUM(rep(1,12),colnames(AL),AL)
#' estimate_SHUM(rep(1,12),colnames(AL),AL, p = 1)
#'
#'
#' estimate_SHUM(1:10 , sample(c(rep("lab1",10),rep("lab2",10),rep("lab3",10))), matrix(rnorm(300), nrow = 10))
#'
#' @name estimate_SHUM
NULL

obj <- function(x, n, p)
{
  if (!p)
    return(1 / (1 + exp(-x * n)))
  else
    return(pnorm(x * n, 0, 1))
}

nys <- function(sock, labels, p, n, vec, hold)
{
  uni <- unique(labels)
  l <- length(uni)

  sum <- 0

  t <- which(labels == l - 1)

  for (v in sock[t])
  {
    vec[l] <- v
    if (l > 1)
      sum <- sum + nys(sock[-t], labels[-t], p, n, vec, hold)
    else

      sum <- sum + prod(obj(vec[-1] - vec[-hold], n, p))

  }
  return(sum)


}

#' @rdname estimate_SHUM
#' @export

estimate_SHUM <- function(beta, labels, x_mat, p = 0)
{
  stopifnot(nrow(x_mat) == length(beta))
  stopifnot(ncol(x_mat) == length(labels))

  uni <- unique(labels)
  l <- length(uni)
  stopifnot(l > 1)

  n <- length(labels)

  sock <- beta %*% x_mat
  #Changing labels to consecutive integers
  for (i in 1:l)
  {
    t <- which(labels == uni[i])
    labels[t] = i - 1
  }
  labels = as.numeric(labels)
  e <- numeric(l)

  for (i in length(labels):1)
    e[labels[i] + 1] <- 1 + e[labels[i] + 1]
  #"vec" stores all combinations from "l" categories
  vec <- numeric(l)
  hold <- l
  #Computing SHUM value by running recursive function of "l" for loops
  return(nys(sock, labels, p, n, vec, hold) / prod(e))

}
