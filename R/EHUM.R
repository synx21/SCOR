#' Empirical Hyper Volume Under Manifolds
#'
#' An estimator of Hyper Volume Under Manifolds
#'
#'
#' @param beta The parameter we measure EHUM based on.
#' @param labels The labels of the Columns of the data matrix.
#' @param x_mat The Data Matrix
#'
#' @return Empirical Hyper-volume Under Maniforlds Estimate
#'
#' @export
#'
#' @examples
#' estimate_EHUM(rep(1,12),colnames(AL),AL)
#'
#'
#' estimate_EHUM(1:10 , sample(c(rep("lab1",10),rep("lab2",10),rep("lab3",10))), matrix(rnorm(300), nrow = 10))
#'

estimate_EHUM <- function(beta, labels, x_mat)
{
  uni <- unique(labels)

  sock <- beta %*% x_mat

  #Changing labels to consecutive integers
  for (i in 1:length(uni))
  {
    t <- which(labels == uni[i])
    labels[t] = i - 1
  }
  labels = as.numeric(labels)

  order <- order(sock)
  labels <- labels[order]
  sock <- sock[order]

  t <- numeric(length(uni) + 1)
  p <- numeric(length(uni))
  t[length(uni) + 1] <- 1
  #Making Graphs where i joins (i+1)and (i-1) and counting such entries.
  for (i in length(labels):1)
  {
    t[labels[i] + 1] <- t[labels[i] + 2] + t[labels[i] + 1]
    p[labels[i] + 1] <- 1 + p[labels[i] + 1]
  }
  return(t[1] / prod(p))
}
