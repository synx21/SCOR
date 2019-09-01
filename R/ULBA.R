#' Upper And Lower Bound Approach
#'
#' `ULBA` is an another approach to Hyper Volume Under Manifold Problem
#'
#' @param beta The parameter we measure ULBA based on.
#' @param labels The labels of the Columns of the data matrix.
#' @param x_mat The Data Matrix
#'
#' @return  Upper and Lower Bound Approach on empirical Hyper-volume Under Manifolds Estimate
#'
#' @export
#'
#' @examples
#' estimate_ULBA(rep(1,12),colnames(AL),AL)
#'
#'
#' estimate_ULBA(1:10 , sample(c(rep("lab1",10),rep("lab2",10),rep("lab3",10))), matrix(rnorm(300), nrow = 10))
#'
#'

estimate_ULBA <- function(beta, labels, x_mat)
{
  stopifnot(nrow(x_mat) == length(beta))
  stopifnot(ncol(x_mat) == length(labels))

  #Finding Unique Labels
  uni <- unique(labels)
  sock <- beta %*% x_mat
  sum <- 0

  for (i in 2:length(uni))
  {
    queue <- 0

    p <- labels == uni[i]

    q <- labels == uni[i - 1]

    #To find number of vectors that are greater than the class preceding when mutiplied by beta
    for (m in sock[p])
      queue <- queue + sum(sock[q] < m)

    sum <- sum + queue / (sum(p) * sum(q))

  }
  return(sum / (length(uni) - 1))
}
