#' Finding Youden Indices
#'
#' A function to find Youden Indices and Cutpoints for number of categories less than equal to 3.
#'
#' @param beta The parameter we do HUM based on
#' @param labels The labels of the Columns of the data matrix.
#' @param x_mat The Data Matrix
#' @param grid_size The size of increment in the grid we check cutpoints against. Default value is 100.
#'
#' @return Youden Indices and Cut Points
#'
#' @examples
#'
#' beta <- c(-0.399,-0.155,-0.265,-0.184,
#'      -0.267,0.666,-0.187,0.273,0.0463,0.167,0.163,0.178)
#'
#' youden_points(beta,colnames(AL),AL)
#'
#'
#' @name youden_points
NULL

youden3 <- function(x_mat, cat, t_minus, t_plus)
{
  F1 <-
    length(which(x_mat[cat[, 1]] < t_minus)) / length(x_mat[, cat[, 1]])
  F2 <-
    length(which(x_mat[cat[, 2]] < t_minus)) / length(x_mat[, cat[, 2]])
  F3 <-
    length(which(x_mat[cat[, 2]] < t_plus)) / length(x_mat[, cat[, 2]])
  F4 <-
    length(which(x_mat[cat[, 3]] < t_plus)) / length(x_mat[, cat[, 3]])

  youden_val <- 0.5 * (F1 - F2 + F3 - F4)

  return(youden_val)
}

youden2 <- function(x_mat, cat, t)
{
  F1 <-
    length(which(x_mat[cat[, 1]] < t)) / length(x_mat[, cat[, 1]])
  F2 <-
    length(which(x_mat[cat[, 2]] < t)) / length(x_mat[, cat[, 2]])

  youden_val <- F1 - F2

  return(youden_val)
}

#' @rdname youden_points
#' @export
youden_points <-
  function(beta, labels, x_mat, grid_size = 100)
  {
    uni <- unique(labels)
    l <- length(uni)
    stopifnot(l < 4)
    stopifnot(nrow(x_mat) == length(beta))
    stopifnot(ncol(x_mat) == length(labels))
    sample_sizes <- numeric(length(uni))
    temp <- matrix(0, dim(x_mat)[1], dim(x_mat)[2])
    temp_count <- 1
    for (i in 1:length(uni))
    {
      sample_sizes[i] <- sum(labels == uni[i])
      temp[, temp_count:(temp_count + sample_sizes[i] - 1)] <-
        x_mat[, labels == uni[i]]
      temp_count = temp_count + sample_sizes[i]
    }
    #arranging x_mat with respect to labels
    x_mat <- temp
    rm(temp)
    rm(temp_count)
    x_mat <- beta %*% x_mat
    categories <- l
    cat <- matrix(FALSE, sum(sample_sizes), l)
    sum = 0
    t <- numeric(l)
    #cat contains TRUE values wherever label is same as row index
    for (i in 1:l)
    {
      cat[sum:(sum + sample_sizes[i]), i] <- TRUE
      sum <- sum + sample_sizes[i]
      t[i] <- median(x_mat[cat[, i]])
    }
    poss <- list()
    #poss is a list of vectors , that are sequence from two median values of clusters.
    for (i in 1:(l - 1))
      poss[[i]] <-
      seq(from = t[i],
          to = t[i + 1],
          by = (t[i + 1] - t[i]) / grid_size)
    if (l == 3)
    {
      youdens <- matrix(0, length(poss[[1]]), length(poss[[2]]))
      for (ii in 1:length(poss[[1]]))
      {
        t_minus <- poss[[1]][ii]
        for (jj in 1:length(poss[[2]]))
        {
          t_plus <- poss[[2]][jj]
          #youden matrix accounts for youden values between clusters
          youdens[ii, jj] <- youden3(x_mat, cat, t_minus, t_plus)
        }
      }
      indices <- which(youdens == max(youdens), arr.ind = TRUE)
      indices <- indices[1,]
      return(list("YI" = youdens[indices[1], indices[2]],
                  "cutpoints" = c(poss[[1]][indices[1]], poss[[2]][indices[2]])))
    }
    else
    {
      youdens <- numeric(length(poss[[1]]))
      for (ii in 1:length(poss[[1]]))
      {
        t <- poss[[1]][ii]
        #youden matrix accounts for youden values between clusters
        youdens[ii] <- youden2(x_mat, cat, t)
      }
      indices <- which(youdens == max(youdens))
      indices <- indices[1]
      return(list("YI" = youdens[indices[1]], "cutpoints" = poss[[1]][indices[1]]))
    }
  }
