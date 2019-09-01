#' Visualization Based On Youden Indices.
#'
#' A Box Plot Visualization Based On Youden Indices for less than equal to 3 categories.
#'
#' @param beta The parameter we do HUM based on
#' @param labels The labels of the Columns of the data matrix
#' @param x_mat The Data Matrix
#' @param cat_names The vector of strings containing category names.
#' @param grid_size The size of increment in the grid we check cutpoints against. Default value is 100.
#'
#' @return Box Plot Visualization Based On Youden Indices
#'
#' @export
#'
#' @examples
#'
#' beta <- c(-0.399,-0.155,-0.265,-0.184,
#'      -0.267,0.666,-0.187,0.273,0.0463,0.167,0.163,0.178)
#'
#' YoupointsBoxPlot(beta,colnames(AL),AL, cat_names = c("Healthy","MCI","AD"))

YoupointsBoxPlot <-
  function(beta,
           labels,
           x_mat,
           cat_names = NULL,
           grid_size = 100)
  {
    uni <- unique(labels)
    l <- length(uni)
    stopifnot(nrow(x_mat) == length(beta))
    stopifnot(ncol(x_mat) == length(labels))
    stopifnot(l < 4)
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
    cat <- matrix(FALSE, sum(sample_sizes), l)
    sum = 0
    for (i in 1:l)
    {
      cat[sum:(sum + sample_sizes[i]), i] <- TRUE
      sum <- sum + sample_sizes[i]
    }
    AA <- youden_points(beta, labels, x_mat, grid_size)
    #Looks out for parameters via function youden_points and uses them to make boxplot.
    x_mat <- temp
    rm(temp)
    rm(temp_count)
    x <- beta %*% x_mat
    if (l == 3)
    {
      boxplot(x[cat[, 1]], x[cat[, 2]], x[cat[, 3]], names = cat_names, ylab = 'Score')
      abline(h = AA$cutpoints[1],
             lwd = 1,
             lty = 2)
      abline(h = AA$cutpoints[2],
             lwd = 1,
             lty = 2)
    }
    else if (l == 2)
    {
      boxplot(x[cat[, 1]], x[cat[, 2]], names = cat_names, ylab = 'Score')
      abline(h = AA$cutpoints,
             lwd = 1,
             lty = 2)
    }
  }
