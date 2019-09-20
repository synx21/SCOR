#' Spherically Constrained Optimization Routine
#'
#' `SCOR` is our optimization algorithm, efficient in estimating maximizing Hyper Volume Under Manifolds Estimators.
#'
#' SCOR is the modified version of RMPS,  Recursive Modified Pattern Search.
#' This is a blackbox algorithm efficienct in optimizing non-differentiable functions.
#' It works great in the shown cases of SHUM, EHUM and ULBA.
#'
#'
#' @param x0 The initial guess by user
#' @param func The function to be optimized
#' @param rho Step Decay Rate with default value 2
#' @param phi Lower Bound Of Global Step Size. Default value is \eqn{10^{-6}}
#' @param max_iter Max Number Of Iterations In each Run. Default Value is 50,000.
#' @param s_init Initial Global Step Size. Default Value is 2.
#' @param tol_fun Termination Tolerance on the function value. Default Value is \eqn{10^{-6}}
#' @param tol_fun_2 Termination Tolerance on the difference of solutions in two consecutive runs. Default Value is \eqn{10^{-6}}
#' @param minimize Binary Command to set SCOR on minimization or maximization. TRUE is for minimization which is set default.
#' @param time Time Alloted for execution of SCOR
#' @param print Binary Command to print optimized value of objective function after each interation. FALSE is set fault
#' @param lambda Sparsity Threshold. Default value is \eqn{10^{-3}}
#' @param parallel Binary Command to ask SCOR to perform parallel computing. Default is set at FALSE.
#'
#' @return The point where the value Of the Function is maximized under a sphere.
#'
#' @examples
#' f <- function(x)
#' return(x[2]^2 + x[3]^3 +x[4]^4)
#'
#' SCOR(rep(1,10),f)
#'
#' SCOR(c(2,4,6,2,1),f, minimize = FALSE, print = TRUE)
#' #Will Print the List and Find the Maximum
#'
#' SCOR(c(1,2,3,4),f, time = 10, lambda = 1e-2)
#' #Will perform no iterations after 10 secs, Sparsity Threshold is 0.01
#'
#' SCOR(c(2,6,2,7,8),f, parallel = TRUE)
#' #Will do Parallel Computing
#'
#' @references
#' \itemize{
#'
#'   \item Das, Priyam and De, Debsurya and Maiti, Raju and Chakraborty, Bibhas and Peterson, Christine B \cr
#'    "Estimating the Optimal Linear Combination of Biomarkers using Spherically Constrained Optimization" \cr
#'          (available at `arXiv \url{http://arxiv.org/abs/1909.04024}).
#' }
#'
#'
#' @name SCOR
NULL

dist <- function(a, b)
{
  return(sqrt(sum((a - b) ^ 2)))
}

syn <-
  function(x0,
           func,
           rho,
           phi,
           max_iter,
           s_init,
           tol_fun,
           minimize,
           time,
           t0,
           R,
           print,
           lambda,
           parallel)
  {
    d <- length(x0)
    j <- 1
    s <- numeric(max_iter)
    s[1] <- s_init
    x <- matrix(x0, nrow = 1)
    s1 <- numeric(2 * d)
    f <- numeric(2 * d)

    #doing several iterations until the distance is very less
    while (!(j > max_iter || s[j] < phi || Sys.time() - t0 > time))
    {
      F1 <- func(x[j, ])
      s2 <- s[j]
      if (parallel)
      {
        eta <- foreach (h = 1:(2 * d), .combine = rbind) %dopar%
          {
            i <- floor((h + 1) / 2)
            beta <- x[j, ]
            q1 <- (abs(beta) < lambda)
            q2 <- !(abs(beta) < lambda)
            q1[i] <- FALSE
            q2[i] <- FALSE
            s1[h] <- (-1) ^ h * s[j]
            D <-
              (2 * sum(beta[q2])) ^ 2 - 4 * length(q2[q2]) * ((2 * beta[i] + s1[h]) * s1[h] - sum(beta[q1] ^ 2))
            while (D < 0 && abs(s1[h]) > phi)
            {
              s1[h] <- s1[h] / rho
              D <-
                (2 * sum(beta[q2])) ^ 2 - 4 * length(q2[q2]) * ((2 * beta[i] + s1[h]) * s1[h] - sum(beta[q1] ^ 2))
            }
            if (!(D < 0))
            {
              t1 <- (sqrt(D) - 2 * sum(beta[q2])) / (2 * length(q2[q2]))
              t2 <-
                (-sqrt(D) - 2 * sum(beta[q2])) / (2 * length(q2[q2]))
              p1 <- numeric(d)
              p2 <- numeric(d)
              p1[i] <- beta[i] + s1[h]
              p2[i] <- p1[i]
              p1[q2] <- beta[q2] + t1
              p2[q2] <- beta[q2] + t2
              if (func(p1) > func(p2))
              {
                if (!minimize)
                  beta <- p1
                else
                  beta <- p2
              }
              else
              {
                if (minimize)
                  beta <- p1
                else
                  beta <- p2
              }
              if (max(abs(beta)) >= 1)
              {
                pos = which.max(abs(beta))
                beta[pos] = sign(beta[pos])
              }
              list(func(beta), beta)
            }
            else
              list(F1, beta)
          }
        f <- unlist(eta[, 1], use.names = FALSE)
        beta1 <-
          matrix(unlist(eta[, 2], use.names = FALSE),
                 nrow = 2 * d ,
                 byrow = TRUE)
      }
      else
      {
        beta1 <- t(replicate(2 * d, x[j, ]))
        for (h in 1:(2 * d))
        {
          i <- floor((h + 1) / 2)
          beta <- x[j, ]
          q1 <- (abs(beta) < lambda)
          q2 <- !(abs(beta) < lambda)
          q1[i] <- FALSE
          q2[i] <- FALSE
          s1[h] <- (-1) ^ h * s[j]
          D <-
            (2 * sum(beta[q2])) ^ 2 - 4 * length(q2[q2]) * ((2 * beta[i] + s1[h]) * s1[h] - sum(beta[q1] ^ 2))
          while (D < 0 && abs(s1[h]) > phi)
          {
            s1[h] <- s1[h] / rho
            D <-
              (2 * sum(beta[q2])) ^ 2 - 4 * length(q2[q2]) * ((2 * beta[i] + s1[h]) * s1[h] - sum(beta[q1] ^ 2))
          }
          if (!(D < 0))
          {
            t1 <- (sqrt(D) - 2 * sum(beta[q2])) / (2 * length(q2[q2]))
            t2 <-
              (-sqrt(D) - 2 * sum(beta[q2])) / (2 * length(q2[q2]))
            p1 <- numeric(d)
            p2 <- numeric(d)
            p1[i] <- beta[i] + s1[h]
            p2[i] <- p1[i]
            p1[q2] <- beta[q2] + t1
            p2[q2] <- beta[q2] + t2
            if (func(p1) > func(p2))
            {
              if (!minimize)
                beta <- p1
              else
                beta <- p2
            }
            else
            {
              if (minimize)
                beta <- p1
              else
                beta <- p2
            }
            if (max(abs(beta)) >= 1)
            {
              pos = which.max(abs(beta))
              beta[pos] = sign(beta[pos])
            }

            f[h] <- func(beta)
            beta1[h, ] <- beta
          }
          else
          {
            f[h] <- func(beta)
            beta1[h, ] <- beta
          }
        }
      }

      if (!minimize)
      {
        pos <- which.max(f)
        B <- beta1[pos, ]
        FF <- f[pos]
        if (FF > F1)
          t <- B
        else
          t <- x[j,]
        x <- rbind(x, t)
        if (dist(max(FF, F1), F1) < tol_fun && s2 > phi)
          s2 <- s2 / rho
      }
      else
      {
        pos = which.min(f)
        B <- beta1[pos, ]
        FF <- f[pos]
        if (FF < F1)
          t <- B
        else
          t <- x[j,]
        x <- rbind(x, t)
        if (dist(min(FF, F1), F1) < tol_fun && s2 > phi)
          s2 <- s2 / rho
      }
      s[j + 1] <- s2
      j <- j + 1
      if (print)
        cat("\t", R, "\t", j - 1, "\t", "\t", func(x[j,]), "\n")

    }
    return(x[j,])
  }

#' @rdname SCOR
#' @export

SCOR <-
  function(x0,
           func,
           rho = 2,
           phi = 10 ^ (-3),
           max_iter = 50000,
           s_init = 2,
           tol_fun = 10 ^ (-6),
           tol_fun_2 = 10 ^ (-6),
           minimize = TRUE,
           time = 600000,
           print = FALSE,
           lambda = 10 ^ (-3),
           parallel = FALSE)
  {
    if (print)
      cat("\t",
          "Run",
          "\t",
          "Iteration",
          "\t",
          "Objective Function",
          "\n")

    t0 = Sys.time()

    #normalize the objective vector
    x0 <- x0 / dist(0, x0)
    z <- matrix(0, 1000, length(x0))

    #initializing parallel threading
    if (parallel)
    {
      cores <- detectCores()
      c1 <- makeCluster(cores[1] - 1)
      registerDoParallel(c1)
    }

    #The Run Count, as in Step 2 Of the Algorithm
    R <- 1

    #Performing step 1 of the Algorithm, as done in the function "syn"
    z[1, ] <-
      syn(
        x0,
        func,
        rho,
        phi,
        max_iter,
        s_init,
        tol_fun,
        minimize,
        time,
        t0,
        R,
        print,
        lambda,
        parallel
      )

    #We'll halt when distance between two runs is not much.
    while (R < 2 || dist(z[R, ], z[R - 1, ]) > tol_fun_2)
    {
      #Time Check
      if (Sys.time() - t0 > time)
        return(z[R, ])

      R <- R + 1
      z[R, ] <-
        syn(
          z[R - 1, ],
          func,
          rho,
          phi,
          max_iter,
          s_init,
          tol_fun,
          minimize,
          time,
          t0,
          R,
          print,
          lambda,
          parallel
        )
    }

    if (parallel)
      stopCluster(c1)

    return(z[R, ])
  }
