#' fitting the model with squared exponential kernel.
#'
#' @details The choice of bounds for the optimization follows the approach used in \pkg{hetGP}.
#' For more details, see the reference below.
#'
#' @references
#' M. Binois and R. B. Gramacy (2021). hetGP: Heteroskedastic Gaussian Process Modeling and Sequential Design in R.
#' \emph{Journal of Statistical Software}, 98(13), 1-44;
#' \doi{doi: 10.18637/jss.v098.i13}
#'
#' @param X vector or matrix of input locations.
#' @param y vector of response values.
#' @param g nugget parameter. Default is 1.490116e-08.
#' @param constant logical indicating for constant mean (constant=TRUE) or zero mean (constant=FALSE). Default is FALSE.
#' @param p quantile on distances. Default is 0.1.
#' @param min_cor minimal correlation between two design points at the defined p quantile distance. Default is 0.2.
#' @param max_cor maximal correlation between two design points at the defined (1-p) quantile distance. Default is 0.5.
#' @param init initial value of optimization. Default is \code{NULL}.
#' @param lower lower bound of optimization. Default is \code{NULL}.
#' @param upper upper bound of optimization. Default is \code{NULL}.
#'
#' @return A list containing hyperparameters, covariance inverse matrix, X, y and logical inputs:
#' \itemize{
#'   \item \code{K}: matrix of covariance.
#'   \item \code{Ki}: matrix of covariance inverse.
#'   \item \code{X}: copy of X.
#'   \item \code{y}: copy of y.
#'   \item \code{theta}: vector of lengthscale hyperparameter.
#'   \item \code{g}: copy of g.
#'   \item \code{mu.hat}: optimized constant mean. If constant=FALSE, 0.
#'   \item \code{tau2hat}: estimated scale hyperparameter.
#'   \item \code{constant}: copy of constant.
#' }
#'
#' @importFrom plgp distance covar.sep
#' @importFrom stats optim quantile
#' @noRd
#' @keywords internal
#' @examples
#' \dontrun{
#' library(lhs)
#' ### synthetic function ###
#' f1 <- function(x) {
#'   sin(8 * pi * x)
#' }
#'
#' ### training data ###
#' n1 <- 15
#'
#' X1 <- maximinLHS(n1, 1)
#' y1 <- f1(X1)
#'
#' GP(X1, y1)
#' }
GP <- function(X, y, g = sqrt(.Machine$double.eps), constant = FALSE, p=0.1, min_cor = 0.2, max_cor = 0.5, init=NULL, lower=NULL, upper=NULL) { # p=0.05 for hetGP
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1) {
    stop("GP: 'p' must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(min_cor) || length(min_cor) != 1 || min_cor <= 0 || min_cor >= 1) {
    stop("GP: 'min_cor' must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(max_cor) || length(max_cor) != 1 || max_cor <= 0 || max_cor >= 1) {
    stop("GP: 'max_cor' must be a single numeric value between 0 and 1.")
  }
  if (min_cor >= max_cor) {
    stop("GP: 'min_cor' must be strictly less than 'max_cor'.")
  }
  if (!is.null(init) && !is.numeric(init)) stop("GP: 'init' must be numeric.")
  if (!is.null(lower) && !is.numeric(lower)) stop("GP: 'lower' must be numeric.")
  if (!is.null(upper) && !is.numeric(upper)) stop("GP: 'upper' must be numeric.")

  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  y <- as.matrix(y)
  n <- nrow(X)

  Xscaled <- crossprod(t(X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)),
                       diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X)))
  if(is.null(lower)) lower <- -quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], p) / log(min_cor) *
    (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2
  if(is.null(upper)) upper <- -quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], 1-p) / log(max_cor) *
    (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2
  if(is.null(init)) init <- sqrt(lower * upper)

  if (constant) {
    one.vec <- matrix(1, ncol = 1, nrow = n)

    nlsep <- function(par, X, Y) {
      theta <- par
      K <- covar.sep(X, d = theta, g = g)
      chol_K <- chol(K)
      ldetK <- 2 * sum(log(diag(chol_K)))

      KinvY <- backsolve(chol_K, forwardsolve(t(chol_K), Y))
      Kinv1 <- backsolve(chol_K, forwardsolve(t(chol_K), one.vec))
      mu.hat <- drop(crossprod(one.vec, KinvY) / crossprod(one.vec, Kinv1))
      tau2hat <- drop(crossprod(Y - mu.hat, KinvY - mu.hat * Kinv1) / n)

      ll <- -(n / 2) * log(tau2hat) - (1 / 2) * ldetK
      return(drop(-ll))
    }

    gradnlsep <- function(par, X, Y) {
      theta <- par
      K <- covar.sep(X, d = theta, g = g)
      chol_K <- chol(K)
      Ki <- chol2inv(chol_K)
      mu.hat <- drop(crossprod(one.vec, crossprod(Ki, Y)) / crossprod(one.vec, crossprod(Ki, one.vec)))

      Ki_resid <- crossprod(Ki, Y - mu.hat)

      dlltheta <- rep(NA, length(theta))
      for (k in 1:length(dlltheta)) {
        dotK <- K * distance(X[, k]) / (theta[k]^2)
        dlltheta[k] <- (n / 2) * drop(crossprod(Ki_resid, crossprod(dotK, Ki_resid))) / drop(crossprod(Y - mu.hat, Ki_resid)) - (1 / 2) * sum(Ki * dotK)
      }
      return(-c(dlltheta))
    }

    outg <- optim(init, nlsep, gradnlsep, method = "L-BFGS-B", lower = lower, upper = upper, X = X, Y = y)

    theta <- outg$par
    K <- covar.sep(X, d = theta, g = g)
    Ki <- chol2inv(chol(K))
    mu.hat <- drop(crossprod(one.vec, crossprod(Ki, y)) / crossprod(one.vec, crossprod(Ki, one.vec)))
    tau2hat <- drop(crossprod((y - mu.hat), crossprod(Ki, (y - mu.hat))) / n)
    names(theta) <- NULL

    return(list(K = K, Ki = Ki, X = X, y = y, theta = theta, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  } else {

    nlsep <- function(par, X, Y) {
      theta <- par
      K <- covar.sep(X, d = theta, g = g)
      chol_K <- chol(K)
      ldetK <- 2 * sum(log(diag(chol_K)))

      KinvY <- backsolve(chol_K, forwardsolve(t(chol_K), Y))

      ll <- -(n / 2) * log(drop(crossprod(Y, KinvY))) - (1 / 2) * ldetK
      return(drop(-ll))
    }

    gradnlsep <- function(par, X, Y) {
      theta <- par
      K <- covar.sep(X, d = theta, g = g)
      chol_K <- chol(K)
      Ki <- chol2inv(chol_K)
      KiY <- crossprod(Ki, Y)

      dlltheta <- rep(NA, length(theta))
      for (k in 1:length(dlltheta)) {
        dotK <- K * distance(X[, k]) / (theta[k]^2)
        dlltheta[k] <- (n / 2) * drop(crossprod(KiY, crossprod(dotK, KiY))) / drop(crossprod(Y, KiY)) - (1 / 2) * sum(Ki * dotK)
      }
      return(-c(dlltheta))
    }

    outg <- optim(init, nlsep, gradnlsep, method = "L-BFGS-B", lower = lower, upper = upper, X = X, Y = y)

    theta <- outg$par
    K <- covar.sep(X, d = theta, g = g)
    Ki <- chol2inv(chol(K))
    mu.hat <- 0
    tau2hat <- drop(crossprod(y, crossprod(Ki, y)) / n)
    names(theta) <- NULL

    return(list(K = K, Ki = Ki, X = X, y = y, theta = theta, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  }
}

#' predictive posterior mean and variance with squared exponential kernel.
#'
#' @param fit an object of class GP.
#' @param xnew vector or matrix of new input locations to predict.
#' @param cov.out logical indicating for returning covariance matrix (cov.out=TRUE) or not (cov.out=FALSE). Default is FALSE.
#'
#' @return A list predictive posterior mean and variance:
#' \itemize{
#'   \item \code{mu}: vector of predictive posterior mean.
#'   \item \code{cov}: matrix of predictive posterior covariance.
#'   \item \code{sig2}: vector of predictive posterior variance.
#' }
#'
#' @importFrom plgp covar.sep
#' @noRd
#' @keywords internal
#' @examples
#' \dontrun{
#' library(lhs)
#' ### synthetic function ###
#' f1 <- function(x) {
#'   sin(8 * pi * x)
#' }
#'
#' ### training data ###
#' n1 <- 15
#'
#' X1 <- maximinLHS(n1, 1)
#' y1 <- f1(X1)
#'
#' fit1 <- GP(X1, y1)
#'
#' ### test data ###
#' x <- seq(0, 1, 0.01)
#' pred.GP(fit1, x)
#' }
pred.GP <- function(fit, xnew, cov.out=FALSE) {
  if(is.null(fit$Ki) || is.null(fit$theta)) stop("Invalid fit object passed to pred.GP")
  xnew <- as.matrix(xnew)

  Ki <- fit$Ki
  theta <- fit$theta
  g <- fit$g
  X <- fit$X
  y <- fit$y
  tau2hat <- fit$tau2hat
  mu.hat <- fit$mu.hat

  KXX <- covar.sep(xnew, d = theta, g = g)
  KX <- covar.sep(xnew, X, d = theta, g = 0)

  mup2 <- mu.hat + crossprod(t(KX), crossprod(Ki, y - mu.hat))
  Sigmap2 <- tau2hat * (KXX - crossprod(t(KX), tcrossprod(Ki, KX)))
  Sigmap2 <- (t(Sigmap2) + Sigmap2) / 2
  if(cov.out){
    return(list(mu = mup2, cov = Sigmap2, sig2 = pmax(0, diag(Sigmap2))))
  }else{
    return(list(mu = mup2, sig2 = pmax(0, diag(Sigmap2))))
  }
}
