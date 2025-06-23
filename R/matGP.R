#' calculating Matern kernel with corresponding smoothness parameter
#'
#'
#' @param r vector or matrix of input.
#' @param nu numerical value of smoothness hyperparameter. It should be 1.5, or 2.5.
#' @noRd
#' @keywords internal
#' @return A value from Matern kernel.

matern.kernel <- function(r, nu) {
  if(nu==1.5){
    out <- (1+r*sqrt(3)) * exp(-r*sqrt(3))
  }else if(nu==2.5){
    out <- (1+r*sqrt(5)+5*r^2/3) * exp(-r*sqrt(5))
  }
  return(out)
}

#' calculating separable Matern kernel
#'
#' @param X vector or matrix of input.
#' @param x vector or matrix of new input. Default is NULL
#' @param theta lengthscale parameter. It should have the length of ncol(X).
#' @param nu numerical value of smoothness hyperparameter. It should be 1.5, or 2.5.
#' @param derivative logical indicating for its first derivative(derivative=1)
#' @noRd
#' @keywords internal
#' @return A covariance matrix of Matern kernel.

cor.mat <- function(X, x = NULL, theta, nu) {
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  d <- NCOL(X)
  n <- NROW(X)
  nu <- rep(nu, d)
  if (is.null(x)) {
    K <- matrix(1, n, n)
    for (i in 1:d) {
      R <- sqrt(distance(X[, i] / theta[i]))
      K <- K * matern.kernel(R, nu = nu[i])
    }
  } else {
    n.new <- NROW(x)
    K <- matrix(1, n, n.new)
    for (i in 1:d) {
      R <- sqrt(distance(X[, i] / theta[i], x[, i] / theta[i]))
      K <- K * matern.kernel(R, nu = nu[i])
    }
  }
  return(K)
}

#' fitting the model with Matern kernel.
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
#' @param nu numerical value of smoothness hyperparameter. It should be 1.5, or 2.5.
#' @param constant logical indicating for constant mean (constant=TRUE) or zero mean (constant=FALSE). Default is FALSE.
#' @param p quantile on distances. Default is 0.1.
#' @param min_cor minimal correlation between two design points at the defined quantile distance. Default is 0.2.
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
#'   \item \code{nu}: copy of nu.
#'   \item \code{g}: copy of g.
#'   \item \code{mu.hat}: optimized constant mean. If constant=FALSE, 0.
#'   \item \code{tau2hat}: estimated scale hyperparameter.
#'   \item \code{constant}: copy of constant.
#' }
#'
#' @importFrom stats optim quantile uniroot
#' @importFrom methods is
#' @noRd
#' @keywords internal
#'
matGP <- function(X, y, nu = 2.5, g = sqrt(.Machine$double.eps), constant = FALSE, p=0.1, min_cor = 0.2, max_cor = 0.5, init=NULL, lower=NULL, upper=NULL) {
  if (constant) {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)
    # hetGP way
    Xscaled <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*%
      diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X))
    tmpfun <- function(theta, repr_dist, value){
      cor.mat(matrix(sqrt(repr_dist/ncol(X)), ncol = ncol(X)), matrix(0, ncol = ncol(X)), theta = rep(theta,ncol(X)), nu = nu) - value
    }
    theta_min <- try(uniroot(tmpfun, interval = c(g, 100), value = min_cor, #nu=nu,
                             repr_dist = quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], p), tol = g)$root)
    if(is(theta_min, "try-error")){
      warning("The automatic selection of lengthscales bounds was not successful. Perhaps provide lower and upper values.")
      theta_min <- 1e-2
    }
    theta_max <- try(uniroot(tmpfun, interval = c(g, 100), value = max_cor,
                             repr_dist = quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], 1-p), tol = g)$root, silent = TRUE)
    if(is(theta_max, "try-error")){
      theta_max <- 5
    }
    if(is.null(lower)) lower <- theta_min * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])
    if(is.null(upper)) upper <- max(1, theta_max) * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])
    if(is.null(init)) init <- sqrt(lower * upper)

    n <- length(y)

    nlsep <- function(par, X, Y) {
      theta <- par # lengthscale
      K <- cor.mat(X, theta = theta, nu = nu)
      Ki <- solve(K + diag(g, n))
      ldetK <- determinant(K, logarithm = TRUE)$modulus

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      tau2hat <- drop(t(Y - mu.hat) %*% Ki %*% (Y - mu.hat) / n)
      ll <- -(n / 2) * log(tau2hat) - (1 / 2) * ldetK
      return(-ll)
    }

    gradnlsep <- function(par, X, Y) {
      theta <- par
      K <- cor.mat(X, theta = theta, nu = nu)
      Ki <- solve(K + diag(g, n))

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      KiY <- Ki %*% (Y - mu.hat)
      ## loop over theta components
      dlltheta <- rep(NA, length(theta))
      for (k in 1:length(dlltheta)) {
        if(nu==1.5){
          u_k <- sqrt(distance(X[, k] / theta[k]))
          dotK <- K * (3 * u_k^2 / theta[k]) / (1 + sqrt(3) * u_k)
        }else{
          u_k <- sqrt(distance(X[, k] / theta[k]))
          dotK <- K * (5/3 * u_k^2 / theta[k]) * (1+sqrt(5)*u_k) / (1 + sqrt(5) * u_k + 5/3*u_k^2)
        }
        dlltheta[k] <- (n / 2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dotK))
      }
      return(-c(dlltheta))
    }

    outg <- optim(init, nlsep, gradnlsep, method = "L-BFGS-B", lower = lower, upper = upper, X = X, Y = y)

    theta <- outg$par
    K <- cor.mat(X, theta = theta, nu = nu)
    Ki <- solve(K + diag(g, n))
    one.vec <- matrix(1, ncol = 1, nrow = n)
    mu.hat <- drop((t(one.vec) %*% Ki %*% y) / (t(one.vec) %*% Ki %*% one.vec))
    tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / nrow(X))
    names(theta) <- NULL

    return(list(K = K, Ki = Ki, X = X, y = y, theta = theta, nu = nu, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  } else {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)
    # hetGP way
    Xscaled <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*%
      diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X))
    tmpfun <- function(theta, repr_dist, value){
      cor.mat(matrix(sqrt(repr_dist/ncol(X)), ncol = ncol(X)), matrix(0, ncol = ncol(X)), theta = rep(theta,ncol(X)), nu = nu) - value
    }
    theta_min <- try(uniroot(tmpfun, interval = c(g, 100), value = min_cor,
                             repr_dist = quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], p), tol = g)$root)
    if(is(theta_min, "try-error")){
      warning("The automatic selection of lengthscales bounds was not successful. Perhaps provide lower and upper values.")
      theta_min <- 1e-2
    }
    theta_max <- try(uniroot(tmpfun, interval = c(g, 100), value = max_cor,
                             repr_dist = quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], 1-p), tol = g)$root, silent = TRUE)
    if(is(theta_max, "try-error")){
      theta_max <- 5
    }
    if(is.null(lower)) lower <- theta_min * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])
    if(is.null(upper)) upper <- max(1, theta_max) * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])
    if(is.null(init)) init <- sqrt(lower * upper)

    n <- length(y)

    nlsep <- function(par, X, Y) {
      theta <- par # lengthscale
      K <- cor.mat(X, theta = theta, nu = nu)
      Ki <- solve(K + diag(g, n))
      ldetK <- determinant(K, logarithm = TRUE)$modulus
      tau2hat <- drop(t(Y) %*% Ki %*% (Y) / n)
      ll <- -(n / 2) * log(tau2hat) - (1 / 2) * ldetK
      return(drop(-ll))
    }

    gradnlsep <- function(par, X, Y) {
      theta <- par
      K <- cor.mat(X, theta = theta, nu = nu)
      Ki <- solve(K + diag(g, n))

      KiY <- Ki %*% Y
      ## loop over theta components
      dlltheta <- rep(NA, length(theta))
      for (k in 1:length(dlltheta)) {
        if(nu==1.5){
          u_k <- sqrt(distance(X[, k] / theta[k]))
          dotK <- K * (3 * u_k^2 / theta[k]) / (1 + sqrt(3) * u_k)
        }else{
          u_k <- sqrt(distance(X[, k] / theta[k]))
          dotK <- K * (5/3 * u_k^2 / theta[k]) * (1+sqrt(5)*u_k) / (1 + sqrt(5) * u_k + 5/3*u_k^2)
        }
        dlltheta[k] <- (n / 2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dotK))
      }
      return(-c(dlltheta))
    }

    outg <- optim(init, nlsep, gradnlsep, method = "L-BFGS-B", lower = lower, upper = upper, X = X, Y = y)

    theta <- outg$par
    K <- cor.mat(X, theta = outg$par, nu = nu)
    Ki <- solve(K + diag(g, n))
    mu.hat <- 0
    tau2hat <- drop(t(y) %*% Ki %*% (y) / n)
    names(theta) <- NULL

    return(list(K = K, Ki = Ki, X = X, y = y, theta = theta, nu = nu, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  }
}

#' predictive posterior mean and variance with Matern kernel.
#'
#' @param fit an object of class matGP.
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
#' @noRd
#' @keywords internal
#'
pred.matGP <- function(fit, xnew, cov.out=FALSE) {
  xnew <- as.matrix(xnew)

  Ki <- fit$Ki
  theta <- fit$theta
  nu <- fit$nu
  g <- fit$g
  X <- fit$X
  y <- fit$y
  tau2hat <- fit$tau2hat
  mu.hat <- fit$mu.hat

  KXX <- cor.mat(xnew, theta = theta, nu = nu)
  KX <- t(cor.mat(X, xnew, theta = theta, nu = nu))

  mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
  Sigmap2 <- tau2hat * (KXX - KX %*% Ki %*% t(KX))
  Sigmap2 <- (t(Sigmap2)+Sigmap2)/2
  if(cov.out){
    return(list(mu = mup2, cov=Sigmap2, sig2 = pmax(0,diag(Sigmap2))))
  }else{
    return(list(mu = mup2, sig2 = pmax(0,diag(Sigmap2))))
  }
}
