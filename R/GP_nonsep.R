#' @title Compute initial theta bounds for nonseparable squared exponential kernel
#'
#' @description The function computes lower and upper bounds for the
#' hyperparameters \code{theta_x} and \code{theta_t}
#' based on distance quantiles in input and tuning parameter space.
#'
#' The choice of bounds for the optimization follows the approach used in \pkg{hetGP}.
#' For more details, see the reference below.
#'
#' @references
#' M. Binois and R. B. Gramacy (2021). hetGP: Heteroskedastic Gaussian Process Modeling and Sequential Design in R.
#'\emph{Journal of Statistical Software}, 98(13), 1-44;
#' \doi{doi:10.18637/jss.v098.i13}
#'
#' @param X A vector or matrix of input locations.
#' @param t A vector of tuning parameters.
#' @param p Quantile on distances. Default is 0.05.
#' @param beta Interaction parameter between 0 and 1. Default is 0.5.
#' @param min_cor minimal correlation between two design points at the defined p quantile distance. Default is 0.1.
#' @param max_cor maximal correlation between two design points at the defined (1-p) quantile distance. Default is 0.9.
#'
#' @return A list with elements \code{lower} and \code{upper} containing estimated bounds.
#' @importFrom stats quantile
#' @importFrom plgp distance
#' @noRd
#'

theta_bounds_sqex <- function(X, t, p = 0.05, beta = 0.5, min_cor = 0.1, max_cor = 0.9) {

  Xscaled <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*%
    diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X))
  tscaled <- t / max(t)

  repr_dist_x_low <- quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], p)
  repr_dist_x_high <- quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], 1-p)
  repr_dist_t <- quantile(distance(tscaled)[lower.tri(distance(tscaled))], 0.5)

  d <- ncol(X)

  tmpfun <- function(theta, repr_dist_t, repr_dist_x, beta, d, value) {
    theta_x <- theta[1]
    theta_t <- theta[2]
    delta <- pmax(1-beta*d/2, 0.5)

    term1 <- (repr_dist_t / theta_t + 1)^(- (beta * d) / 2 - delta)
    term2 <- exp(- (repr_dist_t / theta_t + 1)^(-beta) * (repr_dist_x / theta_x))

    return(term1 * term2 - value)
  }

  theta_min <- tryCatch(
    stats::uniroot(
      function(theta) tmpfun(c(theta, theta), repr_dist_t, repr_dist_x_low, beta, d, min_cor),
      interval = c(sqrt(.Machine$double.eps), 100)
    )$root,
    error = function(e) {
      warning("Lower bound estimation failed. Using default values.")
      return(c(1e-2, 1e-2))
    }
  )
  theta_max <- tryCatch(
    stats::uniroot(
      function(theta) tmpfun(c(theta, theta), repr_dist_t, repr_dist_x_high, beta, d, max_cor),
      interval = c(sqrt(.Machine$double.eps), 100)
    )$root,
    error = function(e) {
      warning("Upper bound estimation failed. Using default values.")
      return(c(5, 5))
    }
  )

  return(list(lower = theta_min, upper = theta_max))
}

#' @title Compute initial theta bounds for nonseparable Matern kernel
#'
#' @description The function computes lower and upper bounds for the
#' hyperparameters \code{theta_x} and \code{theta_t}
#' based on distance quantiles in input and tuning parameter space.
#'
#' The choice of bounds for the optimization follows the approach used in \pkg{hetGP}.
#' For more details, see the reference below.
#'
#' @references
#' M. Binois and R. B. Gramacy (2021). hetGP: Heteroskedastic Gaussian Process Modeling and Sequential Design in R.
#'\emph{Journal of Statistical Software}, 98(13), 1-44;
#' \doi{doi:10.18637/jss.v098.i13}
#'
#' @param X A vector or matrix of input locations.
#' @param t A vector of tuning parameters.
#' @param nu A numerical value of smoothness hyperparameter. It should be 1.5, or 2.5.
#' @param p Quantile on distances. Default is 0.05.
#' @param beta Interaction parameter between 0 and 1. Default is 0.5.
#' @param min_cor minimal correlation between two design points at the defined p quantile distance. Default is 0.1.
#' @param max_cor maximal correlation between two design points at the defined (1-p) quantile distance. Default is 0.9.
#'
#' @return A list with elements \code{lower} and \code{upper} containing estimated bounds.
#' @importFrom stats quantile
#' @importFrom plgp distance
#' @noRd
#'

theta_bounds_matern <- function(X, t, nu, p = 0.05, beta = 0.5, min_cor = 0.1, max_cor = 0.9) {

  Xscaled <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*%
    diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X))
  tscaled <- t / max(t)

  repr_dist_x_low <- quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], p)
  repr_dist_x_high <- quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], 1-p)
  repr_dist_t <- quantile(distance(tscaled)[lower.tri(distance(tscaled))], 0.5)

  d <- ncol(X)

  tmpfun <- function(theta, repr_dist_t, repr_dist_x, beta, d, value) {
    theta_x <- theta[1]
    theta_t <- theta[2]
    delta <- pmax(1-beta*d/2, 0.5)

    term1 <- (repr_dist_t / theta_t + 1)^(- (beta * d) / 2 - delta)
    term2 <- matern.kernel(sqrt(repr_dist_x)/theta_x*(repr_dist_t / theta_t + 1)^(-beta/2), nu = nu)

    return(term1 * term2 - value)
  }

  theta_min <- tryCatch(
    stats::uniroot(
      function(theta) tmpfun(c(theta, theta), repr_dist_t, repr_dist_x_low, beta, d, min_cor),
      interval = c(sqrt(.Machine$double.eps), 100)
    )$root,
    error = function(e) {
      warning("Lower bound estimation failed. Using default values.")
      return(c(1e-2, 1e-2))
    }
  )
  theta_max <- tryCatch(
    stats::uniroot(
      function(theta) tmpfun(c(theta, theta), repr_dist_t, repr_dist_x_high, beta, d, max_cor),
      interval = c(sqrt(.Machine$double.eps), 100)
    )$root,
    error = function(e) {
      warning("Upper bound estimation failed. Using default values.")
      return(c(5, 5))
    }
  )

  return(list(lower = theta_min, upper = theta_max))
}

#' @title Covariance matrix for nonseparable squared exponential kernel
#'
#' @description The function computes the covariance matrix for
#' the nonseparable squared exponential kernel for training or prediction.
#'
#' @param X A vector or matrix of training input locations.
#' @param x A vector or matrix of optional new input locations for prediction.
#' @param t A vector of tuning parameters for training inputs.
#' @param tnew A vector of optional tuning parameters for new inputs.
#' @param param Kernel hyperparameters (\code{theta_x}, \code{theta_t}, \code{beta}, \code{delta}).
#'
#' @return Nonseparable squared exponential covariance matrix \code{K}.
#' @importFrom plgp covar.sep
#' @noRd
#'

cor.sep.sqex <- function(X, x = NULL, t, tnew=NULL, param) {
  d <- NCOL(X)
  n <- NROW(X)
  if (length(t) != n) {
    stop("Length of mesh size should be the number of rows of inputs")
  }

  theta_x <- param[1:d]
  theta_t <- param[d+1]
  beta <- param[d+2]
  delta <- param[d+3]

  if (is.null(x)) {
    if(length(unique(t))==1){
      Tcomponent=1
    }else{
      Tcomponent <- 1/(outer(t, t, "-")^(2)/theta_t + 1)
    }
    K <- Tcomponent^(beta*d/2+delta) * covar.sep(X, d=theta_x, g=0)^(Tcomponent^beta)
  } else {
    n.new <- NROW(x)
    if (length(tnew) != n.new) {
      stop("Length of mesh size should be the number of rows of inputs")
    }
    if(length(unique(t))==1){
      Tcomponent=1
    }else{
      Tcomponent <- 1/(outer(t, tnew, "-")^(2)/theta_t + 1)
    }
    K <- Tcomponent^(beta*d/2+delta) * covar.sep(X, x, d=theta_x, g=0)^(Tcomponent^beta)
  }
  return(K)
}

#' @title Calculate Matern kernel with corresponding smoothness parameter
#'
#' @description The function evaluates the Matern covariance function for given distances.
#'
#' @param r A scalar or vector of distance.
#' @param nu A numerical value of smoothness hyperparameter. It should be 1.5, or 2.5.
#'
#' @return Covariance value(s) from Matern kernel.
#' @noRd
#'

matern.kernel <- function(r, nu){
  if(nu==1.5){
    out <- (1+r*sqrt(3)) * exp(-r*sqrt(3))
  }else if(nu==2.5){
    out <- (1+r*sqrt(5)+5*r^2/3) * exp(-r*sqrt(5))
  }
  return(out)
}

#' @title Covariance matrix for nonseparable Matern kernel
#'
#' @description The function computes the correlation matrix for
#' the nonseparable Matern kernel for training or prediction.
#'
#' @param X A vector or matrix of training input locations.
#' @param x A vector or matrix of optional new input locations for prediction.
#' @param t A vector of tuning parameters for training inputs.
#' @param tnew A vector of optional tuning parameters for new inputs.
#' @param param Kernel hyperparameters (\code{theta_x}, \code{theta_t}, \code{beta}, \code{delta}).
#' @param nu A numerical value of smoothness hyperparameter. It should be 1.5, or 2.5.
#'
#' @return Nonseparable Matern kernel covariance matrix \code{K}.
#' @importFrom plgp distance
#' @noRd
#'

cor.sep.matern <- function(X, x = NULL, t, tnew=NULL, param, nu=2.5) {
  d <- NCOL(X)
  n <- NROW(X)
  if (length(t) != n) {
    stop("Length of mesh size should be the number of rows of inputs")
  }

  theta_x <- param[1:d]
  theta_t <- param[d+1]
  beta <- param[d+2]
  delta <- param[d+3]

  if (is.null(x)) {
    if(length(unique(t))==1){
      Tcomponent=1
    }else{
      Tcomponent <- 1/(outer(t, t, "-")^2/theta_t + 1)
    }
    # Rt <- sqrt(distance(t/sqrt(theta_t)))
    K <- matrix(1, n, n)
    for (i in 1:d) {
      R <- sqrt(distance(X[, i] / theta_x[i]))
      K <- K * matern.kernel(R*Tcomponent^(beta/2), nu = nu)
    }
    K <- Tcomponent^(beta*d/2+delta) * K
  } else {
    n.new <- NROW(x)
    if (length(tnew) != n.new) {
      stop("Length of mesh size should be the number of rows of inputs")
    }
    if(length(unique(t))==1){
      Tcomponent=1
    }else{
      Tcomponent <- 1/(outer(t, tnew, "-")^2/theta_t + 1)
    }
    # Rt <- sqrt(distance(t/sqrt(theta_t), tnew/sqrt(theta_t)))
    K <- matrix(1, n, n.new)
    for (i in 1:d) {
      R <- sqrt(distance(X[, i] / theta_x[i], x[, i] / theta_x[i]))
      K <- K * matern.kernel(R*Tcomponent^(beta/2), nu = nu)
    }
    K <- Tcomponent^(beta*d/2+delta) * K
  }
  return(K)
}

#' @title Fit GP with nonseparable squared exponential kernel
#'
#' @description The function fits the GP model with nonseparable squared exponential kernel for
#' multi-fidelity designs with tuning parameters.
#'
#' @param X A vector or matrix of input locations.
#' @param y A vector or matrix of response values.
#' @param t A vector of tuning parameters.
#' @param g Nugget term for numerical stability. Default is \code{sqrt(.Machine$double.eps)}.
#' @param constant A logical indicating for constant mean of GP (\code{constant=TRUE}) or zero mean (\code{constant=FALSE}). Default is \code{TRUE}.
#' @param p Quantile on distances. Default is 0.1.
#' @param min_cor minimal correlation between two design points at the defined p quantile distance. Default is 0.02.
#' @param max_cor maximal correlation between two design points at the defined (1-p) quantile distance. Default is 0.3.
#' @param init Optional vector of initial parameter values for optimization. Default is \code{NULL}.
#' @param lower Lower bound of optimization. Default is \code{NULL}.
#' @param upper Upper bound of optimization. Default is \code{NULL}.
#' @param multi.start Number of random starting points for optimization. Default is 1.
#'
#' @return A fitted model list containing:
#' \itemize{
#'   \item \code{K}: A matrix of covariance.
#'   \item \code{Ki}: A matrix of covariance inverse.
#'   \item \code{X}: A copy of X.
#'   \item \code{y}: A copy of y.
#'   \item \code{t}: A copy of \code{t}.
#'   \item \code{theta_x}: A vector of optimized lengthscale hyperparameter for \code{X}.
#'   \item \code{theta_t}: A vector of optimized lengthscale hyperparameter for \code{t}.
#'   \item \code{beta}: Optimized interaction hyperparameter \eqn{\beta}.
#'   \item \code{delta}: Optimized hyperparameter \eqn{\delta}.
#'   \item \code{g}: A copy of g.
#'   \item \code{mu.hat}: Optimized constant mean. If \code{constant=FALSE}, 0.
#'   \item \code{tau2hat}: Estimated scale hyperparameter.
#'   \item \code{constant}: A copy of constant.
#' }
#' @importFrom stats optim
#' @importFrom plgp distance
#' @noRd
#'

GP.nonsep.sqex <- function(X, y, t, g = sqrt(.Machine$double.eps), constant = FALSE,
                           p=0.1, min_cor = 0.02, max_cor = 0.3,
                           init=NULL, lower=NULL, upper=NULL, multi.start=1) { # p=0.05 for hetGP
  if (constant) {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)

    theta_bounds <- theta_bounds_sqex(X, t, p = p, min_cor = min_cor, max_cor = max_cor)

    if(is.null(lower)) {
      Xlower <- theta_bounds$lower * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2
      tlower <- 1e-2 * max(t) # theta_bounds_t$lower * (max(max(t)-min(t), g))^2
      lower <- c(Xlower, tlower, g, 0.5)
    }
    if(is.null(upper)) {
      Xupper <- theta_bounds$upper * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2
      tupper <- 1e1 * max(t) # theta_bounds_t$upper * (max(max(t)-min(t), g))^2
      upper <- c(Xupper, tupper, 1-g, 2-g)
    }
    if(is.null(init)) {
      init <- c(sqrt(lower*upper)[-c(length(lower)-1, length(lower))], 0.5, 1)
    }

    lower[init < lower] <- init[init < lower]
    upper[init > upper] <- init[init > upper]

    n <- nrow(X)

    nlsep <- function(par, X, Y, tt) {
      K <- cor.sep.sqex(X, t=tt, param=par)
      Ki <- solve(K + diag(g, n))
      ldetK <- determinant(K + diag(g, n), logarithm = TRUE)$modulus

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      tau2hat <- drop(t(Y - mu.hat) %*% Ki %*% (Y - mu.hat) / n)
      ll <- -(n / 2) * log(tau2hat) - (1 / 2) * ldetK
      return(drop(-ll))
    }

    gradnlsep <- function(par, X, Y, tt) {
      d <-ncol(X)
      dd <- d-1
      n <- nrow(X)
      theta <- par[1:d]
      theta_t <- par[d+1]
      beta <- par[d+2]
      delta <- par[d+3]

      K <- cor.sep.sqex(X, t=tt, param=par)
      Ki <- solve(K + diag(g, n))

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      KiY <- Ki %*% (Y - mu.hat)

      Rt <- sqrt(distance(tt / sqrt(theta_t)))
      u <- (1 + Rt^2)^(-1)
      v <- matrix(rowSums(sapply(1:d, function(i) distance(X[, i])/theta[i])), ncol=n)

      du_dtheta_t <- (Rt^2 / theta_t) * u^2
      dk_dbeta <- K * log(u) * ((dd+1)/2 - u^beta*v)
      dk_ddelta <- K * log(u)

      dlltheta <- rep(NA, length(par))

      dk_du <- ((dd+1)/2 * beta + delta - beta*u^beta*v) * u^((dd+1)/2*beta+delta-1) * exp(-u^beta * v)
      dk_dv <- -u^((dd+3)/2*beta+delta) * exp(-u^beta * v)

      dk_dtheta_t <- dk_du * du_dtheta_t

      for (i in 1:d) {
        dk_dtheta_i <- -distance(X[,i] / theta[i]) * dk_dv
        dlltheta[i] <- (n / 2) * t(KiY) %*% dk_dtheta_i %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_i))
      }

      dlltheta[d+1] <- (n / 2) * t(KiY) %*% dk_dtheta_t %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_t))
      dlltheta[d+2] <- (n / 2) * t(KiY) %*% dk_dbeta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dbeta))
      dlltheta[d+3] <- (n / 2) * t(KiY) %*% dk_ddelta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_ddelta))

      return(-c(dlltheta))
    }

    if(multi.start > 1){
      start <- randomLHS(multi.start - 1, ncol(X)+3)
      start <- t(t(start) * (upper - lower) + lower)
      start <- rbind(init, start)
      for(i in 1:nrow(start)) {
        outi <- optim(start[i,], nlsep, gradnlsep,
                      method = "L-BFGS-B", lower = lower, upper = upper,
                      X = X, Y = y, tt = t)
        if(i == 1) {
          out <- outi
        }else if(outi$value < out$value) {
          out <- outi
        }
      }
    }else{
      out <- optim(init, nlsep, gradnlsep,
                   method = "L-BFGS-B", lower = lower, upper = upper,
                   X = X, Y = y, tt = t)
    }

    theta_x <- out$par[1:ncol(X)]
    theta_t <- out$par[ncol(X)+1]
    beta <- out$par[ncol(X)+2]
    delta <- out$par[ncol(X)+3]
    K <- cor.sep.sqex(X, t=t, param=out$par)
    Ki <- solve(K + diag(g, n))
    one.vec <- matrix(1, ncol = 1, nrow = n)
    mu.hat <- drop((t(one.vec) %*% Ki %*% y) / (t(one.vec) %*% Ki %*% one.vec))
    tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / nrow(X))

    return(list(K = K, Ki = Ki, X = X, y = y, t=t, theta_x=theta_x, theta_t=theta_t, beta=beta, delta=delta, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  } else {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)

    theta_bounds <- theta_bounds_sqex(X, t, p = p, min_cor = min_cor, max_cor = max_cor)

    if(is.null(lower)) {
      Xlower <- theta_bounds$lower * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2
      tlower <- 1e-2 * max(t) # theta_bounds_t$lower * (max(max(t)-min(t), g))^2
      lower <- c(Xlower, tlower, g, 0.5)
    }
    if(is.null(upper)) {
      Xupper <- theta_bounds$upper * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2
      tupper <- 1e1 * max(t) # theta_bounds_t$upper * (max(max(t)-min(t), g))^2
      upper <- c(Xupper, tupper, 1-g, 2-g)
    }
    if(is.null(init)) {
      init <- c(sqrt(lower*upper)[-c(length(lower)-1, length(lower))], 0.5, 1)
    }

    lower[init < lower] <- init[init < lower]
    upper[init > upper] <- init[init > upper]

    n <- nrow(X)

    nlsep <- function(par, X, Y, tt) {
      K <- cor.sep.sqex(X, t=tt, param=par)
      Ki <- solve(K + diag(g, n))
      ldetK <- determinant(K + diag(g, n), logarithm = TRUE)$modulus

      mu.hat <- 0

      tau2hat <- drop(t(Y - mu.hat) %*% Ki %*% (Y - mu.hat) / n)
      ll <- -(n / 2) * log(tau2hat) - (1 / 2) * ldetK
      return(drop(-ll))
    }

    gradnlsep <- function(par, X, Y, tt) {
      d <-ncol(X)
      dd <- d-1
      n <- nrow(X)
      theta <- par[1:d]
      theta_t <- par[d+1]
      beta <- par[d+2]
      delta <- par[d+3]

      K <- cor.sep.sqex(X, t=tt, param=par)
      Ki <- solve(K + diag(g, n))

      mu.hat <- 0

      KiY <- Ki %*% (Y - mu.hat)

      Rt <- sqrt(distance(tt / sqrt(theta_t)))
      u <- (1 + Rt^2)^(-1)
      v <- matrix(rowSums(sapply(1:d, function(i) distance(X[, i])/theta[i])), ncol=n)


      du_dtheta_t <- (Rt^2 / theta_t) * u^2
      dk_dbeta <- K * log(u) * ((dd+1)/2 - u^beta*v)
      dk_ddelta <- K * log(u)

      dlltheta <- rep(NA, length(par))

      dk_du <- ((dd+1)/2 * beta + delta - beta*u^beta*v) * u^((dd+1)/2*beta+delta-1) * exp(-u^beta * v)
      dk_dv <- -u^((dd+3)/2*beta+delta) * exp(-u^beta * v)

      dk_dtheta_t <- dk_du * du_dtheta_t

      for (i in 1:d) {
        dk_dtheta_i <- -distance(X[,i] / theta[i]) * dk_dv
        dlltheta[i] <- (n / 2) * t(KiY) %*% dk_dtheta_i %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_i))
      }

      dlltheta[d+1] <- (n / 2) * t(KiY) %*% dk_dtheta_t %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_t))
      dlltheta[d+2] <- (n / 2) * t(KiY) %*% dk_dbeta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dbeta))
      dlltheta[d+3] <- (n / 2) * t(KiY) %*% dk_ddelta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_ddelta))

      return(-c(dlltheta))
    }

    if(multi.start > 1){
      start <- randomLHS(multi.start - 1, ncol(X)+3)
      start <- t(t(start) * (upper - lower) + lower)
      start <- rbind(init, start)
      for(i in 1:nrow(start)) {
        outi <- optim(start[i,], nlsep, gradnlsep,
                      method = "L-BFGS-B", lower = lower, upper = upper,
                      X = X, Y = y, tt = t)
        if(i == 1) {
          out <- outi
        }else if(outi$value < out$value) {
          out <- outi
        }
      }
    }else{
      out <- optim(init, nlsep, gradnlsep,
                   method = "L-BFGS-B", lower = lower, upper = upper,
                   X = X, Y = y, tt = t)
    }

    theta_x <- out$par[1:ncol(X)]
    theta_t <- out$par[ncol(X)+1]
    beta <- out$par[ncol(X)+2]
    delta <- out$par[ncol(X)+3]
    K <- cor.sep.sqex(X, t=t, param=out$par)
    Ki <- solve(K + diag(g, n))
    mu.hat <- 0
    tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / nrow(X))

    return(list(K = K, Ki = Ki, X = X, y = y, t=t, theta_x=theta_x, theta_t=theta_t, beta=beta, delta=delta, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  }
}

#' @title Fit GP with nonseparable Matern kernel
#'
#' @description The function fits the GP model with nonseparable Matern kernel for
#' multi-fidelity designs with tuning parameters.
#'
#' @param X A vector or matrix of input locations.
#' @param y A vector or matrix of response values.
#' @param nu A numerical value of smoothness hyperparameter. It should be 1.5, or 2.5.
#' @param t A vector of tuning parameters.
#' @param g Nugget term for numerical stability. Default is \code{sqrt(.Machine$double.eps)}.
#' @param constant A logical indicating for constant mean of GP (\code{constant=TRUE}) or zero mean (\code{constant=FALSE}). Default is \code{TRUE}.
#' @param p Quantile on distances. Default is 0.1.
#' @param min_cor minimal correlation between two design points at the defined p quantile distance. Default is 0.02.
#' @param max_cor maximal correlation between two design points at the defined (1-p) quantile distance. Default is 0.3.
#' @param init Optional vector of initial parameter values for optimization. Default is \code{NULL}.
#' @param lower Lower bound of optimization. Default is \code{NULL}.
#' @param upper Upper bound of optimization. Default is \code{NULL}.
#' @param multi.start Number of random starting points for optimization. Default is 1.
#'
#' @return A fitted model list containing:
#' \itemize{
#'   \item \code{K}: A matrix of covariance.
#'   \item \code{Ki}: A matrix of covariance inverse.
#'   \item \code{X}: A copy of X.
#'   \item \code{y}: A copy of y.
#'   \item \code{nu}: A copy of nu.
#'   \item \code{t}: A copy of \code{t}.
#'   \item \code{theta_x}: A vector of optimized lengthscale hyperparameter for \code{X}.
#'   \item \code{theta_t}: A vector of optimized lengthscale hyperparameter for \code{t}.
#'   \item \code{beta}: Optimized interaction hyperparameter \eqn{\beta}.
#'   \item \code{delta}: Optimized hyperparameter \eqn{\delta}.
#'   \item \code{g}: A copy of g.
#'   \item \code{mu.hat}: Optimized constant mean. If \code{constant=FALSE}, 0.
#'   \item \code{tau2hat}: Estimated scale hyperparameter \eqn{\tau^2}.
#'   \item \code{constant}: A copy of constant.
#' }
#' @importFrom stats optim
#' @importFrom plgp distance
#' @importFrom lhs randomLHS
#' @noRd
#'

GP.nonsep.matern <- function(X, y, nu, t, g = sqrt(.Machine$double.eps), constant = FALSE,
                             p=0.1, min_cor = 0.02, max_cor = 0.3,
                             init=NULL, lower=NULL, upper=NULL, multi.start=1) {
  if (constant) {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)

    theta_bounds <- theta_bounds_matern(X, t, nu=nu, p = p, min_cor = min_cor, max_cor = max_cor)

    if(is.null(lower)) {
      Xlower <- theta_bounds$lower * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])
      tlower <- 1e-2 * max(t) # theta_bounds_t$lower * (max(max(t)-min(t), g))^2
      lower <- c(Xlower, tlower, g, 0.5)
    }
    if(is.null(upper)) {
      Xupper <- theta_bounds$upper * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])
      tupper <- 1e1 * max(t) # theta_bounds_t$upper * (max(max(t)-min(t), g))^2
      upper <- c(Xupper, tupper, 0.01-g, 1-g)
    }
    if(is.null(init)) {
      init <- c(sqrt(lower*upper)[-c(length(lower)-1, length(lower))], 0.005, 0.5)
    }

    lower[init < lower] <- init[init < lower]
    upper[init > upper] <- init[init > upper]

    n <- nrow(X)

    nlsep <- function(par, X, Y, tt, nu) {
      K <- cor.sep.matern(X, t=tt, param=par, nu=nu)
      Ki <- solve(K + diag(g, n))
      ldetK <- determinant(K + diag(g, n), logarithm = TRUE)$modulus

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      tau2hat <- drop(t(Y - mu.hat) %*% Ki %*% (Y - mu.hat) / n)
      ll <- -(n / 2) * log(tau2hat) - (1 / 2) * ldetK
      return(drop(-ll))
    }

    gradnlsep <- function(par, X, Y, tt, nu) {
      d <-ncol(X)
      dd <- d-1
      n <- nrow(X)
      theta <- par[1:d]
      theta_t <- par[d+1]
      beta <- par[d+2]
      delta <- par[d+3]

      K <- cor.sep.matern(X, t=tt, param=par, nu=nu)
      Ki <- solve(K + diag(g, n))

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      KiY <- Ki %*% (Y - mu.hat)

      Rt <- sqrt(distance(tt / sqrt(theta_t)))
      u <- (1 + Rt^2)^(-1 / 2)

      if (nu == 1.5) {
        v <- sqrt(3) * sapply(1:d, function(i) distance(X[, i])/theta[i])
        u_beta_2_v <- c(u^(beta/2)) * v

        du_dtheta_t <- (Rt^2 / theta_t) * u^2
        dk_dbeta <- K * 1/2 * log(u) * ((dd + 1) - matrix(rowSums(u_beta_2_v^2/(1+u_beta_2_v)), ncol=n))
        dk_ddelta <- K * log(u)

        dlltheta <- rep(NA, length(par))

        dk_du <- K * (((dd+1)/2*beta+delta)/u - beta/2*u^(-1)*matrix(rowSums(u_beta_2_v^2/(1+u_beta_2_v)), ncol=n))
        dk_dv_without_K <- c(u^(beta/2)) * u_beta_2_v/(1+u_beta_2_v)

        dk_dtheta_t <- dk_du * du_dtheta_t

        for (i in 1:d) {
          dk_dv <- -K * matrix(dk_dv_without_K[,i], ncol=n)
          dk_dtheta_i <- - matrix(v[,i], ncol=n) / theta[i] * dk_dv
          dlltheta[i] <- (n / 2) * t(KiY) %*% dk_dtheta_i %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_i))
        }

        dlltheta[d+1] <- (n / 2) * t(KiY) %*% dk_dtheta_t %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_t))
        dlltheta[d+2] <- (n / 2) * t(KiY) %*% dk_dbeta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dbeta))
        dlltheta[d+3] <- (n / 2) * t(KiY) %*% dk_ddelta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_ddelta))
      } else {
        v <- sqrt(5) * sapply(1:d, function(i) distance(X[, i])/theta[i])
        u_beta_2_v <- c(u^(beta/2)) * v

        du_dtheta_t <- (Rt^2 / theta_t) * u^2
        dk_dbeta <- K * 1/2 * log(u) * ((dd + 1) - matrix(rowSums((u_beta_2_v^2*(1+u_beta_2_v))/(3*(1+u_beta_2_v+u_beta_2_v^2/3))), ncol=n))
        dk_ddelta <- K * log(u)

        dlltheta <- rep(NA, length(par))

        dk_du <- K * (((dd+1)/2*beta+delta)/u - beta/2*u^(-1)*matrix(rowSums((u_beta_2_v^2*(1+u_beta_2_v))/(3*(1+u_beta_2_v+u_beta_2_v^2/3))), ncol=n))
        dk_dv_without_K <- c(u^(beta/2)) * u_beta_2_v * (1+u_beta_2_v) / (3*(1+u_beta_2_v+u_beta_2_v^2/3))

        dk_dtheta_t <- dk_du * du_dtheta_t

        for (i in 1:d) {
          dk_dv <- -K * matrix(dk_dv_without_K[,i], ncol=n)
          dk_dtheta_i <- - matrix(v[,i], ncol=n) / theta[i] * dk_dv
          dlltheta[i] <- (n / 2) * t(KiY) %*% dk_dtheta_i %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_i))
        }

        dlltheta[d+1] <- (n / 2) * t(KiY) %*% dk_dtheta_t %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_t))
        dlltheta[d+2] <- (n / 2) * t(KiY) %*% dk_dbeta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dbeta))
        dlltheta[d+3] <- (n / 2) * t(KiY) %*% dk_ddelta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_ddelta))
      }

      return(-c(dlltheta))
    }

    if(multi.start > 1){
      start <- randomLHS(multi.start - 1, ncol(X)+3)
      start <- t(t(start) * (upper - lower) + lower)
      start <- rbind(init, start)
      for(i in 1:nrow(start)) {
        outi <- optim(start[i,], nlsep, gradnlsep,
                      method = "L-BFGS-B", lower = lower, upper = upper, X = X, Y = y, tt = t, nu=nu)
        if(i == 1) {
          out <- outi
        }else if(outi$value < out$value) {
          out <- outi
        }
      }
    }else{
      out <- optim(init, nlsep, gradnlsep,
                   method = "L-BFGS-B", lower = lower, upper = upper,
                   X = X, Y = y, tt = t, nu=nu)
    }

    theta_x <- out$par[1:ncol(X)]
    theta_t <- out$par[ncol(X)+1]
    beta <- out$par[ncol(X)+2]
    delta <- out$par[ncol(X)+3]
    K <- cor.sep.matern(X, t=t, param=out$par, nu=nu)
    Ki <- solve(K + diag(g, n))
    one.vec <- matrix(1, ncol = 1, nrow = n)
    mu.hat <- drop((t(one.vec) %*% Ki %*% y) / (t(one.vec) %*% Ki %*% one.vec))
    tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / nrow(X))

    return(list(K = K, Ki = Ki, X = X, y = y, nu=nu, t=t, theta_x=theta_x, theta_t=theta_t, beta=beta, delta=delta, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  } else {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)

    theta_bounds <- theta_bounds_matern(X, t, nu=nu, p = p, min_cor = min_cor, max_cor = max_cor)

    if(is.null(lower)) {
      Xlower <- theta_bounds$lower * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])
      tlower <- 1e-2 * max(t) # theta_bounds_t$lower * (max(max(t)-min(t), g))^2
      lower <- c(Xlower, tlower, g, 0.5)
    }
    if(is.null(upper)) {
      Xupper <- theta_bounds$upper * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])
      tupper <- 1e1 * max(t) # theta_bounds_t$upper * (max(max(t)-min(t), g))^2
      upper <- c(Xupper, tupper, 1-g, 2-g)
    }
    if(is.null(init)) {
      init <- c(sqrt(lower*upper)[-c(length(lower)-1, length(lower))], 0.5, 1)
    }

    lower[init < lower] <- init[init < lower]
    upper[init > upper] <- init[init > upper]

    n <- nrow(X)

    nlsep <- function(par, X, Y, tt, nu) {
      K <- cor.sep.matern(X, t=tt, param=par, nu=nu)
      Ki <- solve(K + diag(g, n))
      ldetK <- determinant(K + diag(g, n), logarithm = TRUE)$modulus

      mu.hat <- 0

      tau2hat <- drop(t(Y - mu.hat) %*% Ki %*% (Y - mu.hat) / n)
      ll <- -(n / 2) * log(tau2hat) - (1 / 2) * ldetK
      return(drop(-ll))
    }

    gradnlsep <- function(par, X, Y, tt, nu) {
      d <-ncol(X)
      dd <- d-1
      n <- nrow(X)
      theta <- par[1:d]
      theta_t <- par[d+1]
      beta <- par[d+2]
      delta <- par[d+3]

      K <- cor.sep.matern(X, t=tt, param=par, nu=nu)
      Ki <- solve(K + diag(g, n))

      mu.hat <- 0

      KiY <- Ki %*% (Y - mu.hat)

      Rt <- sqrt(distance(tt / sqrt(theta_t)))
      u <- (1 + Rt^2)^(-1 / 2)

      if (nu == 1.5) {
        v <- sqrt(3) * sapply(1:d, function(i) distance(X[, i])/theta[i])
        u_beta_2_v <- c(u^(beta/2)) * v

        du_dtheta_t <- (Rt^2 / theta_t) * u^2
        dk_dbeta <- K * 1/2 * log(u) * ((dd + 1) - matrix(rowSums(u_beta_2_v^2/(1+u_beta_2_v)), ncol=n))
        dk_ddelta <- K * log(u)

        dlltheta <- rep(NA, length(par))

        dk_du <- K * (((dd+1)/2*beta+delta)/u - beta/2*u^(-1)*matrix(rowSums(u_beta_2_v^2/(1+u_beta_2_v)), ncol=n))
        dk_dv_without_K <- c(u^(beta/2)) * u_beta_2_v/(1+u_beta_2_v)

        dk_dtheta_t <- dk_du * du_dtheta_t

        for (i in 1:d) {
          dk_dv <- -K * matrix(dk_dv_without_K[,i], ncol=n)
          dk_dtheta_i <- - matrix(v[,i], ncol=n) / theta[i] * dk_dv
          dlltheta[i] <- (n / 2) * t(KiY) %*% dk_dtheta_i %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_i))
        }

        dlltheta[d+1] <- (n / 2) * t(KiY) %*% dk_dtheta_t %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_t))
        dlltheta[d+2] <- (n / 2) * t(KiY) %*% dk_dbeta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dbeta))
        dlltheta[d+3] <- (n / 2) * t(KiY) %*% dk_ddelta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_ddelta))
      } else {
        v <- sqrt(5) * sapply(1:d, function(i) distance(X[, i])/theta[i])
        u_beta_2_v <- c(u^(beta/2)) * v

        du_dtheta_t <- (Rt^2 / theta_t) * u^2
        dk_dbeta <- K * 1/2 * log(u) * ((dd + 1) - matrix(rowSums((u_beta_2_v^2*(1+u_beta_2_v))/(3*(1+u_beta_2_v+u_beta_2_v^2/3))), ncol=n))
        dk_ddelta <- K * log(u)

        dlltheta <- rep(NA, length(par))

        dk_du <- K * (((dd+1)/2*beta+delta)/u - beta/2*u^(-1)*matrix(rowSums((u_beta_2_v^2*(1+u_beta_2_v))/(3*(1+u_beta_2_v+u_beta_2_v^2/3))), ncol=n))
        dk_dv_without_K <- c(u^(beta/2)) * u_beta_2_v * (1+u_beta_2_v) / (3*(1+u_beta_2_v+u_beta_2_v^2/3))

        dk_dtheta_t <- dk_du * du_dtheta_t

        for (i in 1:d) {
          dk_dv <- -K * matrix(dk_dv_without_K[,i], ncol=n)
          dk_dtheta_i <- - matrix(v[,i], ncol=n) / theta[i] * dk_dv
          dlltheta[i] <- (n / 2) * t(KiY) %*% dk_dtheta_i %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_i))
        }

        dlltheta[d+1] <- (n / 2) * t(KiY) %*% dk_dtheta_t %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_t))
        dlltheta[d+2] <- (n / 2) * t(KiY) %*% dk_dbeta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dbeta))
        dlltheta[d+3] <- (n / 2) * t(KiY) %*% dk_ddelta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_ddelta))
      }

      return(-c(dlltheta))
    }

    if(multi.start > 1){
      start <- randomLHS(multi.start - 1, ncol(X)+3)
      start <- t(t(start) * (upper - lower) + lower)
      start <- rbind(init, start)
      for(i in 1:nrow(start)) {
        outi <- optim(start[i,], nlsep, gradnlsep,
                      method = "L-BFGS-B", lower = lower, upper = upper, X = X, Y = y, tt = t, nu=nu)
        if(i == 1) {
          out <- outi
        }else if(outi$value < out$value) {
          out <- outi
        }
      }
    }else{
      out <- optim(init, nlsep, gradnlsep,
                   method = "L-BFGS-B", lower = lower, upper = upper,
                   X = X, Y = y, tt = t, nu=nu)
    }

    theta_x <- out$par[1:ncol(X)]
    theta_t <- out$par[ncol(X)+1]
    beta <- out$par[ncol(X)+2]
    delta <- out$par[ncol(X)+3]
    K <- cor.sep.matern(X, t=t, param=out$par, nu=nu)
    Ki <- solve(K + diag(g, n))
    mu.hat <- 0
    tau2hat <- drop(t(y) %*% Ki %*% y / n)

    return(list(K = K, Ki = Ki, X = X, y = y, nu=nu, t=t, theta_x=theta_x, theta_t=theta_t, beta=beta, delta=delta, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  }
}

#' @title Predictive posterior mean and variance with nonseparable squared exponential kernel.
#'
#' @param fit Fitted GP object from \code{GP.nonsep.sqex}.
#' @param xnew A vector or matrix of new input locations to predict.
#' @param tnew A vector or matrix of new tuning parameters to predict.
#'
#' @return A list with predictive posterior containing:
#' \itemize{
#'   \item \code{mu}: A vector of predictive posterior mean.
#'   \item \code{sig2}: A vector of predictive posterior variance.
#' }
#'
#' @noRd
#' @keywords internal
#'

pred.GP.nonsep <- function(fit, xnew, tnew) {
  xnew <- as.matrix(xnew)

  Ki <- fit$Ki
  theta_x <- fit$theta_x
  theta_t <- fit$theta_t
  beta <- fit$beta
  delta <- fit$delta
  g <- fit$g
  X <- fit$X
  y <- fit$y
  t <- fit$t
  tau2hat <- fit$tau2hat
  mu.hat <- fit$mu.hat

  KXX <- cor.sep.sqex(xnew, t=tnew, param=c(theta_x,theta_t,beta,delta))
  KX <- cor.sep.sqex(xnew, X, t=tnew, tnew=t, param=c(theta_x,theta_t,beta,delta))

  mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
  Sigmap2 <- pmax(0, diag(tau2hat * (KXX - KX %*% Ki %*% t(KX))))

  return(list(mu = mup2, sig2 = Sigmap2))
}

#' @title Predictive posterior mean and variance with nonseparable Matern kernel.
#'
#' @param fit Fitted GP object from \code{GP.nonsep.sqex}.
#' @param xnew A vector or matrix of new input locations to predict.
#' @param tnew A vector or matrix of new tuning parameters to predict.
#'
#' @return A list with predictive posterior containing:
#' \itemize{
#'   \item \code{mu}: A vector of predictive posterior mean.
#'   \item \code{sig2}: A vector of predictive posterior variance.
#' }
#'
#' @noRd
#' @keywords internal
#'

pred.matGP.nonsep <- function(fit, xnew, tnew) {
  xnew <- as.matrix(xnew)

  Ki <- fit$Ki
  nu <- fit$nu
  theta_x <- fit$theta_x
  theta_t <- fit$theta_t
  beta <- fit$beta
  delta <- fit$delta
  g <- fit$g
  X <- fit$X
  y <- fit$y
  t <- fit$t
  tau2hat <- fit$tau2hat
  mu.hat <- fit$mu.hat

  KXX <- cor.sep.matern(xnew, t=tnew, param=c(theta_x,theta_t,beta,delta), nu=nu)
  KX <- cor.sep.matern(xnew, X, t=tnew, tnew=t, param=c(theta_x,theta_t,beta,delta), nu=nu)

  mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
  Sigmap2 <- pmax(0, diag(tau2hat * (KXX - KX %*% Ki %*% t(KX))))

  return(list(mu = mup2, sig2 = Sigmap2))
}
