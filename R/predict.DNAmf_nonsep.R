#' @title Closed-form predictive mean for nonseparable squared exponential kernel
#'
#' @description The function computes the closed-form predictive mean
#' at each level for the DNAmf model with the nonseparable squared exponential kernel.
#'
#' @param x A vector or matrix of new input locations to predict.
#' @param mu.hat Estimated GP constant mean of \eqn{f}.
#' @param theta_x A vector of optimized lengthscale hyperparameter for \code{X}.
#' @param X2 A vector or matrix of input locations of \eqn{f}.
#' @param c_i A vector of the effects of tuning parameters.
#' @param a A vector of precomputed weights \eqn{K^{-1}(Y_{-1} - \alpha \mathbf{1}_{N_{-1}})}.
#' @param w1.x2 A vector of outputs from previous level \eqn{c(Y_{-L})}.
#' @param previous_mu Posterior mean from previous level \eqn{\mu^*_{l-1}(\mathbf{x})}.
#' @param previous_sig2 Posterior variance from previous level \eqn{\sigma^{*2}_{l-1}(\mathbf{x})}.
#' @param d Input dimension.
#' @param delta Hyperparameter \eqn{\delta}.
#' @param beta Interaction hyperparameter \eqn{\beta}.
#'
#' @return A computed predictive mean vector.
#'
#' @importFrom plgp distance
#' @keywords internal
#' @noRd
#'

h1_sqex <- function(x, mu.hat, theta_x, X2, c_i, a, w1.x2, previous_mu, previous_sig2, d, delta, beta) {

  dist_matrix <- distance(t(t(x) / sqrt(theta_x[-(d + 1)])), t(t(X2) / sqrt(theta_x[-(d + 1)])))

  exp_dist <- exp(-sweep(dist_matrix, 2, c_i, `*`))

  sqrt_term <- sqrt(theta_x[d + 1] / (theta_x[d + 1] + 2 * previous_sig2))

  diff_mu_w1 <- drop(outer(previous_mu, w1.x2, FUN = "-"))
  denom <- theta_x[d + 1] + 2 * previous_sig2
  exp_mu <- exp(-sweep((diff_mu_w1^2 / denom), 2, c_i, `*`))

  weights <- a * c_i^((d+1)/2 + delta/beta)
  pred <- mu.hat + (exp_dist * sqrt_term * exp_mu) %*% weights

  return(pred)
}

#' @title Closed-form predictive variance for single point for nonseparable squared exponential kernel
#'
#' @description The function computes the closed-form predictive variance
#' at each level for a single prediction location with the nonseparable squared exponential kernel.
#'
#' @param x A vector of new input location to predict.
#' @param mu.hat Estimated GP constant mean of \eqn{f}.
#' @param theta_x A vector of optimized lengthscale hyperparameter for \code{X}.
#' @param X2 A vector or matrix of input locations of \eqn{f}.
#' @param c_i A vector of the effects of tuning parameters.
#' @param a A vector of precomputed weights \eqn{K^{-1}(Y_{-1} - \alpha \mathbf{1}_{N_{-1}})}.
#' @param w1.x2 A vector of outputs from previous level \eqn{c(Y_{-L})}.
#' @param previous_mu Posterior mean from previous level \eqn{\mu^*_{l-1}(\mathbf{x})}.
#' @param previous_sig2 Posterior variance from previous level \eqn{\sigma^{*2}_{l-1}(\mathbf{x})}.
#' @param d Input dimension.
#' @param delta Hyperparameter \eqn{\delta}.
#' @param beta Interaction hyperparameter \eqn{\beta}.
#' @param tau2hat Estimated scale hyperparameter \eqn{\tau^2}.
#' @param Ci Inverse covariance matrix \eqn{K^{-1}}.
#' @param predy Predictive mean at current level from \code{h1_sqex}.
#' @param c_sum An outer sum of \eqn{c_i + c_j}.
#' @param diff_w1_sq A matrix of precomputed squared difference for \code{w1.x2}.
#' @param c_i_outer An outer product of \eqn{c_i c_j}.
#' @param power_term A precomputed term of \eqn{c_i^{\frac{d+1}{2}+\frac{\delta}{\beta}}}.
#'
#' @return A computed predictive variance scalar.
#'
#' @importFrom plgp distance
#' @keywords internal
#' @noRd
#'

h2_sqex_single <- function(x, mu.hat, theta_x, X2, c_i, a, w1.x2, previous_mu, previous_sig2,
                           d, delta, beta, tau2hat, Ci, predy, c_sum, diff_w1_sq, c_i_outer, power_term) {

  dist <- distance(t(x / sqrt(theta_x[-(d + 1)])), t(t(X2) / sqrt(theta_x[-(d + 1)])))
  v <- exp(-dist * t(c_i))
  mat1 <- drop(v %o% v)

  # Zeta components
  denom <- theta_x[d + 1] + 2 * c_sum * previous_sig2
  sqrt_term <- sqrt(theta_x[d + 1] / denom)

  term1 <- drop(outer(c_i * (w1.x2 - previous_mu)^2, c_i * (w1.x2 - previous_mu)^2, "+"))
  term2 <- (2 / theta_x[d + 1]) * c_i_outer * previous_sig2 * diff_w1_sq
  exp_term <- exp(-(term1 + term2) / denom)

  mat <- mat1 * sqrt_term * exp_term * power_term

  # Compute h2
  quad <- drop(t(a) %*% mat %*% a)
  trace <- sum(Ci * mat)
  h2_val <- tau2hat - (predy - mu.hat)^2 + quad - tau2hat * trace
  return(pmax(0, h2_val))
}

#' @title Closed-form predictive variance for nonseparable squared exponential kernel
#'
#' @description The function computes the closed-form predictive variance
#' at each level for the DNAmf model with the nonseparable squared exponential kernel.
#'
#' @param x A vector or matrix of new input locations to predict.
#' @param mu.hat Estimated GP constant mean of \eqn{f}.
#' @param theta_x A vector of optimized lengthscale hyperparameter for \code{X}.
#' @param X2 A vector or matrix of input locations of \eqn{f}.
#' @param c_i A vector of the effects of tuning parameters.
#' @param a A vector of precomputed weights \eqn{K^{-1}(Y_{-1} - \alpha \mathbf{1}_{N_{-1}})}.
#' @param w1.x2 A vector of outputs from previous level \eqn{c(Y_{-L})}.
#' @param previous_mu Posterior mean from previous level \eqn{\mu^*_{l-1}(\mathbf{x})}.
#' @param previous_sig2 Posterior variance from previous level \eqn{\sigma^{*2}_{l-1}(\mathbf{x})}.
#' @param d Input dimension.
#' @param delta Hyperparameter \eqn{\delta}.
#' @param beta Interaction hyperparameter \eqn{\beta}.
#' @param tau2hat Estimated scale hyperparameter \eqn{\tau^2}.
#' @param Ci Inverse covariance matrix \eqn{K^{-1}}.
#' @param predy Predictive mean at current level from \code{h1_sqex}.
#'
#' @return A computed predictive variance vector.
#'
#' @keywords internal
#' @noRd
#'

h2_sqex <- function(x, mu.hat, theta_x, X2, c_i, a, w1.x2, previous_mu, previous_sig2, d, delta, beta, tau2hat, Ci, predy) {
  # Precompute terms outside the loop
  c_sum <- outer(c(c_i), c(c_i), "+")
  diff_w1_sq <- drop(outer(w1.x2, w1.x2, "-"))^2
  c_i_outer <- c_i %*% t(c_i)
  power_term <- (c_i_outer)^((d+1)/2 + delta/beta)

  # Apply h2_single to each test point
  predsig2 <- numeric(nrow(x))
  for (i in 1:nrow(x)) {
    predsig2[i] <- h2_sqex_single(x[i, ], mu.hat, theta_x, X2, c_i, a, w1.x2, previous_mu[i], previous_sig2[i],
                                  d, delta, beta, tau2hat, Ci, predy[i], c_sum, diff_w1_sq, c_i_outer, power_term)
  }
  return(predsig2)
}

#' @title Closed-form predictive mean for nonseparable Matern kernel
#'
#' @description The function computes the closed-form predictive mean
#' at each level for the DNAmf model with the nonseparable Matern kernel.
#'
#' @param x A vector or matrix of new input locations to predict.
#' @param mu.hat Estimated GP constant mean of \eqn{f}.
#' @param theta_x A vector of optimized lengthscale hyperparameter for \code{X}.
#' @param X2 A vector or matrix of input locations of \eqn{f}.
#' @param c_i A vector of the effects of tuning parameters.
#' @param nu A numerical value of smoothness hyperparameter. It should be 1.5, or 2.5.
#' @param a A vector of precomputed weights \eqn{K^{-1}(Y_{-1} - \alpha \mathbf{1}_{N_{-1}})}.
#' @param w1.x2 A vector of outputs from previous level \eqn{c(Y_{-L})}.
#' @param previous_mu Posterior mean from previous level \eqn{\mu^*_{l-1}(\mathbf{x})}.
#' @param previous_sig2 Posterior variance from previous level \eqn{\sigma^{*2}_{l-1}(\mathbf{x})}.
#' @param d Input dimension.
#' @param delta Hyperparameter \eqn{\delta}.
#' @param beta Interaction hyperparameter \eqn{\beta}.
#'
#' @return A computed predictive mean vector.
#'
#' @keywords internal
#' @noRd
#'

h1_matern <- function(x, mu.hat, theta_x, X2, c_i, nu, a, w1.x2, previous_mu, previous_sig2, d, delta, beta) {

  dist_matrix <- cor.mat(x, X2, theta_x[-(d + 1)], nu=nu)

  xi_matrix <- matrix(xifun(c = c_i, w = rep(w1.x2, each = nrow(x)),
                            m = rep(previous_mu, nrow(X2)), s = rep(previous_sig2, nrow(X2)),
                            theta = theta_x[d+1], nu = nu), nrow(x), nrow(X2))

  weights <- a * c_i^((d + 1) + 2 * delta / beta)
  pred <- mu.hat + rowSums((dist_matrix * xi_matrix) %*% weights)

  return(pred)
}

#' @title Closed-form predictive variance for single point for nonseparable Matern kernel
#'
#' @description The function computes the closed-form predictive variance
#' at each level for a single prediction location with the nonseparable Matern kernel.
#'
#' @param x A vector of new input location to predict.
#' @param mu.hat Estimated GP constant mean of \eqn{f}.
#' @param theta_x A vector of optimized lengthscale hyperparameter for \code{X}.
#' @param X2 A vector or matrix of input locations of \eqn{f}.
#' @param c_i A vector of the effects of tuning parameters.
#' @param nu A numerical value of smoothness hyperparameter. It should be 1.5, or 2.5.
#' @param a A vector of precomputed weights \eqn{K^{-1}(Y_{-1} - \alpha \mathbf{1}_{N_{-1}})}.
#' @param w1.x2 A vector of outputs from previous level \eqn{c(Y_{-L})}.
#' @param previous_mu Posterior mean from previous level \eqn{\mu^*_{l-1}(\mathbf{x})}.
#' @param previous_sig2 Posterior variance from previous level \eqn{\sigma^{*2}_{l-1}(\mathbf{x})}.
#' @param d Input dimension.
#' @param delta Hyperparameter \eqn{\delta}.
#' @param beta Interaction hyperparameter \eqn{\beta}.
#' @param tau2hat Estimated scale hyperparameter \eqn{\tau^2}.
#' @param Ci Inverse covariance matrix \eqn{K^{-1}}.
#' @param predy Predictive mean at current level from \code{h1_matern}.
#' @param power_term A precomputed term of \eqn{c_i^{d+1+\frac{2\delta}{\beta}}}.
#'
#' @return A computed predictive variance scalar.
#'
#' @keywords internal
#' @noRd
#'

h2_matern_single <- function(x, mu.hat, theta_x, X2, c_i, nu, a, w1.x2, previous_mu, previous_sig2,
                             d, delta, beta, tau2hat, Ci, predy, power_term) {

  dist <- cor.mat(x, X2, theta_x[-(d + 1)], nu=nu)
  mat1 <- drop(dist %o% dist)

  # Zeta components
  zeta_mat <- matrix(zetafun(c1 = rep(c_i, times = nrow(X2)), c2 = rep(c_i, each = nrow(X2)),
                             w1 = rep(w1.x2, times = nrow(X2)), w2 = rep(w1.x2, each  = nrow(X2)),
                             m = previous_mu, s = previous_sig2,
                             theta = theta_x[d+1], nu = nu), nrow(X2), nrow(X2))

  mat <- mat1 * zeta_mat * power_term

  # Compute h2
  quad <- drop(t(a) %*% mat %*% a)
  trace <- sum(Ci * mat)
  h2_val <- tau2hat - (predy - mu.hat)^2 + quad - tau2hat * trace
  return(pmax(0, h2_val))
}

#' @title Closed-form predictive variance for nonseparable Matern kernel
#'
#' @description The function computes the closed-form predictive variance
#' at each level for the DNAmf model with the nonseparable Matern kernel.
#'
#' @param x A vector or matrix of new input locations to predict.
#' @param mu.hat Estimated GP constant mean of \eqn{f}.
#' @param theta_x A vector of optimized lengthscale hyperparameter for \code{X}.
#' @param X2 A vector or matrix of input locations of \eqn{f}.
#' @param c_i A vector of the effects of tuning parameters.
#' @param nu Smoothness hyperparameter \eqn{\nu}.
#' @param a A vector of precomputed weights \eqn{K^{-1}(Y_{-1} - \alpha \mathbf{1}_{N_{-1}})}.
#' @param w1.x2 A vector of outputs from previous level \eqn{c(Y_{-L})}.
#' @param previous_mu Posterior mean from previous level \eqn{\mu^*_{l-1}(\mathbf{x})}.
#' @param previous_sig2 Posterior variance from previous level \eqn{\sigma^{*2}_{l-1}(\mathbf{x})}.
#' @param d Input dimension.
#' @param delta Hyperparameter \eqn{\delta}.
#' @param beta Interaction hyperparameter \eqn{\beta}.
#' @param tau2hat Estimated scale hyperparameter \eqn{\tau^2}.
#' @param Ci Inverse covariance matrix \eqn{K^{-1}}.
#' @param predy Predictive mean at current level from \code{h1_matern}.
#'
#' @return A computed predictive variance vector.
#'
#' @keywords internal
#' @noRd
#'

h2_matern <- function(x, mu.hat, theta_x, X2, c_i, nu, a, w1.x2, previous_mu, previous_sig2, d, delta, beta, tau2hat, Ci, predy) {
  # Precompute terms outside the loop
  c_i_outer <- c_i %*% t(c_i)
  power_term <- (c_i_outer)^((d+1) + 2*delta/beta)

  # Apply h2_single to each test point
  predsig2 <- numeric(nrow(x))
  for (i in 1:nrow(x)) {
    predsig2[i] <- h2_matern_single(x[i, ], mu.hat, theta_x, X2, c_i, nu, a, w1.x2, previous_mu[i], previous_sig2[i],
                                    d, delta, beta, tau2hat, Ci, predy[i], power_term)
  }
  return(predsig2)
}

#' @title Closed-form prediction for DNAmf model
#'
#' @description The function computes the closed-form posterior mean and variance for the DNAmf model
#' both at the fidelity levels used in model fitting and at any user-specified target fidelity level,
#' using the chosen nonseparable kernel.
#'
#' @param fit1 A fitted GP object for \eqn{f_1}.
#' @param fit2 A fitted GP object for \eqn{f}.
#' @param targett A numeric value of target tuning parameter to predict.
#' @param kernel A character specifying the kernel type to be used. Choices are \code{"sqex"}(nonseparable squared exponential kernel), \code{"matern1.5"}(nonseparable Matern kernel with \eqn{\nu=1.5}), or \code{"matern2.5"}(nonseparable Matern kernel with \eqn{\nu=2.5}). Default is \code{"sqex"}.
#' @param nn A vector specifying the number of design points at each fidelity level.
#' @param tt A vector of tuning parameters for each fidelity level.
#' @param nlevel The number of fidelity levels \eqn{L}.
#' @param x A vector or matrix of new input locations to predict.
#' @param XX A list containing a pseudo-complete inputs \code{X_star}(\eqn{\left\{\mathcal{X}^*_l\right\}_{l=1}^{L}}), an original inputs \code{X_list}(\eqn{\left\{\mathcal{X}_l\right\}_{l=1}^{L}}), and a pseudo inputs \code{X_tilde}(\eqn{\left\{\widetilde{\mathcal{X}}_l\right\}_{l=1}^{L}}) for non-nested design.
#' @param pseudo_yy A list containing a pseudo-complete outputs \code{y_star}(\eqn{\left\{\mathbf{y}^*_l\right\}_{l=1}^{L}}), an original outputs \code{y_list}(\eqn{\left\{\mathbf{y}_l\right\}_{l=1}^{L}}), and a pseudo outputs \code{y_tilde}(\eqn{\left\{\widetilde{\mathbf{y}}_l\right\}_{l=1}^{L}}) imputed by \code{\link{imputer}}.
#'
#' @return A list of predictive posterior mean and variance for each level containing:
#' \itemize{
#'   \item \code{mu_1}, \code{sig2_1}, ..., \code{mu_L}, \code{sig2_L}: A vector of predictive posterior mean and variance at each level.
#'   \item \code{mu}: A vector of predictive posterior mean at target tuning parameter.
#'   \item \code{sig2}: A vector of predictive posterior variance at target tuning parameter.
#' }
#'

closed_form <- function(fit1, fit2, targett, kernel, nn, tt, nlevel, x, XX=NULL, pseudo_yy=NULL){

  tt <- tt
  TT <- rep(tt[-1],nn[-1])

  theta_x <- fit2$theta_x
  theta_t <- fit2$theta_t
  beta <- fit2$beta
  delta <- fit2$delta
  g <- fit2$g
  d <- ncol(fit1$X)

  if(is.null(pseudo_yy)){
    ### prediction ###
    y <- fit2$y
    n <- length(y)
    Ci <- fit2$Ki
    tau2hat <- fit2$tau2hat
    mu.hat <- fit2$mu.hat
    X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
    w1.x2 <- fit2$X[, (d + 1)]
  }else{
    y <- do.call(rbind, pseudo_yy$y_star[-1])
    n <- length(y)

    # update K based on the new pseudo_yy
    L <- length(XX$X_star)
    X_m1 <- do.call(rbind, XX$X_star[-1])
    Y_mL <- do.call(rbind, lapply(2:L, function(l) {
      pseudo_yy$y_star[[l-1]][checkindices(XX$X_star[[l-1]], XX$X_star[[l]]), , drop = FALSE]
    }))
    t_m1 <- unlist(lapply(2:L, function(l) rep(tt[l], nrow(XX$X_star[[l]]))))

    if(kernel=="sqex"){
      K <- cor.sep.sqex(cbind(X_m1,Y_mL), t=t_m1, param=c(theta_x,theta_t,beta,delta))
    }else if(kernel=="matern1.5"){
      K <- cor.sep.matern(cbind(X_m1,Y_mL), t=t_m1, param=c(theta_x,theta_t,beta,delta), nu=1.5)
    }else if(kernel=="matern2.5"){
      K <- cor.sep.matern(cbind(X_m1,Y_mL), t=t_m1, param=c(theta_x,theta_t,beta,delta), nu=2.5)
    }
    Ci <- solve(K + diag(g, n))
    one.vec <- matrix(1, ncol = 1, nrow = n)
    mu.hat <- drop((t(one.vec) %*% Ci %*% y) / (t(one.vec) %*% Ci %*% one.vec))
    tau2hat <- drop(t(y - mu.hat) %*% Ci %*% (y - mu.hat) / n)
    X2 <- X_m1
    w1.x2 <- c(Y_mL)
  }
  a <- Ci %*% (y - mu.hat)
  x <- matrix(x, ncol = d)
  if(kernel=="sqex"){
    pred.fit <- pred.GP(fit1, x)
  }else if(kernel=="matern1.5"){
    pred.fit <- pred.matGP(fit1, x)
  }else if(kernel=="matern2.5"){
    pred.fit <- pred.matGP(fit1, x)
  }

  ### calculate the closed form ###

  if (nlevel > 2){
    pred_result <- list()
    pred_result[["mu_1"]] <- pred.fit$mu
    pred_result[["sig2_1"]] <- pred.fit$sig2

    for (j in 2:nlevel) { ### from level 2 to level of current design
      # mean & var
      if(kernel=="sqex"){
        c_i <- 1/((TT-tt[j])^2/theta_t + 1)^beta # c_i = c(t,t_i)
        predy_temp <- h1_sqex(x, mu.hat, theta_x, X2, c_i, a, w1.x2, pred_result[[paste0("mu_", j-1)]], pred_result[[paste0("sig2_", j-1)]], d, delta, beta)
        predsig2_temp <- h2_sqex(x, mu.hat, theta_x, X2, c_i, a, w1.x2, pred_result[[paste0("mu_", j-1)]], pred_result[[paste0("sig2_", j-1)]], d, delta, beta, tau2hat, Ci, predy_temp)
      }else if(kernel=="matern1.5"){
        c_i <- 1/((TT-tt[j])^2/theta_t + 1)^(beta/2) # c_i = c(t,t_i)
        predy_temp <- h1_matern(x, mu.hat, theta_x, X2, c_i, nu=1.5, a, w1.x2, pred_result[[paste0("mu_", j-1)]], pred_result[[paste0("sig2_", j-1)]], d, delta, beta)
        predsig2_temp <- h2_matern(x, mu.hat, theta_x, X2, c_i, nu=1.5, a, w1.x2, pred_result[[paste0("mu_", j-1)]], pred_result[[paste0("sig2_", j-1)]], d, delta, beta, tau2hat, Ci, predy_temp)
      }else if(kernel=="matern2.5"){
        c_i <- 1/((TT-tt[j])^2/theta_t + 1)^(beta/2) # c_i = c(t,t_i)
        predy_temp <- h1_matern(x, mu.hat, theta_x, X2, c_i, nu=2.5, a, w1.x2, pred_result[[paste0("mu_", j-1)]], pred_result[[paste0("sig2_", j-1)]], d, delta, beta)
        predsig2_temp <- h2_matern(x, mu.hat, theta_x, X2, c_i, nu=2.5, a, w1.x2, pred_result[[paste0("mu_", j-1)]], pred_result[[paste0("sig2_", j-1)]], d, delta, beta, tau2hat, Ci, predy_temp)
      }

      pred_result[[paste0("mu_", j)]] <- predy_temp
      pred_result[[paste0("sig2_", j)]] <- predsig2_temp
    }

    ### for target mesh size ###
    if(any(tt == targett)){ # if target mesh belongs to initial design, use previous outputs with that specific mesh size (for interpolation)
      pred_result[["mu"]] <- pred_result[[paste0("mu_", which(tt == targett))]]
      pred_result[["sig2"]] <- pred_result[[paste0("sig2_", which(tt == targett))]]
    }else{
      # find the smallest mesh which is larger than target t
      tt_place <- findInterval(-targett, -tt)
      # mean & var
      if(kernel=="sqex"){
        c_i <- 1/((TT-targett)^2/theta_t + 1)^beta # c_i = c(t,t_i)
        predy <- h1_sqex(x, mu.hat, theta_x, X2, c_i, a, w1.x2, pred_result[[paste0("mu_", tt_place)]], pred_result[[paste0("sig2_", tt_place)]], d, delta, beta)
        predsig2 <- h2_sqex(x, mu.hat, theta_x, X2, c_i, a, w1.x2, pred_result[[paste0("mu_", tt_place)]], pred_result[[paste0("sig2_", tt_place)]], d, delta, beta, tau2hat, Ci, predy)
      }else if(kernel=="matern1.5"){
        c_i <- 1/((TT-targett)^2/theta_t + 1)^(beta/2) # c_i = c(t,t_i)
        predy <- h1_matern(x, mu.hat, theta_x, X2, c_i, nu=1.5, a, w1.x2, pred_result[[paste0("mu_", tt_place)]], pred_result[[paste0("sig2_", tt_place)]], d, delta, beta)
        predsig2 <- h2_matern(x, mu.hat, theta_x, X2, c_i, nu=1.5, a, w1.x2, pred_result[[paste0("mu_", tt_place)]], pred_result[[paste0("sig2_", tt_place)]], d, delta, beta, tau2hat, Ci, predy)
      }else if(kernel=="matern2.5"){
        c_i <- 1/((TT-targett)^2/theta_t + 1)^(beta/2) # c_i = c(t,t_i)
        predy <- h1_matern(x, mu.hat, theta_x, X2, c_i, nu=2.5, a, w1.x2, pred_result[[paste0("mu_", tt_place)]], pred_result[[paste0("sig2_", tt_place)]], d, delta, beta)
        predsig2 <- h2_matern(x, mu.hat, theta_x, X2, c_i, nu=2.5, a, w1.x2, pred_result[[paste0("mu_", tt_place)]], pred_result[[paste0("sig2_", tt_place)]], d, delta, beta, tau2hat, Ci, predy)
      }

      pred_result[["mu"]] <- predy
      pred_result[["sig2"]] <- predsig2
    }
  } else {
    stop("Level should be larger than 2")
  }
  return(pred_result)
}

#' @title Predictive posterior mean and variance for DNAmf object with nonseparable kernel.
#'
#' @description The function computes the predictive posterior mean and variance
#' for the DNAmf model using closed-form expressions based on the chosen nonseparable kernel
#' at given new input locations.
#'
#' @seealso \code{\link{DNAmf}} for the user-level function.
#'
#' @details The \code{predict.DNAmf} function internally calls \code{\link{closed_form}},
#' which further calls \code{h1_sqex}, \code{h2_sqex}, \code{h2_sqex_single} for \code{kernel="sqex"},
#' or \code{h1_matern}, \code{h2_matern}, \code{h2_matern_single} for \code{kernel="matern1.5"} or\code{kernel="matern2.5"},
#' to recursively compute the closed-form posterior mean and variance at each level.
#'
#' From the fitted model from \code{\link{DNAmf}},
#' the posterior mean and variance are calculated based on the closed-form expression derived by a recursive fashion.
#' The formulas depend on its kernel choices.
#'
#' If the fitted model was constructed with non-nested designs (\code{nested=FALSE}),
#' the function generates \code{nimpute} sets of imputations for pseudo outputs
#' via \code{imputer}.
#'
#' For further details, see Heo, Boutelet, and Sung (2025+, <arXiv:2506.08328>).
#'
#' @param object A fitted DNAmf object.
#' @param x A vector or matrix of new input locations to predict.
#' @param targett A numeric value of target tuning parameter to predict.
#' @param nimpute Number of imputations for non-nested designs. Default is 50.
#' @param ... Additional arguments for compatibility with generic method \code{predict}.
#'
#' @return A list of predictive posterior mean and variance for each level and computation time containing:
#' \itemize{
#'   \item \code{mu_1}, \code{sig2_1}, ..., \code{mu_L}, \code{sig2_L}: A vector of predictive posterior mean and variance at each level.
#'   \item \code{mu}: A vector of predictive posterior mean at target tuning parameter.
#'   \item \code{sig2}: A vector of predictive posterior variance at target tuning parameter \code{targett}.
#'   \item \code{time}: Total computation time in seconds.
#' }
#' @importFrom stats setNames
#'
#' @rdname predict
#' @method predict DNAmf
#' @export
#' @examples
#' ### Non-Additive example ###
#' library(RNAmf)
#'
#' ### Non-Additive Function ###
#' fl <- function(x, t){
#'   term1 <- sin(10 * pi * x / (5+t))
#'   term2 <- 0.2 * sin(8 * pi * x)
#'   term1 + term2
#' }
#'
#' ### training data ###
#' n1 <- 13; n2 <- 10; n3 <- 7; n4 <- 4; n5 <- 1;
#' m1 <- 2.5; m2 <- 2.0; m3 <- 1.5; m4 <- 1.0; m5 <- 0.5;
#' d <- 1
#' eps <- sqrt(.Machine$double.eps)
#' x <- seq(0,1,0.01)
#'
#' ### fix seed to reproduce the result ###
#' set.seed(1)
#'
#' ### generate initial nested design ###
#' NestDesign <- NestedX(c(n1,n2,n3,n4,n5),d)
#'
#' X1 <- NestDesign[[1]]
#' X2 <- NestDesign[[2]]
#' X3 <- NestDesign[[3]]
#' X4 <- NestDesign[[4]]
#' X5 <- NestDesign[[5]]
#'
#' y1 <- fl(X1, t=m1)
#' y2 <- fl(X2, t=m2)
#' y3 <- fl(X3, t=m3)
#' y4 <- fl(X4, t=m4)
#' y5 <- fl(X5, t=m5)
#'
#' ### fit a DNAmf ###
#' fit.DNAmf <- DNAmf(X=list(X1, X2, X3, X4, X5), y=list(y1, y2, y3, y4, y5), kernel="sqex",
#'                    t=c(m1,m2,m3,m4,m5), multi.start=10, constant=TRUE)
#'
#' ### predict ###
#' pred.DNAmf <- predict(fit.DNAmf, x, targett=0)
#' predydiffu <- pred.DNAmf$mu
#' predsig2diffu <- pred.DNAmf$sig2
#'
#' ### RMSE ###
#' print(sqrt(mean((predydiffu-fl(x, t=0))^2))) # 0.1162579
#'
#' ### visualize the emulation performance ###
#' oldpar <- par(mfrow = c(2,3))
#' create_plot_base <- function(i, mesh_size, x, pred_mu, pred_sig2,
#'                              X_points = NULL, y_points = NULL, add_points = TRUE, yylim) {
#'   lower <- pred_mu - qnorm(0.995) * sqrt(pred_sig2)
#'   upper <- pred_mu + qnorm(0.995) * sqrt(pred_sig2)
#'
#'   plot(x, pred_mu, type = "n", ylim = c(-yylim, yylim), xlab = "", ylab = "",
#'        main = paste0("Mesh size = ", mesh_size), axes = FALSE)
#'   box()
#'
#'   polygon(c(x, rev(x)), c(upper, rev(lower)),
#'           col = adjustcolor("blue", alpha.f = 0.2), border = NA)
#'   lines(x, pred_mu, col = "blue", lwd = 2)
#'   lines(x, fl(x, mesh_size), lty = 2, col = "black", lwd = 2)
#'
#'   if (add_points && !is.null(X_points) && !is.null(y_points)) {
#'     points(X_points, y_points, col = "red", pch = 16, cex = 1.3)
#'   }
#' }
#'
#' mesh_sizes <- c(m1, m2, m3, m4, m5, 0)
#' mu_list <- list(pred.DNAmf$mu_1, pred.DNAmf$mu_2, pred.DNAmf$mu_3,
#'                 pred.DNAmf$mu_4, pred.DNAmf$mu_5, pred.DNAmf$mu)
#' sig2_list <- list(pred.DNAmf$sig2_1, pred.DNAmf$sig2_2, pred.DNAmf$sig2_3,
#'                   pred.DNAmf$sig2_4, pred.DNAmf$sig2_5, pred.DNAmf$sig2)
#' X_list <- list(X1, X2, X3, X4, X5, NULL)
#' y_list <- list(y1, y2, y3, y4, y5, NULL)
#'
#' plots <- mapply(function(i, m, mu, sig2, X, y) {
#'   create_plot_base(i, m, x, mu, sig2, X, y, add_points = !is.null(X), yylim=1.5)
#' }, i = 1:6, m = mesh_sizes, mu = mu_list, sig2 = sig2_list,
#' X = X_list, y = y_list, SIMPLIFY = FALSE)
#' par(oldpar)
#'

predict.DNAmf <- function(object, x, targett=0, nimpute=50, ...) {
  t1 <- proc.time()

  ### check the object ###
  if (!inherits(object, "DNAmf")) {
    stop("The object is not of class \"DNAmf\" \n")
  }

  ### prediction ###
  nn <- object$nn
  tt <- object$t
  nlevel  <- object$nlevel
  kernel  <- object$kernel

  fit1 <- object$fit1
  fit2 <- object$fit2

  if(object$nested){
    pred_result <- closed_form(fit1, fit2, targett, kernel=kernel, nn, tt, nlevel, x)
  }else{
    XX <- object$XX
    yy <- object$yy
    pred1 <- object$pred1

    sum_mu <- sum_mu2_plus_sig2 <- sig2_star <- vector("list", nlevel+1)
    names(sum_mu) <- mu_names <- c(paste0("mu_", 1:nlevel), "mu")
    names(sum_mu2_plus_sig2) <- names(sig2_star) <- sig2_names <- c(paste0("sig2_", 1:nlevel), "sig2")

    n <- nrow(matrix(x, ncol = ncol(fit1$X)))
    sum_mu <- setNames(rep(list(rep(0, n)), nlevel + 1), mu_names)
    sum_mu2_plus_sig2 <- setNames(rep(list(rep(0, n)), nlevel + 1), sig2_names)

    for (m in 1:nimpute) {
      # Generate imputed dataset for the m-th imputation
      yy <- imputer(XX, yy, kernel=kernel, tt, pred1, fit2)
      pred <- closed_form(fit1, fit2, targett, kernel=kernel, nn, tt, nlevel, x, XX, yy)

      sum_mu <- mapply(`+`, sum_mu, pred[mu_names], SIMPLIFY = FALSE)
      sum_mu2_plus_sig2 <- mapply(function(old, mu_vec, sig2_vec) old + mu_vec^2 + sig2_vec,
                                  sum_mu2_plus_sig2, pred[mu_names], pred[sig2_names], SIMPLIFY = FALSE)
    }

    # Compute final mu* and sig2* for each component
    mu_star <- lapply(sum_mu, function(s) s / nimpute)
    sig2_star <- mapply(function(mu2_plus_sig2, mu_star_vec) (mu2_plus_sig2 / nimpute) - (mu_star_vec)^2,
                        sum_mu2_plus_sig2, mu_star, SIMPLIFY = FALSE, USE.NAMES = TRUE)

    # Combine results into pred_result
    pred_result <- c(mu_star, sig2_star)
  }

  pred_result[["time"]] <- (proc.time() - t1)[3]

  return(pred_result)
}

