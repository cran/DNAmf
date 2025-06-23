#' @title Imputation step in stochastic EM for the non-nested DNA Model
#'
#' @description The function performs the imputation step of the stochastic EM algorithm for the DNA model when the design is not nested.
#' The function generates pseudo outputs \eqn{\widetilde{\mathbf{y}}_l} at pseudo inputs \eqn{\widetilde{\mathcal{X}}_l}.
#'
#' @details For non-nested designs, pseudo-input locations \eqn{\widetilde{\mathcal{X}}_l}
#' are constructed using the internal \code{makenested} function.
#' The \code{imputer} function then imputes the corresponding pseudo outputs
#' \eqn{\widetilde{\mathbf{y}}_l = f_l(\widetilde{\mathcal{X}}_l)}
#' by drawing samples from the conditional normal distribution,
#' given fixed parameter estimates and previous-level outputs \eqn{Y_{-L}^{*(m-1)}},
#' at the \eqn{m}-th iteration of the EM algorithm.
#'
#' For further details, see Heo, Boutelet, and Sung (2025+, <arXiv:2506.08328>).
#'
#' @param XX A list of design sets for all fidelity levels, containing \code{X_star}, \code{X_list}, and \code{X_tilde}.
#' @param yy A list of current observed and pseudo-responses, containing \code{y_star}, \code{y_list}, and \code{y_tilde}.
#' @param kernel A character specifying the kernel type to be used. Choices are \code{"sqex"}(squared exponential), \code{"matern1.5"}, or \code{"matern2.5"}.
#' @param t A vector of tuning parameters for each fidelity level.
#' @param pred1 Predictive results for the lowest fidelity level \eqn{f_1}. It should include \code{cov} obtained by setting \code{cov.out=TRUE}.
#' @param fit2 A fitted model object for higher fidelity levels \eqn{f} from \eqn{(t_{-1}, X_{-1}, y_{-1})}.
#'
#' @return An updated \code{yy} list containing:
#' \itemize{
#'   \item \code{y_star}: An updated pseudo-complete outputs \eqn{\mathbf{y}^*_l}.
#'   \item \code{y_list}: An original outputs \eqn{\mathbf{y}_l}.
#'   \item \code{y_tilde}: A newly imputed pseudo outputs \eqn{\widetilde{\mathbf{y}}_l}.
#' }
#' @usage imputer(XX, yy, kernel=kernel, t, pred1, fit2)
#' @export
#'

imputer <- function(XX, yy, kernel=kernel, t, pred1, fit2){

  L <- length(XX$X_star)
  X <- fit2$X
  y <- fit2$y
  K <- fit2$K
  g <- fit2$g
  K <- K + diag(g, length(y))

  n_list <- sapply(yy$y_list, length)
  n_star <- sapply(yy$y_star, length)
  ncum_star <- c(0, cumsum(n_star[-1])) # cum sum of n_star
  n_tilde <- sapply(yy$y_tilde, length)
  ncum_tilde <- c(0,cumsum(n_tilde[2:(L-1)])) # cum sum of n_tilde

  X_m1 <- do.call(rbind, XX$X_star[-1])
  Y_m1 <- do.call(rbind, yy$y_star[-1])
  t_m1 <- unlist(lapply(2:L, function(l) rep(t[l], nrow(XX$X_star[[l]]))))
  params <- c(fit2$theta_x, fit2$theta_t, fit2$beta, fit2$delta)
  mu.hat <- fit2$mu.hat
  tau2hat <- fit2$tau2hat


  # Draw from prior distribution
  yy$y_tilde[[1]] <- t(mvtnorm::rmvnorm(1, mean = pred1$mu, sigma = pred1$cov))
  yy$y_star[[1]] <- rbind(yy$y_list[[1]],yy$y_tilde[[1]])

  ### sampling from Y_tilde give Y_list
  list_idx <- unlist(mapply(seq, ncum_star[1:(L-1)]+1, ncum_star[1:(L-1)]+n_list[2:L]))
  K_list <- K[list_idx, list_idx] # K corresponding to X_list
  Ki_list <- solve(K_list)
  K_tilde <- K[-list_idx, -list_idx] # K corresponding to X_tilde
  K_list_tilde <- K[list_idx, -list_idx] # K corresponding to X_list given X_tilde
  y_list <- y[list_idx]

  cond_mean <- mu.hat + t(K_list_tilde) %*% Ki_list %*% (y_list - mu.hat)
  cond_var <- tau2hat * (K_tilde - t(K_list_tilde) %*% Ki_list %*% K_list_tilde)
  cond_var <- (cond_var+t(cond_var))/2

  y_prior <- t(mvtnorm::rmvnorm(1, mean = cond_mean, sigma = cond_var))

  for(l in 2:(L-1)){
    yy$y_tilde[[l]] <- matrix(y_prior[(ncum_tilde[l-1]+1):ncum_tilde[l]],ncol=1)
    yy$y_star[[l]] <- rbind(yy$y_list[[l]], yy$y_tilde[[l]])
  }

  Y_m1 <- do.call(rbind, yy$y_star[-1])
  Y_mL <- do.call(rbind, lapply(2:L, function(l) {
    yy$y_star[[l-1]][checkindices(XX$X_star[[l-1]], XX$X_star[[l]]), , drop = FALSE]
  }))
  y <- Y_m1
  X[,ncol(X)] <- Y_mL
  if(kernel=="sqex"){
    K <- cor.sep.sqex(X, t=t_m1, param=params) + diag(g, length(y))
  }else if(kernel=="matern1.5"){
    K <- cor.sep.matern(X, t=t_m1, param=params, nu=1.5) + diag(g, length(y))
  }else if(kernel=="matern2.5"){
    K <- cor.sep.matern(X, t=t_m1, param=params, nu=2.5) + diag(g, length(y))
  }

  return(yy)
}
