#' @title Internal fitting function for the Diffusion Non-Additive model
#'
#' @description Internal function that performs model fitting for the DNA model.
#' This function handles both nested and non-nested design using distinct fitting procedures.
#' In the nested case, parameters are estimated by maximum likelihood.
#' In the non-nested case, the stochastic EM algorithm is applied using an \code{\link{imputer}} function.
#'
#' @seealso \code{\link{DNAmf}} for the user-level function.
#'
#' @details This function is not intended to be called directly by users.
#' Instead, use \code{\link{DNAmf}} for model fitting.
#'
#' For further details, see Heo, Boutelet, and Sung (2025+, <arXiv:2506.08328>).
#'
#' @param X1 A matrix of input locations for the lowest fidelity level.
#' @param y1 A vector or matrix of response values at the lowest fidelity level.
#' @param X_list A list of matrices containing input locations for all higher fidelity levels \eqn{2,\ldots,L} combined.
#' @param y_list A list of matrices of responses for all higher fidelity levels \eqn{2,\ldots,L} combined.
#' @param kernel A character specifying the kernel type to be used. Choices are \code{"sqex"}(squared exponential), \code{"matern1.5"}, or \code{"matern2.5"}. Default is \code{"sqex"}.
#' @param t A vector of tuning parameters for each fidelity level.
#' @param nn A vector specifying the number of design points at each fidelity level.
#' @param nested A logical indicating whether the design is nested. Default is \code{TRUE}.
#' @param constant A logical indicating for constant mean of GP (\code{constant=TRUE}) or zero mean (\code{constant=FALSE}). Default is \code{TRUE}.
#' @param fitGP1 A logical indicating to fit GP at the lowest fidelity level \eqn{f_1} if \code{TRUE}, or not if \code{FALSE}. Default is \code{TRUE}.
#' @param init Optional vector of initial parameter values for optimization. Default is \code{NULL}.
#' @param n.iter Number of iterations for the stochastic EM algorithm in non-nested designs. Default is 50.
#' @param multi.start Number of random starting points for optimization. Default is 10.
#' @param trace A logical indicating to print progress of iterations if \code{TRUE}, or not if \code{FALSE}. Default is \code{TRUE}.
#' @param g Nugget term for numerical stability.
#' @param burn.ratio Fraction of iterations to discard as burn-in.
#' @param ... Additional arguments for compatibility with internal function \code{GP} or \code{matGP}.
#'
#' @return A fitted model object of class \code{DNAmf}:
#' \itemize{
#'   \item \code{fit1}: A fitted model \eqn{f_1} from \eqn{(X_1, y_1)}.
#'   \item \code{fit2}: A fitted model \eqn{f} from \eqn{(t_{-1}, X_{-1}, y_{-1})}.
#'   \item \code{kernel}: A copy of \code{kernel}.
#'   \item \code{constant}: A copy of \code{constant}.
#'   \item \code{t}: A copy of \code{t}.
#'   \item \code{TT}: A vector of tuning parameters corresponding to each input locations.
#'   \item \code{nn}: A copy of \code{nn}.
#'   \item \code{nlevel}: The number of fidelity levels \eqn{L}.
#'   \item \code{nested}: A copy of \code{nested}.
#'   \item \code{time}: Total computation time in seconds.
#' }
#' @usage DNAmf_internal(X1, y1, X_list, y_list, kernel, t, nn, nested = TRUE,
#' constant = TRUE, fitGP1 = TRUE, init=NULL, multi.start = 10, n.iter = 50,
#' trace = TRUE, g = sqrt(.Machine$double.eps), burn.ratio = 0.75, ...)
#' @noRd
#'

DNAmf_internal <- function(X1, y1, X_list, y_list, kernel, t, nn, nested = TRUE, constant = TRUE, fitGP1 = TRUE, init=NULL,
                           n.iter = 50, multi.start = 10, trace = TRUE, g = g, burn.ratio = burn.ratio, ...) {
  time0 <- proc.time()[3]

  if (nested) { # Nested
    X <- X_list
    y <- matrix(unlist(y_list[-1]))

    if (!checknested(X1, X)) stop("X is not nested by X1")

    tt <- rep(t[-1], times = nn[-1])

    if(kernel != "sqex") nu <- ifelse(kernel == "matern1.5", 1.5, 2.5)

    if(fitGP1){
      if (kernel == "sqex") {
        fit1 <- GP(X1, y1, constant = constant, ...)
      } else if (kernel %in% c("matern1.5", "matern2.5")) {
        fit1 <- matGP(X1, y1, nu = nu, constant = constant, ...)
      } else {
        stop("Unknown kernel")
      }
    }
    if (2 < length(nn)){
      nnn <- c(1, nn[-1])
      yy <- y1[-(1:(nn[1]-nn[2])),, drop=FALSE]
      for(i in seq_len(length(t)-2)){
        yy <- rbind(yy, y[sum(nnn[1:i]):(sum(nnn[1:(i+1)])-1),, drop=FALSE][-(1:(nnn[i+1]-nnn[i+2])),, drop=FALSE])
      }
      if(kernel=="sqex"){
        fit2 <- GP.nonsep.sqex(cbind(X, yy), y, tt, constant = constant, multi.start=multi.start, init=init, ...)
      }else{
        fit2 <- GP.nonsep.matern(cbind(X, yy), y, nu=nu, tt, constant = constant, multi.start=multi.start, init=init, ...)
      }
    } else {
      stop("Level should be larger than 2")
    }

    model <- list()
    if(fitGP1) model$fit1 <- fit1
    model$fit2 <- fit2
    model$kernel <- kernel
    model$constant <- constant
    model$t <- t
    model$nn <- nn
    model$nlevel <- length(t)
    model$nested <- TRUE
    class(model) <- "DNAmf"

  } else { # Non-nested
    XX    <- X_list
    y_all <- y_list
    L     <- length(XX$X_star)
    d     <- ncol(XX$X_star[[1]])

    t_list <- lapply(seq_len(L), function(i) rep(t[i], nrow(XX$X_list[[i]])))
    yy <- list(y_star  = vector("list", L), y_list  = y_all, y_tilde = vector("list", L))

    n.burnin <- ceiling(n.iter * burn.ratio) # number of burn in
    n.param <- n.iter - n.burnin # number of collected estimates
    param_mat <- matrix(NA, nrow=n.param, ncol=d+6) # theta_x, theta_y, theta_t, beta, delta, mu.hat, tau.hat

    if(kernel != "sqex") nu <- ifelse(kernel == "matern1.5", 1.5, 2.5)

    ### Draw initial y tilde ###
    ### sample initial y1_tilde from the imputer ###
    if(kernel=="sqex"){
      fit1 <- GP(XX$X_list[[1]], y_all[[1]], g = g, constant = TRUE, init = NULL)
      pred1 <- pred.GP(fit1, XX$X_tilde[[1]], cov.out = TRUE)
    }else{
      fit1 <- matGP(XX$X_list[[1]], y_all[[1]], nu=nu, g = g, constant = TRUE, init = NULL)
      pred1 <- pred.matGP(fit1, XX$X_tilde[[1]], cov.out = TRUE)
    }
    yy$y_tilde[[1]] <- pred1$mu
    yy$y_star[[1]] <- rbind(yy$y_list[[1]],yy$y_tilde[[1]])

    ### sample initial y_2_tilde using individual GPs ###
    if(kernel=="sqex"){
      fit2 <- GP.nonsep.sqex(X = do.call(rbind,XX$X_list[-1]), y = do.call(rbind,yy$y_list[-1]), t = do.call(c,t_list[-1]), constant = TRUE)
      for (l in 2:(L-1)) {
        pred2 <- pred.GP.nonsep(fit2, XX$X_tilde[[l]], rep(t[l], nrow(XX$X_tilde[[l]])))
        yy$y_tilde[[l]] <- pred2$mu
        yy$y_star[[l]] <- rbind(yy$y_list[[l]],yy$y_tilde[[l]])
      }
    }else{
      fit2 <- GP.nonsep.matern(X = do.call(rbind,XX$X_list[-1]), y = do.call(rbind,yy$y_list[-1]), nu=nu, t = do.call(c,t_list[-1]), constant = TRUE)
      for (l in 2:(L-1)) {
        pred2 <- pred.matGP.nonsep(fit2, XX$X_tilde[[l]], rep(t[l], nrow(XX$X_tilde[[l]])))
        yy$y_tilde[[l]] <- pred2$mu
        yy$y_star[[l]] <- rbind(yy$y_list[[l]],yy$y_tilde[[l]])
      }
    }
    yy$y_star[[L]] <- yy$y_list[[L]]

    ### initial estimates
    fit.DNAmf <- DNAmf(X=XX$X_star, y=yy$y_star, kernel=kernel, t=t, fitGP1=FALSE,
                       constant=constant, init = NULL, multi.start=10, ...)
    fit2 <- fit.DNAmf$fit2
    message(c(fit2$theta_x, fit2$theta_t, fit2$beta, fit2$delta, fit2$mu.hat, fit2$tau2hat))

    for (j in 1:n.iter) { # Imputation and Maximization
      if (trace) cat(j, '\n')

      # Imputation step; impute y tilde using ESS
      yy <- imputer(XX, yy, kernel=kernel, t, pred1, fit2)

      param.init <- c(fit2$theta_x, fit2$theta_t, fit2$beta, fit2$delta)
      # Maximization step; optimize parameters n.iter times
      fit.DNAmf <- DNAmf(X=XX$X_star, y=yy$y_star, kernel=kernel, t=t, fitGP1=FALSE,
                         constant=constant, init = param.init, ...)
      fit2 <- fit.DNAmf$fit2

      if(j > n.burnin){
        param_mat[j-n.burnin,] <- c(fit2$theta_x, fit2$theta_t, fit2$beta, fit2$delta, fit2$mu.hat, fit2$tau2hat)
      }
      if (trace) message(c(fit2$theta_x, fit2$theta_t, fit2$beta, fit2$delta, fit2$mu.hat, fit2$tau2hat))
    } # end of j for loop

    # average with 75% burn-in
    colnames(param_mat) <- c(paste0("theta_x", seq_len(d)), "theta_y", "theta_t", "beta", "delta", "mu_hat", "tau2_hat")
    final_params <- colMeans(param_mat)

    model <- fit.DNAmf
    model$fit1 <- fit1
    model$fit2 <- fit2
    model$kernel <- kernel
    model$fit2$theta_x <- final_params[1:(d+1)]
    model$fit2$theta_t <- final_params[d+2]
    model$fit2$beta <- final_params[d+3]
    model$fit2$delta <- final_params[d+4]
    model$fit2$mu.hat <- final_params[d+5]
    model$fit2$tau2hat <- final_params[d+6]

    model$XX <- XX
    model$yy <- yy
    model$t <- t
    model$pred1 <- pred1
    model$estim <- final_params
    model$nested <- FALSE
  }
  model$time <- proc.time()[3] - time0
  return(model)
}
