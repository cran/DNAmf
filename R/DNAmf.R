#' @method summary DNAmf
#' @export
summary.DNAmf <- function(object,...){
  # Create a shallow copy to avoid modifying the original object in place
  ans <- object
  # Prepend the summary class (safer than overwriting all classes)
  class(ans) <- c("summary.DNAmf", class(object))
  ans
}

#' @export
print.summary.DNAmf <- function(x, ...){
  cat("----------------------------------------------------\n")
  cat(" Diffusion Non-Additive Multi-Fidelity (DNAmf) Model\n")
  cat("----------------------------------------------------\n")
  cat(sprintf("  * Levels:   %d\n", x$level))
  cat(sprintf("  * Kernel:   %s\n", x$kernel))
  cat(sprintf("  * Constant: %s\n", x$constant))
  cat(sprintf("  * nested:   %s\n", x$nested))

  # Calculate sample sizes for all levels
  if (x$nested) {
    # For nested, 'nn' stores the original sample sizes
    sizes <- x$nn
  } else {
    # For non-nested, 'XX$X_list' stores the original design matrices
    sizes <- vapply(x$XX$X_list, nrow, integer(1))
  }
  size_str <- paste(sprintf("Level %d: %d", seq_along(sizes), sizes), collapse = ", ")
  cat(sprintf("  * Samples:  %s\n", size_str))
  tune_str <- paste(sprintf("Level %d: %g", seq_along(x$t), x$t), collapse = ", ")
  cat(sprintf("  * Tuning Parameter:  %s\n", tune_str))

  cat(sprintf("  * Time:     %.3f sec\n", x$time))

  cat("----------------------------------------------------\n")
  invisible(x)
}

#' @method print DNAmf
#' @export
print.DNAmf <- function(x, ...) {
  print(summary(x, ...)) # Delegate printing to the summary method
  invisible(x)           # CRITICAL: Return the original object 'x'
}

#' @method coef DNAmf
#' @export
coef.DNAmf <- function(object, ...) {
  # Level 1 parameters (Standard GP)
  par_L1 <- list(
    theta = object$fit1$theta,
    mu_hat = object$fit1$mu.hat,
    tau2_hat = object$fit1$tau2hat
  )

  # Higher level parameters (Non-separable GP)
  # Structure differs slightly between nested/non-nested, but these are standard
  fit2 <- object$fit2
  par_High <- list(
    theta_x = fit2$theta_x, # or fit2$theta for nested sqex sometimes
    theta_t = fit2$theta_t,
    beta = fit2$beta,
    delta = fit2$delta,
    mu_hat = fit2$mu.hat,
    tau2_hat = fit2$tau2hat
  )

  # Fallback for nested object structure if naming differs
  if (is.null(par_High$theta_x) && !is.null(fit2$theta)) par_High$theta_x <- fit2$theta

  list(Level1 = par_L1, HigherLevels = par_High)
}

#' @method fitted DNAmf
#' @importFrom stats predict
#' @export
fitted.DNAmf <- function(object, ...) {
  # 1. Determine dimensions and setup
  d <- if(object$nested) ncol(object$fit1$X) else ncol(object$XX$X_list[[1]])
  L <- object$nlevel

  if (object$nested) {
    # --- NESTED CASE (Robust Superset Strategy) ---
    # 1. Predict on the full Level 1 design (Superset) to avoid "dims 0" crash
    X1 <- object$fit1$X
    preds_all <- predict(object, X1, ...)

    # 2. Extract the specific rows for each level
    X_rest_combined <- object$fit2$X[, 1:d, drop = FALSE]
    counts <- object$nn[-1]
    fac <- rep(seq_along(counts), times = counts)
    X_rest_list <- split(as.data.frame(X_rest_combined), fac)
    X_rest_list <- lapply(X_rest_list, as.matrix)

    X_list <- c(list(X1), X_rest_list)

    fitted_vals <- lapply(seq_len(L), function(l) {
      target_X <- X_list[[l]]

      # Robust row matching to find indices of this level's points inside X1
      id_X1 <- apply(X1, 1, paste, collapse = "\r")
      id_tgt <- apply(target_X, 1, paste, collapse = "\r")
      idx <- match(id_tgt, id_X1)

      # 3. Extract prediction for level l
      # Try "mu1", "mu2"... first. If not found, try "mu_1", "mu_2"...
      val <- preds_all[[paste0("mu", l)]]
      if (is.null(val)) val <- preds_all[[paste0("mu_", l)]]

      if (is.null(val)) stop(paste("Could not find prediction for level", l, "in predict object."))

      val[idx]
    })
  } else {
    # --- NON-NESTED CASE ---
    # We must predict on specific sets because design points differ.
    X_list <- object$XX$X_list

    fitted_vals <- lapply(seq_len(L), function(l) {
      X_curr <- X_list[[l]][, 1:d, drop = FALSE]

      # Pad single row to prevent crash
      single_row <- nrow(X_curr) == 1
      if (single_row) X_curr <- rbind(X_curr, X_curr)

      preds <- predict(object, X_curr)

      # Extract prediction
      val <- preds[[paste0("mu", l)]]
      if (is.null(val)) val <- preds[[paste0("mu_", l)]]

      if (single_row) return(val[1]) else return(val)
    })
  }
  return(fitted_vals)
}

#' @method residuals DNAmf
#' @importFrom stats fitted
#' @export
residuals.DNAmf <- function(object, ...) {
  # 1. Get Fitted values
  y_hat <- fitted(object)

  # 2. Get Observed values (original y)
  if (object$nested) {
    y1 <- object$fit1$y
    y_rest_combined <- object$fit2$y
    counts <- object$nn[-1]
    y_rest_list <- split(y_rest_combined, rep(seq_along(counts), times = counts))
    y_list <- c(list(y1), y_rest_list)
  } else {
    y_list <- object$yy$y_list
  }

  # 3. Calculate Residuals
  mapply(function(y, y_h) y - y_h, y_list, y_hat, SIMPLIFY = FALSE)
}

#' @method plot DNAmf
#' @importFrom graphics par lines points legend
#' @export
plot.DNAmf <- function(x, ...) {
  # 1. Determine dimensions and number of levels
  # DNAmf uses 'nlevel', not 'level'
  L <- x$nlevel

  if (x$nested) {
    d <- ncol(x$fit1$X)
  } else {
    d <- ncol(x$XX$X_list[[1]])
  }

  # 2. Stop if input is not 1-dimensional
  if (d > 1) {
    stop("The plot method for 'DNAmf' currently only supports 1-dimensional inputs.")
  }

  # 3. Setup plotting layout (1 row, L columns)
  oldpar <- par(mfrow = c(1, L))
  on.exit(par(oldpar))

  # 4. Prepare Data for Plotting (Original Training Data)
  # We need to build a list of X and y for each level to iterate easily
  if (x$nested) {
    # -- Nested Case --
    # Level 1 data
    X1 <- x$fit1$X
    y1 <- x$fit1$y

    # Higher levels data (stored combined in fit2)
    # We take only column 1 for X because we checked d=1
    X_rest <- x$fit2$X[, 1, drop=FALSE]
    y_rest <- x$fit2$y

    # Split them based on sample sizes stored in x$nn
    counts <- x$nn[-1]
    fac <- rep(seq_along(counts), times = counts)

    X_list <- c(list(X1), split(X_rest, fac))
    y_list <- c(list(y1), split(y_rest, fac))

    # Ensure they remain matrices/vectors after split
    X_list <- lapply(X_list, as.matrix)
    y_list <- lapply(y_list, as.matrix)

  } else {
    # -- Non-Nested Case --
    # Data is stored in XX and yy lists
    X_list <- x$XX$X_list
    y_list <- x$yy$y_list
  }

  # 5. Generate Predictions on Grid
  xg <- matrix(seq(0, 1, length.out = 101), ncol = 1)
  preds <- predict(x, xg)

  # 6. Plot each level
  for (l in 1:L) {
    # Extract training data from our prepared lists
    X_train <- X_list[[l]]
    y_train <- y_list[[l]]

    # Extract predictions for this level
    # DNAmf predict returns flat names "mu1", "mu2" etc.
    y_mu  <- preds[[paste0("mu", l)]]
    if(is.null(y_mu)) y_mu <- preds[[paste0("mu_", l)]] # Fallback name check

    # Extract Variance (same flat structure)
    y_sig2 <- preds[[paste0("sig2", l)]]
    if(is.null(y_sig2)) y_sig2 <- preds[[paste0("sig2_", l)]]
    y_sig <- sqrt(y_sig2)

    # Calculate 95% Confidence Interval
    lower <- y_mu - 1.96 * y_sig
    upper <- y_mu + 1.96 * y_sig

    # Determine plot limits
    ylim <- range(c(y_train, lower, upper), na.rm = TRUE)

    # Plot Mean Curve
    plot(xg, y_mu, type = "l", lwd = 2, col = "blue", ylim = ylim,
         main = paste("Level", l, "Fit"),
         xlab = "x", ylab = "y", ...)

    # Add Confidence Interval (Dashed lines)
    lines(xg, lower, col = "blue", lty = 2)
    lines(xg, upper, col = "blue", lty = 2)

    # Add Training Points
    points(X_train, y_train, pch = 16, cex = 1, col = "black")

    # Add Legend (only on the first plot)
    if (l == 1) {
      legend("topleft", legend = c("Mean", "95% CI", "Data"),
             col = c("blue", "blue", "black"),
             lty = c(1, 2, NA), pch = c(NA, NA, 16), bty = "n")
    }
  }

  invisible(x)
}

#' @title Fitting a Diffusion Non-Additive model for multi-fidelity computer experiments with tuning parameters
#'
#' @description The function fits DNA models for multi-fidelity computer experiments with tuning parameters.
#' Available kernel choices include nonseparable squared exponential kernel, and nonseparable Matern kernel with smoothness parameter 1.5 and 2.5.
#' The function returns a fitted model object of class \code{DNAmf}, produced by \code{DNAmf_internal}.
#'
#' @seealso \code{\link{predict.DNAmf}} for prediction.
#'
#' @details The \code{DNAmf} function internally calls \code{DNAmf_internal} to fit the DNA model with nonseparable kernel.
#'
#' The model structure is:
#' \eqn{\begin{cases}
#' & f_1(\bm{x}) = W_1(\bm{x}),\\
#' & f_l(\bm{x}) = W(t_l, \bm{x}, f_{l-1}(\bm{x})),
#' \end{cases}}
#' where \eqn{W(t, \bm{x}, y) \sim GP(\alpha, \tau^2 K((t, \bm{x}, y), (t', \bm{x}', y')))} is a GP model.
#' Hyperparameters \eqn{(\alpha, \tau^2, \bm{\theta})} are estimated by
#' maximizing the log-likelihood via an optimization algorithm "L-BFGS-B".
#' For \code{constant=FALSE}, \eqn{\alpha=0}.
#'
#' The nonseparable covariance kernel is defined as:
#' \deqn{K((t, \bm{x}, y), (t', \bm{x}', y'))=
#' \left(\frac{(t-t')^2}{\theta_t} + 1\right)^{ - \left(\frac{\beta(d+1)}{2}+\delta \right) }
#' \prod^d_{j=1}\phi(x_j,x'_j;\theta_{j})\phi(y,y';\theta_{y}),}
#' where \eqn{\phi(\cdot, \cdot)} depens on the chosen kernel:
#'
#' \itemize{
#' \item For nonseparable squared exponential kernel(\code{kernel = "sqex"}):
#' \deqn{\phi(x, x';\theta) = \exp \left( -\left(\frac{(t-t')^2}{\theta_t} + 1\right)^{-\beta}
#' \frac{ (x-x')^2 }{\theta} \right)}
#'
#' \item For nonseparable Matern kernel with smoothness parameter of \eqn{\nu=1.5} (\code{kernel = "matern1.5"}):
#' \deqn{\phi(x,x';\theta) = \left( 1+\frac{1}{\left(\frac{(t-t')^2}{\theta_t} + 1\right)^{ \frac{\beta}{2} }}\frac{\sqrt{3}|x- x'|}{\theta} \right)
#' \exp \left( -\frac{1}{\left(\frac{(t-t')^2}{\theta_t} + 1\right)^{ \frac{\beta}{2} }}\frac{\sqrt{3}|x- x'|}{\theta} \right)}
#'
#' \item For nonseparable Matern kernel with smoothness parameter of \eqn{\nu=2.5} (\code{kernel = "matern2.5"}):
#' \deqn{\phi(x, x';\theta) = \left( 1+\frac{1}{\left(\frac{(t-t')^2}{\theta_t} + 1\right)^{ \frac{\beta}{2} }}\frac{\sqrt{5}|x- x'|}{\theta}+
#' \frac{1}{3}\left(\frac{1}{\left(\frac{(t-t')^2}{\theta_t} + 1\right)^{ \frac{\beta}{2} }}\frac{\sqrt{5}|x- x'|}{\theta} \right)^2 \right) }
#' \deqn{\times \exp \left( -\frac{1}{\left(\frac{(t-t')^2}{\theta_t} + 1\right)^{ \frac{\beta}{2} }}\frac{\sqrt{5}|x- x'|}{\theta} \right)}
#' }
#'
#' When the input locations are not nested, the internal \code{makenested} function constructs nested designs as
#' \eqn{\mathcal{X}^*_L = \mathcal{X}_L} and
#' \eqn{\mathcal{X}^*_l = \mathcal{X}_l \cup \mathcal{X}^*_{l+1}} for \eqn{l = 1, \dots, L-1}.
#' The function \code{\link{imputer_DNA}} then imputes pseudo outputs \eqn{\widetilde{\mathbf{y}}_l := f_l(\widetilde{\mathcal{X}}_l)}
#' at pseudo inputs \eqn{\widetilde{\mathcal{X}}_l := \mathcal{X}^*_l \setminus \mathcal{X}_l},
#' using a stochastic EM algorithm.
#'
#' For further details, see Heo, Boutelet, and Sung (2025+, <arXiv:2506.08328>).
#'
#' @param X A list of input locations for all fidelity levels \eqn{1,\ldots,L} combined.
#' @param y A list of response values for all fidelity levels \eqn{1,\ldots,L} combined.
#' @param kernel A character specifying the kernel type to be used. Choices are \code{"sqex"}(squared exponential), \code{"matern1.5"}, or \code{"matern2.5"}. Default is \code{"sqex"}.
#' @param t A vector of tuning parameters for each fidelity level.
#' @param constant A logical indicating for constant mean of GP (\code{constant=TRUE}) or zero mean (\code{constant=FALSE}). Default is \code{TRUE}.
#' @param init Optional vector of initial parameter values for optimization. Default is \code{NULL}.
#' @param n.iter Number of iterations for the stochastic EM algorithm for non-nested designs. Default is 50.
#' @param multi.start Number of random starting points for optimization. Default is 10.
#' @param g Nugget term for numerical stability. Default is \code{sqrt(.Machine$double.eps)}.
#' @param burn.ratio Fraction of iterations to discard as burn-in. Default is 0.75.
#' @param ... Additional arguments for compatibility with \code{DNAmf_internal}.
#'
#' @return A fitted model object of class \code{DNAmf}.
#'
#' @usage DNAmf(X, y, kernel = "sqex", t, constant = TRUE, init=NULL,
#' n.iter=50, multi.start=10, g = sqrt(.Machine$double.eps), burn.ratio = 0.75, ...)
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

DNAmf <- function(X, y, kernel="sqex", t, constant = TRUE, init=NULL, n.iter=50, multi.start=10,
                  g = sqrt(.Machine$double.eps), burn.ratio = 0.75, ...) {
  if (!is.list(X)) stop("DNAmf: 'X' must be a list of input matrices.")
  if (!is.list(y)) stop("DNAmf: 'y' must be a list of response matrices")

  L <- length(X)
  if (length(y) != L) stop("DNAmf: Length of 'X' and 'y' must match.")
  if (L < 2) stop("DNAmf: Requires at least two fidelity levels (L >= 2).")
  if (missing(t)) stop("DNAmf: Argument 't' (tuning parameters) is missing.")
  if (!is.numeric(t)) stop("DNAmf: 't' must be a numeric vector.")
  if (length(t) != L) stop(paste0("DNAmf: Length of 't' (", length(t), ") must match the number of levels in 'X' (", L, ")."))
  if (!is.character(kernel) || length(kernel) != 1 || !kernel %in% c("sqex", "matern1.5", "matern2.5")) {
    stop(paste("DNAmf: 'kernel' must be one of:", paste(c("sqex", "matern1.5", "matern2.5"), collapse = ", ")))
  }
  if (!is.numeric(n.iter) || n.iter < 1) stop("DNAmf: 'n.iter' must be a positive integer.")
  if (!is.numeric(multi.start) || multi.start < 1) stop("DNAmf: 'multi.start' must be a positive integer.")
  if (!is.numeric(burn.ratio) || burn.ratio <= 0 || burn.ratio >= 1) stop("DNAmf: 'burn.ratio' must be between 0 and 1.")
  if (is.null(g)) g <- sqrt(.Machine$double.eps)

  d <- NCOL(X[[1]])
  for (i in 1:L) {
    if (is.vector(X[[i]])) X[[i]] <- matrix(X[[i]], ncol = 1)
    if (is.vector(y[[i]])) y[[i]] <- matrix(y[[i]], ncol = 1)
    if (!is.matrix(X[[i]])) stop(paste("DNAmf: Element", i, "of 'X' must be a matrix or numeric vector."))
    if (!is.matrix(y[[i]])) stop(paste("DNAmf: Element", i, "of 'y' must be a matrix or numeric vector."))
    if (ncol(X[[i]]) != d) stop(paste0("DNAmf: Dimension mismatch. Level 1 has ", d, " columns, but Level ", i, " has ", ncol(X[[i]]), "."))
    if (nrow(X[[i]]) != nrow(y[[i]])) stop(paste0("DNAmf: Sample size mismatch at Level ", i, "."))
    if (ncol(y[[i]]) != 1) stop(paste0("DNAmf: Response 'y' at Level ", i, " must be a 1-column matrix (univariate)."))
  }

  X1 <- X[[1]]; y1 <- y[[1]]
  nn <- unlist(lapply(X, nrow))
  X <- do.call(rbind, X[-1])
  y <- do.call(rbind, y[-1])
  lvl <- rep(seq_along(nn[-1]), times = nn[-1])
  idxs  <- split(seq_len(nrow(X)), lvl)
  X_list <- lapply(idxs, function(i) X[i, , drop = FALSE])
  y_list <- lapply(idxs, function(i) y[i, , drop = FALSE])

  # check whether the design is nested #
  nested <- all(unlist(lapply(X_list, checknested, XX1=X1)))

  # make designs nested
  if(!nested) XX <- makenested(c(list(X1), X_list)) else XX <- X
  model <- DNAmf_internal(X1, y1, XX, c(list(y1), y_list), kernel=kernel, t=t, nn=nn, nested = nested, constant = constant,
                          init=init, n.iter = n.iter, multi.start = multi.start, trace = TRUE, g = g, burn.ratio = 0.75, ...)
  return(model)
}
