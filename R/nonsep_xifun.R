#' computing xi component in the closed form posterior mean of nonseparable matern kernel.
#'
#' @param c numerical value of \eqn{c_i(t)}.
#' @param w numerical value of \eqn{(Y_{-L})_i)}.
#' @param m numerical value of \eqn{\mu^*_{l-1}(x)}.
#' @param s numerical value of \eqn{\sigma^{*2}_{l-1}(x)}.
#' @param nu numerical value of smoothness hyperparameter. It should be 1.5 or 2.5.
#' @param theta numerical value of lengthscale hyperparameter \eqn{\theta_y}.
#'
#' @return calculated value of xi component.
#'
#' @importFrom stats pnorm
#' @noRd
#'

xifun <- function(c, w, m, s, theta, nu) {
  L    <- length(w)
  sig    <- sqrt(s)
  t2     <- theta^2
  
  if (nu == 1.5) {
    c3 <- sqrt(3)
    mua <- m - c3 * c * s / theta
    mub <- m + c3 * c * s / theta
    
    e1_1 <- (theta - c3 * c * w) / theta
    e1_2 <- c3 * c / theta
    e2_1 <- (theta + c3 * c * w) / theta
    
    lam11_2 <- mua
    lam21_2 <- -mub
    
    term1 <- exp((3*c^2*s + 2*c3*c*theta*(w - m)) / (2*t2)) *
      ( (e1_1 + e1_2 * lam11_2) * pnorm((mua - w)/sig) +
          e1_2 * sig / sqrt(2 * pi) * exp(-(w - mua)^2/(2*s)) )
    term2 <- exp((3*c^2*s - 2*c3*c*theta*(w - m)) / (2*t2)) *
      ( (e2_1 + e1_2 * lam21_2) * pnorm((w - mub)/sig) +
          e1_2 * sig / sqrt(2 * pi) * exp(-(w - mub)^2/(2*s)) )
  } else if (nu == 2.5) {
    c5 <- sqrt(5)
    mua <- m - c5 * c * s / theta
    mub <- m + c5 * c * s / theta
    
    e1_1 <- 1 - c5*c*w/theta + 5*c^2*w^2/(theta^2*3)
    e1_2 <- c5*c/theta - 10*c^2*w/(theta^2*3)
    e1_3 <- 5*c^2/(theta^2*3)
    e2_1 <- 1 + c5*c*w/theta + 5*c^2*w^2/(theta^2*3)
    e2_2 <- c5*c/theta + 10*c^2*w/(theta^2*3)
    
    lam11_2 <- mua
    lam11_3 <- mua^2 + s
    lam12_3 <- mua + w
    lam21_2 <- -mub
    lam21_3 <- mub^2 + s
    lam22_3 <- -(mub + w)
    
    term1 <- exp((5*c^2*s + 2*c5*c*theta*(w - m)) / (2*t2)) *
      ( (e1_1 + e1_2 * lam11_2 + e1_3 * lam11_3) * pnorm((mua - w)/sig) +
          (e1_2 + e1_3 * lam12_3) * sig / sqrt(2 * pi) * exp(-(w - mua)^2/(2*s)) )
    term2 <- exp((5*c^2*s - 2*c5*c*theta*(w - m)) / (2*t2)) *
      ( (e2_1 + e2_2 * lam21_2 + e1_3 * lam21_3) * pnorm((w - mub)/sig) +
          (e2_2 + e1_3 * lam22_3) * sig / sqrt(2 * pi) * exp(-(w - mub)^2/(2*s)) )
  }
  return(term1 + term2)
}