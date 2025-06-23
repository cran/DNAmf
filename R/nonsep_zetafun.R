#' computing zeta component in the closed form posterior variance of nonseparable matern kernel.
#'
#' @param c1 numerical value of \eqn{c_i(t_l)}.
#' @param c2 numerical value of \eqn{c_k(t_l)}.
#' @param w1 numerical value of \eqn{(Y_{-L})_i)}.
#' @param w2 numerical value of \eqn{(Y_{-L})_k)}.
#' @param m numerical value of \eqn{\mu^*_{l-1}(x)}.
#' @param s numerical value of \eqn{\sigma^{*2}_{l-1}(x)}.
#' @param nu numerical value of smoothness hyperparameter. It should be 1.5 or 2.5.
#' @param theta numerical value of lengthscale hyperparameter \eqn{\theta_y}.
#'
#' @return calculated value of zeta component.
#'
#' @importFrom stats pnorm
#' @noRd
#'

zetafun <- function(c1, c2, w1, w2, m, s, theta, nu) {
  L    <- length(w1)
  sig    <- sqrt(s)
  t2     <- theta^2
  
  wsmall <- pmin(w1, w2)
  wlarge <- pmax(w1, w2)
  csmall <- ifelse(w1 <= w2, c1, c2)
  clarge <- ifelse(w1 <= w2, c2, c1)
  
  if (nu == 1.5) {
    c3    <- sqrt(3)
    muc <- m - c3 * (csmall + clarge) * s / theta
    mud <- m + c3 * (csmall + clarge) * s / theta
    
    lam31_2 <- muc; lam31_3 <- muc^2+s; lam32_3 <- muc + wlarge
    lam41_2 <- m; lam41_3 <- m^2+s; lam42_3 <- m + wsmall; lam43_3 <- m + wlarge
    lam51_2 <- -mud; lam51_3 <- mud^2+s; lam52_3 <- -(mud + wsmall)
    
    e3_1 <- 1 + (3*csmall*clarge*wsmall*wlarge - c3*theta*(csmall*wsmall+clarge*wlarge))/t2
    e3_2 <- (c3*(csmall+clarge)*theta - 3*csmall*clarge*(wsmall+wlarge))/t2
    e3_3 <- 3*csmall*clarge/t2
    
    e4_1 <- 1 + (-3*csmall*clarge*wsmall*wlarge + c3*theta*(clarge*wlarge-csmall*wsmall))/t2
    e4_2 <- (c3*(csmall-clarge)*theta + 3*csmall*clarge*(wsmall+wlarge))/t2
    e4_3 <- -3*csmall*clarge/t2
    
    e5_1 <- 1 + (3*csmall*clarge*wsmall*wlarge + c3*theta*(csmall*wsmall+clarge*wlarge))/t2
    e5_2 <- (c3*(csmall+clarge)*theta + 3*csmall*clarge*(wsmall+wlarge))/t2
    e5_3 <- 3*csmall*clarge/t2
    
    term1 <- exp((3*s*(csmall+clarge)^2 + 2*c3*theta*(csmall*(wsmall - m)+clarge*(wlarge - m))) / (2*t2)) *
      ( (e3_1 + e3_2*lam31_2 + e3_3*lam31_3) * pnorm((muc-wlarge)/sig) +
          (e3_2 + e3_3*lam32_3) * (sig/sqrt(2*pi)) * exp(-(wlarge-muc)^2/(2*s)) )
    term2 <- exp((3 * s * (csmall - clarge)^2 + 2*c3*theta*(csmall*(wsmall - m)-clarge*(wlarge - m))) / (2*t2)) * 
      ( (e4_1 + e4_2*lam41_2 + e4_3*lam41_3) * (pnorm((wlarge-m)/sig) - pnorm((wsmall-m)/sig)) +
          (e4_2 + e4_3*lam42_3) * (sig/sqrt(2*pi)) * exp(-(wsmall-m)^2/(2*s)) -
          (e4_2 + e4_3*lam43_3) * (sig/sqrt(2*pi)) * exp(-(wlarge-m)^2/(2*s)) )
    term3 <- exp((3*s*(csmall+clarge)^2 - 2*c3*theta*(csmall*(wsmall - m)+clarge*(wlarge - m))) / (2*t2)) *
      ( (e5_1 + e5_2*lam51_2 + e5_3*lam51_3) * pnorm((wsmall-mud)/sig) +
          (e5_2 + e5_3*lam52_3) * (sig/sqrt(2*pi)) * exp(-(wsmall-mud)^2/(2*s)) )
  } else if (nu == 2.5) {
    c5  <- sqrt(5)
    muc <- m - c5 * (csmall + clarge) * s/theta
    mud <- m + c5 * (csmall + clarge) * s/theta
    
    lam31_2 <- muc; lam31_3 <- muc^2+s; lam31_4 <- muc^3+3*muc*s; lam31_5 <- muc^4+6*muc^2*s+3*s^2
    lam32_3 <- muc + wlarge; lam32_4 <- muc^2 + 2*s + wlarge^2 + muc*wlarge; lam32_5 <- wlarge^3 + muc^3 + muc^2*wlarge + muc*wlarge^2 + 3*s*wlarge + 5*muc*s
    lam41_2 <- m; lam41_3 <- m^2+s; lam41_4 <- m^3+3*m*s; lam41_5 <- m^4+6*m^2*s+3*s^2
    lam42_3 <- m + wsmall; lam42_4 <- m^2 + 2*s + wsmall^2 + m*wsmall; lam42_5 <- wsmall^3 + wsmall^2*m + wsmall*m^2 + m^3 + 3*s*wsmall + 5*s*m
    lam43_3 <- m + wlarge; lam43_4 <- m^2 + 2*s + wlarge^2 + m*wlarge; lam43_5 <- wlarge^3 + wlarge^2*m + wlarge*m^2 + m^3 + 3*s*wlarge + 5*s*m
    lam51_2 <- -mud; lam51_3 <- mud^2+s; lam51_4 <- -mud^3-3*mud*s; lam51_5 <- mud^4+6*mud^2*s+3*s^2
    lam52_3 <- -(mud + wsmall); lam52_4 <- mud^2 + 2*s + wsmall^2 + mud*wsmall; lam52_5 <- -mud^3 - wsmall^3 - mud^2*wsmall - mud*wsmall^2 - 3*s*wsmall - 5*mud*s
    
    e3_1 <- 1 + (25*csmall^2*clarge^2*wsmall^2*wlarge^2 - 3*c5*theta*(csmall*wsmall+clarge*wlarge)*(3*theta^2+5*csmall*clarge*wsmall*wlarge) + 15*theta^2*(csmall^2*wsmall^2+clarge^2*wlarge^2+3*csmall*clarge*wsmall*wlarge)) / (9*t2^2)
    e3_2 <- (9*c5*(csmall+clarge)*theta^3 + 15*c5*csmall*clarge*theta*(csmall*wsmall^2+clarge*wlarge^2) - 15*theta^2*(csmall*wsmall*(2*csmall+3*clarge)+clarge*wlarge*(2*clarge+3*csmall)) - 50*csmall^2*clarge^2*wsmall*wlarge*(wsmall+wlarge) + 30*c5*csmall*clarge*(csmall+clarge)*theta*wsmall*wlarge) / (9*t2^2)
    e3_3 <- 5*(5*csmall^2*clarge^2*(wsmall^2+wlarge^2+4*wsmall*wlarge)+3*(csmall^2+clarge^2+3*csmall*clarge)*theta^2 - 3*c5*theta*csmall*clarge*(2*csmall*wsmall+2*clarge*wlarge+clarge*wsmall+csmall*wlarge)) / (9*t2^2)
    e3_4 <- 5*(3*c5*csmall*clarge*(csmall+clarge)*theta - 10*csmall^2*clarge^2*(wsmall+wlarge)) / (9*t2^2)
    e3_5 <- 25*csmall^2*clarge^2/(9*t2^2)
    
    e4_1 <- 1 + (25*csmall^2*clarge^2*wsmall^2*wlarge^2 + 3*c5*theta*(clarge*wlarge-csmall*wsmall)*(3*theta^2-5*csmall*clarge*wsmall*wlarge) + 15*theta^2*(csmall^2*wsmall^2+clarge^2*wlarge^2-3*csmall*clarge*wsmall*wlarge)) / (9*t2^2)
    e4_2 <- (9*c5*(csmall-clarge)*theta^2 + 15*c5*csmall*clarge*theta*(clarge*wlarge^2-csmall*wsmall^2) - 15*theta^2*(csmall*wsmall*(2*csmall-3*clarge)+clarge*wlarge*(2*clarge-3*csmall)) - 50*csmall^2*clarge^2*wsmall*wlarge*(wsmall+wlarge) + 30*c5*csmall*clarge*(clarge-csmall)*theta*wsmall*wlarge) / (9*t2^2)
    e4_3 <- 5*(5*csmall^2*clarge^2*(wsmall^2+wlarge^2+4*wsmall*wlarge)+3*(csmall^2+clarge^2-3*csmall*clarge)*theta^2 - 3*c5*theta*csmall*clarge*(2*clarge*wlarge-2*csmall*wsmall+clarge*wsmall-csmall*wlarge)) / (9*t2^2)
    e4_4 <- 5*(3*c5*csmall*clarge*(clarge-csmall)*theta - 10*csmall^2*clarge^2*(wsmall+wlarge)) / (9*t2^2)
    e4_5 <- 25*csmall^2*clarge^2/(9*t2^2)
    
    e5_1 <- 1 + (25*csmall^2*clarge^2*wsmall^2*wlarge^2 + 3*c5*theta*(csmall*wsmall+clarge*wlarge)*(3*theta^2+5*csmall*clarge*wsmall*wlarge) + 15*theta^2*(csmall^2*wsmall^2+clarge^2*wlarge^2+3*csmall*clarge*wsmall*wlarge)) / (9*t2^2)
    e5_2 <- (9*c5*(csmall+clarge)*theta^3 + 15*c5*csmall*clarge*theta*(csmall*wsmall^2+clarge*wlarge^2) + 15*theta^2*(csmall*wsmall*(2*csmall+3*clarge)+clarge*wlarge*(2*clarge+3*csmall)) + 50*csmall^2*clarge^2*wsmall*wlarge*(wsmall+wlarge) + 30*c5*csmall*clarge*(csmall+clarge)*theta*wsmall*wlarge) / (9*t2^2)
    e5_3 <- 5*(5*csmall^2*clarge^2*(wsmall^2+wlarge^2+4*wsmall*wlarge)+3*(csmall^2+clarge^2+3*csmall*clarge)*theta^2 + 3*c5*theta*csmall*clarge*(2*clarge*wlarge+2*csmall*wsmall+clarge*wsmall+csmall*wlarge)) / (9*t2^2)
    e5_4 <- 5*(3*c5*csmall*clarge*(csmall+clarge)*theta + 10*csmall^2*clarge^2*(wsmall+wlarge)) / (9*t2^2)
    e5_5 <- 25*csmall^2*clarge^2/(9*t2^2)
    
    term1 <- exp((5*s*(csmall+clarge)^2 + 2*c5*theta*(csmall*(wsmall-m)+clarge*(wlarge-m))) / (2*t2)) *
      ( (e3_1 + e3_2*lam31_2 + e3_3*lam31_3 + e3_4*lam31_4 + e3_5*lam31_5) * pnorm((muc-wlarge)/sig) +
          (e3_2 + e3_3*lam32_3 + e3_4*lam32_4 + e3_5*lam32_5) * (sig/sqrt(2*pi)) * exp(-(wlarge-muc)^2/(2*s)) )
    term2 <- exp((5*s*(csmall-clarge)^2 + 2*c5*theta*(csmall*(wsmall-m)-clarge*(wlarge-m))) / (2*t2)) *
      ( (e4_1 + e4_2*lam41_2 + e4_3*lam41_3 + e4_4*lam41_4 + e4_5*lam41_5) * (pnorm((wlarge-m)/sig) - pnorm((wsmall-m)/sig)) +
          (e4_2 + e4_3*lam42_3 + e4_4*lam42_4 + e4_5*lam42_5) * (sig/sqrt(2*pi)) * exp(-(wsmall-m)^2/(2*s)) -
          (e4_2 + e4_3*lam43_3 + e4_4*lam43_4 + e4_5*lam43_5) * (sig/sqrt(2*pi)) * exp(-(wlarge-m)^2/(2*s)) )
    term3 <- exp((5*s*(csmall+clarge)^2 - 2*c5*theta*(csmall*(wsmall-m)+clarge*(wlarge-m))) / (2*t2)) *
      ( (e5_1 + e5_2*lam51_2 + e5_3*lam51_3 + e5_4*lam51_4 + e5_5*lam51_5) * pnorm((wsmall-mud)/sig) +
          (e5_2 + e5_3*lam52_3 + e5_4*lam52_4 + e5_5*lam52_5) * (sig/sqrt(2*pi)) * exp(-(wsmall-mud)^2/(2*s)) )
  }
  return(term1 + term2 + term3)
}