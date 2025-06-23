#' @title Check if designs are nested
#'
#' @description The function checks whether the design \code{XX2} (higher fidelity)
#' is nested within the design \code{XX1} (lower fidelity).
#'
#' @param XX1 A vector or matrix of input locations at lower fidelity.
#' @param XX2 A vector or matrix of input locations at higher fidelity.
#'
#' @return A logical indicating if XX2 is nested if \code{TRUE}, or not if \code{FALSE}.
#'
#' @noRd
#'

checknested <- function(XX1, XX2) {
  checknest <- c()
  for (i in 1:nrow(XX2)) {
    checknest <- c(checknest, suppressWarnings(any(apply(XX1, 1, function(xx) {
      all.equal(XX2[i, ], xx, tolerance = sqrt(.Machine$double.eps))
    }))))
  }
  checknest[is.na(checknest)] <- FALSE
  all(checknest)
}

#' @title Construct nested design sets
#'
#' @description The function constructs nested designs for multi-fidelity models as
#' \eqn{\mathcal{X}^*_L = \mathcal{X}_L},
#' \eqn{\mathcal{X}^*_l = \mathcal{X}_l \cup \mathcal{X}^*_{l+1}} for \eqn{l = 1, \dots, L-1},
#' and pseudo inputs \eqn{\widetilde{\mathcal{X}}_l := \mathcal{X}^*_l \setminus \mathcal{X}_l}.
#'
#' @param X_list A list of design sets for all fidelity levels.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{X_star}: A pseudo-complete nested inputs \eqn{\mathcal{X}^*_l}.
#'   \item \code{X_list}: An original inputs \eqn{\mathcal{X}_l}.
#'   \item \code{X_tilde}: A pseudo inputs \eqn{\widetilde{\mathcal{X}}_l}.
#' }
#'
#' @noRd
#'

makenested <- function(X_list){
  # Initialize the list to store nested matrices
  L <- length(X_list)
  X_star <- vector("list",L)

  # Set X^*_L = XL
  X_star[[L]] <- X_list[[L]]

  # Recursively build X^*_l for l = L down to 1
  for (l in (L-1):1) {
    X_star[[l]] <- unique(rbind(X_list[[l]], X_star[[l+1]])) # intersections are on X_list
  }

  # Compute X_tilde where X_tilde[[l]] = X_star[[l]] \ X_list[[l]]
  X_tilde <- vector("list", L)
  for (l in 1:(L-1)) {
    combined <- rbind(X_list[[l]], X_star[[l]])
    dup <- duplicated(combined)
    indices <- (nrow(X_list[[l]]) + 1) : nrow(combined)
    X_tilde[[l]] <- X_star[[l]][!dup[indices], , drop = FALSE]
  }

  return(list(X_star=X_star, X_list=X_list, X_tilde = X_tilde))
}

#' @title Find matching indices in nested designs
#'
#' @description The function identifies matching indices
#' between nested design based on a numerical distance tolerance.
#'
#' @param XX1 A matrix of original design locations.
#' @param XX2 A matrix of nested design locations.
#'
#' @return A vector of matching row indices.
#' @importFrom fields rdist
#'
#' @noRd
#'

checkindices <- function(XX1, XX2) {
  dist_matrix <- fields::rdist(XX2, XX1)
  matches <- which(dist_matrix < sqrt(.Machine$double.eps), arr.ind = TRUE)
  return(matches[, 2])
}
