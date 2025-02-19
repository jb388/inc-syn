#' @author J. Beem-Miller
#' @description various helper functions for running soil R models
#' @details loads functions for calculating initial fraction modern value from k and vice-versa; supplies true half-life of 14C (lambda)
#' @param k turnover rate
#' @param fm fraction modern value
#' @return function 'fm' returns k, function 'k' returns fm; loads lambda into environment

#lambda (true half-life of 14C)
lambda <- 1 / 8267

# get initial fm value from k (pre-bomb)
fm <- function (k){
  k/(k + lambda)
}

# get k from initial fm (pre-bomb)
k <- function (fm) {
  if (fm > 1) {
    warning("fm must be pre-bomb")
  }
  (fm * lambda) / (1 - fm)
}