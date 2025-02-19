#' @author J. Beem-Miller
#' @description Calculates steady-state SOC stocks for 2-pool series type linear compartmental model system from supplied parameters
#' @param pars model parameters for linear compartmental model
#' @param In inputs for system given as single numeric value
#' @param out text string giving desired output for SOC stocks, i.e. by pools or for whole system; default is "pools", alternative is "sum" for total system stocks
#' @param mod_mat should model matrix (A_mat), input vector (in_vector), and parameter names be returned? T/F
#' @return list with components "ss_soc" and optionally "A_mat", "in_vector", "par_names"
soc.fx <- function(pars, In, out = "pools", mod_mat = FALSE) {
  
  # steady-state stock calc fx
  calc.soc <- function(A, in_vector) {
    (-1 * solve(A) %*% in_vector)
  }
  
  # define pool names
  pnms <- c("fast", "slow")
  
  # 2ps mod matrix 
  A <- -diag(pars[1:2])
  A[2, 1] <- pars[3] * pars[1]
  
  # 2ps steady-state C stocks
  in_vector <- c(In, 0)
  ss.cstock <- calc.soc(A, in_vector)
  
  if (out == "sum") {
    soc <- sum(ss.cstock)
    if (mod_mat) {
      list(A_mat = A, in_vector = in_vector, ss_soc = soc, par_names = pnms)
    } else {
      soc
    }
  } else {
    soc <- ss.cstock
    if (mod_mat) {
      list(A_mat = A, in_vector = in_vector, ss_soc = soc, par_names = pnms)
    } else {
      soc
    }
  }
}
