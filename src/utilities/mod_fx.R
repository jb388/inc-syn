#' @author J. Beem-Miller
#' @description Runs SoilR models and returns a data frame for running FME::modFit or returns a data frame for plotting in ggplot2
#' @param pars model parameters for linear compartmental model
#' @param In inputs for system given as single numeric value
#' @param pass should SoilR fit model even if output is unreasonable, i.e. negative SOC stocks/respiration? Default is TRUE
#' @param verbose if out != "modFit", should function print SOC stocks and siteDepth?
#' @param out text string giving desired output: either "modFit" (default) for runnning FME::modFit, or a data frame for plotting in ggplot2 
#' @param siteDepth siteDepth string; used for printing SOC stocks during search for initial parameter values and retrieving site-specific atm14C record
#' @param var_14c name for 14C variable; default "d14c"
#' @param sensR flag for fitting FME::sensRange function, which requires lag = 0; default is FALSE
#' @return if out == "modFit", data frame with columns time, resp, bulkC, and cStock; else data frame with columns years, d14c, pool
#' @import SoilR
#' @import dplyr
modFun <- function(pars, In, atm14C, siteDepth, pass = TRUE, verbose = TRUE, out = "modFit", var_14c = "d14c", sensR = FALSE) {
  
  # soc_fx
  if (!exists("soc.fx")) {
    source('/Users/jbeemmil@umich.edu/LIM_soilCmods/src/utilities/soc_fx copy.R')
  }
  
  # helper fxs
  if (!exists("lambda")) {
    source('/Users/jbeemmil@umich.edu/LIM_soilCmods/src/utilities/mod_help_fxs copy.R')
  }
  
  # run soc.fx to get: mod_mat [[1]], in_vector [[2]], steady-state C [[3]], and pool names [[4]]
  soc.fx_out <- soc.fx(pars, In, mod_mat = TRUE)
  ss.cstock <- soc.fx_out[[3]]
  
  # check for negative stocks
  if (any(ss.cstock <= 0)) {
    cat("pool ", which(ss.cstock <= 0), "< 0\n")
  }
  
  # calculate initial 14C
  F0_Delta14C <- unlist(
    lapply(-diag(soc.fx_out[[1]]), function(x) Delta14C_from_AbsoluteFractionModern(fm(x))))
    
  # multipool model fx
  mod.fx <- function(A,
                     t,
                     in_vector,
                     C0,
                     F0_Delta14C, 
                     xi = 1, # timestep
                     inputFc, 
                     lag = 0,
                     pass = pass) {
    t_start = min(t)
    t_stop = max(t)
    inputFluxes = BoundInFluxes(function(t) {
      matrix(nrow = length(in_vector), ncol = 1, in_vector)
    }, t_start, t_stop)
    if (length(xi) == 1) 
      fX = function(t) {
        xi
      }
    At = BoundLinDecompOp(map = function(t) {
      fX(t) * A
    }, t_start, t_stop)
    Fc = BoundFc(inputFc, lag = lag, format = "Delta14C")
    mod = Model_14(t, At, ivList = C0, initialValF = ConstFc(F0_Delta14C, "Delta14C"), 
                   inputFluxes = inputFluxes, inputFc = Fc, pass = pass)
  }
    
  # run model
  model <- mod.fx(
    A = soc.fx_out[[1]],
    t = atm14C$Year.AD,
    in_vector = soc.fx_out[[2]],
    C0 = as.vector(ss.cstock),
    F0_Delta14C = F0_Delta14C,
    inputFc = atm14C,
    pass = pass) 
  
  # get mod values
  C14m <- getF14C(model)
  C14p <- getF14(model) 
  C14r <- getF14R(model)
  Ctot <- getC(model)
  
  if (var_14c == "fm") {
    dates <- atm14C$Year.AD
    for (i in seq_along(dates)) {
      C14m[i] <- convert_fm_d14c(d14c = C14m[i], obs_date_y = dates[i], verbose = FALSE) 
      C14p[i] <- convert_fm_d14c(d14c = C14p[i], obs_date_y = dates[i], verbose = FALSE) 
      C14r[i] <- convert_fm_d14c(d14c = C14r[i], obs_date_y = dates[i], verbose = FALSE)
    }
  }
  
  if (out == "modFit") {
    # dataframe for modFit fx
    data.frame(
      time = atm14C$Year.AD,
      resp = C14r,
      bulkC = C14m,
      cStock = sum(Ctot[1, ]))
  } else {
    
    # sum c stocks
    ss.cstock <- round(ss.cstock, 2)
    cstock.sum <- ifelse(is.null(dim(ss.cstock)), ss.cstock, colSums(ss.cstock))
    
    # print site and steady-state stocks
    if (verbose) {
      if (!is.null(siteDepth)) cat(paste0(siteDepth, "\n"))
      for (i in seq_along(ss.cstock)) {
        cat(paste(soc.fx_out[[4]][i], ss.cstock[i], "\n"))
      }
      cat(cstock.sum, " (modeled total C stock)\n")
      if (!is.null(siteDepth)) {
        cat(round(obs.soc.ls[[siteDepth]], 1), " (measured total C stock)\n")
        if (cstock.sum < obs.soc.ls[[siteDepth]]) cat("Modeled stocks too low w/ current inputs\n") else
          if (cstock.sum > obs.soc.ls[[siteDepth]]) cat("Modeled stocks too high w/ current inputs\n")
      }
    }
    
    
    # data frame for plotting
    data.frame(
      years = rep(atm14C$Year.AD, (ncol(C14p) + 3)),
      d14C = c(c(C14p),
               C14m,
               C14r,
               atm14C$Delta14C),
      pool = rep(c(soc.fx_out[[4]], "bulkC", "respiration", "atm"), 
                 each = nrow(C14p))) %>%
      distinct
  }
}
