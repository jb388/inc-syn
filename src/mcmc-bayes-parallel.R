# MCMC script for Bayesian parameter optimization
# aut: J. Beem-Miller
# date: 17-Feb-2025

# set date for saving files
date <- Sys.Date()

## Markov Chain Monte Carlo parameter optimization
# Note that the model uses prior variance distribution from 'modFit' optimization

# Model iterations
iter <- 5000 # start with 5000

# for saving script output
save.iter <- paste0(iter, "iter", ".RData")
save.dir <- file.path(paste0("~/inc-syn/dat/derived/bayes-par-fit-", date))
if(!dir.exists(save.dir)) dir.create(file.path(paste0("~/inc-syn/dat/derived/bayes-par-fit-", date)))


## Create list of MCMC fits (delayed rejection option, w/ adaptation & burn-in)

# Parallel implementation
#####
library(doParallel)
library(foreach)

# Detect number of cores on machine
UseCores <- detectCores()

# redefine bayes fit fx w/ parallel syntax
bayes.fit.prll.fx <- function(pars.fit, In.fit, obs.bulk.14c, obs.resp.14c, atm14C, upper, lower, ntrydr = 2,
                              var0, var0_name = "var_ms_unweighted", jump = NULL, burninlength = 0, 
                              niter, updatecov, nms) {
  # start timer
  start <- Sys.time()
  
  # define lambda (true half-life of 14C)
  lambda <- 1 / 8267
  
  # define fxs
  fm <- function (k){
    k/(k + lambda)
  }
  
  # modFun
  source('~/inc-syn/src/utilities/mod_fx.R')
  
  # define cost fx for current iteration 
  # only 14C costs for now...
  mod.Cost <- function(pars) {
    modelOutput <- modFun(pars, In = In.fit, pass = TRUE, atm14C = atm14C)
    cost1 <- FME::modCost(model = modelOutput, obs = obs.bulk.14c, scaleVar = TRUE)
    FME::modCost(model = modelOutput, obs = obs.resp.14c, scaleVar = TRUE, cost = cost1)
  }
  
  # run MCMC
  fit <- tryCatch(
    FME::modMCMC(f = mod.Cost,
            p = pars.fit, 
            var0 = var0[[var0_name]],
            jump = jump,
            upper = upper, 
            lower = lower,
            niter = niter,
            ntrydr = ntrydr,
            burninlength = burninlength,
            updatecov = updatecov),
    error = function (e) {cat("ERROR :", conditionMessage(e), "\n")})
  end <- Sys.time()
  cat(paste0(nms, " time: ", end - start, "\n"))
  return(fit)
}

# Function to calculate system age, pool ages, and transit time for all bayesian parameter combinations
sa.tt.prll.fx2 <- function(bayes_pars, pred_uncert, dlen, soc.meas, res) {
  
  # determine iterations
  iter <- nrow(bayes_pars[["pars"]])
  
  # sample index
  s_ind <- attributes(pred_uncert)$pset
  
  # initialize list
  ls.nms <- c("sysAge", "transT", "ins", "in_vs", "A_mat")
  SA.TT.ls <- lapply(ls.nms, function(ls) {
    ls <- vector(mode = "list", length = length(s_ind))
  })
  names(SA.TT.ls) <- ls.nms
  
  # load soc_fx
  source('~/inc-syn/src/utilities/soc_fx.R')
  
  # run loop
  for (j in seq_along(s_ind)) {
    
    # get pars
    PARS <- bayes_pars[["pars"]][s_ind[j], ]
    
    # run in fit fx
    SOC <- soc.meas
    IN <- 1 
    
    if (length(SOC) > 1) {
      SOC <- mean(SOC, na.rm = TRUE)
    }
    
    if (!is.null(SOC)) {
      if  (SOC < soc.fx(PARS, IN, "sum")) { 
        
        # by step; floor set at .001
        byStep <- (IN - .001) / res 
        
        # in vector
        ins <- seq(.001, IN, byStep)
        
      } else {
        
        # by step; ceiling set at SOC
        byStep <- (SOC - IN) / res 
        
        # in vector
        ins <- seq(IN, 
                   SOC, 
                   byStep)
      }
      
      # modeled stocks
      soc_mod <- lapply(seq_along(ins), function(j) soc.fx(PARS, ins[j], "sum"))
      IN.FIT <- round(ins[which.min(abs(unlist(soc_mod) - SOC))], 3)
      
    } else {
      IN.FIT <- 1
    }
    
    # run soc.fx
    soc.out <- soc.fx(pars = PARS, In = IN.FIT, mod_mat = TRUE)
    
    # calc ages & transit times
    # set index for distributions
    a <- seq(1, dlen)
    sa <- tryCatch(
      systemAge(A = soc.out$A_mat, u = soc.out$in_vector, a = a),
      error = function (e) {cat("ERROR :", conditionMessage(e), "\n")})
    tt <- tryCatch(
      transitTime(soc.out$A_mat, u = soc.out$in_vector),
      error = function (e) {cat("ERROR :", conditionMessage(e), "\n")})
    
    # Append to list
    SA.TT.ls[["sysAge"]][[j]] <- sa
    SA.TT.ls[["transT"]][[j]] <- tt
    SA.TT.ls[["ins"]][[j]] <- IN.FIT
    SA.TT.ls[["in_vs"]][[j]] <- soc.out$in_vector
    SA.TT.ls[["A_mat"]][[j]] <- soc.out$A_mat
  }
  return(SA.TT.ls)
}

## 2ps
# Register the cluster using n - 1 cores, send output to console (outfile = "")
outfile <- file.path(save.dir, "bayes_fit_2ps_ad_5000itr.txt")
cl <- makeCluster(UseCores - 1, outfile = outfile)

# start parallelization
registerDoParallel(cl)

# get req. objects
nms <- names(pars.i.2ps)

# Use foreach loop and %dopar% to compute in parallel
bayes_fit_2ps <- foreach(i = seq_along(pars.i.2ps), .packages = c("SoilR", "ISRaD", "FME")) %dopar%
  (bayes.fit.prll.fx(pars.fit = pars.i.2ps[[i]],
                     In.fit = in.1.ls[[i]],
                     obs.bulk.14c = obs.bulk14C.ls[[i]],
                     obs.resp.14c = obs.resp14C.ls[[i]],
                     atm14C = atm14C.ls[[i]],
                     upper = c(1, 1, 1),
                     lower = c(0, 0, 0),
                     # jump = , # default is par * 0.1
                     var0 = NULL, # best assuming var0 = 1 (i.e., NULL); alternative: modFit.2ps.ls[[i]]
                     niter = 10000,
                     burninlength = 500,
                     updatecov = 50, # algorithm did not even try w/ niter = updatecov [0% acceptance]
                     nms = nms[i]))
names(bayes_fit_2ps) <- names(pars.i.2ps)
# save output 
# *WARNING* will overwrite if file from current date exist!
save(bayes_fit_2ps, file = paste0(save.dir, "/bayes_fit_", "2ps_ad_", save.iter))

#stop the cluster
stopCluster(cl)

# calculate sensitivity
pred_uncert_2ps <- lapply(seq_along(bayes_fit_2ps), function(i) {
  if (!is.null(bayes_fit_2ps[[i]])) {
    cat(paste0("Estimating ", names(bayes_fit_2ps)[i], " sensitivity\n"))
    pars <- bayes_fit_2ps[[i]][["pars"]]
    sensRange(func = modFun, parInput = pars, In = 1, sensvar = c("bulkC", "resp"), sensR = TRUE, atm14C = atm14C.ls[[i]]) 
  }
})
names(pred_uncert_2ps) <- names(bayes_fit_2ps)
save(pred_uncert_2ps, 
     file = paste0(save.dir, "/pred_uncert_", "2ps_ad_", save.iter))

## SA and TT uncertainty
# Register the cluster using n - 1 cores, send output to console (outfile = "")
outfile <- file.path(save.dir, "sa.tt_2ps_10000itr.txt")
cl <- makeCluster(UseCores - 1, outfile = outfile)

# start parallelization
registerDoParallel(cl)

ix <- which(unlist(sapply(bayes_fit_2ps, "[[", "naccepted")) == 0)
# 2ps
SA.TT.2ps.ls <- foreach(i = seq_along(bayes_fit_2ps[-ix]), .packages = c("SoilR")) %dopar%
  (sa.tt.prll.fx2(
    bayes_pars = bayes_fit_2ps[-ix][[i]], 
    pred_uncert = pred_uncert_2ps[-ix][[i]], 
    dlen = 500,
    soc.meas = obs.soc.ls[-ix][[i]],
    res = 500))
names(SA.TT.2ps.ls) <- names(pars.i.2ps[-ix])
save(SA.TT.2ps.ls, file = paste0(save.dir, "/bayes_fit_SA_TT_", "2ps_", save.iter))

#stop the cluster
stopCluster(cl)
