## 2-pool series model initial parameter search script
# Requires the following functions: modFun; soc.fx; con.df.fx; C14.plot.fx
# Requires the following input: obs.bulk14C.ls, obs.resp14C.ls, obs.soc.ls
# These functions and inputs are in Rmd file "mod-strs.RMD" (this directory)

# Initialize lists (5 sites by 3 depths = 27 elements)
pars.i.2ps <- vector(mode = "list", length = length(obs.bulk14C.ls))
names(pars.i.2ps) <- names(obs.bulk14C.ls)

# Hoelstein control
##### 
## 0-5
siteDepth <- names(pars.i.2ps)[1]
k1 <- .2
k2 <- .03
propFast <- .01
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


## 10-20
siteDepth <- names(pars.i.2ps)[2]
k1 <- .2
k2 <- .001
propFast <- .01
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


## 5-10
siteDepth <- names(pars.i.2ps)[3]
k1 <- .2
k2 <- .001
propFast <- .01
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


## LI
siteDepth <- names(pars.i.2ps)[4]
k1 <- .5
k2 <- .1
propFast <- .05
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


#####

# Jaun
##### 
## 0-5
siteDepth <- names(pars.i.2ps)[5]
k1 <- .15
k2 <- .03
propFast <- .03
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


## 10-20
siteDepth <- names(pars.i.2ps)[6]
k1 <- .08
k2 <- .004
propFast <- .3
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


## 20-40
siteDepth <- names(pars.i.2ps)[7]
k1 <- .08
k2 <- .0009
propFast <- .2
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


## 40-90
siteDepth <- names(pars.i.2ps)[8]
k1 <- .05
k2 <- .00023
propFast <- .43
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


## 5-10
siteDepth <- names(pars.i.2ps)[9]
k1 <- .1
k2 <- .01
propFast <- .1
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars
#####

# Reckenholz
#####
## 0-5
siteDepth <- names(pars.i.2ps)[10]
k1 <- .12
k2 <- .0013
propFast <- .15
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


## 10-20
siteDepth <- names(pars.i.2ps)[11]
k1 <- .12
k2 <- .0011
propFast <- .15
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


## 20-40
siteDepth <- names(pars.i.2ps)[12]
k1 <- .1
k2 <- .0008
propFast <- .95
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


## 40-90
siteDepth <- names(pars.i.2ps)[13]
k1 <- .01
k2 <- .00035
propFast <- .95
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars


## 5-10
siteDepth <- names(pars.i.2ps)[14]
k1 <- .01
k2 <- .0011
propFast <- .8
pars <- c(k1, k2, propFast)

# plot
C14.plot.fx(
  plot.df = modFun(pars = pars, In = 1, atm14C = atm14C.ls[[siteDepth]], out = "", verbose = TRUE, siteDepth = siteDepth), 
  con.df = con.df.fx(siteDepth),
  siteDepth = siteDepth)

# save pars
pars.i.2ps[[siteDepth]] <- pars
#####


# save
save(pars.i.2ps, file = "/Users/jbeemmil@umich.edu/LIM_soilCmods/dat/derived/modFit_pars/pars.i.2ps.RData")
