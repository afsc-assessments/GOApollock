## LIBRARIES
library(GOApollock)
library(TMBhelper)
library(ggplot2)
library(TMB)
library(dplyr)

## Converting random walk fishery selectivity deviates to random effects and estimating sd
# - NOTE: 4*sd for inf1_fsh_dev was removed


## READ IN 2022 DATA ----
dat <- readRDS('TMB/data/datfile.RDS')
replist <- readRDS("TMB/data/repfile.RDS")
stdadmb <- readRDS("TMB/data/stdfile.RDS")
pars <- readRDS('TMB/data/pars.RDS')
map <- readRDS('TMB/data/map.RDS')


# MODEL 1 ----
# - Double logistic with random effects on ascending portion
# - Adjust fishery random walk standard deviation to parameter
pars$ln_sel_sd <- log(dat$rwlk_sd[1]) # All the same so reducing to length of 1
map$ln_sel_sd <- as.factor(1)

# - Non-parametric params
pars$selpars_re <- array(rep(0, dat$trmage), dim = c(dat$trmage, 1, 1))
map$selpars_re <- as.factor(rep(NA, dat$trmage))

pars$mean_sel <- 0
map$mean_sel <- as.factor(NA)

# - AR params
pars$sel_rho_y <- 0
pars$sel_rho <- 0

map$sel_rho <- as.factor(NA)
map$sel_rho_y <- as.factor(NA)


# - Ad fishery selectivity switch
dat$seltype <- 1


## COMPILE AND BUILD OBJECT
compile("TMB/src/goa_pk_tmb.cpp", flag='-w')
dyn.load('TMB/src/goa_pk_tmb.dll')
random <- c("slp1_fsh_dev", "inf1_fsh_dev", "slp2_fsh_dev", "inf2_fsh_dev")
obj <- MakeADFun(data=dat, parameters=pars, map=map, random=random, silent=FALSE)


## OPTIMIZE
opt <- with(obj, nlminb(par,fn,gr))
# opt <- TMBhelper::fit_tmb(obj, control=list(trace=50))
sdrep <- sdreport(obj)
quantities <- obj$report(obj$env$last.par.best)


## GET OBJECTS
stdtmb <- with(sdrep, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
params <- obj$env$parList()



# MODEL 2 ----
map2 <- map

# Adjust map
# - AR1 on age from mean sel param
map2$ln_sel_sd <- as.factor(1)

# - Non-parametric params
map2$selpars_re <- as.factor(1:dat$trmage)
map2$mean_sel <- as.factor(NA)

# - AR params
map2$sel_rho <- as.factor(1)
map2$sel_rho_y <- as.factor(NA)

# - Turn off logistic params
map2$slp1_fsh_dev <- as.factor(rep(NA, length(pars$slp1_fsh_dev)))
map2$inf1_fsh_dev <- as.factor(rep(NA, length(pars$inf1_fsh_dev)))
map2$slp2_fsh_dev <- as.factor(rep(NA, length(map2$slp2_fsh_dev)))
map2$inf2_fsh_dev <- as.factor(rep(NA, length(map2$inf2_fsh_dev)))

map2$inf1_fsh_mean <- as.factor(NA)
map2$inf2_fsh_mean <- as.factor(NA)
map2$log_slp1_fsh_mean <- as.factor(NA)
map2$log_slp2_fsh_mean <- as.factor(NA)

# - Ad fishery selectivity switch
dat$seltype <- 2

## COMPILE AND BUILD OBJECT
compile("TMB/src/goa_pk_tmb.cpp", flag='-w')
dyn.load('TMB/src/goa_pk_tmb.dll')
random <- c("selpars_re")
obj2 <- MakeADFun(data=dat, parameters=pars, map=map2, random=random, silent=FALSE)


## OPTIMIZE
opt2 <- with(obj2, nlminb(par,fn,gr))
# opt <- TMBhelper::fit_tmb(obj, control=list(trace=50))
sdrep2 <- sdreport(obj)
quantities2 <- obj$report(obj2$env$last.par.best)


## GET OBJECTS
stdtmb2 <- with(sdrep2, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
params2 <- obj2$env$parList()




## PLOTS ----
ssb <- rbind(cbind(model='TMB-RE',filter(stdtmb, par=='Espawnbio')),
             cbind(model='ADMB',filter(stdadmb, name=='Espawnbio') %>% select(par=name, est,se, year)))

filter(ssb, year>2015) %>% arrange(year)

ggplot(ssb, aes(year, est, color=model, fill=model, ymin=est-1.96*se, ymax=est+1.96*se)) +
  geom_ribbon(alpha=.5) + geom_line(lwd=2)

