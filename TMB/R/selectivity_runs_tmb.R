## LIBRARIES
library(GOApollock)
library(TMBhelper)
library(ggplot2)
library(TMB)
library(dplyr)


## Converting random walk fishery selectivity deviates to random effects and estimating sd
# - NOTE: 4*sd for inf1_fsh_dev was removed


## READ IN 2022 DATA
dat <- readRDS('TMB/data/datfile.RDS')
replist <- readRDS("TMB/data/repfile.RDS")
stdadmb <- readRDS("TMB/data/stdfile.RDS")
pars <- readRDS('TMB/data/pars.RDS')
map <- readRDS('TMB/data/map.RDS')

# - Adjust fishery random walk standard deviation to parameter
pars$ln_rwlk_sd <- log(dat$rwlk_sd)
map$ln_rwlk_sd <- as.factor(rep(1, length(pars$ln_rwlk_sd)))


## SET UP 0 PARAMETERS
pars0 <- lapply(pars, function(x) replace(x, values = rep(0, length(x))))
pars0$natMscalar <- 1
pars0$mean_log_recruit <- 9
pars0$sigmaR <- 1.3


## SET UP MAP
map0 <- lapply(map, function(x) as.numeric(x))
map0 <- lapply(map0, function(x) replace(x, values = 1:length(x)))
map0$natMscalar <- NA
map0$sigmaR <- NA
map0 <- lapply(map0, function(x) as.factor(x))


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


## PLOTS
ssb <- rbind(cbind(model='TMB-RE',filter(stdtmb, par=='Espawnbio')),
             cbind(model='ADBM',filter(stdadmb, name=='Espawnbio') %>% select(par=name, est,se, year)))

filter(ssb, year>2015) %>% arrange(year)

ggplot(ssb, aes(year, est, color=model, fill=model, ymin=est-1.96*se, ymax=est+1.96*se)) +
  geom_ribbon(alpha=.5) + geom_line(lwd=2)

