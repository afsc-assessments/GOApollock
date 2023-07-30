## LIBRARIES
# library(GOApollock)
library(TMBhelper)
library(ggplot2)
library(TMB)
library(dplyr)

## READ IN 2022 DATA ----
dat <- readRDS('TMB/data/datfile.RDS')
replist <- readRDS("TMB/data/repfile.RDS")
stdadmb <- readRDS("TMB/data/stdfile.RDS")
pars <- readRDS('TMB/data/pars.RDS')
map <- readRDS('TMB/data/map.RDS')


## ADD "NEW" PARAMETER OBJECTS
# - Double logistic with random effects on ascending portion
# - Adjust fishery random walk standard deviation to parameter
pars$ln_sel_sd <- log(dat$rwlk_sd[1]) # All the same so reducing to length of 1

# - Non-parametric params
pars$selpars_re <- array(rep(0, dat$trmage), dim = c(dat$trmage, 1, 1))

pars$mean_sel <- 0

# - AR params
pars$sel_rho_y <- 0
pars$sel_rho <- 0

# - Number of years
dat$nyrs <- dat$endyr - dat$styr + 1

# - Ad fishery selectivity switch
dat$seltype <- 1

## sET MAP TO NA FOR ALL PARAMS
map$slp1_fsh_dev <- as.factor(rep(NA,length(pars$slp1_fsh_dev)))
map$slp2_fsh_dev <- as.factor(rep(NA,length(pars$slp2_fsh_dev)))
map$inf1_fsh_dev <- as.factor(rep(NA,length(pars$inf1_fsh_dev)))
map$inf2_fsh_dev <- as.factor(rep(NA,length(pars$inf2_fsh_dev)))

# -- Fixed effect pars
map$ln_sel_sd <- as.factor(NA)
map$inf1_fsh_mean <- as.factor(NA)
map$inf2_fsh_mean <- as.factor(NA)
map$log_slp1_fsh_mean <- as.factor(NA)
map$log_slp2_fsh_mean <- as.factor(NA)


map$ln_sel_sd <- as.factor(NA)

# - Non-parametric params
map$selpars_re <- factor(rep(NA, length(pars$selpars_re)))
map$mean_sel <- as.factor(NA)
map$sel_rho_y <- as.factor(NA)
map$sel_rho <- as.factor(NA)
