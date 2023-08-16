## LIBRARIES
# library(GOApollock)
library(TMBhelper)
library(ggplot2)
library(TMB)
library(dplyr)
library(tidyr)

## READ IN 2022 DATA ----
dat <- readRDS('TMB/data/datfile.RDS')
replist <- readRDS("TMB/data/repfile.RDS")
stdadmb <- readRDS("TMB/data/stdfile.RDS")
pars <- readRDS('TMB/data/pars.RDS')
map <- readRDS('TMB/data/map.RDS')


## Selectivity projection years ----
dat$projfsh_nyrs <- 5
pars$slp1_fsh_dev <- c(pars$slp1_fsh_dev, rep(0, dat$projfsh_nyrs))
pars$slp2_fsh_dev <- c(pars$slp2_fsh_dev, rep(0, dat$projfsh_nyrs))
pars$inf1_fsh_dev <- c(pars$inf1_fsh_dev, rep(0, dat$projfsh_nyrs))
pars$inf2_fsh_dev <- c(pars$inf2_fsh_dev, rep(0, dat$projfsh_nyrs))


## BUILD FULL MAP ----
map_full <- lapply(pars, function(x) factor(1:length(x)))
for(i in 1:length(map)){
  map_full[names(map)[i]] <- map[names(map)[i]]
}
map <- map_full


## ADD "NEW" PARAMETER OBJECTS ----
# - Double logistic with random effects on ascending portion
# - Adjust fishery random walk standard deviation to parameter
pars$ln_sel_sd <- log(dat$rwlk_sd[1]) # All the same so reducing to length of 1

# - Non-parametric params
pars$selpars_re <- array(rep(0, dat$trmage), dim = c(dat$trmage, 1))

pars$mean_sel <- 0

# - AR params
pars$sel_rho_a <- 0
pars$sel_rho_y <- 0
pars$sel_rho_c <- 0


## ADD NEW DATA OBJECTS ----
# - Number of years
dat$nyrs <- (dat$endyr - dat$styr + 1)
dat$nages <- dat$trmage - dat$rcrage + 1;


# - Add fishery selectivity switch
dat$seltype <- 1
dat$sel_vartype <- 0

# Create an index for ages and years to feed into TMB, which helps construct the precision matrix
dat$ay_Index <- as.matrix(expand.grid("age" = seq_len(dat$nages),
                                  "year" = seq_len(dat$nyrs + dat$projfsh_nyrs) ))


## BUILD NEW MAP
# - SET TO NA FOR ALL FISHERIES SELECTIVTY PARAMS
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


# - Non-parametric params
map$selpars_re <- factor(rep(NA, length(pars$selpars_re)))
map$mean_sel <- as.factor(NA)
map$sel_rho_a <- as.factor(NA)
map$sel_rho_y <- as.factor(NA)
map$sel_rho_c <- as.factor(NA)


## Optimization control ----
control <- list(eval.max=10000, iter.max=10000)


## PHASES ----
phases <- list(
  dev_log_initN= 2,
  mean_log_recruit = 1,
  dev_log_recruit = 2,
  sigmaR = 5,
  log_recr_proj = 5,
  log_slp1_fsh_mean = 3,
  inf1_fsh_mean = 3,
  log_slp2_fsh_mean = 3,
  inf2_fsh_mean = 3,
  slp1_fsh_dev = 4,
  inf1_fsh_dev = 4,
  slp2_fsh_dev = 4,
  inf2_fsh_dev = 4,
  log_slp2_srv1 = 3,
  inf2_srv1 = 3,
  log_slp1_srv2 = 3,
  inf1_srv2 = 3,
  log_slp2_srv2 = 3,
  inf2_srv2 = 3,
  log_slp1_srv3 = 3,
  inf1_srv3 = 3,
  log_slp1_srv6 = 3,
  inf1_srv6 = 3,
  log_slp2_srv6 = 3,
  inf2_srv6 = 3,
  mean_log_F = 1,
  dev_log_F = 1,
  log_q1_mean = 3,
  log_q1_dev = 4,
  log_q2_mean = 3,
  log_q2_dev = 4,
  log_q3_mean = 3,
  log_q3_dev = 4,
  log_q4 = 3,
  q4_pow = 3,
  log_q5 = 3,
  q5_pow = 3,
  log_q6 = 3,
  natMscalar = 4,
  ln_sel_sd = 4,
  selpars_re = 4,
  mean_sel = 3,
  sel_rho_y = 4,
  sel_rho_c = 4,
  sel_rho_a = 4
)
