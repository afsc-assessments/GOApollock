## LIBRARIES AND DATA ----
source("TMB/R/prepare_tmb_objects.R")

# map <- lapply(pars, function(x) as.factor(x*NA))


## COMPILE AND BUILD TMB ----
compile("TMB/src/goa_pk_tmb.cpp")
dyn.load('TMB/src/goa_pk_tmb.dll')

## MODEL 1 ----
# - Double logistic with random effects on ascending portion
# - Turn on sel sd and slp and inf parameters
map_mod1 <- map

# -- Random effect pars
map_mod1$slp1_fsh_dev <- as.factor(1:length(pars$slp1_fsh_dev))
map_mod1$inf1_fsh_dev <- as.factor(1:length(pars$inf1_fsh_dev))
#map_mod1$slp2_fsh_dev <- as.factor(1:length(pars$slp2_fsh_dev))
#map_mod1$inf2_fsh_dev <- as.factor(1:length(pars$inf2_fsh_dev))

# -- Fixed effect pars
map_mod1$ln_sel_sd <- as.factor(1)
map_mod1$inf1_fsh_mean <- as.factor(1)
map_mod1$inf2_fsh_mean <- as.factor(1)
map_mod1$log_slp1_fsh_mean <- as.factor(1)
map_mod1$log_slp2_fsh_mean <- as.factor(1)

# -- Data switch
dat$seltype <- 1

# - Build model 1
random <- c("slp1_fsh_dev", "inf1_fsh_dev", "slp2_fsh_dev", "inf2_fsh_dev")
obj_mod1 <- MakeADFun(data=dat, parameters=pars, map=map_mod1, random=random, silent=TRUE)

# - Optimize
opt_mod1 <- with(obj_mod1, nlminb(par,fn,gr))
# opt <- TMBhelper::fit_tmb(obj)
sdrep_mod1 <- sdreport(obj_mod1)
quantities_mod1 <- obj_mod1$report(obj_mod1$env$last.par.best)

# - Get objects
stdtmb_mod1 <- with(sdrep_mod1, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
params_mod1 <- obj_mod1$env$parList()


## MODEL 2 ----
# - Double logistic with AR1 age random effects on function
# - Turn on sel sd and slp and inf parameters and ranef vector and rho
map_mod2 <- map

# -- Random effect pars
map_mod2$selpars_re <- as.factor(1:length(pars$selpars_re))

# -- Fixed effect pars
map_mod2$ln_sel_sd <- as.factor(1)
map_mod2$sel_rho <- as.factor(1)

map_mod2$inf1_fsh_mean <- as.factor(1)
map_mod2$inf2_fsh_mean <- as.factor(1)
map_mod2$log_slp1_fsh_mean <- as.factor(1)
map_mod2$log_slp2_fsh_mean <- as.factor(1)

# -- Data switch
dat$seltype <- 2

# - Build model 1
random <- c("selpars_re")
obj_mod2 <- MakeADFun(data=dat, parameters=pars, map=map_mod2, random=random, silent=TRUE)

# - Optimize
opt_mod2 <- with(obj_mod2, nlminb(par,fn,gr))
# opt <- TMBhelper::fit_tmb(obj)
sdrep_mod2 <- sdreport(obj_mod2)
quantities_mod2 <- obj_mod2$report(obj_mod2$env$last.par.best)

# - Get objects
stdtmb_mod2 <- with(sdrep_mod2, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
params_mod2 <- obj_mod2$env$parList()


## MODEL 3 ----
# - Double logistic with 2D-AR1 age, year random effects on function
# - Turn on sel sd and slp and inf parameters and ranef vector and rho
map_mod3 <- map
pars_mod3 <- pars

# -- Random effect pars
pars_mod3$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs, 1))
map_mod3$selpars_re <- as.factor(1:length(pars_mod3$selpars_re))

# -- Fixed effect pars
map_mod3$ln_sel_sd <- as.factor(1)
map_mod3$sel_rho <- as.factor(1)
map_mod3$sel_rho_y <- as.factor(1)

map_mod3$inf1_fsh_mean <- as.factor(1)
map_mod3$inf2_fsh_mean <- as.factor(1)
map_mod3$log_slp1_fsh_mean <- as.factor(1)
map_mod3$log_slp2_fsh_mean <- as.factor(1)

# -- Data switch
dat$seltype <- 3

# - Build model 1
random <- c("selpars_re")
obj_mod3 <- MakeADFun(data=dat, parameters=pars_mod3, map=map_mod3, random=random, silent=TRUE)

# - Optimize
opt_mod3 <- with(obj_mod3, nlminb(par,fn,gr))
# opt <- TMBhelper::fit_tmb(obj)
sdrep_mod3 <- sdreport(obj_mod3)
quantities_mod3 <- obj_mod3$report(obj_mod3$env$last.par.best)

# - Get objects
stdtmb_mod3 <- with(sdrep_mod3, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
params_mod3 <- obj_mod3$env$parList()


## MODEL 4 ----
# - Non-parametric 1D-AR1 age random effects
# - Turn on mean_sel, rho, sd, and ranef vector
map_mod4 <- map

# -- Random effect pars
map_mod4$selpars_re <- as.factor(1:length(pars$selpars_re))

# -- Fixed effect pars
map_mod4$ln_sel_sd <- as.factor(1)
map_mod4$sel_rho <- as.factor(1)
map_mod4$mean_sel <- as.factor(1)

# -- Data switch
dat$seltype <- 4

# - Build model 1
random <- c("selpars_re")
obj_mod4 <- MakeADFun(data=dat, parameters=pars, map=map_mod4, random=random, silent=TRUE)

# - Optimize
opt_mod4 <- with(obj_mod4, nlminb(par,fn,gr))
# opt <- TMBhelper::fit_tmb(obj)
sdrep_mod4 <- sdreport(obj_mod4)
quantities_mod4 <- obj_mod4$report(obj_mod4$env$last.par.best)

# - Get objects
stdtmb_mod4 <- with(sdrep_mod4, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
params_mod4 <- obj_mod4$env$parList()


## MODEL 5 ----
# - Non-parametric 2D-AR1 age, year random effects
# - Turn on mean_sel, rho, sd, and ranef vector
map_mod5 <- map
pars_mod5 <- pars

# -- Random effect pars
pars_mod5$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs, 1))
map_mod5$selpars_re <- as.factor(1:length(pars_mod5$selpars_re))

# -- Fixed effect pars
map_mod5$ln_sel_sd <- as.factor(1)
map_mod5$sel_rho <- as.factor(1)
map_mod5$sel_rho_y <- as.factor(1)
map_mod5$mean_sel <- as.factor(1)

# -- Data switch
dat$seltype <- 5

# - Build model 1
random <- c("selpars_re")
obj_mod5 <- MakeADFun(data=dat, parameters=pars_mod5, map=map_mod5, random=random, silent=TRUE)

# - Optimize
opt_mod5 <- with(obj_mod5, nlminb(par,fn,gr))
# opt <- TMBhelper::fit_tmb(obj)
sdrep_mod5 <- sdreport(obj_mod5)
quantities_mod5 <- obj_mod5$report(obj_mod5$env$last.par.best)

# - Get objects
stdtmb_mod5 <- with(sdrep_mod5, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
params_mod5 <- obj_mod5$env$parList()


## PLOTS ----
# - Combine
ssb <- rbind(cbind(model='TMB-RE',filter(stdtmb_mod1, par=='Espawnbio')),
             cbind(model='TMB-Log-AR1',filter(stdtmb_mod2, par=='Espawnbio')),
             cbind(model='TMB-Log-2D-AR1',filter(stdtmb_mod3, par=='Espawnbio')),
             cbind(model='TMB-AR1',filter(stdtmb_mod4, par=='Espawnbio')),
             # cbind(model='TMB-2D-AR1',filter(stdtmb_mod5, par=='Espawnbio')),
             # cbind(model='TMB-3D-AR1',filter(stdtmb_mod6, par=='Espawnbio')),
             cbind(model='ADMB',filter(stdadmb, name=='Espawnbio') %>% select(par=name, est,se, year)))

# - Plot
ggplot(ssb, aes(year, est, color=model, fill=model, ymin=est-1.96*se, ymax=est+1.96*se)) +
  geom_ribbon(alpha=.5) + geom_line(lwd=2)

# - Show recent ssb
ssb %>%
  filter(year>2015) %>%
  arrange(year)



