## LIBRARIES AND DATA ----
source("TMB/R/prepare_tmb_objects.R")
source("TMB/R/phaser.R")

map <- lapply(map, function(x) as.factor(as.numeric(x)*NA))


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

# -- Phase
# phased_pars_mod1 <- TMBphase(
#   data = dat,
#   parameters = pars,
#   map = map_mod1,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
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
map_mod2$sel_rho_a <- as.factor(1)

map_mod2$inf1_fsh_mean <- as.factor(1)
map_mod2$inf2_fsh_mean <- as.factor(1)
map_mod2$log_slp1_fsh_mean <- as.factor(1)
map_mod2$log_slp2_fsh_mean <- as.factor(1)

# -- Data switch
dat$seltype <- 2

# - Build model 2
random <- c("selpars_re")

# -- Phase
# phased_pars_mod2 <- TMBphase(
#   data = dat,
#   parameters = pars,
#   map = map_mod2,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
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
pars_mod3$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs))
map_mod3$selpars_re <- as.factor(1:length(pars_mod3$selpars_re))

# -- Fixed effect pars
map_mod3$ln_sel_sd <- as.factor(1)
map_mod3$sel_rho_a <- as.factor(1)
map_mod3$sel_rho_y <- as.factor(1)

map_mod3$inf1_fsh_mean <- as.factor(1)
map_mod3$inf2_fsh_mean <- as.factor(1)
map_mod3$log_slp1_fsh_mean <- as.factor(1)
map_mod3$log_slp2_fsh_mean <- as.factor(1)

# -- Data switch
dat$seltype <- 3

# - Build model 3
random <- c("selpars_re")

# -- Phase
# phased_pars_mod3 <- TMBphase(
#   data = dat,
#   parameters = pars_mod3,
#   map = map_mod3,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
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
map_mod4$sel_rho_a <- as.factor(1)
map_mod4$mean_sel <- as.factor(1)

# -- Data switch
dat$seltype <- 4

# - Build model 4
random <- c("selpars_re")

# -- Phase
# phased_pars_mod4 <- TMBphase(
#   data = dat,
#   parameters = pars,
#   map = map_mod4,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
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
pars_mod5$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs))
map_mod5$selpars_re <- as.factor(1:length(pars_mod5$selpars_re))

# -- Fixed effect pars
map_mod5$ln_sel_sd <- as.factor(1)
map_mod5$sel_rho_a <- as.factor(1)
map_mod5$sel_rho_y <- as.factor(1)
map_mod5$mean_sel <- as.factor(1)

# -- Data switch
dat$seltype <- 5

# - Build model 5
random <- c("selpars_re")

# -- Phase
# phased_pars_mod5 <- TMBphase(
#   data = dat,
#   parameters = pars_mod5,
#   map = map_mod5,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod5 <- MakeADFun(data=dat, parameters=pars_mod5, map=map_mod5, random=random, silent=TRUE)

# - Optimize
opt_mod5 <- with(obj_mod5, nlminb(par,fn,gr))
# opt <- TMBhelper::fit_tmb(obj_mod5)
sdrep_mod5 <- sdreport(obj_mod5)
quantities_mod5 <- obj_mod5$report(obj_mod5$env$last.par.best)

# - Get objects
stdtmb_mod5 <- with(sdrep_mod5, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
params_mod5 <- obj_mod5$env$parList()


## MODEL 6 ----
# - Non-parametric 3D-AR1 age, year, cohort random effects using conditional var
# - Turn on mean_sel, rho (a,y,c), sd, and ranef vector
map_mod6 <- map
pars_mod6 <- pars

# -- Random effect pars
pars_mod6$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs))
map_mod6$selpars_re <- as.factor(1:length(pars_mod6$selpars_re))

# -- Fixed effect pars
map_mod6$ln_sel_sd <- as.factor(1)
map_mod6$sel_rho_a <- as.factor(1)
map_mod6$sel_rho_y <- as.factor(1)
map_mod6$sel_rho_c <- as.factor(1)
map_mod6$mean_sel <- as.factor(1)

# -- Data switch
dat$seltype <- 6
dat$sel_vartype <- 0

# - Build model 6
random <- c("selpars_re")

# -- Phase
# phased_pars_mod6 <- TMBphase(
#   data = dat,
#   parameters = pars_mod6,
#   map = map_mod6,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod6 <- MakeADFun(data=dat, parameters=pars_mod6, map=map_mod6, random=random, silent=TRUE)

# - Optimize
opt_mod6 <- with(obj_mod6, nlminb(par,fn,gr))
# opt <- TMBhelper::fit_tmb(obj_mod6)
sdrep_mod6 <- sdreport(obj_mod6)
quantities_mod6 <- obj_mod6$report(obj_mod6$env$last.par.best)

# - Get objects
stdtmb_mod6 <- with(sdrep_mod6, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
params_mod6 <- obj_mod6$env$parList()


## MODEL 7 ----
# - Non-parametric 3D-AR1 age, year, cohort random effects using marginal var
# - Turn on mean_sel, rho (a,y,c), sd, and ranef vector
map_mod7 <- map
pars_mod7 <- pars

# -- Random effect pars
pars_mod7$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs))
map_mod7$selpars_re <- as.factor(1:length(pars_mod7$selpars_re))

# -- Fixed effect pars
map_mod7$ln_sel_sd <- as.factor(1)
map_mod7$sel_rho_a <- as.factor(1)
map_mod7$sel_rho_y <- as.factor(1)
map_mod7$sel_rho_c <- as.factor(1)
map_mod7$mean_sel <- as.factor(1)

# -- Data switch
dat$seltype <- 6
dat$sel_vartype <- 1

# - Build model 7
random <- c("selpars_re")

# -- Phase
# phased_pars_mod7 <- TMBphase(
#   data = dat,
#   parameters = pars_mod7,
#   map = map_mod7,
#   random = random,
#   phases = phases,
#   model_name = "goa_pk_tmb",
#   silent = TRUE,
#   use_gradient = TRUE)

# -- Build model obj
obj_mod7 <- MakeADFun(data=dat, parameters=pars_mod7, map=map_mod7, random=random, silent=TRUE)

# - Optimize
opt_mod7 <- with(obj_mod7, nlminb(par,fn,gr))
# opt <- TMBhelper::fit_tmb(obj_mod7)
sdrep_mod7 <- sdreport(obj_mod7)
quantities_mod7 <- obj_mod7$report(obj_mod7$env$last.par.best)

# - Get objects
stdtmb_mod7 <- with(sdrep_mod7, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
params_mod7 <- obj_mod7$env$parList()


## PLOTS ----
# - Combine
ssb <- rbind(cbind(model='TMB-Mod1-RE',filter(stdtmb_mod1, par=='Espawnbio')),
             cbind(model='TMB-Mod2-Log-AR1',filter(stdtmb_mod2, par=='Espawnbio')),
             cbind(model='TMB-Mod3-Log-2D-AR1',filter(stdtmb_mod3, par=='Espawnbio')),
             cbind(model='TMB-Mod4-AR1',filter(stdtmb_mod4, par=='Espawnbio')),
             cbind(model='TMB-Mod5-2D-AR1',filter(stdtmb_mod5, par=='Espawnbio')),
             cbind(model='TMB-Mod6-3D-AR1cond',filter(stdtmb_mod6, par=='Espawnbio')),
             cbind(model='TMB-Mod7-3D-AR1mar',filter(stdtmb_mod7, par=='Espawnbio')),
             cbind(model='ADMB',filter(stdadmb, name=='Espawnbio') %>% select(par=name, est,se, year)))

# - Plot
ggplot(ssb, aes(year, est, color=model, fill=model, ymin=est-1.96*se, ymax=est+1.96*se)) +
  geom_ribbon(alpha=.5) + geom_line(lwd=2)

write.csv(ssb, file = "TMB/Output/Initial_runs.csv")

# - Show recent ssb
ssb %>%
  filter(year>2015) %>%
  arrange(year) %>%
  pivot_wider(names_from = model, values_from = c(est, se)) %>%
  select(year, paste0(c("est_","se_"),rep(unique(ssb$model), each = 2))) %>%
  as.data.frame()



