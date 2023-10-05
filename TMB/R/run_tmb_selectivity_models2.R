### Run a series of fisheries selex models to compare. Modified
### from earlier work by Grant Adams.

## LIBRARIES AND DATA ----
source("TMB/R/prepare_tmb_objects.R")
source("TMB/R/Functions/phaser.R")
source("TMB/R/Functions/create_bounds.R")
library(tmbstan)
library(shinystan)

## Quick plot of selex
plot_fsh_selex <- function(fit){
  modelname <- paste('Model', fit$modelnum)[1]
  years <- 1970:(2022+5)
  ## df <- fit$std %>% filter(par=='slctfsh') %>%
  ##   mutate(age=rep(1:10, each=length(years)),
  ##          year=rep(years, times=10), model=model)
  sel <- t(fit$report$slctfsh)
  persp(y = years, x =  1:10,
        z = sel,
        col = "white",
        xlab = "Age", ylab = "\n\nYear",
        zlab = "\n\nSelectivity", expand = 0.5,
        box = TRUE, ticktype = "detailed", phi = 35,
        theta = -19, main = modelname)
}

## COMPILE AND BUILD TMB ----
compile("TMB/src/goa_pk_tmb.cpp")
dyn.load('TMB/src/goa_pk_tmb.dll')

## MODEL 0 ----
dat$seltype <- 0
## - Double logistic constant in time
map0 <- map
random <- NULL
obj0 <- MakeADFun(data=dat, parameters=pars, map=map0, random=random, silent=TRUE)
lwr <- get_bounds(obj0)$lwr
upr <- get_bounds(obj0)$upr
fit0 <- fit_tmb(obj0, upper=upr, lower=lwr,
                control=control, newtonsteps=1)
fit0$report <- obj0$report(obj0$env$last.par.best)
fit0$params <- obj0$env$parList(obj0$env$last.par.best)
fit0$std <- with(fit0$SD, data.frame(par=names(value), est=value, se=sd)) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
fit0$parList <- obj0$env$parList()
fit0$modelnum <- 0
plot_fsh_selex(fit0)

## MODEL 1 ----
dat$seltype <- 1
## - Double logistic with random effects on ascending portion
## - Turn on sel sd and slp and inf parameters
map1 <- map
map1$slp1_fsh_dev <- as.factor(1:length(pars$slp1_fsh_dev))
map1$inf1_fsh_dev <- as.factor(1:length(pars$inf1_fsh_dev))
map1$ln_sel_sd <- as.factor(1)
### cannot estimate these as they are perfectly confounded with
### the devs b/c not using a dev_vector approach, so we have to
### drop a degree of freedom. Arbitrarily setting the means to
### the previous MLE.
# map1$inf1_fsh_mean <- as.factor(1)
map1$inf2_fsh_mean <- as.factor(1)
# map1$log_slp1_fsh_mean <- as.factor(1)
map1$log_slp2_fsh_mean <- as.factor(1)
random <- c("slp1_fsh_dev", "inf1_fsh_dev", "slp2_fsh_dev", "inf2_fsh_dev")[1:2]
obj1 <- MakeADFun(data=dat, parameters=pars, map=map1, random=random, silent=TRUE)
lwr <- get_bounds(obj1)$lwr
upr <- get_bounds(obj1)$upr
fit1 <- fit_tmb(obj1, upper=upr, lower=lwr,
                control=control, newtonsteps=1)
fit1$report <- obj1$report(obj1$env$last.par.best)
fit1$std <- with(fit1$SD, data.frame(par=names(value), est=value, se=sd)) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
fit1$parList <- obj1$env$parList()
fit1$modelnum <- 1
plot_fsh_selex(fit1)


## MODEL 2 ----
## - Semiparametric double logistic with AR1 **age** random effects
## - Turn on sel sd and slp and inf parameters and ranef vector and rho
dat$seltype <- 2
map2 <- map
# -- Random effect pars
map2$selpars_re <- as.factor(1:length(pars$selpars_re))
# -- Fixed effect pars
map2$ln_sel_sd <- as.factor(1)
map2$sel_rho_a <- as.factor(1)
map2$inf1_fsh_mean <- as.factor(1)
map2$inf2_fsh_mean <- as.factor(1)
map2$log_slp1_fsh_mean <- as.factor(1)
map2$log_slp2_fsh_mean <- as.factor(1)
random <- c("selpars_re")
obj2 <- MakeADFun(data=dat, parameters=pars, map=map2, random=random, silent=TRUE)
lwr <- get_bounds(obj2)$lwr
upr <- get_bounds(obj2)$upr
fit2 <- fit_tmb(obj2, upper=upr, lower=lwr,
                control=control, newtonsteps=1)
fit2$report <- obj2$report(obj2$env$last.par.best)
fit2$std <- with(fit2$SD, data.frame(par=names(value), est=value, se=sd)) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
fit2$parList <- obj2$env$parList()
fit2$modelnum <- 2
plot_fsh_selex(fit2)

## MODEL 3 ----
## - Semiparametric double logistic with AR1 **year** random effects on function
## - Turn on sel sd and slp and inf parameters and ranef vector and rho y
dat$seltype <- 3
map3 <- map
pars3 <- pars
pars3$ln_sel_sd <- 0
pars3$selpars_re <- matrix(0, nrow = 1, ncol = dat$nyrs + dat$projfsh_nyrs)
map3$selpars_re <- as.factor(1:length(pars3$selpars_re))
map3$ln_sel_sd <- as.factor(1)
map3$sel_rho_y <- as.factor(1)
map3$inf1_fsh_mean <- as.factor(1)
map3$inf2_fsh_mean <- as.factor(1)
map3$log_slp1_fsh_mean <- as.factor(1)
map3$log_slp2_fsh_mean <- as.factor(1)
random <- c("selpars_re")
obj3 <- MakeADFun(data=dat, parameters=pars3, map=map3, random=random, silent=TRUE)
lwr <- get_bounds(obj3)$lwr
upr <- get_bounds(obj3)$upr
fit3 <- fit_tmb(obj3, upper=upr, lower=lwr, loopnum=3,
                control=control, newtonsteps=2)
fit3$report <- obj3$report(obj3$env$last.par.best)
fit3$std <- with(fit3$SD, data.frame(par=names(value), est=value, se=sd)) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
fit3$parList <- obj3$env$parList()
fit3$modelnum <- 3
plot_fsh_selex(fit3)


## MODEL 4 ----
## - Semiparametric double logistic with 2D-AR1 age, year random effects
## - Turn on sel sd and slp and inf parameters and ranef vector and rho
dat$seltype <- 4
map4 <- map
pars4 <- pars
# -- Random effect pars
pars4$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs + dat$projfsh_nyrs))
map4$selpars_re <- as.factor(1:length(pars4$selpars_re))
# -- Fixed effect pars
map4$ln_sel_sd <- as.factor(1)
map4$sel_rho_a <- as.factor(1)
map4$sel_rho_y <- as.factor(1)
map4$inf1_fsh_mean <- as.factor(1)
map4$inf2_fsh_mean <- as.factor(1)
map4$log_slp1_fsh_mean <- as.factor(1)
map4$log_slp2_fsh_mean <- as.factor(1)
random <- c("selpars_re")
obj4 <- MakeADFun(data=dat, parameters=pars4, map=map4, random=random, silent=TRUE)
lwr <- get_bounds(obj4)$lwr
upr <- get_bounds(obj4)$upr
fit4 <- fit_tmb(obj4, upper=upr, lower=lwr, loopnum=3,
                control=control, newtonsteps=2)
fit4$report <- obj4$report(obj4$env$last.par.best)
fit4$std <- with(fit4$SD, data.frame(par=names(value), est=value, se=sd)) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
fit4$parList <- obj4$env$parList()
fit4$modelnum <- 4
plot_fsh_selex(fit4)


## MODEL 5 ----
## - Age-specific fixed effects (no random effects)
## - Turn on mean_sel
dat$seltype <- 5
map5 <- map
map5$mean_sel <- as.factor(1:dat$nages)
obj5 <- MakeADFun(data=dat, parameters=pars, map=map5, silent=TRUE)
lwr <- get_bounds(obj5)$lwr
upr <- get_bounds(obj5)$upr
fit5 <- fit_tmb(obj5, upper=upr, lower=lwr, loopnum=3,
                control=control, newtonsteps=2)
fit5$report <- obj5$report(obj5$env$last.par.best)
fit5$std <- with(fit5$SD, data.frame(par=names(value), est=value, se=sd)) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
fit5$parList <- obj5$env$parList()
fit5$modelnum <- 5
plot_fsh_selex(fit5)


## MODEL 6 ----
## - Non-parametric 1D-AR1 year random effects
## - Turn on mean_sel, rho_y, sd, and ranef vector
dat$seltype <- 6
map6 <- map
pars6 <- pars
pars6$selpars_re <- matrix(0, nrow = 1, ncol = dat$nyrs + dat$projfsh_nyrs)
map6$selpars_re <- as.factor(1:length(pars6$selpars_re))
map6$ln_sel_sd <- as.factor(1)
map6$sel_rho_y <- as.factor(1)
map6$mean_sel <- as.factor(1:dat$nages)
random <- c("selpars_re")
obj6 <- MakeADFun(data=dat, parameters=pars6, map=map6, random=random, silent=TRUE)
lwr <- get_bounds(obj6)$lwr
upr <- get_bounds(obj6)$upr
fit6 <- fit_tmb(obj6, upper=upr, lower=lwr, loopnum=3,
                control=control, newtonsteps=3)
fit6$report <- obj6$report(obj6$env$last.par.best)
fit6$std <- with(fit6$SD, data.frame(par=names(value), est=value, se=sd)) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
fit6$parList <- obj6$env$parList()
fit6$modelnum <- 6
plot_fsh_selex(fit6)


## MODEL 7 ----
## - Non-parametric 2D-AR1 age, year random effects
## - Turn on mean_sel, rho, sd, and ranef vector
dat$seltype <- 7
map7 <- map
pars7 <- pars
pars7$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs + dat$projfsh_nyrs))
map7$selpars_re <- as.factor(1:length(pars7$selpars_re))
map7$ln_sel_sd <- as.factor(1)
map7$sel_rho_a <- as.factor(1)
map7$sel_rho_y <- as.factor(1)
map7$mean_sel <- as.factor(1:dat$nages)
random <- c("selpars_re")
obj7 <- MakeADFun(data=dat, parameters=pars7, map=map7, random=random, silent=TRUE)
lwr <- get_bounds(obj7)$lwr
upr <- get_bounds(obj7)$upr
fit7 <- fit_tmb(obj7, upper=upr, lower=lwr, loopnum=3,
                control=control, newtonsteps=2)
fit7$report <- obj7$report(obj7$env$last.par.best)
fit7$std <- with(fit7$SD, data.frame(par=names(value), est=value, se=sd)) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
fit7$parList <- obj7$env$parList()
fit7$modelnum <- 7
plot_fsh_selex(fit7)


## MODEL 8 ----
## - Non-parametric 3D-AR1 age, year, cohort random effects using **conditional** var
## - Turn on mean_sel, rho (a,y,c), sd, and ranef vector
dat$seltype <- 8
dat$sel_vartype <- 0
map8 <- map
pars8 <- pars
pars8$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs + dat$projfsh_nyrs))
map8$selpars_re <- as.factor(1:length(pars8$selpars_re))
map8$ln_sel_sd <- as.factor(1)
map8$sel_rho_a <- as.factor(1)
map8$sel_rho_y <- as.factor(1)
map8$sel_rho_c <- as.factor(1)
map8$mean_sel <- as.factor(1:dat$nages)
random <- c("selpars_re")
obj8 <- MakeADFun(data=dat, parameters=pars8, map=map8, random=random, silent=TRUE)
lwr <- get_bounds(obj8)$lwr
upr <- get_bounds(obj8)$upr
fit8 <- fit_tmb(obj8, upper=upr, lower=lwr, loopnum=3,
                control=control, newtonsteps=2)
fit8$report <- obj8$report(obj8$env$last.par.best)
fit8$std <- with(fit8$SD, data.frame(par=names(value), est=value, se=sd)) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
fit8$parList <- obj8$env$parList()
fit8$modelnum <- 8
plot_fsh_selex(fit8)



## MODEL 9 ----
## - Non-parametric 3D-AR1 age, year, cohort random effects using **marginal** var
## - Turn on mean_sel, rho (a,y,c), sd, and ranef vector
dat$seltype <- 8 # not 9! see cpp Q_sparse calcs
dat$sel_vartype <- 1
map9 <- map
pars9 <- pars
## pars9$sel_rho_a <- pars9$sel_rho_y <- pars9$sel_rho_c <- .5 # better inits?
pars9$mean_sel <- c(seq(-1,5,len=8),2,1)
pars9$ln_sel_sd <- .1
pars9$selpars_re <- array(0, dim = c(dat$trmage, dat$nyrs + dat$projfsh_nyrs))
map9$selpars_re <- as.factor(1:length(pars9$selpars_re))
map9$ln_sel_sd <- as.factor(1)
map9$sel_rho_a <- as.factor(1)
map9$sel_rho_y <- as.factor(1)
map9$sel_rho_c <- as.factor(1)
map9$mean_sel <- as.factor(1:dat$nages)
random <- c("selpars_re")
obj9 <- MakeADFun(data=dat, parameters=pars9, map=map9, random=random, silent=TRUE)
obj9$fn()
lwr <- get_bounds(obj9)$lwr
upr <- get_bounds(obj9)$upr
fit9 <- fit_tmb(obj9, upper=upr, lower=lwr, loopnum=3,
                control=control, newtonsteps=3)
fit9$report <- obj9$report(obj9$env$last.par.best)
## fit9$std <- with(fit9$SD, data.frame(par=names(value), est=value, se=sd)) %>%
##   group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup
fit9$parList <- obj9$env$parList()
fit9$modelnum <- 9
plot_fsh_selex(fit9)


## Save Image ----
fits <- list(fit0=fit0,fit1=fit1,fit2=fit2,fit3=fit3,fit4=fit4,fit5=fit5,fit6=fit6,fit7=fit7,fit8=fit8,fit9=fit9)
saveRDS(fits, 'TMB/Output/fits.RDS')
## save.image(file='TMB/Selectivity_runs.RData')


## q <- fit9$report$Q_sparse %>% as.matrix
## q[1:20, 1:20]
## sum(q)

### quick checks on the 3D models
cor7 <- cov2cor(solve(fits$fit7$report$Q_sparse))
cor8 <- cov2cor(solve(fits$fit8$report$Q_sparse))
corrplot::corrplot(cor7[1:20,1:20])
corrplot::corrplot(cor8[1:20,1:20])

with(fit8$SD, data.frame(par=names(par.fixed), est=par.fixed, se=sqrt(diag(cov.fixed)))) %>%
  filter(grepl('rho|ln', par)) %>% mutate(lwr=est-1.96*se, upr=est+1.96*se)%>% print(digits=2)

fit9$opt$diagnostics %>% filter(grepl('rho|ln', Param))  %>%
  print(digits=2)


sm8 <- fit8$parList$mean_sel
sm9 <- fit9$parList$mean_sel

plot(sm9, type='b', ylim=c(-10,10))
lines(sm8)



par(mfrow=c(3,4), mar = c(3.2, 3.2 , 1 , 0.5) ,
    oma = c(0 , 0 , 0 , 0), tcl = -0.8, mgp = c(10, 0.6, 0))
trash <- lapply(fits, plot_fsh_selex)

sp8 <- fits$fit8$parList$selpars_re
sp7 <- fits$fit7$parList$selpars_re
zlim <- range(c(sp7, sp8))
par(mfrow=c(1,2))
image(t(sp7), y=1:10, x=years, ylab='age',  zlim=zlim, main='2D AR(1)')
image(t(sp8), y=1:10, x=years, ylab='age', zlim=zlim, main='3D AR(1) conditional')
