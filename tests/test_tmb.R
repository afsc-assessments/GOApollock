
## devtools::install()
library(GOApollock)


dyn.unload('model_23/goa_pk_tmb.dll')
file.copy('../source/goa_pk_tmb.cpp', to='model_23/goa_pk_tmb.cpp', overwrite = TRUE)

## make a few changes to be able to explore better
input <- prepare_pk_input('model_23', datfile='goa_pk.dat')
tmb <- fit_pk(input, do.fit=TRUE)
input1 <- input
input1$dat$Ftarget <- c(input1$dat$Ftarget,0*input1$dat$Ftarget)
input1$version <- 'estimate sigmaR'
input1$map$sigmaR <- factor(1)
input1$random <- 'dev_log_recruit'
tmb1 <- fit_pk(input1, do.fit=TRUE)

z <- list(tmb,tmb1)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=TRUE, uselog=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=TRUE, uselog=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE, uselog=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=FALSE, uselog=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=TRUE, uselog=FALSE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=TRUE, uselog=FALSE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE, uselog=FALSE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=FALSE, uselog=FALSE)

plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=TRUE, uselog=TRUE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=TRUE, uselog=TRUE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE, uselog=TRUE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=FALSE, uselog=TRUE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=TRUE, uselog=FALSE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=TRUE, uselog=FALSE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE, uselog=FALSE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=FALSE, uselog=FALSE, addproj=TRUE)

ssb <- plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=FALSE,
            uselog=FALSE, plot=FALSE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE,
            uselog=FALSE, alpha1=.3)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE,
            uselog=FALSE, alpha2=.3)


x <- get_rep(z)
(test <- get_rep(z, 'q5'))
(test <- get_rep(z, 'ssb'))
(test <- get_rep(z, 'Espawnbio'))
(test <- get_rep(z, 'Espawnbio_proj'))
(test <- get_rep(z, 'slctfsh'))
(test <- get_rep(z))


## restrospectives
fits <- fit_pk_retros(fit=tmb)
fits <- fit_pk_retros(fit=tmb, peels=0:15, parallel=TRUE)
plot_pk_ssb(fits)

plot_pk_ssb(fits[[1]])

## Bridge to dnorm, dmultinom, and SIMULATE statements, make sure
## everything matches
library(GOApollock)
input0 <- prepare_pk_input(path='model_23', datfile='goa_pk.dat',
                           modfile='goa_pk_tmb_old', version='original')
fit0 <- fit_pk(input0, do.fit=TRUE, newtonsteps=5)
## fit0 <- readRDS("../data/2023/fit.RDS")
## saveRDS(fit0, file='fit_before_dnorm.RDS')
##fit0 <- readRDS('fit_before_dnorm.RDS')

## this is the old file before modifying the likelihood
## statements to have constants
input1 <- prepare_pk_input(path='model_23', datfile='goa_pk.dat', version='original')
input1$pars <- fit0$obj$env$parList(fit0$opt$par)
input1$version <- 'updated'
file.copy('../source/goa_pk_tmb.cpp', to='model_23/goa_pk_tmb.cpp', overwrite = TRUE)
rm(fit1)
fit1 <- fit_pk(input1, do.fit=TRUE, newtonsteps=5)

## difference between NLL to cancel out the constants
p1 <- fit0$opt$par
p2 <- p1+.1
r0 <- fit0$obj$rep(p1)$loglik-fit0$obj$rep(p2)$loglik
r1 <- fit1$obj$rep(p1)$loglik-fit1$obj$rep(p2)$loglik
round(r0-r1,5)
r0 <- fit0$obj$rep(p1)$llcatp-fit0$obj$rep(p2)$llcatp
r1 <- fit1$obj$rep(p1)$llcatp-fit1$obj$rep(p2)$llcatp
round(r0-r1,5)
r0 <- fit0$obj$rep(p1)$llsrvp1-fit0$obj$rep(p2)$llsrvp1
r1 <- fit1$obj$rep(p1)$llsrvp1-fit1$obj$rep(p2)$llsrvp1
round(r0-r1,5)

fit0$obj$rep(p1)$catp-fit1$obj$rep(p1)$catp
fit0$obj$rep(p1)$Ecatp-fit1$obj$rep(p1)$Ecatp

all.equal(fit0$rep, fit1$rep)
## drop proj stuff since didn't work in fit0
all.equal(filter(fit0$sd, !grepl('_proj',name)), filter(fit1$sd, !grepl('_proj',name)))
all.equal(fit0$opt$par,fit1$opt$par)
cbind(fit0$opt$par,fit1$opt$par)
cbind(fit0$rep$loglik,fit1$rep$loglik)

all.equal(fit0$rep$res_fish,fit1$rep$res_fish)
all.equal(fit0$rep$res_srv1,fit1$rep$res_srv1)
all.equal(fit0$rep$res_srv2,fit1$rep$res_srv2)
all.equal(fit0$rep$res_srv2,fit1$rep$res_srv2)
all.equal(fit0$rep$res_srv6,fit1$rep$res_srv6)

all.equal(fit0$rep$pearson_fish,fit1$rep$pearson_fish)
all.equal(fit0$rep$pearson_srv1,fit1$rep$pearson_srv1)
(fit0$rep$pearson_srv1-fit1$rep$pearson_srv1)
all.equal(fit0$rep$pearson_srv2,fit1$rep$pearson_srv2)
all.equal(fit0$rep$pearson_srv2,fit1$rep$pearson_srv2)
all.equal(fit0$rep$pearson_srv6,fit1$rep$pearson_srv6)

test <- get_rep(list(fit0,fit1), 'pearson_fish')
test <- test %>% pivot_wider(names_from=version, values_from=value) %>%
  mutate(diff=original-updated)


sdtest <- rbind(fit0$sd, fit1$sd) %>% filter(!grepl('_proj',name))
sdtest %>% group_by(name,year) %>%
  summarize(estdiff=est[version=='updated']-est[version!='updated']) %>%
  arrange(desc(abs(estdiff)))

xx <- fit1$obj$simulate()
xx$catp
xx$catp %>% rowSums
xx$srvp1
xx$srvp1 %>% rowSums



fit0$obj$fn()
fit1$obj$fn()
