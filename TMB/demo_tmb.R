library(GOApollock)
library(TMBhelper)
library(ggplot2)
library(TMB)
library(compResidual) ## https://github.com/fishfollower/compResidual
source("functions.R")

## The dat and rep files from the 2022 final model
dat <- readRDS('datfile.RDS')
replist <- readRDS("repfile.RDS")
stdadmb <- readRDS("stdfile.RDS")
pars <- readRDS('pars.RDS')
map <- readRDS('map.RDS')
years <- 1970:2022

## Fit the model and compare to ADMB. Should be identical
compile("../source/goa_pk_tmb.cpp")
dyn.load('../source/goa_pk_tmb.dll')
obj <- MakeADFun(data=dat, parameters=pars, map=map, random=NULL, silent=TRUE)

## optimize and compare
opt <- with(obj, nlminb(par,fn,gr))
##opt <- TMBhelper::fit_tmb(obj, control=list(trace=50))
rep <- obj$report()
sdrep <- sdreport(obj)
stdtmb <- with(sdrep, data.frame(par=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(par) %>% mutate(year=1969+1:n()) %>% ungroup

ssb <- rbind(cbind(model='TMB',filter(stdtmb, par=='Espawnbio')),
             cbind(model='ADBM',filter(stdadmb, name=='Espawnbio') %>% select(par=name, est,se, year)))
filter(ssb, year>2015) %>% arrange(year)
ggplot(ssb, aes(year, est, color=model, fill=model, ymin=est-1.96*se, ymax=est+1.96*se)) +
  geom_ribbon(alpha=.5) + geom_line(lwd=2)

r1 <- plot_osa_comps(obs=replist$Fishery_observed_and_expected_age_comp[,1:10],
               exp=replist$Fishery_observed_and_expected_age_comp[,11:20],
               ages=1:10, years=dat$fshyrs, Neff=dat$multN_fsh,
               model='ADMB final 2022')
g1 <- plot_osa_comps(obs=replist$Fishery_observed_and_expected_age_comp[,1:10],
               exp=replist$Fishery_observed_and_expected_age_comp[,11:20],
               ages=1:10, years=dat$fshyrs, Neff=dat$multN_fsh,
               model='ADMB final 2022', plot.type='ggplot')

rep$res_fish[,c(1,11)] <- 0 # need to setZero in template

r2 <- plot_osa_comps(obs=rep$res_fish[,1:10], exp=rep$res_fish[,11:20],
               ages=1:10, years=dat$fshyrs, Neff=dat$multN_fsh,
               model='TMB port', plot.type='default')
g2 <- plot_osa_comps(obs=rep$res_fish[,1:10], exp=rep$res_fish[,11:20],
               ages=1:10, years=dat$fshyrs, Neff=dat$multN_fsh,
               model='TMB port', plot.type='ggplot')

## I'm really surprised models with such similar expectations
## have such different OSA residuals. Not sure what to think
## about that
r1-r2
re <- (replist$Fishery_observed_and_expected_age_comp[,11:20] - rep$res_fish[,11:20])/rep$res_fish[,11:20]
hist(re)

cowplot::plot_grid(g1,g2)
