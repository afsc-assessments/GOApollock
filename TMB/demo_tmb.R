library(TMBhelper)
library(ggplot2)
library(TMB)
library(compResidual) ## https://github.com/fishfollower/compResidual
## source("functions.R")
## devtools::load_all('C:/Users/cole.monnahan/GOApollock/')
## devtools::install('C:/Users/cole.monnahan/GOApollock/')
library(GOApollock)

## The dat and rep files from the 2022 final model
setwd("Data")
stdadmb <- readRDS("stdfile.RDS")
dat <- readRDS('datfile.RDS')
replist <- readRDS("repfile.RDS")
## pars <- readRDS('pars.RDS')
## map <- readRDS('map.RDS')
## years <- 1970:2022
setwd("..")
input <- prepare_pk_input(dat)


## Fit the model and compare to ADMB. Should be identical
fit <- fit_pk(input)

## optimize and compare
stdtmb <- fit$sd
stds <- bind_rows(cbind(model='TMB',stdtmb),
              cbind(model='ADBM',stdadmb))
## saveRDS(stds, 'admb_tmb_stds.RDS') # used in Sep PT presentation folder
ssb <- bind_rows(cbind(model='TMB',filter(stdtmb, name=='Espawnbio')),
             cbind(model='ADMB',filter(stdadmb, name=='Espawnbio'))) %>%
  select(model, par=name, est,se, year)
filter(ssb, year==2015) %>% arrange(year)
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

r2 <- plot_osa_comps(obs=rep$res_fish[,1:10], exp=rep$res_fish[,11:20],
               ages=1:10, years=dat$fshyrs, Neff=dat$multN_fsh,
               model='TMB port', plot.type='default')
g2 <- plot_osa_comps(obs=rep$res_fish[,1:10], exp=rep$res_fish[,11:20],
               ages=1:10, years=dat$fshyrs, Neff=dat$multN_fsh,
               model='TMB port', plot.type='ggplot')
cowplot::plot_grid(g1,g2)

