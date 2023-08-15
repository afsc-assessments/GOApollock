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
obj <- MakeADFun(data=dat, parameters=pars, map=map, random=NULL,
                 silent=TRUE, DLL='goa_pk_tmb')

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


### Try to get the model to be more stable
control <- list(eval.max=10000, iter.max=10000)
obj <- MakeADFun(data=dat, parameters=pars, map=map, random=NULL,
                 silent=TRUE, DLL='goa_pk_tmb')
lwr <- get_bounds(obj)$lwr
upr <- get_bounds(obj)$upr
opt0 <- with(obj, nlminb(start=par,fn,gr, control=control,
                         lower=lwr, upper=upr))

## Run for some random starting values (with some
## exceptions). Still fails sometimes for unknown reasons.
df <- list()
for(i in 1:10){
  set.seed(i)
  init1 <- obj$par*runif(length(lwr), .8,1.2)
  init1[grep('dev', names(init1))] <- 0 ## zero out the devs
  init1['inf2_fsh_mean'] <- runif(1, 8,9)           ## if >12 ish gradient lost
  init1['inf2_srv1'] <- runif(1, 8,9)
  init1['inf2_srv6'] <- runif(1, 8,9)
  ## opt1 <- with(obj, nlminb(start=init1,fn,gr, control=control,
  ##                          lower=lwr, upper=upr))
  obj$par <- init1
  if(tail(obj$report(obj$par)$Espawnbio,1)<.01){
    warning('Initial pop crash for i=', i)
    next
  }
  opt1 <- TMBhelper::fit_tmb(obj, control=control, lower=lwr, loopnum=3,
                             upper=upr, getsd=FALSE, newtonsteps=0)
  df[[i]] <- data.frame(i=i, maxgrad=opt1$max_gradient, nll=opt1$objective, par=names(obj$par), est0=opt0$par, init1=init1, est1=opt1$par)
  print(opt0$objective-opt1$objective)
}
df <- bind_rows(df)
df <- df %>% mutate(pardiff=est1-est0)
s <- filter(df, !grepl('dev', par))

ggplot(s, aes(i, pardiff, color=nll)) + geom_point() +
  facet_wrap('par')
ggplot(s, aes(init1, pardiff, color=log(abs(maxgrad)))) + geom_point() +
  facet_wrap('par', scales='free_x')



