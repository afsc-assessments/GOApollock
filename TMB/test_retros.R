library(GOApollock)
library(TMB)
devtools::load_all('C:/Users/cole.monnahan/GOApollock/')

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


retros <- lapply(0:7, function(i) fit_peel(obj, i, getsd=0))

get_std <- function(x, version) {
  ss <- x$SD
  with(ss, data.frame(version=version, name=names(value), est=value, se=sqrt(diag(cov)))) %>%
  group_by(name) %>% mutate(year=1969+1:n(), lwr=est-1.96*se, upr=est+1.96*se) %>% ungroup
}


test <- lapply(0:7, function(x) get_std(retros[[x+1]], version=paste0('peel',x)))

ssbs <- lapply(0:7, function(x){
  data.frame(version=paste('peel',x),
             ssb=retros[[x+1]]$rep$Espawnbio) %>%
    mutate(year=1969+1:n())}) %>% bind_rows

ggplot(ssbs, aes(year, ssb, color=version)) + geom_line()

reps <- lapply(retros, \(x) x$rep)
calculate_rho(reps)

mymelt(reps[[1]], 'Espawnbio')
