## This file runs the GOA pollock assessement in WHAM and
## compares it to the ADMB estimates. Because of some fundamental
## differences I was not able to exactly bridge the original
## model. However, given some documented changes in pk_wham.tpl I
## can match the two extremely close given the same inputs (NAA
## and likelihoods for all components). The only change to the
## WHAM cpp file is to manually add a prior on the NMFS BT
## catchability term


##!!! this has the modification of the Q prior that needs to be
## put into the .cpp and then installed to run teh following
## code. YOu will have to do that on your own for now. !!!
## devtools::install('C:/Users/cole.monnahan/wham')

library(dplyr)
library(tidyr)
library(GOApollock)
library(wham)
library(ggplot2)
theme_set(theme_bw())
library(purrr)
source("functions.R")


### --------------------------------------------------
## Phase 1 bridging: modifying original admb model
## Original ADMB model, final from 2021
arep0 <- read_pk_rep(version='ADMB original', endyr=2021)

## Compile the modified ADMB model, pk_wham, and refit so we can
## compare it to the WHAM version
file.copy('pk_wham.tpl', to='goa_pk_wham/pk_wham.tpl', overwrite=TRUE)
setwd("goa_pk_wham/")
system('admb_new pk_wham ')
## Start from MLEs of original model and do not estimate to see
## differences
system("pk_wham -ind pk21_9.txt -maxfn 0 -nohess -ainp goa_pk.par")
arep1 <- read_pk_rep(model='pk_wham', version='pkwham initial', endyr=2021)
## Same thing but optimize it to see how it moves
system("pk_wham -ind pk21_9.txt -ainp goa_pk.par -nohess -iprint 0")
arep2 <- read_pk_rep( model='pk_wham', version='pkwham optimized', endyr=2021)
clean_pk_dir()
setwd("..")

## Check the difference in estimates from modifications to ADMB
g <- mymelt(list(arep0,arep1,arep2), 'Expected_spawning_biomass') %>%
  ggplot(aes(year,value, color=model)) + geom_line(alpha=.5, lwd=1.2) +
 labs(y='SSB (M mt)') + ylim(0,NA) + theme(legend.position='top')
ggsave('plots/phase1_comparisons.png', g, width=7, height=4, dpi=300)

### --------------------------------------------------
### Phase 2 bridging: matching the model and likelihoods
## Get WHAM initial values setup to be close by using the output
## from ADMB.  To compare NLL take two points and subtract,
## Should match for both WHAM and admb Mnaually modified the .par
## file to add 0.5 to recruit mean in log space and some
## arbitrary devs to create difference in age structure which
## should propagate through all data types
rm(list=ls())
source("functions.R")
file.copy('pk_wham.tpl', to='goa_pk_wham_nll/pk_wham.tpl', overwrite=TRUE)
setwd("goa_pk_wham_nll/")
trash <- file.remove(c("pk_wham.rep", "pk_wham2.rep"))
system("admb_new pk_wham")
file.copy('pk_wham.exe', 'pk_wham2.exe', overwrite=TRUE)
system("pk_wham -ind pk21_9.txt -ainp init.pin -maxfn 0 -nohess")
system('pk_wham2 -ind pk21_9.txt -ainp init2.pin -maxfn 0 -nohess')
clean_pk_dir()
setwd("..")
asap3 <- read_asap3_dat("goa_pk_asap3.txt")
arep <- read_pk_rep('pk_wham', 'goa_pk_wham_nll', version='admb', endyr=2021)
arep2 <- read_pk_rep('pk_wham2', 'goa_pk_wham_nll', version='admb2', endyr=2021)

input <- match_input_nll(arep, asap3)
input2 <- match_input_nll(arep2, asap3)
m <- fit_wham(input, do.osa=FALSE, do.fit=FALSE, do.retro=FALSE,
              do.sdrep=FALSE,, MakeADFun.silent=TRUE )
m2 <- fit_wham(input2, do.osa=FALSE, do.fit=FALSE,
               do.retro=FALSE, do.sdrep=FALSE,
               MakeADFun.silent=TRUE)
wrep <- m$report(); wrep2 <- m2$report()
## Check that the two models have nearly identical output
plot_checks(arep, wrep)
plot_checks(arep2, wrep2)
## yep the dynamics match

## Now check that the difference in NLL is the same (accounts for
## diff in constants)
nllnames <- c("total", "total catch", "fsh age comps", "fsh len comps", "surv1 index",
  "surv1 age comps", "surv1 len comps", "surv2 index",
  "surv2 age comps", "surv2 len comps", "unused", "surv3 index",
  "surv3 age comps", "surv3 len comps", "surv4 index",
  "surv5 index", "surv6 index", "surv6 age + len comps",
  "recruit penalties", "TV selex penalties", "proj recruits",
  "TV Q penalties", "unused", "prior trawl Q", "unused")
nllnames <- gsub(" ","_", nllnames)
nllout <- data.frame(name=nllnames,
                     ## admb is LL not NLL
                     admb1=c(arep$Objective_function, -arep$Likelihood_components),
                     admb2=c(arep2$Objective_function, -arep2$Likelihood_components),
                     wham1=get_wham_nll(wrep),
                     wham2=get_wham_nll(wrep2)) %>%
  mutate(admb=round(admb1-admb2,3), wham=round(wham1-wham2,3),
         diff=admb-wham) %>%
    filter(abs(admb2)+abs(wham2)> 1e-1)  #%>% select(name, admb,wham)
nllout %>% print(digits=2)
nllout %>% select(name, diff)
nllout %>% select(name, diff) %>% filter(abs(diff)>0)

### --------------------------------------------------
### Phase 3 bridging: Estimation

## First optimize the pkwham model, in a separate folder from the
## _nll ones above
source("functions.R")
file.copy('pk_wham.tpl', to='goa_pk_wham/pk_wham.tpl', overwrite=TRUE)
setwd("goa_pk_wham/")
trash <- file.remove(c("pk_wham.rep", "pk_wham.cor"))
system("admb_new pk_wham")
system("pk_wham -ind pk21_9.txt -iprint 50")
clean_pk_dir()
setwd("..")
asap3 <- read_asap3_dat("goa_pk_asap3.txt")
arep <- read_pk_rep('pk_wham', 'goa_pk_wham', version='pkwham', endyr=2021)

## Now match WHAM and fit it
asap3 <- read_asap3_dat("goa_pk_asap3.txt")
input <- match_input(arep, asap3)
## map process variances off to start to match ADMB (penalized
## likelihood)
mapoff('sel_repars')                    # selex variance
mapoff('log_NAA_sigma')                 # recruit variance
mapoff('q_repars')                      # q variance

## Compare penalized likelihood models
fit <- fit_wham(input, do.osa=FALSE, do.fit=TRUE, do.retro=FALSE,
                do.sdrep=TRUE, MakeADFun.silent=TRUE)
g0 <- plot_ssb(fit, asdrep, log=FALSE)
ggsave('plots/phase3_bridging_SSB.png', g0, width=7, height=4)
plot_checks(arep, fit$rep)
