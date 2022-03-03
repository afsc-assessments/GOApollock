library(GOApollock)
library(wham)
## To compare NLL take two points and subtract, Should match for
## both WHAM and admb
## Mnaually modified the .par file to add 0.5 to recruit mean in
## log space and some arbitrary devs to create difference in age
## structure which should propagate through all data types
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
  mutate(admb=round(admb1-admb2,2), wham=round(wham1-wham2,2),
         diff=admb-wham) %>%
    filter(abs(admb1)+abs(wham1)> 1e-1)  #%>% select(name, admb,wham)
nllout
nllout %>% select(name, diff)
nllout %>% select(name, diff) %>% filter(abs(diff)>0)

