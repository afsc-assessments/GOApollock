### need to make sure the TMB model retro approach can replicate
### the old ADMB model

library(GOApollock)
library(TMB)
devtools::load_all('C:/Users/cole.monnahan/GOApollock/')

## The dat and rep files from the 2022 final model
setwd("Data")
stdadmb <- readRDS("stdfile.RDS")
dat <- readRDS('datfile.RDS')
replist <- readRDS("repfile.RDS")
setwd("..")
input <- prepare_pk_input(dat)
## Fit the full model
fit <- fit_pk(input, getsd=FALSE)
obj <- fit$obj

retros <- fit_pk_retros(obj=obj)



stds <- get_std(retros)
ssbtmb <- filter(stds, name=='Espawnbio')
ssbadmb <-
  readRDS('C:/Users/cole.monnahan/Work/assessments/GOA_pollock/2022/results/retro_stds.RDS') %>%
  filter(name=='Espawnbio')


ssb <- bind_rows(cbind(model='TMB', ssbtmb),
                 cbind(model='ADMB',ssbadmb)) %>%
  mutate(peel=as.numeric(gsub('peel','', version))) %>%
  filter(peel<=7) %>% arrange(peel, year) %>%
  select(model, est, se, year, peel)
library(tidyr)
estdiff <- ssb %>% select(-se) %>%
  pivot_wider(names_from='model', values_from='est')%>%
  mutate(estdiff=(TMB-ADMB)/ADMB)
sediff <- ssb %>% select(-est) %>%
  pivot_wider(names_from='model', values_from='se')%>%
  mutate(sediff=(TMB-ADMB)/ADMB)

## yep matches
sediff$sediff %>% abs %>% max
estdiff$estdiff %>% abs %>% max

diffs <- cbind(estdiff, sediff=sediff$sediff)
ggplot(diffs, aes(year, estdiff, color=factor(peel))) +
  geom_line()
ggplot(diffs, aes(year, sediff, color=factor(peel))) +
  geom_line()
