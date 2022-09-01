## This file tests whether the GOA pollock model can run in WHAM

library(dplyr)
library(tidyr)
library(GOApollock)
## devtools::install_github('timjmiller/wham', ref='master')
library(wham)
library(ggplot2)
theme_set(theme_bw())
library(purrr)
source("functions.R")

## Original ADMB model, final from 2021
arep <- read_pk_rep(path='2021_final', version='ADMB original', endyr=2021)
asdrep <- read_pk_cor(path='2021_final', version='ADMB original', endyr=2021)
## Get WHAM initial values setup to be close by using the output
## from ADMB
asap3 <- read_asap3_dat("goa_pk_asap3.txt")
input <- match_input(arep, asap3)
## map process variances off to start to match ADMB (penalized
## likelihood)
mapoff('sel_repars')                    # selex variance
mapoff('log_NAA_sigma')                 # recruit variance
mapoff('q_repars')                      # q variance

## Compare initial values to ADMB
fit0 <- fit_wham(input, do.osa=FALSE, do.fit=FALSE, do.retro=FALSE,
                do.sdrep=TRUE, MakeADFun.silent=TRUE)
fit0$rep[grep('nll',names(fit0$rep))] %>% lapply(sum) %>% unlist
plot_checks(arep, fit0$rep) # special plot function to explore inits

## Compare platforms after optimizing.
fit1 <- fit_wham(input, do.osa=FALSE, do.fit=TRUE, do.retro=FALSE,
                do.sdrep=TRUE, MakeADFun.silent=TRUE)
fit1$rep[grep('nll',names(fit1$rep))] %>% lapply(sum) %>% unlist
plot_ssb(fit1, asdrep)


