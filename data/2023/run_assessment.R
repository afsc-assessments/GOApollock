## A quick script to reproduce the final 2023 assessment for GOA pollock

## devtools::install_github('afsc-assessments/GOApollock', ref='v0.1.2')
library(GOApollock)
library(TMB)

input <- prepare_pk_input(path=getwd(), datfile='pk23_10.txt', version='2023 final')
str(input$dat)
str(input$par)

fit <- fit_pk(input)

str(fit$opt)
str(fit$rep)
str(fit$sd)

