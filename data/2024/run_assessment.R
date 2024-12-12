## A quick script to reproduce the final 2023 assessment for GOA pollock

## devtools::install_github('afsc-assessments/GOApollock', ref='v0.1.3')
library(GOApollock)

## assume working directory is data/2024

input <- prepare_pk_input(path=getwd(), datfile='pk24_12.txt',
                            version='23d: 2024 final',
                            complike = 'D-M')
# age 1 and 2 Shelikof indices turned off so don't estimate catchabilities
input$map$log_q4 <- input$map$log_q5 <- factor(NA)
# bump up so log_DM pars doesn't hit bounds. Needs to be moved into dat file later
input$dat$multN_srv1 <- input$dat$multN_srv1*2
input$dat$multN_srv3 <- input$dat$multN_srv3*2
input$dat$multN_srv6 <- input$dat$multN_srv6*2

str(input$dat)
str(input$par)

# fit the model with defaults
fit <- fit_pk(input)

str(fit$opt)
str(fit$rep)
str(fit$sd)

plot_ssb(fit)
get_std(fit, 'log_recruit') %>%
  ggplot(aes(year, est, ymin=lwr, ymax=upr)) + geom_pointrange() +
  labs(y='log recruits (billions)')
