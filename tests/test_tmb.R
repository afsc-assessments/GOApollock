
## devtools::install()
library(GOApollock)
setwd('tests')

## dyn.unload('model_23/goa_pk_tmb.dll')
## file.copy('../source/goa_pk_tmb.cpp', to='model_23/goa_pk_tmb.cpp', overwrite = TRUE)

## make a few changes to be able to explore better
input <- prepare_pk_input('model_23', datfile='goa_pk.dat')
tmb <- fit_pk(input, do.fit=TRUE, filename=NULL)
input1 <- input
input1$dat$Ftarget <- c(input1$dat$Ftarget,0*input1$dat$Ftarget)
input1$version <- 'estimate sigmaR'
input1$map$sigmaR <- factor(1)
input1$random <- 'dev_log_recruit'
tmb1 <- fit_pk(input1, do.fit=TRUE, filename=NULL)

z <- list(tmb,tmb1)
stopifnot(is.pkfits(z))
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=TRUE, uselog=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=TRUE, uselog=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE, uselog=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=FALSE, uselog=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=TRUE, uselog=FALSE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=TRUE, uselog=FALSE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE, uselog=FALSE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=FALSE, uselog=FALSE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=TRUE, uselog=TRUE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=TRUE, uselog=TRUE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE, uselog=TRUE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=FALSE, uselog=TRUE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=TRUE, uselog=FALSE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=TRUE, uselog=FALSE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE, uselog=FALSE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=FALSE, uselog=FALSE, addproj=TRUE)

ssb <- plot_pk_ssb(z, add_uncertainty=FALSE, plotlog=FALSE,
            uselog=FALSE, plot=FALSE, addproj=TRUE)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE,
            uselog=FALSE, alpha1=.3)
plot_pk_ssb(z, add_uncertainty=TRUE, plotlog=FALSE,
            uselog=FALSE, alpha2=.3)


x <- get_rep(z)
(test <- get_rep(z, 'q5'))
(test <- get_rep(z, 'ssb'))
(test <- get_rep(z, 'Espawnbio'))
(test <- get_rep(z, 'Espawnbio_proj'))
(test <- get_rep(z, 'slctfsh'))
(test <- get_rep(z))


## restrospectives
fits <- fit_pk_retros(fit=tmb)
fits <- fit_pk_retros(fit=tmb, peels=0:6, parallel=TRUE)
plot_pk_ssb(fits)

plot_pk_ssb(fits[[1]])

jitters <- run_jitter(tmb, njitter=3, parallel=FALSE)
jitters <- run_jitter(tmb, njitter=10, parallel=TRUE)
jitters <- run_jitter(tmb, njitter=10, parallel=TRUE, scalar=.5, type='mle')

