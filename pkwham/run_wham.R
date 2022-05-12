## This file runs the GOA pollock assessement in WHAM under some
## different configurations. Use the devel branch of wham not the
## goapk_bridge one like with the bridging exercise.

## an original version of this script assumed that the
## 'goapk_bridge' version of wham was installed. It's probably
## best to break that now so these results won't match previous
## ones. Specifically, the wham models won't match the admb as
## well.

devtools::install_github('timjmiller/wham', ref='devel')
devtools::install_github('Cole-Monnahan-NOAA/wham', ref='goapk_bridge')
library(dplyr)
library(tidyr)
library(GOApollock)
library(wham)
packageVersion('wham') #  '1.0.6.9000'
library(ggplot2)
theme_set(theme_bw())
library(purrr)
source("functions.R")
mapoff <- function(name) input$map[[name]] <<- as.factor(input$par[[name]]*NA)


## Get WHAM initial values setup to be close by using the output
## from ADMB
arep <- read_pk_rep('pk_wham', 'goa_pk_wham', version='pkwham', endyr=2021)
asdrep <- read_pk_cor('pk_wham', 'goa_pk_wham', version='pkwham', endyr=2021)
asap3 <- read_asap3_dat("goa_pk_asap3.txt")
input <- match_input(arep, asap3)
saveRDS(input, 'akwham_input_2021.RDS')
## map process variances off to start to match ADMB (penalized
## likelihood)
mapoff('sel_repars')                    # selex variance
mapoff('log_NAA_sigma')                 # recruit variance
mapoff('q_repars')                      # q variance

## SSB and uncertainty for the two for base case model
fit0 <- fit_wham(input, do.osa=FALSE, do.fit=TRUE, do.retro=FALSE,
                do.sdrep=TRUE, MakeADFun.silent=TRUE)
g0 <- plot_ssb(fit0, asdrep)
ggsave('plots/SSB0.png', g0, width=5, height=3.5)
(tab0 <- par_table(fit0))
plot_checks(arep, fit0$rep)

## Try turning on estimation of the recdevs as RE
input$map$log_NAA_sigma <- factor(1)    # variance estimated
input$random <- c("log_NAA")            # recdevs are RE
fit1 <- fit_wham(input, do.osa=FALSE, do.fit=TRUE, do.retro=FALSE,
                do.sdrep=TRUE, MakeADFun.silent=TRUE)
g1 <- plot_ssb(fit1, asdrep)
ggsave('plots/SSB1.png', g1, width=5, height=3.5)
(tab1 <- par_table(fit1))
est1 <- data.frame(par=names(fit1$sdrep$value),
                   est=fit1$sdrep$value,
                   sd=fit1$sdrep$sd)
sig <- filter(est1, par=='NAA_sigma')
x <- seq(.8,1.75, len=1000);y=dnorm(x,sig$est, sig$sd)
png('plots/sigmaR.png', width=3, height=2, units='in', res=200)
par(mar=c(2.5,2,.5,.5), mgp=c(1.5,.5,0), tck=-.02)
plot(x,y, type='l', xlab='SigmaR', ylab=''); abline(v=1, lty=3)
dev.off()

## Try turning on estimation of the catchability variances
tmp <- as.vector(input$par$q_repars)*NA
tmp[c(1,3)] <- 1:2
input$map$q_repars <- factor(tmp)
input$random <- c(input$random,"q_re")
fit2 <- fit_wham(input, do.osa=FALSE, do.fit=TRUE, do.retro=FALSE,
                do.sdrep=TRUE, MakeADFun.silent=TRUE)
g2 <- plot_ssb(fit2, asdrep)
(tab2 <- par_table(fit2))
ggsave('plots/SSB2.png', g2, width=5, height=3.5)
est2 <- summary(fit2$sdrep) %>% data.frame
sigs <- est2[grepl(x=rownames(est2), 'q_repars'),]
x <- seq(-3.3,0, len=1000);
y1 <- dnorm(x,sigs[1,1], sigs[1,2])
y2 <- dnorm(x,sigs[2,1], sigs[2,2])
png('plots/log_q_sd.png', width=3, height=2, units='in', res=200)
par(mar=c(2.5,2,.5,.5), mgp=c(1.5,.5,0), tck=-.02)
plot(0,0, xlim=range(x), ylim=range(c(y1,y2)), type='n',
     xlab='log q SD', ylab=''); abline(v=1, lty=3)
lines(x,y1); abline(v=-3.27)
lines(x,y2, lty=3); abline(v=-2.9957, lty=3)
dev.off()
## Plot CI of actual q. First get it in logit space
est2 <- data.frame(par=names(fit2$sdrep$value),
                   est=fit2$sdrep$value,
                   sd=fit2$sdrep$sd) %>% filter(par=='logit_q_mat')
qfn <- function(x,a=0,b=1000) a+(b-a)/(1+exp(-x))
q1 <- data.frame(est=matrix(est2$est, ncol=6)[,1],
                 sd=matrix(est2$sd, ncol=6)[,1]) %>%
  mutate(lwr=est-1.96*sd, upr=est+1.96*sd) %>%
  mutate(est=qfn(est), lwr=qfn(lwr), upr=qfn(upr), q='q1', platform='wham')
q3 <- data.frame(est=matrix(est2$est, ncol=6)[,3],
                 sd=matrix(est2$sd, ncol=6)[,3]) %>%
  mutate(lwr=est-1.96*sd, upr=est+1.96*sd) %>%
  mutate(est=qfn(est), lwr=qfn(lwr), upr=qfn(upr), q='q3', platform='wham')
wq <- bind_rows(q1,q3) %>% cbind(year=1970:2021)
aq <- asdrep %>% filter(name %in% c('log_q1', 'log_q3')) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se, q=gsub('log_','', name)) %>%
  mutate(platform='admb', est=exp(est), lwr=exp(lwr), upr=exp(upr)) %>%
  select(est, lwr, upr, platform, year,q)
q <- bind_rows(wq,aq) %>%
  select(est, lwr, upr, platform, year,q) %>%
  filter(year>1990)
gq <- ggplot(q, aes(year, est, ymin=lwr, ymax=upr, fill=platform, color=platform))+
  geom_ribbon(alpha=.5) + geom_line() +
  facet_wrap('q', scales='free_y') +
  theme(legend.position='top') +
  geom_hline(yintercept=1) + labs(y='Catchability', x=NULL)
ggsave('plots/q_comparison.png', gq, width=5, height=3.5)



### Try turning on estimation of the selectivity variances. This
### didn't work when also estimating recdevs and catchabilities
input <- match_input(arep, asap3)
## map process variances off to start to match ADMB (penalized
## likelihood)
mapoff('sel_repars')                    # selex variance
mapoff('log_NAA_sigma')                 # recruit variance
mapoff('q_repars')                      # q variance
## Turn on RE estimation of selex
## tmp <- as.vector(input$par$sel_repars)*NA
## tmp[1:2] <- 1
## input$map$sel_repars <- factor(tmp)
input$map$sel_repars ## for the hypervariance terms
input$map$selpars_re ## map for the RE vectors
input$par$selpars_re
input$par$sel_repars
#input$random <- c("selpars_re")
rm(fit3)
fit3 <- fit_wham(input, do.osa=FALSE, do.fit=TRUE, do.retro=FALSE,
                 do.sdrep=FALSE, MakeADFun.silent=FALSE)
selfsh <- fit3$rep$selAA[[1]]
matplot(t(selfsh), type='l')
fit3$par
table(names(fit3$par))
par_table(fit3)
fit3$par
g3 <- plot_ssb(fit3, asdrep)
g3
plot_checks(arep, fit3$rep)


### Explore full state space capabilities, holding other RE
### estimation off.
input <- match_input(arep, asap3,
                     NAA_re=list(sigma="rec+1", cor='iid'))
mapoff('sel_repars')                    # selex variance
mapoff('q_repars')                      # q variance
input$random <- c("log_NAA")            # recdevs are RE
fit4 <- fit_wham(input, do.osa=FALSE, do.fit=TRUE, do.retro=FALSE,
                do.sdrep=TRUE, MakeADFun.silent=TRUE)
plot_checks(arep, fit4$rep)
g4 <- plot_ssb(fit4, asdrep)
ggsave('plots/SSB4.png', g4, width=5, height=3.5)
(tab4 <- par_table(fit4))
est4 <- summary(fit4$sdrep) %>% data.frame
sigs <- est4[grepl(x=rownames(est4), 'log_NAA_sigma'),]
sigs %>% round(2)


### Explore state space capabilities, inflated Neff to see if
### that affects process error estimation
input <- match_input(arep, asap3, NAA_re=list(sigma="rec+1", cor='iid'))
input$data$index_Neff <- 10*input$data$index_Neff
mapoff('sel_repars')                    # selex variance
mapoff('q_repars')                      # q variance
input$random <- c("log_NAA")            # recdevs are RE
fit5 <- fit_wham(input, do.osa=FALSE, do.fit=TRUE, do.retro=FALSE,
                do.sdrep=TRUE, MakeADFun.silent=TRUE)
plot_checks(arep, fit5$rep)
g5 <- plot_ssb(fit5, asdrep)
ggsave('plots/SSB5.png', g5, width=5, height=3.5)
(tab5 <- par_table(fit5))
est5 <- summary(fit5$sdrep) %>% data.frame
est5 <- est5[grepl(x=rownames(est5), 'log_NAA_sigma|rho'),]
sigs %>% round(2)


### Explore state space capabilities
input <- match_input(arep, asap3, NAA_re=list(sigma="rec+1", cor='ar1_y'))
input$data$index_Neff <- 10*input$data$index_Neff
mapoff('sel_repars')                    # selex variance
mapoff('q_repars')                      # q variance
input$random <- c("log_NAA")            # recdevs are RE
fit6 <- fit_wham(input, do.osa=FALSE, do.fit=TRUE, do.retro=FALSE,
                do.sdrep=TRUE, MakeADFun.silent=TRUE)
plot_checks(arep, fit6$rep)
g6 <- plot_ssb(fit6, asdrep)
ggsave('plots/SSB6.png', g6, width=5, height=3.5)
(tab6 <- par_table(fit6))
est6 <- summary(fit6$sdrep) %>% data.frame
est6 <- est6[grepl(x=rownames(est6), 'log_NAA_sigma|rho'),]


### Explore state space capabilities
input <- match_input(arep, asap3, NAA_re=list(sigma="rec+1", cor='ar1_a'))
input$data$index_Neff <- 10*input$data$index_Neff
mapoff('sel_repars')                    # selex variance
mapoff('q_repars')                      # q variance
input$random <- c("log_NAA")            # recdevs are RE
fit7 <- fit_wham(input, do.osa=FALSE, do.fit=TRUE, do.retro=FALSE,
                do.sdrep=TRUE, MakeADFun.silent=TRUE)
plot_checks(arep, fit7$rep)
g7 <- plot_ssb(fit7, asdrep)
ggsave('plots/SSB7.png', g7, width=5, height=3.5)
(tab7 <- par_table(fit7))
est7 <- summary(fit7$sdrep) %>% data.frame
est7 <- est7[grepl(x=rownames(est7), 'log_NAA_sigma|rho'),]

### Explore state space capabilities
input <- match_input(arep, asap3, NAA_re=list(sigma="rec+1", cor='2dar1'))
input$data$index_Neff <- 10*input$data$index_Neff
mapoff('sel_repars')                    # selex variance
mapoff('q_repars')                      # q variance
input$random <- c("log_NAA")            # recdevs are RE
fit8 <- fit_wham(input, do.osa=FALSE, do.fit=TRUE, do.retro=FALSE,
                do.sdrep=TRUE, MakeADFun.silent=TRUE)
plot_checks(arep, fit8$rep)
g8 <- plot_ssb(fit8, asdrep)
ggsave('plots/SSB8.png', g8, width=5, height=3.5)
(tab8 <- par_table(fit8))
est8 <- summary(fit8$sdrep) %>% data.frame
est8 <- est8[grepl(x=rownames(est8), 'log_NAA_sigma|rho'),]



## Compare different state space options
ssb0 <- get_wham_ssb(fit0, 'penalized ML')
ssb1 <- get_wham_ssb(fit1, 'rec iid')
ssb4 <- get_wham_ssb(fit4, 'iid')
ssb5 <- get_wham_ssb(fit5, 'iid (10*Neff)')
ssb6 <- get_wham_ssb(fit6, 'ar1_y (10*Neff)')
ssb7 <- get_wham_ssb(fit7, 'ar1_a (10*Neff)')
ssb8 <- get_wham_ssb(fit8, '2dar1 (10*Neff)')
ssb.all <- bind_rows(ssb0,ssb1, ssb4, ssb5, ssb6, ssb7)
ssb.all <- bind_rows(ssb5, ssb6, ssb7, ssb8)
g <- ggplot(ssb.all, aes(year, est, ymin=est-1.96*sd, ymax=est+1.96*sd,
                     fill=version, color=version)) +
  geom_ribbon(alpha=.25) + geom_line(lwd=1) +
 ## theme(legend.position='top')+
  labs(y='log SSB',x=NULL, color=NULL, fill=NULL)
g
ggsave('plots/full_ss.png', g, width=7, height=3)

ssb.all <- bind_rows(ssb4, ssb5)
g <- ggplot(ssb.all, aes(year, est, ymin=est-1.96*sd, ymax=est+1.96*sd,
                     fill=version, color=version)) +
  geom_ribbon(alpha=.5) + geom_line(lwd=1) +
  theme(legend.position='top')+
  labs(y='log SSB',x=NULL, color=NULL, fill=NULL)
g
ggsave('plots/full_ss_iid.png', g, width=5, height=3.5)



## ## Look at NAA. Need to fix this. What output to use?
## get_naa <- function(x, v=NULL){
##   x <- x$rep$pred_NAA %>% reshape2::melt()
##   names(x) <- c('year', 'age', 'NAA')
##   x <- x %>% mutate(year=year+1970)
##   if(!is.null(v)) x <- cbind(x, version=v)
##   x
## }
## naa.all <- rbind(get_naa(fit4, 'iid'),
##                  get_naa(fit5, 'iid (10*Neff)'),
##                  get_naa(fit6, 'ar1_y (10*Neff)'),
##                  get_naa(fit7, 'ar1_a (10*Neff)'),
##                  get_naa(fit8, '2dar1 (10*Neff)'))
## ggplot(naa.all, aes(age,year,fill=log(NAA))) + geom_raster() +
##   scale_y_reverse() + facet_wrap('version')
