library(tidyverse)
theme_set(theme_bw())
setwd(here::here())
source("TMB/R/Functions/osa_comp_plots.R")
source("TMB/R/Functions/construct_Q.R")
source("TMB/R/Functions/plot_selectivity.R")
fits <- readRDS("TMB/Output/fits.RDS")
years <- 1970:(2022+5)
rho_trans <- function(x)2/(1+exp(-2*x))-1
ilogit <- function(x) 1/(1+exp(-x))
## Quick plot of selex
plot_selex_persp <- function(fit){
  modelname <- model_names[fit$modelnum+1]
  years <- 1970:(2022+5)
  ## df <- fit$std %>% filter(par=='slctfsh') %>%
  ##   mutate(age=rep(1:10, each=length(years)),
  ##          year=rep(years, times=10), model=model)
  sel <- t(fit$report$slctfsh)
  ## ind <- 1:ncol(sel)
  ## sel <- sel[,rev(ind)]
  ## years <- rev(years)
  persp(y = years, x =  1:10,
        z = sel,
        border=gray(.5),
        col = "white",
        ## xlab = "Age", ylab = "\n\nYear",
        xlab=NA, ylab=NA, zlab=NA,
        ## zlab = "\n\nSelectivity",
        expand = 0.35,
        box = TRUE, ticktype = "detailed", phi = 35,
        theta = -19, main = modelname)
}


## run time (mins)
lapply(fits, function(x) as.numeric(x$time_for_run)/60)
### for making tables
## ParDevs ests
m1 <- with(fits$fit1$SD,
           data.frame(
             par=names(par.fixed),
             est=par.fixed,
             se=sqrt(diag(cov.fixed)))) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se) %>%
  filter(par %in% c("sel_rho_a", "sel_rho_y", "ln_sel_sd"))
m1
exp(-3.0695260599)
exp(-3.0695260599+ 1.96*c(-1,1)*1.702977e-01) %>% round(2)
## get AR1 par estimates out
m7 <- with(fits$fit7$SD,
           data.frame(
             par=names(par.fixed),
             est=par.fixed,
             se=sqrt(diag(cov.fixed)))) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se) %>%
  filter(par %in% c("sel_rho_a", "sel_rho_y", "ln_sel_sd")) %>%
  select(-se)
m7[1:2,-1] <- rho_trans(m7[1:2,-1])
m7[3,-1] <- exp(m7[3,-1])
m7 %>% print(digits=2)
m8 <- with(fits$fit8$SD,
           data.frame(
             par=names(par.fixed),
             est=par.fixed,
             se=sqrt(diag(cov.fixed)))) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se) %>%
  filter(par %in% c("sel_rho_a", "sel_rho_c", "sel_rho_y",
           "ln_sel_sd")) %>% select(-se)
m8[4,-1] <- exp(m8[4,-1])
m8 %>% print(digits=2)


partable <- lapply(fits, \(x) c(x$modelnum, x$number_of_coefficients)) %>%
bind_rows %>% t %>% data.frame %>%
  setNames(c('model', 'total', 'fixed', 'random')) %>%
  mutate(fixed_selex=fixed-115)
partable

## Combine objects ----
model_names <- c(
  'TMB-Mod0: Constant',
  'TMB-Mod1: RE',
  'TMB-Mod2: Log-AR1-Age',
  'TMB-Mod3: Log-AR1-Yr',
  'TMB-Mod4: Log-2D-AR1',
  'TMB-Mod5: Age-specific',
  'TMB-Mod6: AR1-Yr',
  'TMB-Mod7: 2D-AR1',
  "TMB-Mod8: 3D-AR1cond",
  'TMB-Mod9: 3D-AR1mar'
)
model_names <- c(
'Constant',
'ParDevs',
'Log-AR1-Age',
'Log-AR1-Yr',
'Log-2D-AR1',
'AR1-Age',
'AR1-Yr',
'2D-AR1',
'3D-AR1cond',
'3D-AR1mar'
)

## didnt' converge so have to do some hacks to get code below to
## work
fits$fit9$objective <- fits$fit9$opt$objective
fits$fit9$par <- fits$fit9$opt$par
fits$fit9$max_gradient <- fits$fit9$opt$max_gradient
fits$fit9$iterations <- fits$fit9$opt$iterations
fits$fit9 <- c(fits$fit9, fits$fit9$opt)

aic_table <- data.frame(
  model = model_names,
  total_nll = sapply(fits, function(x) round(x$objective,1)),
  fish_age_nll= sapply(fits, function(x) round(-x$report$loglik[2],1)),
  k = sapply(fits, function(x) length(x$par)),
  maxgrad=round(sapply(fits, function(x) x$max_gradient),2),
  #convergence = sapply(opt_list, function(x) x$message),
  iterations = sapply(fits, function(x) x$iterations),
  AIC = sapply(fits, function(x) round(x$AIC,1))
) %>%
  mutate(dAIC = round(AIC -  min(AIC),1))
write.csv(aic_table, file = "TMB/Output/Selectivity_runs_AIC_full_estimation2.csv")



# - ssb
ssb <- rbind(cbind(model=model_names[1],filter(fits$fit0$std, par=='Espawnbio')),
             cbind(model=model_names[2],filter(fits$fit1$std, par=='Espawnbio')),
             #cbind(model=model_names[3],filter(fits$fit3$std, par=='Espawnbio')),
             #cbind(model=model_names[4],filter(fits$fit4$std, par=='Espawnbio')),
             #cbind(model=model_names[5],filter(fits$fit5$std, par=='Espawnbio')),
             #cbind(model=model_names[6],filter(fits$fit6$std, par=='Espawnbio')),
             ## cbind(model=model_names[7],filter(fits$fit7$std, par=='Espawnbio')),
             cbind(model=model_names[8],filter(fits$fit7$std, par=='Espawnbio')),
              cbind(model=model_names[9],filter(fits$fit8$std, par=='Espawnbio')),
             cbind(model='ADMB',filter(stdadmb, name=='Espawnbio') %>%
                                select(par=name, est,se,
                                       year)))%>%
  mutate(CV=se/est) %>% pivot_longer(c('CV','est')) %>%
  mutate(name=factor(name, levels=c('est', 'CV'),
                        labels=c('SSB', 'CV(SSB)')))
g <- ggplot(ssb, aes(year, value, color=model)) + geom_line()  +
  ylim(0,NA) + facet_grid(name~.) + labs(x=NULL, y=NULL, color='Fishery\nSelectivity')
ggsave("TMB/Output/ssb_compare.png", g, width=7, height=5,dpi=300)
ggsave("TMB/Output/ssb_compare2.png", g, width=5, height=3,dpi=300)


png('TMB/Output/selex_persp.png', width=7, height=5, units='in', res=300)
par(mfrow=c(2,2), mar = c(.5, .5 , 1 , 0.5) ,
    oma = c(0 , 0 , 0 , 0), tcl = -0.8, mgp = c(10, 0.6, 0))
trash <- lapply(fits[c(2,8:10)], plot_selex_persp)
dev.off()
png('TMB/Output/selex_persp2.png', width=7, height=4, units='in', res=300)
par(mfrow=c(2,2), mar = c(.5, .5 , 1 , 0.5) ,
    oma = c(0 , 0 , 0 , 0), tcl = -0.8, mgp = c(10, 0.6, 0))
trash <- lapply(fits[c(2,8:10)], plot_selex_persp)
dev.off()



get_selexF <- function(x, model){
  years <- 1970:2022
 ## df <- data.frame(name=names(x$value), est=x$value, se=sqrt(diag(x$cov)))
  df <- filter(df, name=='slctfsh_times_F') %>%
    mutate(age=rep(1:10, each=length(years)),
           year=rep(years, times=10), model=model)
  df
}

get_fsh_selex <- function(x, model){
  ## df <- data.frame(name=names(x$value), est=x$value, se=sqrt(diag(x$cov)))
  df <- filter(x, par=='slctfsh') %>%
    mutate(age=rep(1:10, each=length(years)),
           year=rep(years, times=10), model=model)
  df
}
out <- lapply(1:length(fits[-9]), function(x)
  get_fsh_selex(fits[[x]]$std, model_names[x])) %>%
  bind_rows %>%
  filter(model %in% model_names[c(2,8,9)]) %>%
  mutate(lwr=pmax(0, est-1*se), upr=pmin(1, est+1*se))
  ## mutate(lwr=est-1*se, upr=est+1*se)
  ##  mutate(est=ilogit(est), upr=ilogit(upr), lwr=ilogit(lwr))

group_by(out, model) %>% summarize(maxsel=max(est), maxage=age[which.max(est)])


g <- ggplot(out, aes(year, est, ymin=lwr, ymax=upr, fill=model, color=model)) + facet_wrap('age', scales='fixed') +
  geom_line() + geom_ribbon(alpha=.3) +
   geom_vline(xintercept=2021) + theme(legend.position='top') +
  ylim(0,1) + labs(x=NULL,y='Selectivity')
ggsave("TMB/Output/selex_by_age.png", g, width=9, height=6, dpi=300)
g2 <- g + theme(legend.position=c(.7,.15)) +
  labs(fill=NULL, color=NULL) +
  scale_x_continuous(breaks=seq(1970,2022, by=20))
g2
ggsave("TMB/Output/selex_by_age2.png", g2, width=7, height=4, dpi=300)

g <- ggplot(filter(out, age %in% c(2:5, 9:10)), aes(year, est, ymin=lwr, ymax=upr, fill=model, color=model)) + facet_wrap('age', scales='fixed') +
  geom_line() + geom_ribbon(alpha=.3) +
   geom_vline(xintercept=2021) + theme(legend.position='top') +
  ylim(0,1) + labs(x=NULL,y='Selectivity')  +
  scale_x_continuous(breaks=c(2020,2022,2024), limits=c(2018,2025))
ggsave("TMB/Output/selex_by_age_proj.png", g, width=8, height=5, dpi=300)
ggsave("TMB/Output/selex_by_age_proj2.png", g, width=7, height=3.5, dpi=300)
g$data %>% filter(age==4 & year ==2023)

g <- ggplot(filter(out, year>2021 & year < 2024), aes(age, est,
  ymin=lwr, ymax=upr, fill=model, color=model)) +
  facet_wrap('year', scales='free') +
  geom_line() + geom_ribbon(alpha=.2) + ylim(0,1) +
  theme(legend.position='top') + labs(y='Selectivity')
ggsave("TMB/Output/selex_by_year_proj.png", g, width=7, height=3.5, dpi=300)

xx <- filter(out, year>2022 & year < 2024) %>%
  select(est, age, model)
## add on average of last 5 years from ParDevs
yy <- filter(out, year <2022 & year > 2016) %>%
  group_by(age ,model) %>% summarize(est=mean(est), .groups='drop') %>%
  filter(model=='ParDevs') %>%  select(est, age, model) %>%
  mutate(model='ParDevs: 5 year avg')
xx <- bind_rows(xx,yy)

g <- ggplot(xx, aes(age, est, color=model))  +
  geom_line() +  ylim(0,1) +
  theme(legend.position='top') + labs(color=NULL, x=NULL, y='Selectivity')
ggsave("TMB/Output/selex_by_year_proj2.png", g, width=5, height=3, dpi=300)


## x <- 1:20
## y <- rnorm(20, mean=1:20)
## z <- numeric(5)
## for(i in 1:10) {
##   ## print(c(i, i+10-5, i+10-1, mean(y[(i+10-5):(i+10-1)])))
##   z[i] <- mean(y[(i+10-5):(i+10-1)])
## }
## plot(x,y)
## lines(10+1:10, z)

# * Plot OSA residuals ----
## quantities_list <- list(quantities_mod1, quantities_mod2, quantities_mod3, quantities_mod4, quantities_mod5, quantities_mod6, quantities_mod7, quantities_mod8, quantities_mod9)
r_osa_list <- list()
g_osa_list <- list()
# - ADMB
r_osa_list[[1]] <- plot_osa_comps(obs=replist$Fishery_observed_and_expected_age_comp[,1:10],
                                  exp=replist$Fishery_observed_and_expected_age_comp[,11:20],
                                  ages=1:10, years=dat$fshyrs, Neff=dat$multN_fsh,
                                  model='ADMB final 2022')
## - TMB
for(i in 1:length(fits)){
  x <- fits[[i]]$report
  r_osa_list[[i]] <- try(plot_osa_comps(obs=x$res_fish[,1:10],
                                          exp=x$res_fish[,11:20],
                                          ages=1:10,
                                          years=dat$fshyrs,
                                          Neff=dat$multN_fsh,
                                          model=model_names[i],
                                          plot.type='default'), silent = TRUE)
  g_osa_list[[i]] <- try(plot_osa_comps(obs=x$res_fish[,1:10],
                                          exp=x$res_fish[,11:20],
                                          ages=1:10,
                                          years=dat$fshyrs,
                                          Neff=dat$multN_fsh,
                                          model=model_names[i],
                                          plot.type='ggplot'), silent = TRUE)
}

sapply(r_osa_list, class) # 1 doesnt work
sapply(g_osa_list, class) # 1 and 5 doesnt work

# Only show models that make sense (2D ar1 and 3d ar)
## cowplot::plot_grid(g_osa_list[[1]], g_osa_list[[2]], g_osa_list[[8]], g_osa_list[[9]], g_osa_list[[10]])

tmp <- list(g_osa_list[[1]], g_osa_list[[2]], g_osa_list[[8]], g_osa_list[[9]])
xx <- lapply(tmp, \(x) cbind(model=x$labels$title, x$data)) %>% bind_rows
xx$model <- factor(xx$model, levels=model_names[c(1,2,8,9)])
g <- ggplot(xx, aes(year,age, size=abs(resid), color=resid>0, shape=abs(resid)>3)) +
  geom_point(alpha=.5) + facet_wrap('model') + scale_size(range=c(0,6))
ggsave('TMB/Output/osa_resids.png', g, width=9, height=5, dpi=500)
g2 <- g+labs(x=NULL)
g2
ggsave('TMB/Output/osa_resids2.png', g2, width=7, height=3.5, dpi=500)

## rebuild admb one to better match
g <- plot_osa_comps(obs=replist$Fishery_observed_and_expected_age_comp[,1:10],
                                  exp=replist$Fishery_observed_and_expected_age_comp[,11:20],
                                  ages=1:10, years=dat$fshyrs, Neff=dat$multN_fsh,
                                  model='ADMB final 2022', plot.type='ggplot')
g <- ggplot(g$data, aes(year,age, size=abs(resid), color=resid>0, shape=abs(resid)>3)) +
  geom_point(alpha=.3) + scale_size(range=c(0,6)) + ggtitle("2022 Final (ADMB)")
ggsave('TMB/Output/admb_osa.png', g, width=7, height=4, dpi=300)


xx <- lapply(c(1,2,8,9), \(x) data.frame(model=model_names[x],
                                         resid=as.numeric(r_osa_list[[x]]))) %>%
   bind_rows()
g <- ggplot(xx, aes(sample=resid, color=model)) + stat_qq() + geom_abline(intercept=0, slope=1)
ggsave('TMB/Output/qqnorm_resids.png', g, width=7, height=3.5, dpi=500)


# * Plot cor ----
mod = 9
plot_cor(n_years = 3, n_ages = dat$nages,
         rho_y = par_list[[mod]]$sel_rho_y,
         rho_a = par_list[[mod]]$sel_rho_a,
         rho_c = par_list[[mod]]$sel_rho_c, log_sigma = par_list[[mod]]$ln_sel_sd)





## ## Catch is almost identical but F varies, presumably due to
## differneces in max sel??
## ecat <- lapply(1:length(fits), function(x)
##   data.frame(model=model_names[x], year=1970:2022,
##              catch=fits[[x]]$report$Ecattot)) %>%
##   bind_rows
## ggplot(ecat, aes(year, catch, color=model)) + geom_line()
## fishF <- lapply(1:length(fits), function(x)
##   data.frame(model=model_names[x], year=1970:2022,
##              catch=fits[[x]]$report$F)) %>%
##   bind_rows
## ggplot(fishF, aes(year, catch, color=model)) + geom_line()
## q1 <- lapply(1:length(fits[-9]), function(x)
##   filter(fits[[x]]$std, par=='log_q1') %>%
##   cbind(model=model_names[x])) %>%
##   bind_rows() %>% mutate(lwr=est-1.96*se, upr=est+1.96*se)
## ggplot(q1, aes(year, est, ymin=lwr, ymax=upr, fill=model, color=model)) +
##   geom_line() + geom_ribbon(alpha=.1)
## q3 <- lapply(1:length(fits[-9]), function(x)
##   filter(fits[[x]]$std, par=='log_q3') %>%
##   cbind(model=model_names[x])) %>%
##   bind_rows() %>% mutate(lwr=est-1.96*se, upr=est+1.96*se)
## ggplot(q3, aes(year, est, ymin=lwr, ymax=upr, fill=model, color=model)) +
##   geom_line() + geom_ribbon(alpha=.1)
