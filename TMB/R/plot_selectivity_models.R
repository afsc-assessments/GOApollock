library(tidyverse)
theme_set(theme_bw())
source("TMB/R/Functions/osa_comp_plots.R")
source("TMB/R/Functions/construct_Q.R")
source("TMB/R/Functions/plot_selectivity.R")
load("TMB/Selectivity_runs.RData")
years <- 1970:(2022+5)
rho_trans <- function(x)2/(1+exp(-2*x))-1

m7 <- data.frame(
  par=names(sdrep_mod7$par.fixed),
  est=sdrep_mod7$par.fixed,
           se=sqrt(diag(sdrep_mod7$cov.fixed))) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se) %>%
  filter(par %in% c("sel_rho_a", "sel_rho_y", "ln_sel_sd")) %>%
  mutate(est=ifelse(par=='ln_sel_sd', exp(est), rho_trans(est)),
         lwr=ifelse(par=='ln_sel_sd', exp(lwr), rho_trans(lwr)),
         upr=ifelse(par=='ln_sel_sd', exp(upr), rho_trans(upr)))
m7 %>% print(digits=3)

m9 <- data.frame(
  par=names(sdrep_mod9$par.fixed),
  est=sdrep_mod9$par.fixed,
           se=sqrt(diag(sdrep_mod9$cov.fixed))) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se) %>%
  filter(par %in% c("sel_rho_a", "sel_rho_y", "ln_sel_sd")) %>%
  mutate(est=ifelse(par=='ln_sel_sd', exp(est), est),
         lwr=ifelse(par=='ln_sel_sd', exp(lwr), lwr),
         upr=ifelse(par=='ln_sel_sd', exp(upr), upr))
m9 %>% print(digits=3)

## Combine objects ----
model_names <- c(
  'TMB-Mod1: RE',
  'TMB-Mod2: Log-AR1-Age',
  'TMB-Mod3: Log-AR1-Yr',
  'TMB-Mod4: Log-2D-AR1',
  'TMB-Mod5: AR1-Age',
  'TMB-Mod6: AR1-Yr',
  'TMB-Mod7: 2D-AR1',
  "TMB-Mod8: 3D-AR1cond",
  'TMB-Mod9: 3D-AR1mar'
)


# - Sdrep objects
mod_list <- list(sdrep_mod1, sdrep_mod2, sdrep_mod3, sdrep_mod4, sdrep_mod5, sdrep_mod6, sdrep_mod7, sdrep_mod8, sdrep_mod9)

# - Opt objects
opt_list <- list(opt_mod1, opt_mod2, opt_mod3, opt_mod4, opt_mod5, opt_mod6, opt_mod7, opt_mod8, opt_mod9)
lapply(opt_list, function(x) x$evaluations)

obj_list <- list(obj_mod1, obj_mod2, obj_mod3, obj_mod4, obj_mod5, obj_mod6, obj_mod7, obj_mod8, obj_mod9)

par_list <- lapply(obj_list, function(x) x$env$parList())

# - ssb
ssb <- rbind(cbind(model=model_names[1],filter(stdtmb_mod1, par=='Espawnbio')),
             cbind(model=model_names[2],filter(stdtmb_mod2, par=='Espawnbio')),
             cbind(model=model_names[3],filter(stdtmb_mod3, par=='Espawnbio')),
             cbind(model=model_names[4],filter(stdtmb_mod4, par=='Espawnbio')),
             cbind(model=model_names[5],filter(stdtmb_mod5, par=='Espawnbio')),
             cbind(model=model_names[6],filter(stdtmb_mod6, par=='Espawnbio')),
             cbind(model=model_names[7],filter(stdtmb_mod7, par=='Espawnbio')),
             cbind(model=model_names[8],filter(stdtmb_mod8, par=='Espawnbio')),
             cbind(model=model_names[9],filter(stdtmb_mod9, par=='Espawnbio')),
             cbind(model='ADMB',filter(stdadmb, name=='Espawnbio') %>% select(par=name, est,se, year)))


## TABLES ----
write.csv(ssb, file = "TMB/Output/Selectivity_runs_ssb_full_estimation.csv")

# - Show recent ssb
save_ssb <- ssb %>%
  filter(year>2015) %>%
  arrange(year) %>%
  pivot_wider(names_from = model, values_from = c(est, se)) %>%
  select(year, paste0(c("est_","se_"),rep(unique(ssb$model), each = 2))) %>%
  as.data.frame()
write.csv(save_ssb, file = "TMB/Output/Selectivity_runs_terminal_ssb_full_estimation.csv")

# - AIC
aic_table <- data.frame(
  model = model_names,
  nll = sapply(opt_list, function(x) x$objective),
  k = sapply(opt_list, function(x) length(x$par)),
  AIC = sapply(opt_list, function(x) TMBAIC(x))
) %>%
  mutate(dAIC = AIC - min(AIC))

write.csv(aic_table, file = "TMB/Output/Selectivity_runs_AIC_full_estimation.csv")


## PLOTS ----
# * Plot SSB ----
ggplot(ssb, aes(year, est, color=model, fill=model, ymin=est-1.96*se, ymax=est+1.96*se)) +
  geom_ribbon(alpha=.5) + geom_line(lwd=2)

# * Plot selectivity ----
plot_selectivity(sdrep = mod_list,
                 model_names = model_names)

get_fsh_selex <- function(x, model){
  df <- data.frame(name=names(x$value), est=x$value, se=sqrt(diag(x$cov)))
  df <- filter(df, name=='slctfsh') %>%
    mutate(age=rep(1:10, each=length(years)),
           year=rep(years, times=10), model=model)
  df
}
out <- lapply(1:length(mod_list), function(x)
  get_fsh_selex(mod_list[[x]], model_names[x])) %>%
  bind_rows %>% filter(grepl(x=model, '4|7|9')) %>%
  mutate(lwr=pmax(0, est-1*se), upr=pmin(1.5, est+1*se))
 ## mutate(lwr=pmax(0, est-1.96*se), upr=pmin(1, est+1.96*se))

ggplot(out, aes(year, est, color=model)) + facet_wrap('age', scales='free') +
  geom_line()

ggplot(out, aes(year, est, ymin=lwr, ymax=upr, fill=model, color=model)) + facet_wrap('age', scales='free') +
  geom_line() + geom_ribbon(alpha=.3)
ggplot(filter(out, age %in% 2:5), aes(year, est, ymin=lwr, ymax=upr, fill=model, color=model)) + facet_wrap('age', scales='free') +
  geom_line() +## geom_ribbon(alpha=.3) +
  xlim(2020,2027)

ggplot(filter(out, year>2020), aes(age, est, ymin=lwr, ymax=upr, fill=model, color=model)) + facet_wrap('year', scales='free') +
  geom_line() + geom_ribbon(alpha=.3)


# * Plot OSA residuals ----
quantities_list <- list(quantities_mod1, quantities_mod2, quantities_mod3, quantities_mod4, quantities_mod5, quantities_mod6, quantities_mod7, quantities_mod8, quantities_mod9)

r_osa_list <- list()
g_osa_list <- list()

# - ADMB
r_osa_list[[1]] <- plot_osa_comps(obs=replist$Fishery_observed_and_expected_age_comp[,1:10],
                                  exp=replist$Fishery_observed_and_expected_age_comp[,11:20],
                                  ages=1:10, years=dat$fshyrs, Neff=dat$multN_fsh,
                                  model='ADMB final 2022')
g_osa_list[[1]] <- plot_osa_comps(obs=replist$Fishery_observed_and_expected_age_comp[,1:10],
                                  exp=replist$Fishery_observed_and_expected_age_comp[,11:20],
                                  ages=1:10, years=dat$fshyrs, Neff=dat$multN_fsh,
                                  model='ADMB final 2022', plot.type='ggplot')

# - TMB
for(i in 1:length(quantities_list)){
  r_osa_list[[i+1]] <- try(plot_osa_comps(obs=quantities_list[[i]]$res_fish[,1:10],
                                          exp=quantities_list[[i]]$res_fish[,11:20],
                                          ages=1:10,
                                          years=dat$fshyrs,
                                          Neff=dat$multN_fsh,
                                          model=model_names[i],
                                          plot.type='default'), silent = TRUE)
  g_osa_list[[i+1]] <- try(plot_osa_comps(obs=quantities_list[[i]]$res_fish[,1:10],
                                          exp=quantities_list[[i]]$res_fish[,11:20],
                                          ages=1:10,
                                          years=dat$fshyrs,
                                          Neff=dat$multN_fsh,
                                          model=model_names[i],
                                          plot.type='ggplot'), silent = TRUE)
}

sapply(r_osa_list, class) # 1 doesnt work
sapply(g_osa_list, class) # 1 and 5 doesnt work

# Only show models that make sense (2D ar1 and 3d ar)
cowplot::plot_grid(g_osa_list[[1]], g_osa_list[[2]], g_osa_list[[5]], g_osa_list[[8]], g_osa_list[[9]], g_osa_list[[10]])


# * Plot cor ----
mod = 9
plot_cor(n_years = 3, n_ages = dat$nages,
         rho_y = par_list[[mod]]$sel_rho_y,
         rho_a = par_list[[mod]]$sel_rho_a,
         rho_c = par_list[[mod]]$sel_rho_c, log_sigma = par_list[[mod]]$ln_sel_sd)



