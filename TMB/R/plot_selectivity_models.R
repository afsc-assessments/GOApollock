## Load functions and packages ----
library(tidyverse)
theme_set(theme_bw())
source("TMB/R/prepare_tmb_objects.R")
source("TMB/R/Functions/osa_comp_plots.R")
source("TMB/R/Functions/construct_Q.R")
source("TMB/R/Functions/plot_selectivity.R")
source("TMB/R/proj_calcs/run_projections.R")

## Load data/runs ----
# load("TMB/Selectivity_runs.RData")
fits <- readRDS("TMB/Output/fits.RDS")


years <- 1970:(2022+5)
rho_trans <- function(x)2/(1+exp(-2*x))-1

m7 <- data.frame(
  par=names(fits$fit7$SD$par.fixed),
  est=fits$fit7$SD$par.fixed,
  se=sqrt(diag(fits$fit7$SD$cov.fixed))) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se) %>%
  filter(par %in% c("sel_rho_a", "sel_rho_y", "ln_sel_sd")) %>%
  mutate(est=ifelse(par=='ln_sel_sd', exp(est), rho_trans(est)),
         lwr=ifelse(par=='ln_sel_sd', exp(lwr), rho_trans(lwr)),
         upr=ifelse(par=='ln_sel_sd', exp(upr), rho_trans(upr)))
m7 %>% print(digits=3)

# m9 <- data.frame(
#   par=names(fits$fit9$SD$par.fixed),
#   est=fits$fit9$SD$par.fixed,
#   se=sqrt(diag(fits$fit9$SD$cov.fixed))) %>%
#   mutate(lwr=est-1.96*se, upr=est+1.96*se) %>%
#   filter(par %in% c("sel_rho_a", "sel_rho_y", "ln_sel_sd")) %>%
#   mutate(est=ifelse(par=='ln_sel_sd', exp(est), est),
#          lwr=ifelse(par=='ln_sel_sd', exp(lwr), lwr),
#          upr=ifelse(par=='ln_sel_sd', exp(upr), upr))
# m9 %>% print(digits=3)

## Combine objects ----
model_names <- c(
  'Mod 0: Constant',
  'Mod 1: ParDevs',
  'Mod 2: Log-AR1-Age',
  'Mod 3: Log-AR1-Yr',
  'Mod 4: Log-2D-AR1',
  'Mod 5: Age-specific',
  'Mod 6: AR1-Yr',
  'Mod 7: 2D-AR1',
  "Mod 8: 3D-AR1cond"
  # , 'Mod 9: 3D-AR1mar'
)


# - Sdrep objects
# mod_list <- list(sdrep_mod0, sdrep_mod1, sdrep_mod2, sdrep_mod3, sdrep_mod4, sdrep_mod5, sdrep_mod6, sdrep_mod7, sdrep_mod8, sdrep_mod9)
mod_list <- lapply(fits[1:9], function(x) x$SD)

# - Opt objects
# opt_list <- list(opt_mod0, opt_mod1, opt_mod2, opt_mod3, opt_mod4, opt_mod5, opt_mod6, opt_mod7, opt_mod8, opt_mod9)
opt_list <- lapply(fits[1:9], function(x) x)

# - Obj objects
# obj_list <- list(obj_mod0, obj_mod1, obj_mod2, obj_mod3, obj_mod4, obj_mod5, obj_mod6, obj_mod7, obj_mod8, obj_mod9)

# - Reported quantities
# quantities_list <- list(quantities_mod0, quantities_mod1, quantities_mod2, quantities_mod3, quantities_mod4, quantities_mod5, quantities_mod6, quantities_mod7, quantities_mod8, quantities_mod9)
quantities_list <- lapply(fits[1:9], function(x) x$report)

# - Parameters
par_list <- lapply(fits[1:9], function(x) x$parList)
parse_list <- lapply(fits[1:9], function(x) data.frame(name = names(x$SD$par.fixed), par = x$SD$par.fixed, se = diag(x$SD$cov.fixed)))
rho_list <- lapply(parse_list, function(x) x %>%
                     mutate(lwr=par-1.96*se, upr=par+1.96*se) %>%
                     filter(name %in% c("sel_rho_a",  "sel_rho_y", "sel_rho_c", "ln_sel_sd")) %>%
                     mutate(est=ifelse(name=='ln_sel_sd', exp(par), rho_trans(par)),
                            lwr=ifelse(name=='ln_sel_sd', exp(lwr), rho_trans(lwr)),
                            upr=ifelse(name=='ln_sel_sd', exp(upr), rho_trans(upr)))
)

## ABC ----
abc_list <- lapply(quantities_list, function(x) abc_calc_tmb(x, datlist = dat))


# - Check number of evaluations
for(i in 1:length(opt_list)){
  print(opt_list[[i]]$evaluations)
}


# - ssb
ssb <- rbind(cbind(model=model_names[1],filter(fits[[1]]$std, par=='Espawnbio')),
             cbind(model=model_names[2],filter(fits[[2]]$std, par=='Espawnbio')),
             cbind(model=model_names[3],filter(fits[[3]]$std, par=='Espawnbio')),
             cbind(model=model_names[4],filter(fits[[4]]$std, par=='Espawnbio')),
             cbind(model=model_names[5],filter(fits[[5]]$std, par=='Espawnbio')),
             cbind(model=model_names[6],filter(fits[[6]]$std, par=='Espawnbio')),
             cbind(model=model_names[7],filter(fits[[7]]$std, par=='Espawnbio')),
             cbind(model=model_names[8],filter(fits[[8]]$std, par=='Espawnbio')),
             cbind(model=model_names[9],filter(fits[[9]]$std, par=='Espawnbio')),
             # cbind(model=model_names[10],filter(fits[[10]]$std, par=='Espawnbio')),
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
  nll = unlist(sapply(opt_list, function(x) as.numeric(x$objective))),
  fish_age_nll= sapply(quantities_list, function(x) round(-x$loglik[2],1)),
  k = sapply(opt_list, function(x) length(x$par)),
  AIC = unlist(sapply(opt_list, function(x) as.numeric(x$AIC)))
) %>%
  mutate(dAIC = AIC - min(AIC),
         SSB2022 = sapply(quantities_list, function(x) x$Espawnbio[length(x$Espawnbio)] * 1000000),
         SSB2023 = sapply(abc_list, function(x) x$exec_tableF[3,2]),
         B0 = sapply(abc_list, function(x) x$exec_tableF[4,2]),
         B40 = sapply(abc_list, function(x) x$exec_tableF[5,2]),
         B35 = sapply(abc_list, function(x) x$exec_tableF[6,2]),
         abc2023 = sapply(abc_list, function(x) x$exec_tableF[12,2]),
         ofl2023 = sapply(abc_list, function(x) x$exec_tableF[10,2]))

write.csv(aic_table, file = "TMB/Output/Selectivity_runs_AIC_full_estimation.csv")


par_table <- data.frame(
  model = model_names,
  sigma = sapply(quantities_list, function(x) x$rho_a),
  rho_a = sapply(quantities_list, function(x) x$rho_a),
  rho_y = sapply(quantities_list, function(x) x$rho_y),
  rho_c = sapply(quantities_list, function(x) x$rho_c))


## PLOTS ----
# * Plot SSB ----
ssb %>%
  filter(grepl(x=model, paste0("Mod ", c(0,1,7,8), collapse = "|")) | model == "ADMB") %>%
ggplot(aes(year, est, color=model, fill=model, ymin=est-1.96*se, ymax=est+1.96*se)) +
  #geom_ribbon(alpha=.5) +
  geom_line(lwd=2)

# * Plot selectivity ----
plot_selectivity(sdrep = mod_list,
                 model_names = model_names)

get_fsh_selex <- function(x, model){
  df <- data.frame(name=names(x$value), est=x$value, se=x$sd)
  df <- filter(df, name=='slctfsh') %>%
    mutate(age=rep(1:10, each=length(years)),
           year=rep(years, times=10), model=model)
  df
}
out <- lapply(1:length(mod_list), function(x)
  get_fsh_selex(mod_list[[x]], model_names[x])) %>%
  bind_rows %>% filter(grepl(x=model, paste0("Mod ", c(0,1,7,8), collapse = "|"))) %>%
  mutate(lwr=pmax(0, est-1*se), upr=pmin(1.5, est+1*se))
## mutate(lwr=pmax(0, est-1.96*se), upr=pmin(1, est+1.96*se))

ggplot(out, aes(year, est, color=model)) + facet_wrap('age', scales='free') +
  geom_line()

# * Plot Ribbon ----
out$agetext <- paste("Age", out$age)
ggplot(out, aes(year, est, ymin=lwr, ymax=upr, fill=model, color=model)) + facet_wrap('agetext') +
  geom_line() + geom_ribbon(alpha=.3) + ylab("Selectivity") + xlab("Year") + theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.background = element_blank(),
    strip.text = element_text(face="bold", size=9),
    legend.position = "top")


ggplot(filter(out, age %in% 2:5), aes(year, est, ymin=lwr, ymax=upr, fill=model, color=model)) + facet_wrap('age', scales='free') +
  geom_line() +## geom_ribbon(alpha=.3) +
  xlim(2020,2027)

ggplot(filter(out, year>2020), aes(age, est, ymin=lwr, ymax=upr, fill=model, color=model)) + facet_wrap('year', scales='free') +
  geom_line() + geom_ribbon(alpha=.3)


# * Plot OSA residuals ----
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
prow <- cowplot::plot_grid(g_osa_list[[2]] + theme(legend.position="none"),
                   g_osa_list[[3]] + theme(legend.position="none"),
                   g_osa_list[[9]] + theme(legend.position="none"),
                   g_osa_list[[10]] + theme(legend.position="none"))

legend <- get_legend(
  # create some space to the left of the legend
  g_osa_list[[2]] + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of
# the width of one plot (via rel_widths).
plot_grid(prow, legend, rel_widths = c(3, .4))


# * Plot cor ----
mod = 9
plot_cor(n_years = 3, n_ages = dat$nages,
         rho_y = par_list[[mod]]$sel_rho_y,
         rho_a = par_list[[mod]]$sel_rho_a,
         rho_c = par_list[[mod]]$sel_rho_c, log_sigma = par_list[[mod]]$ln_sel_sd)



