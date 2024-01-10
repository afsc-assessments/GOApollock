
library(GOApollock)
library(TMB)

file.copy('../source/goa_pk.tpl', to='model_23/goa_pk.tpl', overwrite = TRUE)
setwd('model_23')
clean_pk_dir(full=TRUE)
system("admb goa_pk", ignore.stdout = TRUE)
system("goa_pk -iprint 1000")
clean_pk_dir(full=FALSE)
setwd('..')


library(GOApollock)
library(TMB)
dyn.unload('model_23/goa_pk_tmb.dll')
file.copy('../source/goa_pk_tmb.cpp', to='model_23/goa_pk_tmb.cpp', overwrite = TRUE)
input <- prepare_pk_input('model_23', datfile='goa_pk.dat')
tmb <- fit_pk(input, do.fit=TRUE)
trep <- get_rep(tmb)
arep <- read_pk_rep(version='ADMB', endyr=2023, path='model_23')

## check the two match
trep$N_proj-arep$Projection_numbers_at_age
cbind(arep$Projection_spawning_biomass, trep$Espawnbio_proj)
cbind(arep$Projection_total_catches,trep$Ecattot_proj)
cbind(arep$Projection_fishery_selectivity,trep$slctfsh_proj)
cbind(arep$Projection_summary_biomass,trep$Esumbio_proj)
trep$Z_proj-arep$Projection_total_mortality


## Try TMB model with some new Ftargets
input <- prepare_pk_input('model_23', datfile='goa_pk.dat')
nyrs <- 20

input1 <- input
input1$dat$Ftarget <- rep(0,nyrs)
input1$version <- 'No fishing'
tmb1 <- fit_pk(input1, do.fit=TRUE)
sd1 <- get_std(tmb1)
ssb1 <- filter(sd1, name %in% c('Espawnbio', 'Espawnbio_proj')) %>%
  mutate(year=ifelse(name=='Espawnbio', year, year+54))

input2 <- input
input2$dat$Ftarget <- rep(.5,nyrs)
input2$version <- 'F=.5'
tmb2 <- fit_pk(input2, do.fit=TRUE)
sd2 <- get_std(tmb2)
ssb2 <- filter(sd2, name %in% c('Espawnbio', 'Espawnbio_proj')) %>%
  mutate(year=ifelse(name=='Espawnbio', year, year+54))

input3 <- input
input3$dat$Ftarget <- rep(1,nyrs)
input3$version <- 'F=1'
tmb3 <- fit_pk(input3, do.fit=TRUE)
sd3 <- get_std(tmb3)
ssb3 <- filter(sd3, name %in% c('Espawnbio', 'Espawnbio_proj')) %>%
  mutate(year=ifelse(name=='Espawnbio', year, year+54))

input4 <- input
input4$dat$Ftarget <- rep(input$dat$Ftarget[1],nyrs)
input4$version <- 'F=F40'
tmb4 <- fit_pk(input4, do.fit=TRUE)
sd4 <- get_std(tmb4)
ssb4 <- filter(sd4, name %in% c('Espawnbio', 'Espawnbio_proj')) %>%
  mutate(year=ifelse(name=='Espawnbio', year, year+54))


ssb <- bind_rows(ssb1,ssb2,ssb3,ssb4) %>% filter(year>2010)
ggplot(ssb, aes(year, est, ymin=lwr, ymax=upr, color=version,
                fill=version))+ geom_vline(xintercept=2023) +
  geom_hline(yintercept=input$dat$B40) +
  geom_ribbon(alpha=.5)+ geom_line() + scale_y_log10()



get_rep(list(tmb1, tmb2, tmb3,tmb4)) %>% mymelt('F_proj')

sbio <- seq(0,.6, len=100)
d <- sbio/input$dat$B40
F <- .2
F2 <- ifelse(d>1, F, F*((d-0.05)/(1-0.05)))
plot(d,F2, xlim=c(0,2), ylim=c(0,input$dat$B40))
