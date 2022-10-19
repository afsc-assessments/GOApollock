
## super dangerous to do this but oh well
mapoff <- function(name) input$map[[name]] <<- as.factor(input$par[[name]]*NA)

par_table <- function(fit, v=NULL){
  x <- table(names(fit$opt$par)) %>% as.data.frame %>%
    setNames(c("par","length"))
  x <- rbind(x, data.frame(par='total', length=sum(x$length)))
  if(!is.null(v)) x <- cbind(x, version=v)
  x
}

plot_ssb <- function(fit, asdrep, log=TRUE){
  ssb <- get_ssb(asdrep, fit)
  if(log){
  g <- ggplot(ssb, aes(year, est, ymin=est-1.96*sd, ymax=est+1.96*sd,
                       fill=platform, color=platform))+
    labs(y='log SSB')
  } else {
  g <- ggplot(ssb, aes(year, exp(est), ymin=exp(est-1.96*sd), ymax=exp(est+1.96*sd),
                       fill=platform, color=platform)) +
  ylim(0,NA) + labs(y='SSB')
  }
  g <- g+
    geom_ribbon(alpha=.5) + geom_line(lwd=1) +
    theme(legend.position='top')+
    labs( color=NULL, fill=NULL)
  g
}

get_wham_ssb <- function(fit, v=NULL){
  if(is.null(fit$sdrep)){
   x <- data.frame(name='log_SSB', est=log(fit$rep$SSB), sd=NA)

  } else {
    x <- data.frame(name=names(fit$sdrep$value),
                    est=fit$sdrep$value, sd=fit$sdrep$sd)
  }
  out <- cbind(platform='wham', year=1970:2021,
               filter(x, name=='log_SSB') %>% select(-name))
  if(!is.null(v)) out <- cbind(out, version=v)
  out
}
get_admb_ssb <- function(x){
  if(!'log_ssb' %in% x$name) stop("old model version without log_ssb")
  cbind(platform='admb',
        filter(x, name=='log_ssb') %>% select(year, est, sd=se))

}
get_ssb <- function(admb, wham){ rbind(get_wham_ssb(wham), get_admb_ssb(admb))}


get_wham_nll <- function(tmp) {
  nlls <- sapply(tmp[grep('nll', names(tmp))], sum)[-4] %>% round(3)
  sum(nlls)-  tmp$nll
  wham <- with(tmp,c(nll, sum(nll_agg_catch), sum(nll_catch_acomp), 0,
                     sum(nll_agg_indices[,1]), sum(nll_index_acomp[,1]), 0,
                     sum(nll_agg_indices[,2]), sum(nll_index_acomp[,2]), 0, 0,
                     sum(nll_agg_indices[,3]), sum(nll_index_acomp[,3]), 0,
                     sum(nll_agg_indices[,4]), sum(nll_agg_indices[,5]),
                     sum(nll_agg_indices[,6]),
                     sum(nll_index_acomp[,6]),
                     nll_NAA, nll_sel, 0, sum(nll_q), 0, sum(nll_cole),0))
  wham
}


match_input <- function(aa, asap3, NAA_re=list(sigma="rec", cor="iid")){
  a <- 0;b <- 1
  invsel <- function(y) -log((b-a)/(y-a)-1)
  sel <- function(x) a+(b-a)/(1+exp(-x))
  tmp <- aa$Fishery_selectivity
  tmp[tmp>=1] <- 1-1e-15
  tmp <- invsel(tmp)
  fshselmean <- sel(colMeans(tmp))
  ## fshselmean <- log(y/(1-y))
  y <- colMeans(aa$Fishery_selectivity)
  test <- (exp(y)/(1+exp(y)))
  ## y2 <- sel(colMeans(invsel(aa$Fishery_selectivity)))
  ## test2 <- (exp(y2)/(1+exp(y2)))
  fshseldevs <-  t(apply(tmp,1, FUN=function(x) x-y))
  ## This oen is in selex space (0,1) and is trasnformed as log(sel/(1-sel))
  selmods <- rep('double-logistic',7)
  ## survey 4 and 5 age based to be able to fix 0s and 1s
  ## appropriately
  selmods[c(2,5:6)] <- 'age-specific'
  selres <- c('iid',rep('none',6))
  selinits <- list()
  selinits[[1]] <- test#sel(fshselmean)
  selinits[[2]] <- aa$Survey_1_selectivity
  selinits[[3]] <- aa$Survey_2_selectivity
  selinits[[4]] <- aa$Survey_3_selectivity
  selinits[[5]] <- c(1,rep(0,9))
  selinits[[6]] <- c(0, 1, rep(0,8))
  selinits[[7]] <- aa$Survey_6_selectivity
  ## shelikof and adf&g are timevarying
  qres <- rep("none",6); qres[c(1,3)] <- 'ar1'
  NAA_re <- c(NAA_re,
              list(N1_model=1,
              N1_pars=c(exp(13.5),0)))
  input <- prepare_wham_input(asap3, recruit_model=2,
                              model_name="GOA pollock",
                              selectivity=
                                list(model=selmods,
                                     selinits=selinits,
                                     re=selres),
                                NAA_re = NAA_re,
                              catchability=list(re=qres))
  input$random <- NULL
  input$par$log_NAA <- log(1e6*aa$Numbers_at_age)[-1,]
  input$par$log_NAA[,1] <- log(1e6*aa$Numbers_at_age)[-1,1]
  ## Get initial age structure setup the same
  ## Only first column (recruits) used b/ cmapped off due to the
  ## NAA_re setting (age 2+ deterministic). To perfectly match
  ## ADMB need to set N1_model=0 and uncomment these as inits
  ## input$par$log_N1_pars <- log(1e6*c(aa$Numbers_at_age[1,1], aa$Initial_age_comp))
  input$map$log_N1_pars <- factor(c(1,NA)) # assume F=0 at beginning
  input$par$mean_rec_pars <- mean(log(1e6*aa$Recruits))
  input$par$log_F1 <- log(aa$Fishing_mortalities[1]) ## mean(log(aa$Fishing_mortalities))
  F_devs <- log(aa$Fishing_mortalities)[-1]- log(aa$Fishing_mortalities)[-52]
  input$par$F_devs[,1] <- F_devs# matrix(ncol=1,F_devs[-52])
  input$data$mature <- input$data$mature*.5
  ## catchabilities
  invq <- function(y,a=0,b=1000) -log((b-a)/(y-a)-1)
  q <- function(x,a=0,b=1000) a+(b-a)/(1+exp(-x))
  ## set the means
  input$par$logit_q <- c(mean(invq(aa$Survey_1_q[-1])),
                         invq(aa$Survey_2_q[1]),
                         mean(invq(aa$Survey_3_q)),
                         invq(aa$Survey_4_q),
                         invq(aa$Survey_5_q),
                         invq(aa$Survey_6_q))
  ### tried this to get prior on BT to match but didn't work
  input$data$use_q_prior[2] <- 1
  input$data$logit_q_prior_sigma[2] <- .1
  input$par$q_prior_re[2] <- invq(.85)

  input$par$q_re[,1] <- invq(aa$Survey_1_q[-1])-mean(invq(aa$Survey_1_q[-1]))
  input$par$q_re[,3] <- invq(aa$Survey_3_q)-mean(invq(aa$Survey_3_q))
  input$data$use_q_re
  ## first column log SDs, second rhos
  input$par$q_repars[1,1] <- log(0.038)
  input$par$q_repars[3,1] <- log(0.05)
  input$par$q_repars[1,2] <- 10
  input$par$q_repars[3,2] <- 10
  ## currently assumes first year is constatn but is not for some
  ## surveys, this is hacked in the TPL for bridging
  input$data$fracyr_indices <- matrix(c(0.209, .543, 0.60989, 0,0,.519), nrow=52, ncol=6, byrow=TRUE)
  ## WAA matrices are rounded in ADMB but not WHAM, so recreate
  ## that here.
  input$data$waa[1,,] <- round(input$data$waa[1,,],3)
  input$data$waa[2,,] <- round(input$data$waa[2,,],5)
  input$data$waa[3,,] <- round(input$data$waa[3,,],3)
  input$data$waa[4,,] <- round(input$data$waa[4,,],3)
  input$data$bias_correct_oe <- 0
  input$data$bias_correct_pe <- 0
  ## Get selex to match
  invsel <- function(y, a=-10,b=20) -log((b-a)/(y-a)-1)
  sel <- function(x,a=-10,b=20) a+(b-a)/(1+exp(-x))
  ## srv1
  input$data$selpars_lower[,13:16] <- -10
  input$data$selpars_upper[,13:16] <- 20
  ## fishery, see wham.xlsx
  input$par$logit_selpars[1,13:16] <-
    c(-0.16209842,-0.62409404,0.65841736,-0.63755954)
  ##  invsel(c(3.78792789566,0.462458026600594,9.76714339228,0.373952507323135))
  ## input$par$logit_selpars[2,13:16] <-
  ##   invsel(c(-10,0.367879441171442,9.32456556011,0.444118022712467))
  input$par$logit_selpars[3,13:16] <-
    invsel(c(3.47046773474,1.28332083148355,20,0.367879441171442))
  input$par$logit_selpars[4,13:16] <-
    invsel(c(5.24912512874,0.810378025411032,20,0.367879441171442))
### these need to be age specified becuase they're fixed at 0 or 1 for
### different ages
  ## input$par$logit_selpars[5,13:16] <-
  ##   invsel(c(-5,0.367879441171442,15,0.135335283236613))
  ## input$par$logit_selpars[6,13:16] <-
  ##   invsel(c(-5,0.367879441171442,15,0.135335283236613))
  ## Need to be careful b/c invsel returns an Inf sometimes and
  ## that breaks things downstream
  input$par$logit_selpars[2,1:10] <- c(-Inf, -Inf, 10, invsel(aa$Survey_1_selectivity,a=0,b=1)[-(1:3)])
  input$par$logit_selpars[5,1:10] <- c(Inf, rep(-Inf,9))
  input$par$logit_selpars[6,1:10] <- c(-Inf, Inf, rep(-Inf,8))
  input$par$logit_selpars[7,13:16] <-
    invsel(c(0.5,0.00744658307092434,20,0.367879441171442))
  ## Build RE devs to be equivalent in WHAM. See sheet selex in
  ## wham.xlsx.
  devs1 <- c(0.00351, 0.00352, 0.00351, 0.00340, 0.00308, 0.00222, 0.00171, 0.00366, -0.00043, -0.00369, 0.01951, 0.01731, -0.00835, -0.02644, -0.02785, 0.00637, -0.00262, 0.01242, 0.01798, 0.03623, 0.08130, 0.10906, 0.10846, 0.09502, 0.09394, 0.09269, 0.09209, 0.08701, 0.06782, 0.06086, 0.04796, 0.00008, -0.04083, -0.07066, -0.10437, -0.09489, -0.09661, -0.10785, -0.09021, -0.09329, -0.08678, -0.06532, -0.02936, -0.00898, 0.02366, 0.02585, 0.01483, 0.00145, -0.02066, -0.03768, -0.05984, -0.05984)
  devs2 <- c(-0.00927, -0.00927, -0.00927, -0.00927, -0.00926, -0.00925, -0.00926, -0.00926, -0.00910, -0.00890, -0.00688, -0.00475, -0.00159, 0.00038, 0.00207, 0.00356, 0.00554, 0.00666, 0.00778, 0.00767, 0.00716, 0.00773, 0.00850, 0.01009, 0.01060, 0.01057, 0.01079, 0.01099, 0.01106, 0.01020, 0.00967, 0.01057, 0.01088, 0.00989, 0.00866, 0.00574, 0.00370, 0.00198, -0.00130, -0.00242, -0.00372, -0.00553, -0.00761, -0.00892, -0.00995, -0.00928, -0.00838, -0.00749, -0.00651, -0.00587, -0.00507, -0.00507)
  input$par$selpars_re[1:104] <- c(devs1,devs2)
  input$map$selpars_re <-
    factor(c(1:104, rep(NA,length(input$map$selpars_re)-104)))
  ## match estimation of selex
  tmp <- matrix(input$map$logit_selpars, nrow=7)
  tmp[,11:12] <- NA ## turn off all logistics parameters
  ##tmp[2,13:14] <- NA # survey one has first two pars fixed
  ##(decreasing)
  tmp[2,c(1:2,11:16)] <- NA # survey one has first two ages fixed at 0
  tmp[3:4,15:16] <- NA # survey two and three last two pars fixed (increasing)
  tmp[7,] <- NA ## survey 6 assumed constant 1
  input$map$logit_selpars <- factor(tmp)
  input
}

match_input_nll <- function(aa, asap3, NAA_re=list(sigma="rec", cor="iid")){
  a <- 0;b <- 1
  invsel <- function(y) -log((b-a)/(y-a)-1)
  sel <- function(x) a+(b-a)/(1+exp(-x))
  tmp <- aa$Fishery_selectivity
  tmp[tmp>=1] <- 1-1e-15
  tmp <- invsel(tmp)
  fshselmean <- sel(colMeans(tmp))
  ## fshselmean <- log(y/(1-y))
  y <- colMeans(aa$Fishery_selectivity)
  test <- (exp(y)/(1+exp(y)))
  ## y2 <- sel(colMeans(invsel(aa$Fishery_selectivity)))
  ## test2 <- (exp(y2)/(1+exp(y2)))
  fshseldevs <-  t(apply(tmp,1, FUN=function(x) x-y))
  ## This oen is in selex space (0,1) and is trasnformed as log(sel/(1-sel))
  selmods <- rep('double-logistic',7)
  ## survey 4 and 5 age based to be able to fix 0s and 1s
  ## appropriately
  selmods[c(2,5:6)] <- 'age-specific'
  selres <- c('iid',rep('none',6))
  selinits <- list()
  selinits[[1]] <- test#sel(fshselmean)
  selinits[[2]] <- aa$Survey_1_selectivity
  selinits[[3]] <- aa$Survey_2_selectivity
  selinits[[4]] <- aa$Survey_3_selectivity
  selinits[[5]] <- c(1,rep(0,9))
  selinits[[6]] <- c(0, 1, rep(0,8))
  selinits[[7]] <- aa$Survey_6_selectivity
  ## shelikof and adf&g are timevarying
  qres <- rep("none",6); qres[c(1,3)] <- 'ar1'
  NAA_re <- c(NAA_re, list(N1_model=0))
  input <- prepare_wham_input(asap3, recruit_model=2,
                              model_name="GOA pollock",
                              selectivity=
                                list(model=selmods,
                                     selinits=selinits,
                                     re=selres),
                                NAA_re = NAA_re,
                              catchability=list(re=qres))
  input$random <- NULL
  input$par$log_NAA <- log(1e6*aa$Numbers_at_age)[-1,]
  input$par$log_NAA[,1] <- log(1e6*aa$Numbers_at_age)[-1,1]
  ## Get initial age structure setup the same
  ## Only first column (recruits) used b/ cmapped off due to the
  ## NAA_re setting (age 2+ deterministic).
  ## To perfectly match ADMB need to set N1_model=0 and uncomment
  ## these as inits
  input$par$log_N1_pars <-
    log(1e6*c(aa$Numbers_at_age[1,1], aa$Initial_age_comp))
  ## input$map$log_N1_pars <- factor(c(1,NA)) # assume F=0 at beginning
  input$par$mean_rec_pars <- mean(log(1e6*aa$Recruits))
  input$par$log_F1 <- log(aa$Fishing_mortalities[1]) ## mean(log(aa$Fishing_mortalities))
  F_devs <- log(aa$Fishing_mortalities)[-1]- log(aa$Fishing_mortalities)[-52]
  input$par$F_devs[,1] <- F_devs# matrix(ncol=1,F_devs[-52])
  input$data$mature <- input$data$mature*.5
  ## catchabilities
  invq <- function(y,a=0,b=1000) -log((b-a)/(y-a)-1)
  q <- function(x,a=0,b=1000) a+(b-a)/(1+exp(-x))
  ## set the means
  input$par$logit_q <- c(mean(invq(aa$Survey_1_q[-1])),
                         invq(aa$Survey_2_q[1]),
                         mean(invq(aa$Survey_3_q)),
                         invq(aa$Survey_4_q),
                         invq(aa$Survey_5_q),
                         invq(aa$Survey_6_q))
  ### tried this to get prior on BT to match but didn't work
  ## input$data$use_q_prior[2] <- 1
  ## input$data$logit_q_prior_sigma[2] <- .1
  ## input$par$q_prior_re[2] <- invq(.85)

  input$par$q_re[,1] <- invq(aa$Survey_1_q[-1])-mean(invq(aa$Survey_1_q[-1]))
  input$par$q_re[,3] <- invq(aa$Survey_3_q)-mean(invq(aa$Survey_3_q))
  input$data$use_q_re
  ## first column log SDs, second rhos
  input$par$q_repars[1,1] <- log(0.038)
  input$par$q_repars[3,1] <- log(0.05)
  input$par$q_repars[1,2] <- 10
  input$par$q_repars[3,2] <- 10
  ## currently assumes first year is constatn but is not for some
  ## surveys, this is hacked in the TPL for bridging
  input$data$fracyr_indices <- matrix(c(0.209, .543, 0.60989, 0.209,0.209,.519), nrow=52, ncol=6, byrow=TRUE)
  ## WAA matrices are rounded in ADMB but not WHAM, so recreate
  ## that here.
  input$data$waa[1,,] <- round(input$data$waa[1,,],3)
  input$data$waa[2,,] <- round(input$data$waa[2,,],5)
  input$data$waa[3,,] <- round(input$data$waa[3,,],3)
  input$data$waa[4,,] <- round(input$data$waa[4,,],3)
  input$data$bias_correct_oe <- 0
  input$data$bias_correct_pe <- 0
  ## Get selex to match
  invsel <- function(y, a=-10,b=20) -log((b-a)/(y-a)-1)
  sel <- function(x,a=-10,b=20) a+(b-a)/(1+exp(-x))
  ## srv1
  input$data$selpars_lower[,13:16] <- -10
  input$data$selpars_upper[,13:16] <- 20
  ## fishery, see wham.xlsx
  input$par$logit_selpars[1,13:16] <-
    c(-0.16209842,-0.62409404,0.65841736,-0.63755954)
  ##  invsel(c(3.78792789566,0.462458026600594,9.76714339228,0.373952507323135))
  ## input$par$logit_selpars[2,13:16] <-
  ##   invsel(c(-10,0.367879441171442,9.32456556011,0.444118022712467))
  input$par$logit_selpars[3,13:16] <-
    invsel(c(3.47046773474,1.28332083148355,20,0.367879441171442))
  input$par$logit_selpars[4,13:16] <-
    invsel(c(5.24912512874,0.810378025411032,20,0.367879441171442))
### these need to be age specified becuase they're fixed at 0 or 1 for
### different ages
  ## input$par$logit_selpars[5,13:16] <-
  ##   invsel(c(-5,0.367879441171442,15,0.135335283236613))
  ## input$par$logit_selpars[6,13:16] <-
  ##   invsel(c(-5,0.367879441171442,15,0.135335283236613))
  input$par$logit_selpars[2,1:10] <- c(-Inf, -Inf, invsel(aa$Survey_1_selectivity,a=0,b=1)[-(1:2)])
  input$par$logit_selpars[5,1:10] <- c(Inf, rep(-Inf,9))
  input$par$logit_selpars[6,1:10] <- c(-Inf, Inf, rep(-Inf,8))
  input$par$logit_selpars[7,13:16] <-
    invsel(c(0.5,0.00744658307092434,20,0.367879441171442))
  ## Build RE devs to be equivalent in WHAM. See sheet selex in
  ## wham.xlsx.
  devs1 <- c(0.00351, 0.00352, 0.00351, 0.00340, 0.00308, 0.00222, 0.00171, 0.00366, -0.00043, -0.00369, 0.01951, 0.01731, -0.00835, -0.02644, -0.02785, 0.00637, -0.00262, 0.01242, 0.01798, 0.03623, 0.08130, 0.10906, 0.10846, 0.09502, 0.09394, 0.09269, 0.09209, 0.08701, 0.06782, 0.06086, 0.04796, 0.00008, -0.04083, -0.07066, -0.10437, -0.09489, -0.09661, -0.10785, -0.09021, -0.09329, -0.08678, -0.06532, -0.02936, -0.00898, 0.02366, 0.02585, 0.01483, 0.00145, -0.02066, -0.03768, -0.05984, -0.05984)
  devs2 <- c(-0.00927, -0.00927, -0.00927, -0.00927, -0.00926, -0.00925, -0.00926, -0.00926, -0.00910, -0.00890, -0.00688, -0.00475, -0.00159, 0.00038, 0.00207, 0.00356, 0.00554, 0.00666, 0.00778, 0.00767, 0.00716, 0.00773, 0.00850, 0.01009, 0.01060, 0.01057, 0.01079, 0.01099, 0.01106, 0.01020, 0.00967, 0.01057, 0.01088, 0.00989, 0.00866, 0.00574, 0.00370, 0.00198, -0.00130, -0.00242, -0.00372, -0.00553, -0.00761, -0.00892, -0.00995, -0.00928, -0.00838, -0.00749, -0.00651, -0.00587, -0.00507, -0.00507)
  input$par$selpars_re[1:104] <- c(devs1,devs2)
  input$map$selpars_re <-
    factor(c(1:104, rep(NA,length(input$map$selpars_re)-104)))
  ## match estimation of selex
  tmp <- matrix(input$map$logit_selpars, nrow=7)
  tmp[,11:12] <- NA ## turn off all logistics parameters
  ##tmp[2,13:14] <- NA # survey one has first two pars fixed
  ##(decreasing)
  tmp[2,c(1:2,11:16)] <- NA # survey one has first two ages fixed at 0
  tmp[3:4,15:16] <- NA # survey two and three last two pars fixed (increasing)
  tmp[7,] <- NA ## survey 6 assumed constant 1
  input$map$logit_selpars <- factor(tmp)
  input
}






plot_checks <- function(aa, ww){
  ## Check NAA
  par(mfrow=c(2,3), mar=c(3,3,2,.5), mgp=c(1.5,.5,0), tck=-.01)
  ## NAA initial and final years
  x1 <- (-1e6*aa$Numbers_at_age[1,]+ww$NAA[1,])/(1e6*aa$Numbers_at_age[1,])
  x2 <- (-1e6*aa$Numbers_at_age[25,]+ww$NAA[25,])/(1e6*aa$Numbers_at_age[25,])
  x3 <- (-1e6*aa$Numbers_at_age[52,]+ww$NAA[52,])/(1e6*aa$Numbers_at_age[52,])
  plot(1:10, x1, ylim=range(c(x1,x2,x3)), type='l',
       xlab='age', ylab='Relative difference', main='Numbers at age')
  lines(1:10, x2, col=2)
  lines(1:10, x3, col=3)
  abline(h=0)
  legend('bottomleft', legend=c('initial','middle', 'final'), col=1:3, lty=1)
  ## Recruits and F
  x1 <- (ww$F[,1]-aa$Fishing_mortalities)/aa$Fishing_mortalities
  plot(x1, ylim=range(c(0,x1)), xlab='year',
       ylab='Relative difference', type='l', main='Fishing effort');
  abline(h=0)
  ## Check selex
  x1 <- colMeans(ww$selAA[[1]]-aa$Fishery_selectivity)
  x2 <- (ww$selAA[[2]][1,]-aa$Survey_1_selectivity)
  x3 <- (ww$selAA[[3]][1,]-aa$Survey_2_selectivity)
  x4 <- (ww$selAA[[4]][1,]-aa$Survey_3_selectivity)
  x5 <- (ww$selAA[[5]][1,]-c(1,rep(0,9)))
  x6 <- (ww$selAA[[6]][1,]-c(0,1,rep(0,8)))
  x7 <- (ww$selAA[[7]][1,]-aa$Survey_6_selectivity)
  x8 <- (ww$selAA[[1]][1,]-aa$Fishery_selectivity[1,])
  x9 <- (ww$selAA[[1]][52,]-aa$Fishery_selectivity[52,])
  plot(1:10, x1, ylim=range(c(x1,x2,x3,x4,x5,x6,x7,x8,x9)),
       ylab='Absolute difference', type='l', xlab='age', main='Selectivity');
  abline(h=0)
  lines(1:10, x2, col=2)
  lines(1:10, x3, col=3)
  lines(1:10, x4, col=4)
  lines(1:10, x5, col=5)
  lines(1:10, x6, col=6)
  lines(1:10, x7, col=7)
  lines(1:10, x8, col=8)
  lines(1:10, x9, col=9)
  legend('topright', legend=c('fsh mean', 'svy1', 'svy2', 'svy3',
                            'svy4', 'svy5', 'svy6', 'fsh y=1',
                            'fsh y=52'), col=1:9, ncol=2, lty=1)
  x1 <- (ww$SSB-aa$Expected_spawning_biomass*1e6)/(aa$Expected_spawning_biomass*1e6)
  x2 <- (ww$pred_catch[,1]-aa$Expected_total_catch)/(aa$Expected_total_catch)
  abline(h=0)
  plot(1:52, x1, ylim=range(c(x1,x2,0)),ylab='Relative difference',
       type='l', xlab='year', main='SSB and catch');
  lines(1:52, x2,  col=2)
  legend('topleft', legend=c('SSB', 'Catch'), col=1:2, lty=1)
  abline(h=0)
  ## Get the indices right
  x1 <- (ww$pred_indices[,1]-1e6*aa$Expected_survey_1_index[2,])/(1e6*aa$Expected_survey_1_index[2,])
  x2 <- (ww$pred_indices[,2]-1e6*aa$Expected_survey_2_index)/(1e6*aa$Expected_survey_2_index)
  x3 <- (ww$pred_indices[,3]-1e6*aa$Expected_survey_3_index)/(1e6*aa$Expected_survey_3_index)
  x4 <- (ww$pred_indices[,4]-1e6*aa$Expected_survey_4_index)/(1e6*aa$Expected_survey_4_index)
  x5 <- (ww$pred_indices[,5]-1e6*aa$Expected_survey_5_index)/(1e6*aa$Expected_survey_5_index)
  x6 <-
    (ww$pred_indices[,6]-1e6*aa$Expected_survey_6_index)/(1e6*aa$Expected_survey_6_index)
  abline(h=0)
  matplot(cbind(x1,x2,x3,x4,x5,x6), ylab='Relative difference',
          xlab='year', main='Expected survey indices')
  x1 <- ww$q[,1]-aa$Survey_1_q[-1]
  x2 <- ww$q[,2]-aa$Survey_2_q
  x3 <- ww$q[,3]-aa$Survey_3_q
  x4 <- ww$q[,4]-aa$Survey_4_q
  x5 <- ww$q[,5]-aa$Survey_5_q
  x6 <- ww$q[,6]-aa$Survey_6_q
  matplot(cbind(x1,x2,x3,x4,x5,x6), ylab='Absolute differneces',
          main='Catchabilities',
          xlab='year')
  abline(h=0)
}

