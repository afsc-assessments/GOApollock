
#' Fit a sequence of models to peeled data set for retrospective
#' calculations
#'
#' @param obj A fitted model from \code{fit_pk}.
#' @param peels Vector of peels to fit, with 0 being the original
#'   data set
#' @param getsd Whether to run sdreport
#' @param parallel Whether to run peels in parallel or not (default).
#' @param ... Additional arguments to pass to \code{\link{fit_pk}}
#' @return A list of model fits of class 'pkfit'
#' @details This function fits a series of models to peels based
#'   off an original fit. The data, parameter, and map lists are
#'   modified based on the peel in a way that the original MLEs
#'   are preserved so that it runs faster. Pass an unfitted obj
#'   if this is undesired behavior.
#' @export
fit_retros <- function(fit, peels=0:7, getsd=TRUE, parallel=FALSE,...){
  if(class(fit)[1]!='pkfit')
    stop("fit argument is not a fitted model")
  if(parallel){
    message("Preparing parallel session..")
    if(!require(snowfall))
      stop("snowfall package required for parallel execution")
    sfInit(parallel=TRUE, cpus=parallel::detectCores()-1)
    sfLibrary(GOApollock)
    retros <- sfLapply(peels, function(i) GOApollock:::fit_peel(fit, peel=i, getsd=getsd, ...))
    sfStop()
  } else {
    retros <- lapply(peels, function(i) GOApollock:::fit_peel(fit, peel=i, getsd=getsd, ...))
  }
  return(retros)
}

#' Modify parameter list to match peel
peel_pars <- function(pars, dat,  peel){
  p <- pars
  stopifnot(peel>=0)
 ## if(any(sapply(p, function(x) !is.null(dim(x))))) stop("matrix
  ## or arrays break retro code")
  endyr <- dat$endyr-peel
  yrs <- dat$styr:endyr
  nyrs <- length(p$dev_log_recruit)
  nages <- length(dat$mat)
  ind <- 1:(nyrs-peel)
  x <- c("dev_log_recruit", "slp1_fsh_dev", "inf1_fsh_dev",
         "slp2_fsh_dev", "inf2_fsh_dev", "dev_log_F", "log_q1_dev",
         "log_q2_dev", "log_q3_dev", 'Ecov_exp')
  for(i in x){
    if(!is.null(p[[i]])) p[[i]] <- p[[i]][ind]
  }
  x <- 'selpars_re'
  for(i in x){
    if(!is.null(p[[i]])) p[[i]] <- p[[i]][ind,]
  }
  trash <- sapply(1:length(p), function(i) NROW(p[[i]]))
  trash <- which(trash >length(ind))
  if(length(trash)>0)
    stop("Some pars too long in peel ",names(p)[trash])
  if(peel==0) stopifnot(all.equal(p,pars))
  return(p)
}

#' Modify map to match peel
peel_map <- function(map, pars){
  ## tricky part is map elements may or may not be specified so
  ## do this generically
  m <- map
  p <- pars
  stopifnot(is.list(map))
  ## ind <- 1:length(yrs)
  ## nyrs <- length(ind)
  ## for(i in 1:length(m)){
  ##   if(length(m[[i]])>nyrs) m[[i]] <- m[[i]][ind]
  ## }
  ## if(any(sapply(m, NROW)>nyrs))
  ##   stop("Some pars too long in peel ",peel)
  ##  if(peel==0) stopifnot(all.equal(m,map))
  for(i in names(m)){
    #print(i)
    #if(i=='log_q2_dev') browser()
    m[[i]] <- m[[i]][1:length(p[[i]])]
  }
  return(m)
}

#' Modify data list to match peel
#'
peel_data <- function(dat, peel){
  stopifnot(peel>=0)
  stopifnot(is.list(dat))
  endyr <- as.integer(dat$endyr-peel)
  yrs <- dat$styr:endyr
  nyrs <- length(yrs)
  i0 <- yrs-dat$styr+1
  d <- dat
  d$endyr <- endyr
  d$cattot <- d$cattot[i0]
  d$cattot_log_sd <- d$cattot_log_sd[i0]
  if(!is.null(d$Ecov_obs)) d$Ecov_obs <- d$Ecov_obs[i0]
  i1 <- which(d$fshyrs<=endyr)
  d$fshyrs <- d$fshyrs[i1]                # age comp yrs
  d$nyrs_fsh <- length(i1)
  d$multN_fsh <- d$multN_fsh[i1]
  d$ac_yng_fsh <- d$ac_yng_fsh[i1]
  d$ac_old_fsh <- d$ac_old_fsh[i1]
  i2 <- which(d$fshlenyrs<=endyr)
  d$fshlenyrs <- d$fshlenyrs[i2]
  d$multNlen_fsh <- d$multNlen_fsh[i2]
  d$rwlk_sd <- d$rwlk_sd[i0[-nyrs]]
  d$catp <- d$catp[i1,]
  d$lenp <- d$lenp[i2,]
  d$wt_fsh <- d$wt_fsh[i0,]
  i3 <- which(d$srvyrs1 <=endyr)
  d$srvyrs1 <- d$srvyrs1[i3]
  d$nyrs_srv1 <- length(d$srvyrs1)
  d$indxsurv1 <- d$indxsurv1[i3]
  d$indxsurv_log_sd1 <- d$indxsurv_log_sd1[i3]
  d$q1_rwlk_sd <- d$q1_rwlk_sd[i0[-nyrs]]
  d$yrfrct_srv1 <- d$yrfrct_srv1[i0]
  i4 <- which(d$srv_acyrs1<=endyr)
  d$srv_acyrs1 <- d$srv_acyrs1[i4]
  d$nyrsac_srv1 <- length(d$srv_acyrs1)
  d$multN_srv1 <- d$multN_srv1[i4]
  d$ac_yng_srv1 <- d$ac_yng_srv1[i4]
  d$ac_old_srv1 <- d$ac_old_srv1[i4]
  i5 <- which(d$srv_lenyrs1<=endyr)
  d$srv_lenyrs1 <- d$srv_lenyrs1[i5]
  d$nyrslen_srv1 <- length(d$srv_lenyrs1)
  d$multNlen_srv1 <- d$multNlen_srv1[i5]
  d$srvp1 <- d$srvp1[i4,]
  d$srvlenp1 <- d$srvlenp1[i5,]
  d$wt_srv1 <- d$wt_srv1[i0,]
  i3 <- which(d$srvyrs2 <=endyr)
  d$srvyrs2 <- d$srvyrs2[i3]
  d$nyrs_srv2 <- length(d$srvyrs2)
  d$indxsurv2 <- d$indxsurv2[i3]
  d$indxsurv_log_sd2 <- d$indxsurv_log_sd2[i3]
  d$q2_rwlk_sd <- d$q2_rwlk_sd[i0[-nyrs]]
  d$yrfrct_srv2 <- d$yrfrct_srv2[i0]
  i4 <- which(d$srv_acyrs2<=endyr)
  d$srv_acyrs2 <- d$srv_acyrs2[i4]
  d$nyrsac_srv2 <- length(d$srv_acyrs2)
  d$multN_srv2 <- d$multN_srv2[i4]
  d$ac_yng_srv2 <- d$ac_yng_srv2[i4]
  d$ac_old_srv2 <- d$ac_old_srv2[i4]
  i5 <- which(d$srv_lenyrs2<=endyr)
  d$srv_lenyrs2 <- d$srv_lenyrs2[i5]
  d$nyrslen_srv2 <- length(d$srv_lenyrs2)
  d$multNlen_srv2 <- d$multNlen_srv2[i5]
  d$srvp2 <- d$srvp2[i4,]
  d$srvlenp2 <- d$srvlenp2[i5,]
  d$wt_srv2 <- d$wt_srv2[i0,]
  i3 <- which(d$srvyrs3 <=endyr)
  d$srvyrs3 <- d$srvyrs3[i3]
  d$nyrs_srv3 <- length(d$srvyrs3)
  d$indxsurv3 <- d$indxsurv3[i3]
  d$indxsurv_log_sd3 <- d$indxsurv_log_sd3[i3]
  d$q3_rwlk_sd <- d$q3_rwlk_sd[i0[-nyrs]]
  d$yrfrct_srv3 <- d$yrfrct_srv3[i0]
  i4 <- which(d$srv_acyrs3<=endyr)
  d$srv_acyrs3 <- d$srv_acyrs3[i4]
  d$nyrsac_srv3 <- length(d$srv_acyrs3)
  d$multN_srv3 <- d$multN_srv3[i4]
  d$ac_yng_srv3 <- d$ac_yng_srv3[i4]
  d$ac_old_srv3 <- d$ac_old_srv3[i4]
  i5 <- which(d$srv_lenyrs3<=endyr)
  d$srv_lenyrs3 <- d$srv_lenyrs3[i5]
  d$nyrslen_srv3 <- length(d$srv_lenyrs3)
  d$multNlen_srv3 <- d$multNlen_srv3[i5]
  d$srvp3 <- d$srvp3[i4,]
  d$srvlenp3 <- d$srvlenp3[i5,]
  d$wt_srv3 <- d$wt_srv3[i0,]
  i3 <- which(d$srvyrs6 <=endyr)
  if(d$srvyrs6[1] > endyr){
    warning('peeling to endyr=', endyr,
            ' leaves no survey 6 index data, so dummy values used w/ CV=0')
    d$srvyrs6 <- dat$styr
    d$nyrs_srv6 <- 1
    d$indxsurv6 <- 1
    d$indxsurv_log_sd6 <- 0
    d$srv_acyrs6 <- dat$styr
    d$nyrsac_srv6 <- 1
    d$multN_srv6 <- 0
    d$srvp6 <- d$srvp6[1,,drop=FALSE]
    d$ac_yng_srv6 <- 1
    d$ac_old_srv6 <- d$trmage
  } else {
    d$srvyrs6 <- d$srvyrs6[i3]
    d$nyrs_srv6 <- length(d$srvyrs6)
    d$indxsurv6 <- d$indxsurv6[i3]
    d$indxsurv_log_sd6 <- d$indxsurv_log_sd6[i3]
    i4 <- which(d$srv_acyrs6<=endyr)
    d$srv_acyrs6 <- d$srv_acyrs6[i4]
    d$srvp6 <- d$srvp6[i4,, drop=FALSE]
    d$nyrsac_srv6 <- length(d$srv_acyrs6)
    d$multN_srv6 <- d$multN_srv6[i4]
    d$ac_yng_srv6 <- d$ac_yng_srv6[i4]
    d$ac_old_srv6 <- d$ac_old_srv6[i4]
  }
  d$yrfrct_srv6 <- d$yrfrct_srv6[i0]
  if(d$srv_lenyrs6[1] > endyr){
    warning('peeling to endyr=', endyr,
            ' leaves no survey 6 length data, so dummy values used w/ Nsamp=0')
    d$srv_lenyrs6 <- dat$styr
    d$nyrslen_srv6 <- 1
    d$multNlen_srv6 <- 0
    d$srvlenp6 <- d$srvlenp6[1,, drop=FALSE]
  } else {
    i5 <- which(d$srv_lenyrs6<=endyr)
    d$srv_lenyrs6 <- d$srv_lenyrs6[i5]
    d$nyrslen_srv6 <- length(d$srv_lenyrs6)
    d$multNlen_srv6 <- d$multNlen_srv6[i5]
    d$srvp6 <- d$srvp6[i4,, drop=FALSE]
    d$srvlenp6 <- d$srvlenp6[i5,, drop=FALSE]
  }
  d$wt_srv6 <- d$wt_srv6[i0,, drop=FALSE]
  i3 <- which(d$srvyrs4 <=endyr)
  d$srvyrs4 <- d$srvyrs4[i3]
  d$nyrs_srv4 <- length(d$srvyrs4)
  d$indxsurv4 <- d$indxsurv4[i3]
  d$indxsurv_log_sd4 <- d$indxsurv_log_sd4[i3]
  i3 <- which(d$srvyrs5 <=endyr)
  d$srvyrs5 <- d$srvyrs5[i3]
  d$nyrs_srv5 <- length(d$srvyrs5)
  d$indxsurv5 <- d$indxsurv5[i3]
  d$indxsurv_log_sd5 <- d$indxsurv_log_sd5[i3]
  i4 <- which(d$Ecov_obs_year <= endyr)
  d$Ecov_obs_year <- d$Ecov_obs_year[i4]
  d$Ecov_obs <- d$Ecov_obs[i4]
  if(peel==0) stopifnot(all.equal(d,dat))
  return(d)
}

#' Internal wrapper function to fit a single peel
#'
fit_peel <- function(fit, peel, getsd=TRUE, ...){
  control <- list(eval.max=10000, iter.max=10000, trace=0)
  stopifnot(peel>=0)
  dat2 <- peel_data(dat=fit$obj$env$data, peel)
  attributes(dat2) <- attributes(fit$obj$env$data)
  attributes(dat2)$check.passed <- NULL
  pars2 <- peel_pars(pars=fit$obj$env$parList(), dat=dat2, peel)
  map2 <- peel_map(map=fit$obj$env$map, pars2)
  if(dat2$endyr<2013){
    warning("No survey 6 data so mapping off log_q6")
    map2$log_q6 <- factor(NA)
  }
  input2 <- list(version=paste0('peel',peel), path=fit$path,
                 modfile=fit$modfile,
                dat=dat2, pars=pars2, map=map2, random=fit$input$random)
  message("Starting optimization for peel=",peel)
  fit <- fit_pk(input=input2, getsd=getsd, control=control, ...)
  return(fit)
}


#' Plot estimates of recruits from retro lists
#' @param fits List of retro fits including peel=0
#' @param thisyear The current year
#' @param plot Whether to plot or return
#' @export
plot_retro_recruits <- function(fits, thisyear, plot=TRUE){
  stopifnot(is.pkfits(fits))
  stds <- get_std(fits)
  g <- stds %>% filter(name=='log_recruit' & year <= thisyear) %>%
    mutate(peel=as.numeric(gsub('peel','',version)),
           ayear=thisyear-peel, cohort=year-1,
           cohortf=paste(cohort,'cohort'),
           years_data=ayear-cohort,
           est=exp(est), lwr=exp(lwr), upr=exp(upr)) %>%
    filter(cohort>=thisyear-11)
  if(!plot) return(g)
  g <- g %>% # only makes sense so far back
    ggplot(aes(years_data, est, ymin=lwr, ymax=upr))+
    geom_pointrange(fatten=2) + facet_wrap('cohortf', scales='free_y')+
    scale_x_continuous(breaks=0:10)+
    labs(x='Years of data', y='Recruits (billions)')
}
