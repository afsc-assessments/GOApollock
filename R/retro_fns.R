
peel_pars <- function(pars, peel){
  p <- pars
  stopifnot(peel>=0)
  if(any(sapply(p, function(x) !is.null(dim(x))))) stop("matrix or arrays break retro code")
  endyr <- dat$endyr-peel
  yrs <- dat$styr:endyr
  nyrs <- length(p$dev_log_recruit)
  ind <- 1:(nyrs-peel)
  x <- c("dev_log_recruit", "slp1_fsh_dev", "inf1_fsh_dev",
         "slp2_fsh_dev", "inf2_fsh_dev", "dev_log_F", "log_q1_dev",
         "log_q2_dev", "log_q3_dev")
  for(i in x) p[[i]] <- p[[i]][ind]
  if(any(sapply(p, NROW)>length(ind)))
    stop("Some pars too long in peel ",peel)
  if(peel==0) stopifnot(all.equal(p,pars))
  return(p)
}

peel_map <- function(map, pars, yrs){
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
  d$srvyrs6 <- d$srvyrs6[i3]
  d$nyrs_srv6 <- length(d$srvyrs6)
  d$indxsurv6 <- d$indxsurv6[i3]
  d$indxsurv_log_sd6 <- d$indxsurv_log_sd6[i3]
  d$q6_rwlk_sd <- d$q6_rwlk_sd[i0]
  d$yrfrct_srv6 <- d$yrfrct_srv6[i0]
  i4 <- which(d$srv_acyrs6<=endyr)
  d$srv_acyrs6 <- d$srv_acyrs6[i4]
  d$nyrsac_srv6 <- length(d$srv_acyrs6)
  d$multN_srv6 <- d$multN_srv6[i4]
  d$ac_yng_srv6 <- d$ac_yng_srv6[i4]
  d$ac_old_srv6 <- d$ac_old_srv6[i4]
  i5 <- which(d$srv_lenyrs6<=endyr)
  d$srv_lenyrs6 <- d$srv_lenyrs6[i5]
  d$nyrslen_srv6 <- length(d$srv_lenyrs6)
  d$multNlen_srv6 <- d$multNlen_srv6[i5]
  d$srvp6 <- d$srvp6[i4,]
  d$srvlenp6 <- d$srvlenp6[i5,, drop=FALSE]
  d$wt_srv6 <- d$wt_srv6[i0,]
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
  if(peel==0) stopifnot(all.equal(d,dat))
  ## this breaks makeadfun for some reason if not done??
##   for(i in 1:length(dat)){
## ##     browser()
##     print(cbind(i=i, names(dat)[i], class(dat[[i]])[1], class(d[[i]])[1]))
##     ## if(class(dat[[i]])[1]=='integer'){
##     ##   message("converting ", names(dat)[i])
##     ##   class(d[[i]]) <- 'integer'
##     ## }
##     class(d[[i]])[1] <- class(dat[[i]])[1]
##     print(cbind(i=i, names(dat)[i], class(dat[[i]])[1], class(d[[i]])[1]))
##   }
  ## attributes(d) <- attributes(dat)
  ## attributes(d) <- attributes(dat)
  ## attributes(d)$checks.passed <- NULL
  return(d)
}

fit_peel <- function(obj, peel, getsd=FALSE, ...){
  control <- list(eval.max=10000, iter.max=10000, trace=0)
  stopifnot(peel>=0)
  dat2 <- peel_data(dat=obj$env$data, peel)
  attributes(dat2) <- attributes(obj$env$data)
  attributes(dat2)$check.passed <- NULL
  yrs <- dat$styr:dat$endyr
  ## for(ii in names(map2)){
  ##   if(length(pars2[[ii]]) !=    length(map2[[ii]]))
  ##     stop("wrong length in map for ",ii)
  ## }
  pars2 <- peel_pars(pars=obj$env$parList(), peel)
  map2 <- peel_map(map=obj$env$map, pars2, yrs)
  obj2 <- TMB::MakeADFun(data=dat2, parameters=pars2, map=map2,
                         random=obj$random, DLL='goa_pk_tmb', silent=TRUE)
  ## print(obj$fn())
  ## print(obj2$fn())
  ## print(obj$gr())
  ## print(obj2$gr())
  lwr <- get_bounds(obj2)$lwr
  upr <- get_bounds(obj2)$upr
  message("Starting optimization for peel=",peel)
  opt <- TMBhelper::fit_tmb(obj2, lower=lwr, loopnum=3,
                            upper=upr, getsd=getsd, control=control,
                            newtonsteps=1)
  opt$rep <- obj2$report()
  return(opt)
}
