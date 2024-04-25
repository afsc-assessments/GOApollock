
#' Simulate a data set from a given rep file and input data
#'
#' @param datlist A data list as returned by
#'   \link{read_dat}. Observed data overwritten with simulated
#'   data.
#' @param replist A replist as returned by
#'   \link{read_pk_rep}. Uses expected values as truth.
#' @param fileout The new dat file name to be written. If NULL no
#'   file is written.
#' @param path Optional path for fileout
#' @param type Either "model" which uses the model expectation as
#'   the means, or "data" which uses the observed data means.
#' @return An invisible list is returned, and optionally a file
#'   written. The age compositions are returned as proportions
#'   because the model assumes this.
#' @export
#'
#' @details The expected values from either the dat or report
#'   file are used as the "truth" depending on \code{type}. The
#'   "model" routine is similar to the approach used in SS.Data
#'   simulated either as lognormal (indices) or multinomial
#'   (compositions) given the CV and sample sizes given in the
#'   datlist.  This implicitly assumes they have been weighted as
#'   desired. The "data" approach is akin to a non-parametric
#'   bootstrap. Sample sizes of 0 result in no observations, and
#'   fractions will be rounded internally by rmultinom. Currently
#'   the total fishery catch is not resampled, but the age
#'   compositions are.
sim_dat <- function(datlist, replist, fileout=NULL, path=NULL,
                    type=c('model','data')){
  cv2se <- function(cv)  sqrt(log(cv^2+1))
  myrmultinom <- function(N, true){
    stopifnot(length(N)==nrow(true))
    if(any(N<1 & N>0)) message("0<N<1 detected, using N=1")
    N[N<1 & N>0] <- 1
    sim <- t(sapply(1:length(N), function(i){
      tmp <- rmultinom(1, size=N[i], prob=true[i,])
      ## sum(tmp) not always N b/c rmultinom rounds it internally
      ## when its a real number
      if(sum(tmp)>0) tmp <- tmp/sum(tmp)
      tmp
    }))
    sim
  }
  testind <- function(x,y,z)   stopifnot(all(sum(x)==length(y), sum(x)==length(z)))
  testac <- function(x,y,z)   stopifnot(all(sum(x)==length(y), sum(x)==nrow(z)))
  ## If the retro option is used then the expected value
  ## vectors/matrices will be shorter than the data ones. So keep
  ## track of that and only resample those with
  type <- match.arg(type)
  if(type=='model') model <- TRUE else model <- FALSE
  d <- datlist; r <- replist
  endyr <- tail(r$years,1)

### fishery
  ## catch
  se <- cv2se(d$cattot_log_sd)
  if(model){
    true <- r$Expected_total_catch
  } else {
    true <- d$cattot
  }
  testind(length(d$styr:d$endyr), se, true)
  d$cattot <- rlnorm(n=length(se), meanlog=log(true), sdlog=se)

  ## age comps
  indr <- (r$years %in% d$fshyrs)
  ind <- indr & (r$years <=endyr)
  N <- d$multN_fsh
  if(model){
    true <- r$Fishery_expected_age_composition[ind,]
  } else {
    true <- d$catp
  }
  testac(ind, N, true)
  d$catp <- myrmultinom(N, true)


### Survey 1
  ## index
  ind <- r$years %in% d$srvyrs1
  se <- cv2se(d$indxsurv_log_sd1)
  if(model){
    true <- r$Survey_1_expected_index[ind] ## BS survey
  } else {
    true <- d$indxsurv1
  }
  testind(ind, se, true)
  d$indxsurv1 <- rlnorm(n=length(se), meanlog=log(true), sdlog=se)
  ## age comps
  ind <- r$years %in% d$srv_acyrs1
  N <- d$multN_srv1
  if(model){
    true <- r$Survey_1_expected_age_composition[ind,]
  } else {
    true <- d$srvp1
  }
  testac(ind,N,true)
  d$srvp1 <- myrmultinom(N, true)
  ## len comps
  ind <- r$years %in% d$srv_lenyrs1
  N <- d$multNlen_srv1
  if(model){
    true <- r$Survey_1_expected_length_composition[ind,]
  } else {
    true <- d$srvlenp1
  }
  testac(ind,N,true)
  d$srvlenp1 <- myrmultinom(N, true)

### Survey 2
  ## index
  ind <- r$years %in% d$srvyrs2
  se <- cv2se(d$indxsurv_log_sd2)
  if(model){
    true <- r$Survey_2_expected_index[ind] ## BS survey
  } else {
    true <- d$indxsurv2
  }
  testind(ind, se, true)
  d$indxsurv2 <- rlnorm(n=length(se), meanlog=log(true), sdlog=se)
  ## age comps
  ind <- r$years %in% d$srv_acyrs2
  N <- d$multN_srv2
  if(model){
    true <- r$Survey_2_expected_age_composition[ind,]
  } else {
    true <- d$srvp2
  }
  testac(ind,N,true)
  d$srvp2 <- myrmultinom(N, true)
  ## len comps
  ind <- r$years %in% d$srv_lenyrs2
  N <- d$multNlen_srv2
  if(model){
    true <- r$Survey_2_expected_length_composition[ind,]
  } else {
    true <- d$srvlenp2
  }
  testac(ind,N,true)
  d$srvlenp2 <- myrmultinom(N, true)

### Survey 3
  ## index
  ind <- r$years %in% d$srvyrs3
  se <- cv2se(d$indxsurv_log_sd3)
  if(model){
    true <- r$Survey_3_expected_index[ind] ## BS survey
  } else {
    true <- d$indxsurv3
  }
  testind(ind, se, true)
  d$indxsurv3 <- rlnorm(n=length(se), meanlog=log(true), sdlog=se)
  ## age comps
  ind <- r$years %in% d$srv_acyrs3
  N <- d$multN_srv3
  if(model){
    true <- r$Survey_3_expected_age_composition[ind,]
  } else {
    true <- d$srvp3
  }
  testac(ind,N,true)
  d$srvp3 <- myrmultinom(N, true)
  ## len comps
  ind <- r$years %in% d$srv_lenyrs3
  N <- d$multNlen_srv3
  if(model){
    true <- r$Survey_3_expected_length_composition[ind,]
  } else {
    true <- d$srvlenp3
  }
  testac(ind,N,true)
  d$srvlenp3 <- myrmultinom(N, true)
### Survey 4
  ## index
  ind <- r$years %in% d$srvyrs4
  se <- cv2se(d$indxsurv_log_sd4)
  if(model){
    true <- r$Survey_4_expected_index[ind] ## BS survey
  } else {
    true <- d$indxsurv4
  }
  testind(ind, se, true)
  d$indxsurv4 <- rlnorm(n=length(se), meanlog=log(true), sdlog=se)
### Survey 5
  ## index
  ind <- r$years %in% d$srvyrs5
  se <- cv2se(d$indxsurv_log_sd5)
  if(model){
    true <- r$Survey_5_expected_index[ind] ## BS survey
  } else {
    true <- d$indxsurv5
  }
  testind(ind, se, true)
  d$indxsurv5 <- rlnorm(n=length(se), meanlog=log(true), sdlog=se)

### Survey 6
  ## index
  ind <- r$years %in% d$srvyrs6
  se <- cv2se(d$indxsurv_log_sd6)
  if(model){
    true <- r$Survey_6_expected_index[ind] ## BS survey
  } else {
    true <- d$indxsurv6
  }
  testind(ind, se, true)
  d$indxsurv6 <- rlnorm(n=length(se), meanlog=log(true), sdlog=se)
  ## age comps
  ind <- r$years %in% d$srv_acyrs6
  N <- d$multN_srv6
  if(model){
    true <- r$Survey_6_expected_age_composition[ind,]
  } else {
    true <- d$srvp6
  }
  testac(ind,N,true)
  d$srvp6 <- myrmultinom(N, true)
  ## len comps
  ind <- r$years %in% d$srv_lenyrs6
  N <- d$multNlen_srv6
  stopifnot(sum(ind)==length(N))
  true <- r$Survey_6_expected_length_composition[ind,]
  d$srvlenp6 <- myrmultinom(N, true)

  ## plot(r$year[ind], true, cex=2, ylim=c(0,.5))
  ## trash <- sapply(1:500, function(i){
  ##   s <- rlnorm(n=length(se), meanlog=log(true), sdlog=se)
  ##   points(r$year[ind], s, col=rgb(1,0,0,.5))})
  ## points(r$year[ind], true, cex=2, pch=16)

  if(!is.null(fileout)) write_dat(d, fileout=fileout, path=path)
  return(invisible(d))
}

## ## some ugly code to test the simulator seems to work
## replist <- readRDS('C:/Users/cole.monnahan/Work/assessments/GOA_pollock/2022/results/repfile.RDS')
## datlist <- readRDS('C:/Users/cole.monnahan/Work/assessments/GOA_pollock/2022/results/datfile.RDS')
## library(ggplot2)
## library(dplyr)
## Nreps <- 1000
## x <- sim_dat(datlist, replist, type='model')
## x <- sim_dat(datlist, replist, type='data')
## rowSums(x$catp)
## rowSums(x$srvlenp6)
## write_dat(x, fileout='test.dat', path=getwd())
## out1 <- lapply(1:Nreps, function(x) sim_dat(datlist, replist, type='model'))
## out2 <- lapply(1:Nreps, function(x) sim_dat(datlist, replist,
##                                            type='data'))
## etac <- function(x,i){
##   d0 <- data.frame(rep=i, source='fishery', age=1:10, fish=as.numeric(tail(x$catp,1)))
##   d1 <- data.frame(rep=i, source='survey1', age=1:10, fish=as.numeric(tail(x$srvp1,1)))
##   d2 <- data.frame(rep=i, source='survey2', age=1:10, fish=as.numeric(tail(x$srvp2,1)))
##   d3 <- data.frame(rep=i, source='survey3', age=1:10, fish=as.numeric(tail(x$srvp3,1)))
##   d6 <- data.frame(rep=i, source='survey6', age=1:10, fish=as.numeric(tail(x$srvp6,1)))
##   rbind(d0,d1,d2,d3,d6) %>%
##     group_by(source) %>% mutate(proportion=fish/sum(fish)) %>% ungroup
## }
## ac2 <- lapply(1:Nreps, function(i) getac(out2[[i]], i)) %>%
##   bind_rows %>% mutate(type='data')
## acmean2 <- group_by(ac2, source, age) %>%
##   summarize(proportion=mean(proportion), .groups='drop') %>%
##  mutate(type='data')
## ac1 <- lapply(1:Nreps, function(i) getac(out1[[i]], i)) %>%
##   bind_rows %>% mutate(type='model')
## acmean1 <- group_by(ac1, source, age) %>%
##   summarize(proportion=mean(proportion), .groups='drop') %>%
##  mutate(type='model')
## acmean <- bind_rows(acmean1, acmean2)
## ac <- bind_rows(ac1, ac2)
## acobs <- rbind(getac(datlist, 0) %>% mutate(type='data'),
##                getac(datlist, 0) %>% mutate(type='model'))
## g1 <- ggplot(ac, aes(age, proportion)) +
##   geom_jitter(width=.25, height=.01, alpha=.1, size=.1) +
##   facet_grid(source~type) +
##   geom_line(data=acobs, color=2, lwd=2)+
##   geom_line(data=acmean, color=3, lwd=2)
## getind <- function(x,i){
##   d0 <- data.frame(rep=i, source='fishery', year=x$styr:x$endyr,index=x$cattot/1e6)
##   d1 <- data.frame(rep=i, source='survey1', year=x$srvyrs1, index=x$indxsurv1)
##   d2 <- data.frame(rep=i, source='survey2', year=x$srvyrs2, index=x$indxsurv2)
##   d3 <- data.frame(rep=i, source='survey3', year=x$srvyrs3, index=x$indxsurv3)
##   d4 <- data.frame(rep=i, source='survey4', year=x$srvyrs4, index=x$indxsurv4)
##   d5 <- data.frame(rep=i, source='survey5', year=x$srvyrs5, index=x$indxsurv5)
##   d6 <- data.frame(rep=i, source='survey6', year=x$srvyrs6, index=x$indxsurv6)
##   rbind(d0,d1,d2,d3,d4,d5,d6)
## }
## ind2 <- lapply(1:Nreps, function(i) getind(out2[[i]], i)) %>%
##   bind_rows %>% mutate(type='data')
## indmean2 <- group_by(ind2, source, year) %>%
##   summarize(index=mean(index), .groups='drop') %>%
##  mutate(type='data')
## ind1 <- lapply(1:Nreps, function(i) getind(out1[[i]], i)) %>%
##   bind_rows %>% mutate(type='model')
## indmean1 <- group_by(ind1, source, year) %>%
##   summarize(index=mean(index), .groups='drop') %>%
##  mutate(type='model')
## indmean <- bind_rows(indmean1, indmean2)
## ind <- bind_rows(ind1, ind2)
## indobs <- rbind(getind(datlist, 0) %>% mutate(type='data'),
##                 getind(datlist, 0) %>% mutate(type='model'))
## g2 <- ggplot(ind, aes(year, index)) +
##   geom_jitter(width=.25, height=.01, alpha=.1, size=.1) +
##   facet_grid(source~type, scales='free_y') +
##   geom_line(data=indobs, color=2, lwd=2)+
##   geom_line(data=indmean, color=3, lwd=2) + scale_y_log10()
## g1
## g2
## ## cowplot::plot_grid(g1,g2, ncol=2)


#' Write an R data list to file for reading into ADMB
#'
#' @param datlist An object as read in by \link{read_dat}
#' @param fileout The new filename to be written
#' @param path Optional path if not in the working directory
#' @export
#' @return Nothing, new file written to disk
write_dat <- function(datlist, fileout, path=NULL){
  if(!is.null(path)){
    oldwd <- getwd()
    on.exit(setwd(oldwd))
    setwd(path)
  }
  trash <- lapply(1:length(datlist), function(x){
    ## write name then data
    write.table(paste0("# ",names(datlist)[x]), fileout, append=x>TRUE,
                row.names=FALSE, col.names=FALSE, quote=FALSE)
    if(!is.null(nrow(datlist[[x]]))){
      ## matrices we don't want to transpose
      y <- datlist[[x]]
    } else {
      y <- t(as.data.frame(datlist[x]))
    }
    write.table(y, fileout, append= T,
                sep=' ', quote = F, col.names=FALSE, row.names=FALSE)
  })
}

#' Read text file input for ADMB model
#'
#' @param filename The filename to be read in
#' @param path Path to the file if not in working directory
#' @export
read_dat <- function(filename,path=NULL){
  .Deprecated('read_pk_dat')
  read_pk_dat(filename,path)
}



#' Read text file input for ADMB model
#'
#' @param filename The filename to be read in
#' @param path Path to the file if not in working directory
#' @export
#' @return A named list of all elements with corresponding names
#'   to the ADMB model
read_pk_dat <- function(filename, path=NULL, writedat=FALSE){
  if(!is.null(path)){
    oldwd <- getwd()
    on.exit(setwd(oldwd))
    setwd(path)
  }
  if(!file.exists(filename)) stop("File ", filename, " does not exist")
  char.lines <- suppressMessages(readLines(filename, warn=FALSE))
  com.ind <- which(substring(char.lines, 1, 1) %in% c(" ", "#"))
  blank.ind <- which(substring(char.lines, 1, 1) == "\t")
  com.ind <- sort(c(com.ind,blank.ind))
  dat.start <- com.ind[c(which(diff(com.ind) > 1), length(com.ind))]
  comments <- char.lines[dat.start]
  ## dat <- char.lines[-com.ind]
  ## ## drop anything after # on that line
  ## dat <- lapply(strsplit(dat, split='#'), function(x) x[1])
  d <- list()
  ind <- 0
  si <- function(ind, n, nrow=NULL, ncol=NULL){
    s <- tryCatch(scan(filename, quiet = T, what = integer(), comment.char="#",
              skip=dat.start[ind <<- ind+1], n=n),
              error=function(e) 'error')
    if(is.character(s)) {print(d); stop("failed to read integer")}
    if(is.null(nrow)) return(s)
    return(matrix(s, nrow=nrow, ncol=ncol, byrow=TRUE))
  }
  sn <- function(ind, n, nrow=NULL, ncol=NULL){
    s <- scan(filename, quiet = T, what = numeric(), comment.char="#",
              skip=dat.start[ind <<- ind+1], n=n)
    if(is.null(nrow)) return(s)
    return(matrix(s, nrow=nrow, ncol=ncol, byrow=TRUE))
  }
  if(writedat){
    d$styr <- si(ind, n=1)
    d$endyr <- si(ind,n=1)
  } else {
    tmp <- si(ind , n = 2)
    d$styr <- tmp[1]; d$endyr <- tmp[2]
  }
  nyrs <- d$endyr-d$styr+1
  nyrs2 <- nyrs-1
  if(writedat){
    d$rcrage <- si(ind,n=1)
    d$trmage <- si(ind,n=1)
  } else {
    tmp <- si(ind=ind, n = 2)
    d$rcrage <- tmp[1]; d$trmage <- tmp[2]
  }
  maxage <- d$trmage
  if(writedat){
    d$nbins1 <- si(ind,n=1)
    d$nbins2 <- si(ind,n=1)
    d$nbins3 <- si(ind,n=1)
  } else {
    tmp <- si(ind=ind, n = 3)
    d$nbins1 <- tmp[1]; d$nbins2 <- tmp[2]; d$nbins3 <- tmp[3]
  }
  d$cattot <- sn( ind=ind, n = nyrs)
  d$cattot_log_sd <- sn( ind=ind, n = nyrs)
  d$nyrs_fsh <- si(ind=ind, n =1)
  d$fshyrs <- sn( ind=ind, n = d$nyrs_fsh)
  d$multN_fsh <- sn( ind=ind, n = d$nyrs_fsh)
  d$ac_yng_fsh <- si(ind=ind, n = d$nyrs_fsh)
  d$ac_old_fsh <- si(ind=ind, n = d$nyrs_fsh)
  d$nyrslen_fsh <- si(ind=ind, n =1)
  d$fshlenyrs <- sn( ind=ind, n = d$nyrslen_fsh)
  d$multNlen_fsh <- sn( ind=ind, n = d$nyrslen_fsh)
  d$rwlk_sd <- sn( ind=ind, n = nyrs2)
  d$catp <- sn( ind=ind, n = d$nyrs_fsh*maxage, nrow=d$nyrs_fsh, ncol=maxage)
  d$lenp <- sn( ind=ind, n = d$nyrslen_fsh*d$nbins1, nrow=d$nyrslen_fsh, ncol=d$nbins1)
  d$wt_fsh <- sn( ind=ind, n = nyrs*maxage, nrow=nyrs, ncol=maxage)
  ## The used survey 1
  d$nyrs_srv1 <- si(ind=ind, n =1)
  d$srvyrs1 <- si(ind=ind, n =d$nyrs_srv1)
  d$indxsurv1 <- sn( ind=ind, n =d$nyrs_srv1)
  d$indxsurv_log_sd1 <- sn( ind=ind, n =d$nyrs_srv1)
  d$q1_rwlk_sd <- sn( ind=ind, n = nyrs2)
  d$yrfrct_srv1 <- sn( ind=ind, n = nyrs)
  d$nyrsac_srv1 <- si(ind=ind, n =1)
  d$srv_acyrs1 <- sn( ind=ind, n = d$nyrsac_srv1)
  d$multN_srv1 <- sn( ind=ind, n =d$nyrsac_srv1)
  d$ac_yng_srv1 <- tmp <- si(ind=ind, n = d$nyrsac_srv1)
  d$ac_old_srv1 <- tmp <- si(ind=ind, n = d$nyrsac_srv1)
  d$nyrslen_srv1 <- si(ind=ind, n =1)
  d$srv_lenyrs1 <- sn( ind=ind, n = d$nyrslen_srv1)
  d$multNlen_srv1 <- sn( ind=ind, n =d$nyrslen_srv1)
  d$srvp1 <- sn(ind=ind, n =d$nyrsac_srv1*maxage,
                nrow=d$nyrsac_srv1, ncol=maxage)
  d$srvlenp1 <- sn(ind=ind, n =d$nyrslen_srv1*d$nbins3,
                   nrow=d$nyrslen_srv1, ncol=d$nbins3)
  d$wt_srv1 <- sn( ind=ind, n = nyrs*maxage, nrow=nyrs, ncol=maxage)
  ## Survey 2
  d$nyrs_srv2 <- si(ind=ind, n =1)
  d$srvyrs2 <- si(ind=ind, n =d$nyrs_srv2)
  d$indxsurv2 <- sn( ind=ind, n =d$nyrs_srv2)
  d$indxsurv_log_sd2 <- sn( ind=ind, n =d$nyrs_srv2)
  d$q2_rwlk_sd <- sn( ind=ind, n = nyrs2)
  d$yrfrct_srv2 <- sn( ind=ind, n = nyrs)
  d$nyrsac_srv2 <- si(ind=ind, n =1)
  d$srv_acyrs2 <- sn( ind=ind, n = d$nyrsac_srv2)
  d$multN_srv2 <- sn( ind=ind, n =d$nyrsac_srv2)
  d$ac_yng_srv2 <- si(ind=ind, n = d$nyrsac_srv2)
  d$ac_old_srv2 <- si(ind=ind, n = d$nyrsac_srv2)
  d$nyrslen_srv2 <- si(ind=ind, n =1)
  d$srv_lenyrs2 <- sn( ind=ind, n = d$nyrslen_srv2)
  d$multNlen_srv2 <- sn( ind=ind, n =d$nyrslen_srv2)
  d$srvp2 <- sn(ind=ind, n =d$nyrsac_srv2*maxage,
                nrow=d$nyrsac_srv2, ncol=maxage)
  d$srvlenp2 <- sn(ind=ind, n =d$nyrslen_srv2*d$nbins3,
                   nrow=d$nyrslen_srv2, ncol=d$nbins3)
  d$wt_srv2 <- sn( ind=ind, n = nyrs*maxage, nrow=nyrs, ncol=maxage)
  ## Survey 3
  d$nyrs_srv3 <- si(ind=ind, n =1)
  d$srvyrs3 <- si(ind=ind, n =d$nyrs_srv3)
  d$indxsurv3 <- sn( ind=ind, n =d$nyrs_srv3)
  d$indxsurv_log_sd3 <- sn( ind=ind, n =d$nyrs_srv3)
  d$q3_rwlk_sd <- sn( ind=ind, n = nyrs2)
  d$yrfrct_srv3 <- sn( ind=ind, n = nyrs)
  d$nyrsac_srv3 <- si(ind=ind, n =1)
  d$srv_acyrs3 <- sn( ind=ind, n = d$nyrsac_srv3)
  d$multN_srv3 <- sn( ind=ind, n =d$nyrsac_srv3)
  ## tmp <- si(##                    ind=ind, n = 2*d$nyrsac_srv3)
  ## d$ac_yng_srv3 <- tmp[1:d$nyrs_srv3]; d$ac_old_srv3 <- tmp[-(1:d$nyrs_srv3)]
  d$nyrslen_srv3 <- si(ind=ind, n =1)
  d$srv_lenyrs3 <- sn( ind=ind, n = d$nyrslen_srv3)
  d$multNlen_srv3 <- sn( ind=ind, n =d$nyrslen_srv3)
  d$srvp3 <- sn(ind=ind, n =d$nyrsac_srv3*maxage, nrow=d$nyrsac_srv3, ncol=maxage)
  d$srvlenp3 <- sn(ind=ind, n =d$nyrslen_srv3*d$nbins3, nrow=d$nyrslen_srv3, ncol=d$nbins3)
  d$wt_srv3 <- sn( ind=ind, n = nyrs*maxage, nrow=nyrs, ncol=maxage)
  ## survey 4
  d$nyrs_srv4 <- si(ind=ind, n =1)
  d$srvyrs4 <- si(ind=ind, n =d$nyrs_srv4)
  d$indxsurv4 <- sn( ind=ind, n =d$nyrs_srv4)
  d$indxsurv_log_sd4 <- sn( ind=ind, n =d$nyrs_srv4)
  ## survey 5
  d$nyrs_srv5 <- si(ind=ind, n =1)
  d$srvyrs5 <- si(ind=ind, n =d$nyrs_srv5)
  d$indxsurv5 <- sn( ind=ind, n =d$nyrs_srv5)
  d$indxsurv_log_sd5 <- sn( ind=ind, n =d$nyrs_srv5)
  ## Survey 6
  d$nyrs_srv6 <- si(ind=ind, n =1)
  d$srvyrs6 <- si(ind=ind, n =d$nyrs_srv6)
  d$indxsurv6 <- sn( ind=ind, n =d$nyrs_srv6)
  d$indxsurv_log_sd6 <- sn( ind=ind, n =d$nyrs_srv6)
  ## d$q6_rwlk_sd <- sn( ##                    ind=ind, n = nyrs2)
  d$yrfrct_srv6 <- sn( ind=ind, n = nyrs)
  d$nyrsac_srv6 <- si(ind=ind, n =1)
  d$srv_acyrs6 <- sn( ind=ind, n = d$nyrsac_srv6)
  d$multN_srv6 <- sn( ind=ind, n =d$nyrsac_srv6)
  d$ac_yng_srv6 <-si(ind=ind, n = d$nyrsac_srv6)
  d$ac_old_srv6 <- si(ind=ind, n = d$nyrsac_srv6)
  d$nyrslen_srv6 <- si(ind=ind, n =1)
  d$srv_lenyrs6 <- sn( ind=ind, n = d$nyrslen_srv6)
  d$multNlen_srv6 <- sn( ind=ind, n =d$nyrslen_srv6)
  d$srvp6 <- sn(ind=ind, n =d$nyrsac_srv6*maxage,
                nrow=d$nyrsac_srv6, ncol=maxage)
  d$srvlenp6 <- sn(ind=ind, n =d$nyrslen_srv6*d$nbins2,
                   nrow=d$nyrslen_srv6, ncol=d$nbins2)
  d$wt_srv6 <- sn( ind=ind, n = nyrs*maxage, nrow=nyrs, ncol=maxage)
  d$age_trans <- sn( ind=ind, n = maxage*maxage, nrow=maxage, ncol=maxage)
  d$len_trans1 <- sn( ind=ind, n = d$nbins1*maxage, nrow=maxage, ncol=d$nbins1)
  d$len_trans2 <- sn( ind=ind, n = d$nbins2*maxage, nrow=maxage, ncol=d$nbins2)
  d$len_trans3 <- sn( ind=ind, n = d$nbins3*maxage, nrow=maxage, ncol=d$nbins3)
  ## the projection module inputs
  d$mat <- sn( ind=ind, n =maxage)
  d$Ftarget <- sn( ind=ind, n =5)
  d$B40 <- sn( ind=ind, n =1)
  d$log_mean_recr_proj <- sn( ind=ind, n =1)
  d$sigmasq_recr <- sn( ind=ind, n =1)
  d$check <- sn(ind=ind, n=1)
  if(d$check != -999) stop("Failed to read in dat file, final value=",d$check)
  return(d)
}

