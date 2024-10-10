#' Read text file input for model. Excel sheet is saved as .txt and read into R with this function.
#'
#' @param filename The filename to be read in
#' @param path Path to the file if not in working directory
#' @return A named list of all elements with corresponding names
#'   to the ADMB model
read_dat <- function(filename, path=NULL, writedat=FALSE){
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
    if(is.character(s)) stop("Error reading integer on line ", ind, " of length ", n, ".\nList=", d)
    if(is.null(nrow)) return(s)
    return(matrix(s, nrow=nrow, ncol=ncol, byrow=TRUE))
  }
  sn <- function(ind, n, nrow=NULL, ncol=NULL){
    s <- tryCatch(scan(filename, quiet = T, what = numeric(), comment.char="#",
              skip=dat.start[ind <<- ind+1], n=n),
     error=function(e) 'error')
    if(is.character(s)) stop("Error reading numeric on line ", ind, " of length ", n, ".\nList=", d)
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


#' Extract sdreport data.frames from fitted models
#' @param fits A single model or list of models as returned by
#'   \code{\link{fit_pk}}
#' @param slot A variable to filter on, if NULL will return all
#'   variables
#' @param year.scale A scalar to get year2 right for offsetting
#'   plots downstream
#' @return A cleanly formatted data frame of sdreport objects,
#'   with columns 'name', 'est', 'se', 'year', 'lwr', 'upr', and
#'   'version'
#' @export
get_std <- function(fits, slot=NULL, year.scale=1) {
  ## single fit
  if(class(fits[[1]])[1]!='pkfit'){
    out <- fits$sd
    out$version <- factor(out$version)
  } else {
    ## list of fits
    labs <- lapply(fits, function(x) x$version)
    out <- lapply(fits, function(x) x$sd) %>% bind_rows
    ## preserve order in labels for downstream plots
    out$version <- factor(out$version, levels=labs)
  }
  if(!is.null(slot)) {
    out <- filter(out, name==slot)
    if(NROW(out)==0)
      stop("No rows found in sdreport for variable ", slot)
  }
  fs <- as.numeric(out$version)
  fs <- (fs-mean(fs))/length(unique(fs))
  out <- mutate(out, year2=year+fs/year.scale)
  return(out)
}

#' Extract report list from fitted models
#' @param fits A single model or list of models as returned by
#'   \code{\link{fit_pk}}
#' @param slot Which slot to extract and format into data.frame
#'   for easy plotting. If NULL (default) the whole report file
#'   is returned.
#' @return Either a list of reports or a formatted data.frame
#'   containing the specified slot variable.
#' @export
get_rep <- function(fits, slot=NULL) {
  ## single fit
  if(class(fits[[1]])[1]!='pkfit'){
    out <- fits$rep
  } else {
    ## list of fits
    out <- lapply(fits, function(x) x$rep)
  }
  if(!is.null(slot)){
    out <- mymelt(out, slot=slot)
    if(nrow(out)==0)
      warning("Slot", slot, "was not found in report file")
  }
  return(out)
}


#' Fit a GOA pollock model (BETA)
#' @param input Input list as returned by
#'   \code{prepare_pk_input}.
#' @param newtonsteps,control Optimizer settings
#' @param getsd Calculate and return sdreport?
#' @param do.fit Optimize or return obj? Used for testing.
#' @param use_bounds Whether to use bounds, slows optimization
#'   down but may be necessary. See \code{\link{get_bounds}}
#'   which is called internally.
#' @param save.sdrep Whether to return the sdreport object in the
#'   fitted model. This is rarely used and large so turned off by
#'   default. When returned it is named `sdrep`.
#' @param filename Character string giving a file name to save the fitted object
#'   as an RDS object. Defaults to NULL which indicates not to save it. If
#'   specified, it must end in .RDS. The file is written to folder given by
#'   \code{input$path}.
#' @param verbose Whether to print output (default) or suppress
#'   as much as possible.
#' @return A list object of class 'pkfit' which contains a
#'   "version" model name, rep, parList (MLE in list format), opt
#'   as returned by \code{TMBHelper::fit_tmb} but without the SD
#'   slot, std (formatted data frame) and sdrep if
#'   \code{getsd=TRUE}, and the obj.
#' @details This function is a beta version still and subject to change
#'   without warning.
#' @export
fit_pk <- function(input, getsd=TRUE, newtonsteps=1,
                   control=NULL, do.fit=TRUE,
                   use_bounds=FALSE, save.sdrep=FALSE,
                   filename=NULL, verbose=TRUE){
  cpp <- paste0(file.path(input$path, input$modfile),'.cpp')
  if(!file.exists(cpp)) stop("file does not exist: ", cpp)
  tryCatch(dyn.unload(dynlib(file.path(input$path, input$modfile))), error=function(e) 'error')
  suppressWarnings(suppressMessages(compile(cpp)))
  suppressWarnings(suppressMessages(dyn.load(dynlib(file.path(input$path, input$modfile)))))
  obj <- MakeADFun(data=input$dat, parameters=input$pars,
                   map=input$map, random=input$random,
                   DLL=input$modfile, silent=TRUE)
  if(!do.fit) return(obj)
  if(use_bounds){
    lwr <- get_bounds(obj)$lwr
    upr <- get_bounds(obj)$upr
  } else {
    lwr <- -Inf; upr <- Inf
  }
  if(is.null(control))
    control <- list(eval.max=10000, iter.max=10000, trace=0)
  if(!verbose) control$trace <- 0
  ## optimize and compare
  if(!require(TMBhelper)){
    stop('TMBhelper package required to fit models, install using:\ndevtools::install_github("kaskr/TMB_contrib_R/TMBhelper")')
  }
  if(verbose){
    opt <- TMBhelper::fit_tmb(obj, control=control,
                              newtonsteps=newtonsteps, getsd=FALSE,
                              lower=lwr, upper=upr)
  } else {
    suppressWarnings(suppressMessages(opt <- TMBhelper::fit_tmb(obj, control=control,
                              newtonsteps=newtonsteps, getsd=FALSE,
                              lower=lwr, upper=upr)))
  }
  rep <- c(version=input$version, obj$report())
  sdrep <- std <- NULL
  if(verbose) message("Finished optimization")
  if(getsd){
    sdrep <- sdreport(obj)
    std <- summary(sdrep)
    std <- data.frame(dimnames(std)[[1]], std)
    names(std) <- c('name', 'est', 'se')
    row.names(std) <- NULL
    std <- group_by(std, name) %>%
      mutate(year=1969+1:n(), lwr=est-1.96*se,
             upr=est+1.96*se, version=input$version) %>%
      ungroup
    if(verbose) message("Finished sdreport")
  }
  parList <- obj$env$parList()
  parnames <- names(obj$par)
  parnames <- as.vector((unlist(sapply(unique(parnames), function(x){
    temp <- parnames[parnames==x]
    if(length(temp)>1) paste0(temp,'[',1:length(temp),']') else temp
  }))))
  fit <- list(version=input$version, path=input$path,
              modfile=input$modfile, rep=rep, opt=opt, sd=std,
              obj=obj, parList=parList, input=input, parnames=parnames)
    if(save.sdrep) fit <- c(fit,sdrep=sdrep)
  class(fit) <- c('pkfit', 'list')
  if(verbose) print(fit)
  if(!is.null(filename)) {
    saveRDS(fit, file=paste0(input$path,'/', filename))
  }
  return(invisible(fit))
}

#' Update a data list with a new maximum age, summing age compositions, but truncating WAA, ageing error, length transition matrices, and maturity.
#'
#' @param dat A data list as returned by \code{read_dat}
#' @param maxage A new maximum age that is less than the one in the dat file.
#' @return An updated data input list that is reformatted to the new maximum age.
update_dat_maxage <- function(dat, maxage){
  message("Reducing maxage age from ", dat$trmage, " to ", maxage)
  i0 <- maxage:dat$trmage
  i1 <- 1:maxage
  d0 <- d <- dat
  d$trmage <- maxage
  ## sum across dropped columns
  d$catp[,maxage] <-rowSums(d$catp[,i0])
  d$catp <- d$catp[,i1]
#  stopifnot(0==sum(rowSums(d0$catp) - rowSums(d$catp)))
  d$wt_fsh <- d$wt_fsh[,i1]
  d$ac_old_fsh <- d$ac_old_fsh*0 + maxage
  d$srvp1[,maxage] <-rowSums(d$srvp1[,i0])
  d$srvp1 <- d$srvp1[,i1]
 # stopifnot(0==sum(rowSums(d0$srvp1) - rowSums(d$srvp1)))
  d$wt_srv1 <- d$wt_srv1[,i1]
  d$ac_old_srv1 <- d$ac_old_srv1*0 + maxage
  d$srvp2[,maxage] <-rowSums(d$srvp2[,i0])
  d$srvp2 <- d$srvp2[,i1]
#  stopifnot(0==sum(rowSums(d0$srvp2) - rowSums(d$srvp2)))
  d$wt_srv2 <- d$wt_srv2[,i1]
  d$ac_old_srv2 <- d$ac_old_srv2*0 + maxage
  d$srvp3[,maxage] <-rowSums(d$srvp3[,i0])
  d$srvp3 <- d$srvp3[,i1]
 # stopifnot(0==sum(rowSums(d0$srvp3) - rowSums(d$srvp3)))
  d$wt_srv3 <- d$wt_srv3[,i1]
  d$ac_old_srv3 <- d$ac_old_srv3*0 + maxage
  d$srvp6[,maxage] <-rowSums(d$srvp6[,i0])
  d$srvp6 <- d$srvp6[,i1]
 # stopifnot(0==sum(rowSums(d0$srvp6) - rowSums(d$srvp6)))
  d$wt_srv6 <- d$wt_srv6[,i1]
  d$ac_old_srv6 <- d$ac_old_srv6*0 + maxage

  d$age_trans <- d$age_trans[i1,i1]
  d$len_trans1 <- d$len_trans1[i1,]
  d$len_trans2 <- d$len_trans2[i1,]
  d$len_trans3 <- d$len_trans3[i1,]
  d$mat <- d$mat[i1]
  return(d)
}

#' Prepare inputs for the TMB GOA pollock model.
#' @param path Directory containing the model and dat
#'   files. Passed to \code{fit_pk}.
#' @param datfile Name of dat file to be read in.
#' @param version A character string for the model name, used for
#'   plotting downstream
#' @param random A character vector declaring random effect
#'   vectors to be integrated. Defaults to NULL (fully penalized
#'   ML)
#' @param modfile Model name assumed to be 'goa_pk' unless
#'   specified
#' @param maxage Specify a lower maximum age than in the input dat file if desired.
#' @return A list with the version, dat, pars, map and random
#'   which are used to build a TMB 'obj' in \code{\link{fit_pk}}.
#' @export
prepare_pk_input <- function(path, datfile, version='none',
                             random=NULL, modfile='goa_pk',
                             maxage=NULL){
  if(!dir.exists(path)) stop("directory does not exist: ",path)
  dat <- read_dat(filename=datfile, path=path)
  if(!is.null(maxage) && maxage<dat$trmage){
    dat <- update_dat_maxage(dat,maxage)
  }

  ## Prepare the parameter list based on the data
  pars <- prepare_par(dat)
  map <- prepare_map(pars)
  out <- list(version=version, path=path, modfile=modfile,
              dat=dat, pars=pars, map=map, random=random)
  return(out)
}

#' Internal fucntion to prepare parameter list from dat
#' list. Uses initial values from the 2022 model for scalars, and
#' vectors are all 0.
prepare_par <- function(dat){
  ## build a default starting list, based off MLE in 2022 but
  ## truncated a bit
  nyrs <- length(dat$styr:dat$endyr)
  nages <- length(dat$rcrage:dat$trmage)
  stopifnot(nyrs>30)
  stopifnot(nages>9)
  pars <- list(dev_log_initN = rep(0,nages),
               mean_log_recruit = 1.10,
               dev_log_recruit = rep(0,nyrs),
               sigmaR = 1.3,
               ## not sure what to do about this one for now,
               ## presumably modify it later when doing proj stuff?
               log_recr_proj = c(1.1032, 1.1032, 1.1032, 1.1032, 1.1032),
               log_slp1_fsh_mean = 0.77,
               inf1_fsh_mean = 3.74,
               log_slp2_fsh_mean = 0.93,
               inf2_fsh_mean = 9.70,
               slp1_fsh_dev = rep(0,nyrs),
               inf1_fsh_dev = rep(0,nyrs),
               slp2_fsh_dev = rep(0,nyrs),
               inf2_fsh_dev = rep(0,nyrs),
               log_slp2_srv1 = 0.53,
               inf2_srv1 = 9.80,
               log_slp1_srv2 = -0.46,
               inf1_srv2 = 4.07,
               log_slp2_srv2 = 1,
               inf2_srv2 = 20,
               log_slp1_srv3 = 0.46,
               inf1_srv3 = 4.37,
               log_slp1_srv6 = 4.9, inf1_srv6 = 0.5, log_slp2_srv6 = 0.24,
               inf2_srv6 = 7.87, mean_log_F = -1.97,
               dev_log_F = rep(0,nyrs),
               log_q1_mean = -0.53, log_q1_dev = rep(0,nyrs),
               log_q2_mean = -0.20, log_q2_dev = rep(0,nyrs), log_q3_mean = -1.54,
               log_q3_dev = rep(0,nyrs), log_q4 = -1.20,
               q4_pow = 0, log_q5 = -1.08, q5_pow = 0, log_q6 = -0.26,
               natMscalar = 1)
  return(pars)
}

#' Internal function to prepare map list from parameter list
prepare_map <- function(pars){
  parsoff <- c("dev_log_initN", "q4_pow", "q5_pow", "log_q2_dev", "slp2_fsh_dev",
               "inf2_fsh_dev", "log_slp2_srv2", "inf2_srv2", "log_slp1_srv6",
               "inf1_srv6", "sigmaR", "natMscalar", "log_recr_proj", "log_q1_mean",
               "log_q3_mean", "mean_log_F", "inf1_fsh_mean",
               "log_slp1_fsh_mean")
  map <- list()
  for(x in parsoff) map[[x]] <- factor(pars[[x]]*NA)
  return(map)
}




#' Calculate bounds for parameters
#' @param obj Object
#' @return list of lower (lwr) and upper (upr) vectors matching
#'   the object
#' @export
get_bounds <- function(obj){
  lwr <- obj$par-Inf
  upr <- obj$par+Inf
  upr['log_slp2_fsh_mean'] <- 5
  lwr['log_slp2_fsh_mean'] <- -5
  upr['inf2_fsh_mean'] <- 15
  lwr['inf2_fsh_mean'] <- 7
  upr['log_slp1_srv1'] <- 5
  lwr['log_slp1_srv1'] <- -5
  upr['inf1_srv1'] <- 10
  lwr['inf1_srv1'] <- 1
  upr['log_slp2_srv1'] <- 5
  lwr['log_slp2_srv1'] <- -2
  upr['inf2_srv1'] <- 12
  lwr['inf2_srv1'] <- 5
  upr['log_slp1_srv2'] <- 5
  lwr['log_slp1_srv2'] <- -5
  upr['inf1_srv2'] <- 50
  lwr['inf1_srv2'] <- 1
  upr['log_q2_mean'] <- 10
  lwr['log_q2_mean'] <- -10

  # - Non-parametric fish pars
  upr['sel_rho_c'] <- 10
  lwr['sel_rho_c'] <- -10

  upr['sel_rho_a'] <- 10
  lwr['sel_rho_a'] <- -10

  upr['sel_rho_y'] <- 10
  lwr['sel_rho_y'] <- -10

  ind <- which(names(upr)=='mean_sel')
  upr[ind] <- 10
  lwr[ind] <- -10

  ## upr['selpars_re'] <- 20
  ## lwr['selpars_re'] <- -20

  upr['ln_sel_sd'] <- 5
  lwr['ln_sel_sd'] <- -5

  # - CEATTLE bounds
  ## upr['inf1_fsh_dev'] <- 5
  ## upr['slp1_fsh_dev'] <- 5
  ## upr['dev_log_recruit'] <- 15
  ## upr['dev_log_initN'] <- 23
  ##  upr['dev_log_F'] <- 10


  ## lwr['inf1_fsh_dev'] <- -5
  ## lwr['slp1_fsh_dev'] <- -5
  ## lwr['dev_log_recruit'] <- -1000
  ## lwr['dev_log_initN'] <- -1000
  ## lwr['dev_log_F'] <- -1000
  lwr <- lwr[names(lwr ) %in% names(obj$par)]
  upr <- upr[names(upr ) %in% names(obj$par)]
  ## TODO need to ensure order is right here I think
  return(list(lwr=lwr, upr=upr))
}

#' Run a jitter analysis on a fitted model
#' @param fit A fitted model object from \code{fit_pk}.
#' @param njitter The number of jitter runs to do.
#' @param scalar How much to scale the values, where a new init value is generated by multiplying by U(1-scalar,1+scalar)
#' @param parallel Whether to run in parallel (default and recommended).
#' @param type Either "mle" which jitters around the MLE value or "par" where the initial values in obj$par are used.
#' @return A data.frame with columns for the jitter number, terminal SSB, parameter inits and estimations, and maximum gradient. It also prints to console
#' @export
run_jitter <- function(fit, njitter=20, scalar=.1, parallel=TRUE, type='par'){
  stopifnot(is.pkfit(fit))
  type <- match.arg(type, c("mle", 'par'))
  input <- fit$input
  if(type=='mle') par0 <- fit$opt$par
  if(type=='par') par0 <- fit$obj$par
  fit_jitter <- function(i){
    set.seed(i)
     dyn.load(dynlib(file.path(input$path, input$modfile))  )
     obj <- MakeADFun(data=input$dat, parameters=input$pars,
                      map=input$map, random=input$random,
                      DLL=input$modfile, silent=TRUE)
    if(i==0) newpar <- par0 else
    newpar <- par0*runif(n=length(par0), min=1-scalar, max=1+scalar)
    fit <- TMBhelper::fit_tmb(obj, startpar=newpar, loopnum=1, newtonsteps=0, getsd = FALSE,
                              control=list(trace=0))
    ssb <- tail(obj$rep()$Espawnbio,1)
    fit <- c(fit, jitter=i, ssb=ssb)
    fit$init <- newpar
    return(fit)
  }
  if(parallel){
    library(snowfall)
    cores <- parallel::detectCores()-1
    sfInit(cpus=cores, parallel=TRUE)
    sfLibrary(TMB)
    sfExportAll()
    message('Starting jitters in parallel...')
    jitterfits <- sfLapply(0:njitter, fit_jitter)
  } else {
    jitterfits <- lapply(0:njitter, fit_jitter)
  }
  pars <- lapply(jitterfits, function(x)
    data.frame(jitter=x$jitter, nll=as.numeric(x$objective),
               log_maxgrad=log10(x$max_gradient), ssb=x$ssb,
               init=x$init, version=fit$version,
               parnum=1:length(par0), parname=fit$parnames, value=x$par)) %>%
    bind_rows %>% as_tibble
  pars
  if(sum(pars$jitter==0 & pars$parnum==1)==1){
    ssbre <- (pars$ssb - pars$ssb[pars$jitter==0])/pars$ssb[pars$jitter==0]
    ssbre <- ssbre[pars$parnum==1]
    message(paste("Max absolute terminal SSB error relative to initial fit=", sprintf("%.3g", max(abs(ssbre)))))
  } else {
    warning("No jitter=0 found so RE calculations not done")
  }
  message(paste("Max absolute gradient=", sprintf("%.3g", max(10^pars$log_maxgrad))))
  return(pars)
}

#' Perform a self-test estimation for a fitted model.
#'
#' @param fit A fitted object
#' @param reps A vector of integers of reps to run
#' @param parallel Whether to run in parallel (default)
#'
#' @export
#' @return A data.frame containing parameter true, estimates,  absolute and relative errors, and maxgrad.
fit_self_test <- function(fit, reps, parallel=TRUE){
  fit_test <- function(fit, rep){
    set.seed(rep)
    input <- siminput <- fit$input
    dyn.load(dynlib(file.path(input$path, input$modfile))  )
    obj <- MakeADFun(data=input$dat, parameters=fit$parList,
                     map=input$map, random=input$random,
                     DLL=input$modfile, silent=TRUE)
    siminput$dat <- obj$simulate(complete = TRUE)
    attributes(siminput$dat)$check.passed <- NULL
    siminput$dat$catp[1,] <- 0
    tmp <- fit_pk(siminput, newtonsteps = 0, getsd=FALSE, do.fit = TRUE, verbose=FALSE, filename=FALSE)
    pars <- data.frame(rep=rep, year=NA, true=fit$opt$par, est=tmp$opt$par,
                       par=tmp$parnames, max_gradient=tmp$opt$max_gradient)
    ssb <- data.frame(rep=rep, year=fit$rep$years, true=fit$rep$Espawnbio, est=tmp$rep$Espawnbio,
                       par='ssb', max_gradient=tmp$opt$max_gradient)
    rec <- data.frame(rep=rep, year=fit$rep$years, true=fit$rep$recruit, est=tmp$rep$recruit,
                      par='recruit', max_gradient=tmp$opt$max_gradient)
    return(bind_rows(pars,ssb,rec))
      }
  if(parallel){
    message("Preparing parallel session..")
    if(!require(snowfall))
      stop("snowfall package required for parallel execution")
    library(snowfall)
    cores <- parallel::detectCores()-1
    sfInit(cpus=cores, parallel=TRUE)
    sfLibrary(TMB)
    sfLibrary(GOApollock)
    sfExportAll()
    message('Starting self-test runs in parallel...')
    out <- sfLapply(reps, function(i) fit_test(fit, i))
  } else {
    out <- lapply(reps, function(i) fit_test(fit, i))
  }
  out <- do.call(rbind, out) %>%
    mutate(version=fit$version, relerror=(est-true)/true, abserror=est-true)
}



#' Run and plot likelihood profiles
#'
#' @param fits A list of fits
#' @param xseq The vector over which to run the profile
#' @param par Which parameter, currently supports "M", "q2" and "q6".
#' @param plot Whether to plot it
#' @param maxNLL The maximum NLL below which to filter out individual components, for decluttering the plots.
#' @param ylim Y limit for plots
#' @param estSigR Whether to estimate sigmaR
#' @param return.fits Whether to return the fits (not implemented yet)
#' @return A ggplot object of the profiles
#' @export
fit_profile <- function(fits, xseq, par=c('M','q2', 'q6'), plot=TRUE, maxNLL=1,
                        ylim=c(0,10),
                        estSigR=FALSE, return.fits=FALSE){
  stopifnot(!return.fits)
  if(is.pkfits(fits)){
    nmods <- length(fits)
  } else if(is.pkfit(fits)){
    nmods <- 1
  } else { stop("invalid fit object")}
  par <- match.arg(par)
  isM <- isq2 <- isq6 <- FALSE
  if(par=='M'){
    xlab <- 'M at age 6'
    isM <- TRUE
  } else if(par=='q2') {
    xlab <- 'Catchability (q) for NMFS BT'
    xseq <- log(xseq)
    isq2 <- TRUE
  } else if(par=='q6') {
  xlab <- 'Catchability (q) for Summer AT'
  xseq <- log(xseq)
  isq6 <- TRUE
}
  xx <- c("Total Catch", "Fishery Ages", "Fishery Lengths", "Shelikof 3+ Index",
          "Shelikof 3+ Ages", "Shelikof Lengths", "NMFS BT Index",
          "NMFS BT Ages", "NMFS BT Lengths", "Unused", "ADF&G Index",
          "ADF&G Ages", "ADF&G Lengths", "Shelikof Age 1",
          "Shelikof Age 2", "Summer AT Index", "Summer AT Ages",
          "Recruit Penalties", "TV Selex Penalties", "Proj Recruits",
          "TV Q penalties", "Unused", "Priors",
          "Unused", 'Ecov_exp', 'Ecov_obs')
  components <- data.frame(component=1:26, name=xx)#, ctype=ctype, csurv=csurv, name2=paste(csurv, ctype))

  tmp <- list()
  for(jj in 1:nmods){
    if(nmods==1) fit <- fits else fit <- fits[[jj]]
    k <- 1
    tmpfits <- list()
    message("Starting fit ", fit$version)
    for(ii in 1:length(xseq)){
      inputtmp <- fit$input; inputtmp$pars <- fit$parList
      inputtmp$version <- xseq[ii]
      if(isM){
        inputtmp$pars$natMscalar <- xseq[ii]
        inputtmp$map$natMscalar <- factor(NA)
      }
      if(isq2){
        inputtmp$pars$log_q2_mean <- xseq[ii]
        inputtmp$map$log_q2_mean <- factor(NA)
      }
      if(isq6){
        inputtmp$pars$log_q6 <- xseq[ii]
        inputtmp$map$log_q6 <- factor(NA)
      }
      if(estSigR) inputtmp$map$sigmaR <- factor(1)
      if(estSigR) inputtmp$random <- 'dev_log_recruit'
      tmpfits[[k]] <- fit_pk(inputtmp, newtonsteps=0, getsd=FALSE, control=list(trace=0), verbose = FALSE)
      k <- k+1
    }
    tmp[[jj]] <-
      get_rep(tmpfits, 'loglik')    %>% select(-name)  %>%
      group_by(version) %>%
      mutate(component=1:26, par=version, nll=-value) %>%
      mutate(version=fit$version) %>% select(-value)
  }
  out <- bind_rows(tmp) %>% merge(components, by='component')
  if(isM) out$par <- out$par*.3
  if(isq2 | isq6) out$par <- exp(out$par)
  nlls <- group_by(out, version, component) %>%
    mutate(deltaNLL=nll-min(nll)) %>%
    filter(max(deltaNLL)>maxNLL) %>% ungroup
  totals <- nlls %>% group_by(par,version) %>%
    summarize(nll=sum(nll), .groups='drop') %>%
    group_by(version) %>% mutate(name='total', deltaNLL=nll-min(nll))
  mins <- group_by(nlls, name, version) %>% filter(deltaNLL==min(deltaNLL))
  g <- ggplot(nlls, aes(par, deltaNLL, color=name)) +
    geom_line()+facet_wrap('version') +
    geom_line(data=totals, mapping=aes(color=NULL), lwd=1) +
    geom_point(data=mins, size=2) +
    geom_point(data=filter(totals, deltaNLL==min(deltaNLL)), pch=16, size=2, col=1)+
    theme(legend.position='top') + labs(color=NULL, x=xlab, y='Change in NLL') +
    coord_cartesian(ylim=ylim)
  if(plot) print(g)
  return(g)
  if(return.fits) return(tmpfit)
}

#' Turn off data components one at a time to see the effect
#' @param fits A fitted object. Currently does not support a list of fits.
#' @param plot Whether to plot it
#' @param Nmult How much to multiply the ISS by (<1 to downweight it)
#' @param CVmult How much to multiple the index CVs by (>1 to downweight it)
#' @param return.fits Whether to return the fits, if FALSE it returns the ggplot
#' @param ... Further arguments passed to \code{plot_ssb}.
#' @return Either a ggplot object of a list of fits.
#' @details Needs to be updated when the model accepts CVs and multN of zero for all surveys.
#' @export
fit_drop_surveys <- function(fits, plot=TRUE, Nmult=0, CVmult=10, return.fits=FALSE, ...){

  if(is.pkfits(fits)){
    stop("this actually doesn't make sense since version is overwritten, do them individually")
    fitsout <- sapply(fits, function(fit)
      run_drop_surveys(fit, plot=FALSE, return.fits=TRUE, Nmult=Nmult, CVmult=CVmult ))
  } else {
    stopifnot(is.pkfit(fits))
    fit <- fits
    fit$input$pars <- fit$parList
    mapoff <- function(x, slots) {
      for(slot in slots){
        if(length(x$pars[[slot]])==0){
          warning("slot ", slot, " not found in pars")
        } else {
          x$map[[slot]] <- factor(x$pars[[slot]]*NA)
        }
      }
      return(x)
    }

    x <- fit$input
    isDM <- length(x$pars$log_DM_pars)>0
    x$version <- 'Drop Shelikof'
    x$dat$multN_srv1 <- x$dat$multN_srv1*Nmult
    x$dat$multNlen_srv1 <- x$dat$multNlen_srv1*0
    x$dat$indxsurv_log_sd1 <- x$dat$indxsurv_log_sd1*CVmult
    x$dat$indxsurv_log_sd4 <- x$dat$indxsurv_log_sd4*CVmult
    x$dat$indxsurv_log_sd5 <- x$dat$indxsurv_log_sd5*CVmult
    x <- mapoff(x, c('log_q1_dev', 'log_q4', 'log_q5', 'log_q1_mean'))
    x <- mapoff(x, c('Ecov_beta'))
    y <- 1:5; y[2] <- NA
    if(isDM) x$map$log_DM_pars <- factor(y)
    drop1 <- fit_pk(x, filename=NULL, verbose=FALSE)

    x <- fit$input
    x$version <- 'Drop NMFS BT'
    x$dat$multN_srv2 <- x$dat$multN_srv2*Nmult
    x$dat$multNlen_srv2 <- x$dat$multNlen_srv2*0
    x$dat$indxsurv_log_sd2 <- x$dat$indxsurv_log_sd2*CVmult
    x <- mapoff(x, 'log_q2_mean')
    y <- 1:5; y[3] <- NA
    if(isDM) x$map$log_DM_pars <- factor(y)
    drop2 <- fit_pk(x, filename=NULL, verbose=FALSE)

    x <- fit$input
    x$version <- 'Drop ADF&G'
    x$dat$multN_srv3 <- x$dat$multN_srv3*Nmult
    x$dat$multNlen_srv3 <- x$dat$multNlen_srv3*0
    x$dat$indxsurv_log_sd3 <- x$dat$indxsurv_log_sd3*CVmult
    x <- mapoff(x, c('log_q3_mean', 'log_q3_dev'))
    y <- 1:5; y[4] <- NA
    if(isDM) x$map$log_DM_pars <- factor(y)
    drop3 <- fit_pk(x, filename=NULL, verbose=FALSE)

    x <- fit$input
    x$version <- 'Drop Summer AT'
    x$dat$multN_srv6 <- x$dat$multN_srv6*Nmult
    x$dat$multNlen_srv6 <- x$dat$multNlen_srv6*0
    x$dat$indxsurv_log_sd6 <- x$dat$indxsurv_log_sd6*CVmult
    x <- mapoff(x, 'log_q6')
    y <- 1:5; y[5] <- NA
    if(isDM) x$map$log_DM_pars <- factor(y)
    drop6 <- fit_pk(x, filename=NULL, verbose=FALSE)

    fitsout <- list(fit, drop1, drop2, drop3, drop6)
  }
  if(return.fits) return(fitsout)
  if(plot) plot_ssb(fitsout, ...)
}
