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
#' @return A named list of all elements with corresponding names
#'   to the ADMB model
read_dat <- function(filename, path=NULL){
  if(!is.null(path)){
    oldwd <- getwd()
    on.exit(setwd(oldwd))
    setwd(path)
  }
  char.lines <- readLines(filename)
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
    s <- scan(filename, quiet = T, what = integer(), comment.char="#",
              skip=dat.start[ind <<- ind+1], n=n)
    if(is.null(nrow)) return(s)
    return(matrix(s, nrow=nrow, ncol=ncol, byrow=TRUE))
  }
  sn <- function(ind, n, nrow=NULL, ncol=NULL){
    s <- scan(filename, quiet = T, what = numeric(), comment.char="#",
              skip=dat.start[ind <<- ind+1], n=n)
    if(is.null(nrow)) return(s)
    return(matrix(s, nrow=nrow, ncol=ncol, byrow=TRUE))
  }
  tmp <- si(ind , n = 2)
  d$styr <- tmp[1]; d$endyr <- tmp[2]
  nyrs <- tmp[2]-tmp[1]+1
  nyrs2 <- nyrs-1
  tmp <- si(ind=ind, n = 2)
  d$rcrage <- tmp[1]; d$trmage <- tmp[2]
  tmp <- si(ind=ind, n = 3)
  d$nbins1 <- tmp[1]; d$nbins2 <- tmp[2]; d$nbins3 <- tmp[3]
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
  d$catp <- sn( ind=ind, n = d$nyrs_fsh*10, nrow=d$nyrs_fsh, ncol=10)
  d$lenp <- sn( ind=ind, n = d$nyrslen_fsh*d$nbins1, nrow=d$nyrslen_fsh, ncol=d$nbins1)
  d$wt_fsh <- sn( ind=ind, n = nyrs*10, nrow=nyrs, ncol=10)
  ## This BS survey is not longer used but the inputs are
  ## here. Will be dropped at some point.
  d$nyrs_srv1_bs <- si(ind=ind, n =1)
  d$srvyrs1_bs <- si(ind=ind, n =d$nyrs_srv1_bs)
  d$indxsurv1_bs <- sn( ind=ind, n =d$nyrs_srv1_bs)
  d$indxsurv_log_sd1_bs <- sn( ind=ind, n =d$nyrs_srv1_bs)
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
  d$srvp1 <- sn(ind=ind, n =d$nyrsac_srv1*10,
                nrow=d$nyrsac_srv1, ncol=10)
  d$srvlenp1 <- sn(ind=ind, n =d$nyrslen_srv1*d$nbins3,
                   nrow=d$nyrslen_srv1, ncol=d$nbins3)
  d$wt_srv1 <- sn( ind=ind, n = nyrs*10, nrow=nyrs, ncol=10)
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
  d$srvp2 <- sn(ind=ind, n =d$nyrsac_srv2*10,
                nrow=d$nyrsac_srv2, ncol=10)
  d$srvlenp2 <- sn(ind=ind, n =d$nyrslen_srv2*d$nbins3,
                   nrow=d$nyrslen_srv2, ncol=d$nbins3)
  d$wt_srv2 <- sn( ind=ind, n = nyrs*10, nrow=nyrs, ncol=10)
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
  d$srvp3 <- sn(ind=ind, n =d$nyrsac_srv3*10, nrow=d$nyrsac_srv3, ncol=10)
  d$srvlenp3 <- sn(ind=ind, n =d$nyrslen_srv3*d$nbins3, nrow=d$nyrslen_srv3, ncol=d$nbins3)
  d$wt_srv3 <- sn( ind=ind, n = nyrs*10, nrow=nyrs, ncol=10)
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
  d$srvp6 <- sn(ind=ind, n =d$nyrsac_srv6*10,
                nrow=d$nyrsac_srv6, ncol=10)
  d$srvlenp6 <- sn(ind=ind, n =d$nyrslen_srv6*d$nbins2,
                   nrow=d$nyrslen_srv6, ncol=d$nbins2)
  d$wt_srv6 <- sn( ind=ind, n = nyrs*10, nrow=nyrs, ncol=10)
  d$age_trans <- sn( ind=ind, n = 10*10, nrow=10, ncol=10)
  d$len_trans1 <- sn( ind=ind, n = d$nbins1*10, nrow=10, ncol=d$nbins1)
  d$len_trans2 <- sn( ind=ind, n = d$nbins2*10, nrow=10, ncol=d$nbins2)
  d$len_trans3 <- sn( ind=ind, n = d$nbins3*10, nrow=10, ncol=d$nbins3)
  ## the projection module inputs
  d$wt_pop <- sn( ind=ind, n = nyrs*10, nrow=nyrs, ncol=10)
  d$wt_spawn <- sn( ind=ind, n = nyrs*10, nrow=nyrs, ncol=10)
  d$mat_old <- sn( ind=ind, n=10)
  d$mat <- sn( ind=ind, n =10)
  d$wt_pop_proj <- sn( ind=ind, n =10)
  d$wt_spawn_proj <- sn( ind=ind, n =10)
  d$wt_fsh_proj <- sn( ind=ind, n =10)
  d$wt_srv_proj <- sn( ind=ind, n =10)
  d$Ftarget <- sn( ind=ind, n =5)
  d$B40 <- sn( ind=ind, n =1)
  d$log_mean_recr_proj <- sn( ind=ind, n =1)
  d$sigmasq_recr <- sn( ind=ind, n =1)
  return(d)
}
