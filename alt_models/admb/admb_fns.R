
#' Clean temporary files from a directory containing a model run
#'
#' @param path Path to directory to clean
#' @param full Whether to delete everything except .dat and .tpl files
#' @return Nothing
clean_pk_dir <- function(path=getwd(), full=FALSE){
  stopifnot(is.logical(full))
  if(full){
    r <- list.files(path)
    r <- r[-grep(pattern='.tpl|.dat', r)]
    if(length(r)>0) trash <- file.remove(r)
  } else {
    r <- list.files(path, pattern='r0|b0|p0')
    r <- c(r,list.files(path, pattern='\\.cpp|\\.obj|\\.log|\\.eva|\\.bar|\\.htp|\\.dep'))
    r <- r[-grep('tmb',r)]
    if(length(r)>0) trash <- file.remove(file.path(path, r))
    s <- file.size(r <- file.path(path, 'mceval.dat'))
    if(!is.na(s)){
      if(s<.0001){
        trash <- file.remove(r)
      } else {
        message("Not removing mceval.dat becuase it appears to have results")
      }
    }
    r <- list.files(path, '\\.psv')
    s <- file.size(file.path(path,r))
    if(length(s)>0){
      if(s<.0001){
        trash <- file.remove(r)
      } else {
        message("Not removing .psv file becuase it appears to have results")
      }
    }
  }
}


#' Copy and build ADMB source from package
#'
#' @param path Path to folder
#' @param name Executable name (default goa_pk)
#' @param compile Whether to compile (default is TRUE)
setup_exe <- function(path=getwd(), name='goa_pk', compile=TRUE){
  tpl <- 'C:/Users/cole.monnahan/GOApollock/source/goa_pk.tpl'
  stopifnot( file.exists(tpl))
  dir.exists(path)
  test <- file.copy(from=tpl, to=file.path(path,paste0(name, '.tpl')), overwrite=TRUE)
  if(!test) warning("Failed to copy file")
  if(compile){
    message("Compiling model..")
    old.wd <- getwd(); on.exit(setwd(old.wd))
    setwd(path)
    system(paste('admb', name))
  }
}



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
read_dat <- function(filename,path=NULL){
  .Deprecated('read_pk_dat')
  read_pk_dat(filename,path)
}



#' Read the GOA pollock model output report
#' @param model model name
#' @param path path to folder
#' @param endyr,styr The start and end years in the model
#' @param version A version name which is added, e.g., 'change_selex'
#' @return A list of outputs
read_pk_rep <- function(model='goa_pk', path=getwd(), version='none', endyr,
                       styr=1970){
  ## named vectors
  fyrs <- styr:endyr
  fages <- 1:10
  add.names <- function(x){
    if(is.matrix(x)){
      if(nrow(x)==length(fyrs) & ncol(x)==length(fages))
        dimnames(x) <- list(year=fyrs, age=fages)
    } else {
      #if(length(x)==length(fyrs)) names(x) <- list('year'=fyrs)
      #if(length(x)==length(fages)) names(x) <- fages
    }
    x
  }

  file <- file.path(path,model)
  if(length(grep('.rep', file))==0) file <- paste0(file,'.rep')
  if(!file.exists(file)) stop("File not found=", file)

  x <- scan(file, what="", sep="\n")
  ## drop initial spaces
  y <- unname(sapply(x, trimws))
  ## Separate elements by one or more whitepace
  z <- strsplit(y, "[[:space:]]+")
  ## Drop empty lines
  z <- purrr::discard(lapply(z, function(x) if(length(x)>0) x), is.null)
  ind <-  suppressWarnings(which(is.na(as.numeric(sapply(z, `[[`, 1)))))
  k <- 1; dummy <- NULL
  mynames <- NA
  myvals <- list()
  for(i in 1:(length(ind)-1)){
    if(i==64) next # ragged shape, need to fix
    if(ind[i]==ind[i+1]-1){
      ## two chars in a row
      dummy <- z[ind[i]]
      next
    } else { dummy <- NULL }
    tmp <- z[ind[i]:(ind[i+1]-1)]
    name <- paste(c(dummy, tmp[[1]]), collapse='_')
    ## special case
    if(tmp[[1]][1]=='Ftarget'){
      next
      ## myvals[[k]] <- list(Ftarget=as.numeric(tmp[[2]]),
      ##                     B40=as.numeric(tmp[[3]]))
      ## k <- k+1
    }
    if(length(tmp)==2){
      vals <- tmp[[-1]]
      val <- strsplit(vals, "[[:space:]]+")
      if(length(val)==1){
        ## scalar
        val <- as.numeric(val[[1]])
      } else {
        ## vector
        val <- as.numeric(unlist(val))
      }
    } else {
      ## matrix of some sort
      vals <- tmp[-1]
      val <- lapply(vals, function(x) strsplit(x, "[[:space:]]"))
      if(length(val[[1]])==1){
        ## column vector
        val <- as.numeric(unlist(val))
      } else {
        ## matrices
        val <- do.call(rbind, lapply(val, as.numeric))
        stopifnot(length(vals)==nrow(val))
      }
    }
    val <- add.names(val)
    myvals[[k]] <- val
    ##  names(mylist[[k]]) <- name
    mynames[k] <- name
    k <- k+1
  }
  names(myvals) <- mynames
  myvals <- c( version=version, ages=list(fages), years=list(fyrs),  myvals)
  return(myvals)
}


## ## tests that it works on all the data types I want
## x1 <- read_pk_rep('model_runs/m01_2020_final_original/pk20_8.rep',
##                   version='pk20_8', endyr=2020)
## x2 <- x1; x2$model <- 'dummy'
## replists <- list(x1,x2)

## mymelt(x1, 'Expected_spawning_biomass') %>% str
## mymelt(x1, 'Numbers_at_age') %>% str
## mymelt(x1, 'Natural_mortality') %>% str
## mymelt(replists, 'Expected_spawning_biomass') %>% str
## mymelt(replists, 'Numbers_at_age') %>% str
## mymelt(replists, 'Natural_mortality') %>% str


#' Read correlation file
#' @param model model name
#' @param path path to folder
#' @param version,endyr,styr see rep function
read_pk_cor <- function(model='goa_pk', path=getwd(), version='none', endyr, styr=1970){
  if(!file.exists(ff <- file.path(path, paste0(model, '.cor'))))
    stop("file does not exists: ",ff)
  yrs <- styr:endyr
  oldwd <- getwd(); on.exit(setwd(oldwd))
  setwd(path)
  sdrep <- suppressWarnings(R2admb:::read_admb(model))
  ## Parse these out into scalars and vectors
  xx <- strsplit(names(sdrep$coefficients), '\\.') %>% sapply(function(x) x[1])
  df <- data.frame(version=version, name=xx, est=as.numeric(sdrep$coefficients),
    se=as.numeric(sdrep$se)) %>%
    group_by(name) %>% mutate(i=1:n(), year=yrs[1:n()]) %>% ungroup() %>%
    mutate(lwr=est-1.96*se, upr=est+1.96*se)
  df                                       # }
}


#' Read standard deviation file
#' @param model model name
#' @param path path to folder
#' @param version,endyr,styr see rep function
read_pk_std <- function(model='goa_pk', path=getwd(), version='none', endyr, styr=1970){
  yrs <- styr:endyr
  ff <- file.path(path, paste0(model,'.std'))
  if(!file.exists(ff))
    stop("file ",ff, " does not exist")
  df <- read.table(ff, header=TRUE) %>% cbind(version=version) %>%
    rename(se=std.dev, est=value) %>% select(-index) %>%
    group_by(name) %>%
    mutate(year=yrs[1:n()], i=1:n()) %>% ungroup %>%
    mutate(lwr=est-1.96*se, upr=est+1.96*se)
  df
}
