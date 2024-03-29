

#' Plot OSA composition residuals as bubble plots
#' @param res A matrix of obs and exp results cbinded together,
#'   as returned by the model
#' @param years Vector of years for plotting
#' @param Neff Effective sample size vector. Non-integer values
#'   are rounded currently.
#' @param ind Which age to drop, modified as need when 'drop' is
#'   used.
#' @param survey Character name for plotting
#' @param drop Optional vector of ages to drop. Alters meaning of
#'   'ind'
#' @param plot Whether to plot or return the data instead
#' @param minyear Filter out years before this
#' @return A ggplot object if plot is TRUE or the melted data if
#'   FALSE
#' @export
plot_osa_comps <- function(res, years, Neff, ind=1, survey,
                           drop=NULL, plot=TRUE, minyear=NULL){
  obs <- res[,1:10]
  exp <- res[,11:20]
  index <- 1:10
  if(!is.null(minyear)){
    indminyear <- which(years==minyear-1)
    obs <- obs[-(1:indminyear),]
    exp <- exp[-(1:indminyear),]
    years <- years[-(1:indminyear)]
    Neff <- Neff[-(1:indminyear)]
  }
  if(!is.null(drop)){
    obs <- obs[,-drop]
    exp <- exp[,-drop]
    index <- index[-drop]
  }
  stopifnot(all.equal(nrow(obs), nrow(exp), length(years)))
  stopifnot(all.equal(ncol(obs), ncol(exp), length(index)))
  o <- round(Neff*obs/rowSums(obs),0); p=exp/rowSums(exp)
  Neff2 <- apply(o,1, sum);
  if(any(Neff2==0)) warning("Some Neff were rounded to 0")
  index2 <- index[-ind]                 #move one to drop to end
  o2 <- cbind(o[,-ind], o[,ind])
  p2 <- cbind(p[,-ind], p[,ind])
  resid <- compResidual::resMulti(t(o2), t(p2))
  ## not sure why these fail sometimes?
  if(!all(is.finite(resid))) {warning('failed when ind=',ind)}
  mat <- t(matrix(resid, nrow=nrow(resid), ncol=ncol(resid)))
  dimnames(mat) <- list(year=years, index=index2)
  reslong <- reshape2::melt(mat, value.name='resid') #%>% filter(is.finite(resid))
  g <- ggplot(reslong, aes(year, index, size=abs(resid),
                           color=resid>0)) +
    geom_point()+
    ggtitle(paste0('OSA Residuals for the ', survey, ' w/o age=', index[ind])) +
    labs(y='Age', x=NULL) + scale_y_continuous(breaks=1:10,
                                               limits=c(1,10)) +
    scale_size(range=c(0,3))
  if(plot){
    print(g)
    return(invisible(g))
  } else {
    return(reslong)
  }
}


#' Calculate SDNR for indices
#' @param fit A fit or list of fits as returned by
#'   \code{fit_pk}. Not currently implemented for the ADMB model.
#' @return A data.frame of SDNR for all surveys with a column for
#'   version
#' @details SDNR=standard deviation of normalized residuals sd((o-e)/std)
get_sdnr <- function(fit){
  calc_sdnr <- function(obs,exp,CV) sd((obs-exp)/sqrt(log(CV^2+1)))
  if(class(fit)[1]=='pkfit'){
    dat <- fit$obj$env$data
    rep <- fit$rep
    yrs <- rep$years
    out <- lapply(c(1:6), function(i){
      obs <- dat[[paste0('indxsurv',i)]]
      exp <- rep[[paste0('Eindxsurv',i)]][yrs %in% dat[[paste0('srvyrs',i)]]]
      CV <- dat[[paste0('indxsurv_log_sd',i)]]
      stopifnot(all.equal(length(obs), length(exp), length(CV)))
      data.frame(survey=i, sdnr=calc_sdnr(obs, exp,CV),
                 version=fit$version)
    }
    ) %>% do.call(rbind,.)
  } else {
    ## list of fits so use recursion
    out <- lapply(fit, get_sdnr) %>% do.call(rbind,.)
  }
  return(out)
}


#' Calculate Mohn's rho
#' @param reps A list of replists for each peel as returned by
#'   \code{read_pk_rep}
#' @param max_peels Optional numeric to control the max number of
#'   peels to use in the calculation. E.g. can calculate rho with
#'   4,5,6, etc. peels and compare
#' @return The rho value rounded to 3 decimal places, and prints
#'   to console
#' @export
calculate_rho <- function(reps, max_peels=NULL){
  if(is.null(reps[[1]]$Expected_spawning_biomass)){
    ssb <- mymelt(reps, 'Espawnbio')
  } else {
    ssb <- mymelt(reps, 'Expected_spawning_biomass')
  }
  ssb <- ssb %>%
    mutate(peel=as.numeric(gsub('peel','', model))) %>%
    rename(SSB=value) %>%
    group_by(year) %>%
    mutate(SSB_Pct_Diff=100*(SSB-SSB[peel==0])/SSB[peel==0]) %>%
    ungroup()
  if(!is.null(max_peels)) ssb <- filter(ssb, peel<=max_peels)
  ## Calculate Mohn's rho. For each terminal point in the peel SSB,
  ## compare to that same year in the full run and calculate
  ## relative differences. Already did this above so just grab the
  ## right years to average
  rho <- ssb%>% filter(max(year)-peel==year & year!=max(year)) %>%
    pull(SSB_Pct_Diff)
  message("Percentage range of annual differences:",round(min(rho),1)," to ", round(max(rho),1))
  message("Median absolute error of rho:",round(median(abs(rho))))
  rho <- mean(rho)/100
  rho.lab <- paste0("Mohn's rho= ", round(rho,3))
  print(rho.lab)
  return(round(rho,3))
}


#' Read the GOA pollock model output report
#' @param model model name
#' @param path path to folder
#' @param endyr,styr The start and end years in the model
#' @param version A version name which is added, e.g., 'change_selex'
#' @return A list of outputs
#' @export
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


#' Extract and melt a named object from a replist (beta)
#' @param replist Output from \code{read_pk_rep}. Can be a single
#'   report object or a list of them.
#' @param slot The slot name to extract
#' @export
#' @return A data frame with named arguments sometimes
#'
mymelt <- function(replist, slot){

  ## multiple runs together or single?
  multi <- ifelse(length(replist)>100, FALSE,TRUE)

  if(!multi){
    ## Matrix already has dimnames for ages and years as appropriate
    y <- replist[[slot]]
    if(is.matrix(y)){
      ## onlyh need this for TMB output since doens't have
      ## dimnames set by read_pk_rep
      if(nrow(y)==length(replist$years) & ncol(y)==length(replist$ages))
        dimnames(y) <- list(year=replist$years, age=replist$ages)
      temp <- data.frame(version=replist$version, reshape2::melt(y))
      if(nrow(temp)==length(replist$ages)) temp <- cbind(temp, age=replist$ages)
      if(nrow(temp)==length(replist$years)) temp <- cbind(temp, year=replist$years)
      temp
    } else {
      temp <- data.frame(model=replist$version, value=y)
      if(nrow(temp)==length(replist$ages)) temp <- cbind(temp, age=replist$ages)
      if(nrow(temp)==length(replist$years)) temp <- cbind(temp, year=replist$years)
      temp
    }
    x <- temp
  }  else {
    ## use recursion to get individual ones
    x <- do.call(rbind, lapply(replist, function(z) mymelt(z, slot=slot)))
  }
  return(x)
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
#' @export
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
#' @export
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



