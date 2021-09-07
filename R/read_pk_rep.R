## Function to read the GOA pollock model output report into R
library(dplyr)
library(ggplot2)

read_pk_rep <- function(file, endyr=2020, styr=1970, model.name='none'){
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

  x <- scan(file, what="", sep="\n")
  ## drop initial spaces
  y <- sapply(x, trimws) %>% unname
  ## Separate elements by one or more whitepace
  z <- strsplit(y, "[[:space:]]+")
  ## Drop empty lines
  z <- lapply(z, function(x) if(length(x)>0) x) %>% purrr::discard(is.null)
  ind <-  suppressWarnings(which(is.na(as.numeric(sapply(z, `[[`, 1)))))
  k <- 1; dummy <- NULL
  mynames <- NA
  myvals <- list()
  for(i in 1:(length(ind)-1)){
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
        val <- unlist(val) %>% as.numeric()
      } else {
        ## matrices
        val <- lapply(val, as.numeric) %>% do.call(rbind,.)
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
  myvals <- c( model=model.name, ages=list(fages), years=list(fyrs),  myvals)
  return(myvals)
}

mymelt <- function(replist, slot){

  ## multiple runs together or single?
  multi <- ifelse(length(replist)>100, FALSE,TRUE)

  if(!multi){
    ## Matrix already has dimnames for ages and years as appropriate
    y <- replist[[slot]]
    if(is.matrix(y)){
      temp <- data.frame(model=replist$model, reshape2::melt(y))
      if(nrow(temp)==length(replist$ages)) temp <- cbind(temp, age=replist$ages)
      if(nrow(temp)==length(replist$years)) temp <- cbind(temp, year=replist$years)
      temp
    } else {
      temp <- data.frame(model=replist$model, value=y)
      if(nrow(temp)==length(replist$ages)) temp <- cbind(temp, age=replist$ages)
      if(nrow(temp)==length(replist$years)) temp <- cbind(temp, year=replist$years)
      temp
    }
    x <- temp
  } else {
    y <- replist[[1]][[slot]]
    if(is.matrix(y)){
      x <- lapply(replist, function(z){
        temp <- data.frame(model=z$model,reshape2::melt(z[[slot]]))
        if(nrow(temp)==length(z$ages)) temp <- cbind(temp, age=z$ages)
        if(nrow(temp)==length(z$years)) temp <- cbind(temp, year=z$years)
        temp
      }) %>% bind_rows
    } else {
      x <- lapply(replist, function(z){
        temp <- data.frame(model=z$model,value=z[[slot]])
        if(nrow(temp)==length(z$ages)) temp <- cbind(temp, age=z$ages)
        if(nrow(temp)==length(z$years)) temp <- cbind(temp, year=z$years)
        temp
      }) %>% bind_rows
    }
  }
  return(x)
}

## ## tests that it works on all the data types I want
## x1 <- read_pk_rep('model_runs/m01_2020_final_original/pk20_8.rep',
##                   model.name='pk20_8', endyr=2020)
## x2 <- x1; x2$model <- 'dummy'
## replists <- list(x1,x2)

## mymelt(x1, 'Expected_spawning_biomass') %>% str
## mymelt(x1, 'Numbers_at_age') %>% str
## mymelt(x1, 'Natural_mortality') %>% str
## mymelt(replists, 'Expected_spawning_biomass') %>% str
## mymelt(replists, 'Numbers_at_age') %>% str
## mymelt(replists, 'Natural_mortality') %>% str



read_pk_cor <- function(model, path, version, endyr, styr=1970){
  if(!file.exists(ff <- file.path(path, paste0(model, '.cor'))))
    stop("file does not exists: ",ff)
  yrs <- styr:endyr
  oldwd <- getwd(); on.exit(setwd(oldwd))
  setwd(path)
  sdrep <- R2admb:::read_admb(model)
  ## Parse these out into scalars and vectors
  xx <- strsplit(names(sdrep$coefficients), '\\.') %>% sapply(function(x) x[1])
  df <- data.frame(version=version, name=xx, est=as.numeric(sdrep$coefficients),
    se=as.numeric(sdrep$se)) %>%
    group_by(name) %>% mutate(i=1:n(), year=yrs[1:n()]) %>% ungroup() %>%
    mutate(lwr=est-1.96*se, upr=est+1.96*se)
  df                                       # }
}
