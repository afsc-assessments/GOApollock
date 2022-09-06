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
  ssb <- mymelt(reps, 'Expected_spawning_biomass') %>%
    mutate(model=as.numeric(gsub('peel','', model))) %>%
    rename(peel=model, SSB=value) %>%
    group_by(year) %>%
    mutate(SSB_Pct_Diff=100*(SSB-SSB[peel==0])/SSB[peel==0]) %>%
    ungroup()
  if(!is.null(max_peels)) ssb <- filter(ssb, peel<=max_peels)
  ## Calculate Mohn's rho. For each terminal point in the peel SSB,
  ## compare to that same year in the full run and calculate
  ## relative differences. Already did this above so just grab the
  ## right years to average
  rho <- ssb%>% filter(2021-peel==year & year!=2021) %>%
    pull(SSB_Pct_Diff)
  message("Percentage range of annual differences:",round(min(rho),1)," to ", round(max(rho),1))
  message("Median absolute error of rho:",round(median(abs(rho))))
  rho <- rho %>% mean %>% '/'(100) # ugly syntax to /100
  rho.lab <- paste0("Mohn's rho= ", round(rho,3))
  print(rho.lab)
  round(rho,3)
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
      temp <- data.frame(model=replist$version, reshape2::melt(y))
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
