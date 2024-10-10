#' Create factor from numeric survey vector. Useful for
#' informative plots and tables
#' @param x Vector of survey numbers, 1,2, etc. which are
#' transformed into a named factor
#' @return Vector of named surveys
#' @export
surveyf <- function(x){
  factor(x, levels=0:6,
         labels=c('Fishery', 'Shelikof', 'NMFS BT',
                  'ADF&G BT',
                  'Shelikof age 1',
                  'Shelikof age 2',
                  'Summer AT'))
}

#' Calculate ADSB (avg diff in SSB) for use in determining model
#' names. If ADSB>0.1 then it's considered a major version.
#'
#' @param base The base model std df or rep list
#' @param new The propsed model std df or rep list. Can have more
#'   years than base model
#' @return ADSB value and message
#' @export
calc_adsb <- function(base, new){
  if(is_tibble(base)){
    ## tibble from read_pk_cor
    yr0 <- filter(base, year >= 1977 & name=='Espawnbio') %>% pull('year')
    sb0 <- filter(base, year >= 1977 & name=='Espawnbio') %>% pull('est')
  } else {
    ## replist from read_pk_rep
    yr0 <- base$years
    sb0 <- base$Expected_spawning_biomass
    stopifnot(length(yr0)==length(sb0))
    ind <- which(yr0>=1977)
    yr0 <- yr0[ind]; sb0 <- sb0[ind]
  }
  if(is_tibble(new)){
    sb1 <- filter(new, year >= 1977 & name=='Espawnbio') %>% pull('est')
  } else {
    yr1 <- new$years
    sb1 <- new$Expected_spawning_biomass
    stopifnot(length(yr1)==length(sb1))
    ind <- which(yr1>=1977)
    yr1 <- yr1[ind]; sb1 <- sb1[ind]
  }
  stopifnot(length(sb0)>=length(sb1))
  adsb <- sqrt(mean((sb1[1:length(sb0)]/sb0-1)^2))
  message('This is a ', ifelse(adsb>.1, '**major**',
                               '**minor**'), ' model version (ADSB=',round(adsb,2),')')
  return(adsb)
}



#' Constructor for the "pkfit" (A-D fit) class
#' @param x Fitted object from \code{\link{fit_tmb}}
#' @return An object of class "pkfit"
#' @export
pkfit <- function(x){
  if(!is.list(x)){
    warning("Object passed to pkfit is not a list -- something went wrong in fitting?")
    return(x)
  }
  if(is.null(x$version)) stop("No version found, something went wrong")
  class(x) <- c('pkfit', 'list')
  x
}

#' Check if an object is of class pkfit
#' @param x Returned list from \code{\link{fit_tmb}}
#' @export
is.pkfit <- function(x) inherits(x, "pkfit")

#' Check if an object is a list of pkfit objects
#' @param x List of fits returned from \code{\link{fit_tmb}}
#' @export
is.pkfits <- function(x){
  if(!is.list(x)) {
    warning("Object passed to is.pkfits is not a list -- something went wrong")
    return(FALSE)
  }
  all(sapply(x, function(i) inherits(i, "pkfit")))
}



#' Print summary of pkfit object
#' @param fit Fitted object from \code{\link{fit_pk}}
#' @param ... Ignored for now
#' @return Summary printed to console
#' @method print pkfit
#' @export
print.pkfit <- function(fit, ...){
  cat("GOA pollock model version: ", fit$version, "\n")
  rt <- as.numeric(fit$opt$time_for_run, units='secs')
  ru <- 'seconds'
  if(rt>60*60*24) {
    rt <- rt/(60*60*24); ru <- 'days'
  } else if(rt>60*60) {
    rt <- rt/(60*60); ru <- 'hours'
  } else if(rt>60){
    rt <- rt/60; ru <- 'minutes'
  }
  cat("Total run time was", round(rt,2),  ru, '\n')
  cat("Number of parameters:", paste(names(fit$opt$number_of_coefficients), fit$opt$number_of_coefficients, sep='='),"\n")
  cat("Final maximum gradient=",
      sprintf("%.3g", fit$opt$max_gradient), "\n")
  cat("Marginal NLL=",  round(fit$opt$objective,5), "\n")
}


#' Run Francis weighting on a fitted model
#'
#' @param fit A fitted object returned by \code{fit_pk}
#' @param iter The number of iterations to do, defaults to 10
#' @param print.check Whether to print a plot showing the weights
#'   by iteration for each survey
#' @return A fitted object using from the last iteration
#' @details Implements Francis weighting for age composition data
#'   except for survey 3 which needs to be fixed still.
#' @export
run_francis_weighting <- function(fit, iter=10, print.check=TRUE){
  if(class(fit)[1]!='pkfit') stop("fit is not the correct type of object")
  calc_francis_weight <- function(obsexp, N0){
    o <- obsexp[,1:10]
    e <- obsexp[,11:20]
    ages <- 1:10
    n <- length(N0)
    resids <- numeric(n)
    for(i in 1:n){
      meano <- sum(o[i,] *ages)
      meane <- sum(e[i,] *ages)
      diff <- meano-meane
      x1 <- sum(ages^2*e[i,])-meane^2
      resids[i] <- diff/(sqrt(x1/N0[i]))
    }
    return(1/var(resids))
  }
  get_francis_weights <- function(ri, r0, iteration){
    wfsh <- calc_francis_weight(ri$res_fish, r0$multN_fsh)
    wsrv1 <- calc_francis_weight(ri$res_srv1, r0$multN_srv1)
    wsrv2 <- calc_francis_weight(ri$res_srv2, r0$multN_srv2)
    wsrv3 <- calc_francis_weight(ri$res_srv3, r0$multN_srv3)
    wsrv6 <- calc_francis_weight(ri$res_srv6, r0$multN_srv6)
    data.frame(iteration=iteration, fsh=wfsh, srv1=wsrv1,
               srv2=wsrv2, srv3=wsrv3, srv6=wsrv6)
  }
  input0 <- fit$input
  input0$pars <- fit$parList            # start from MLE of last run
  r0 <- fit$rep
  w0 <- wnew <- get_francis_weights(r0,r0, 0)
  weights <- w0
  ## Modify the original dat file with the latest estiamtes of weights
  for(ii in 1:iter){
    inputnew <- input0
    inputnew$dat$multN_fsh  <- wnew$fsh*inputnew$dat$multN_fsh
    inputnew$dat$multN_srv1 <- wnew$srv1*inputnew$dat$multN_srv1
    inputnew$dat$multN_srv2 <- wnew$srv2*inputnew$dat$multN_srv2
    ## inputnew$dat$multN_srv3 <- wnew$srv3*inputnew$dat$multN_srv3
    inputnew$dat$multN_srv6 <- wnew$srv6*inputnew$dat$multN_srv6
    if(ii==iter){
      fitnew <- fit_pk(inputnew)
    } else {
      fitnew <- fit_pk(inputnew, newtonsteps=0, getsd=FALSE, control=list(trace=0))
    }
    wnew <- get_francis_weights(fitnew$rep,r0, ii)
    weights <- rbind( weights, wnew)
  }
  ## check it converged
  check <- pivot_longer(weights, -iteration) %>%
    filter(!is.na(resid)) %>%
    ggplot(aes(iteration, value, color=name)) + geom_line() +
    geom_point() + ylim(0,NA)
  if(print.check) print(check)
  return(fitnew)
}

#' Map off a parameter in a given input list
#' @param x Input list
#' @param slots A character vector to map off
#' @return The modified input list
#' @export
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
