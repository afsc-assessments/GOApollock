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

#' Check object of class pkfit
#' @param x Returned list from \code{\link{fit_tmb}}
#' @export
is.pkfit <- function(x) inherits(x, "pkfit")


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
