

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
#' @param model Currently supports 'multinomial' and
#'   'Dirichlet-multinomial', when using hte latter the parameter
#'   theta must be passed to recalculate the "alpha" term
#'   internally.
#' @param theta A scalar parameter estimated as exp(log_DM_pars)
#'   corresponding to this survey
#' @return A ggplot object if plot is TRUE or the melted data if
#'   FALSE
#' @export
plot_osa_comps <- function(res, years, Neff, ind=1, survey,
                           drop=NULL, plot=TRUE, minyear=NULL,
                           model=c('multinomial','Dirichlet-multinomial'),
                           theta){
  model <- match.arg(model)
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
  if(model=='Dirichlet-multinomial'){
    ## recreate the alpha matrix, must match the cpp side, here
    ## it is linear form
    alpha <- rowSums(o)*p*theta
    o2 <- cbind(o[,-ind], o[,ind])
    alpha <- cbind(alpha[,-ind], alpha[,ind])
    resid <- compResidual::resDirM(t(o2), t(alpha))
  } else {
    o2 <- cbind(o[,-ind], o[,ind])
    p2 <- cbind(p[,-ind], p[,ind])
    resid <- compResidual::resMulti(t(o2), t(p2))
  }
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
    scale_size_continuous(##range  = c(0, 3),
      ## limits = c(0, 4),
                          breaks = c(0, 1,2,3))
    ##scale_size_area(max_size=2)
  if(plot){
    print(g)
    return(invisible(g))
  } else {
    return(cbind(reslong, survey=survey))
  }
}


#' Calculate SDNR for indices
#' @param fit A fit or list of fits as returned by
#'   \code{fit_pk}. Not currently implemented for the ADMB model.
#' @return A data.frame of SDNR for all surveys with a column for
#'   version
#' @export
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
      df <- length(obs)-1
      data.frame(survey=i, sdnr=calc_sdnr(obs, exp,CV),
                 version=fit$version, upper=sqrt(qchisq(.95, df)/df))
    }
    ) %>% do.call(rbind,.)
  } else {
    ## list of fits so use recursion
    out <- lapply(fit, get_sdnr) %>% do.call(rbind,.)
  }
  return(out)
}


#' Calculate Mohn's rho
#' @param fits A list of pkfits for each peel as returned by
#'   \code{fit_pk_retro}.
#' @param type Which metric to calculate: 'ssb', 'F' or 'recruit'
#' @param max_peels Optional numeric to control the max number of
#'   peels to use in the calculation. E.g. can calculate rho with
#'   4,5,6, etc. peels and compare
#' @return The rho value rounded to 3 decimal places, and prints
#'   to console
#' @export
calculate_rho <- function(fits, type=c('ssb', 'F', 'recruit'),
                          max_peels=NULL){
  ## if(is.null(reps[[1]]$Expected_spawning_biomass)){
  ##   ssb <- mymelt(reps, 'Espawnbio')
  ## } else {
  ##   ssb <- mymelt(reps, 'Expected_spawning_biomass')
  ## }
  stopifnot(all(sapply(fits, is.pkfit)))
  stopifnot(length(fits)>5)
  type <- match.arg(type)
  if(type=='ssb'){
    x <- get_rep(fits, 'Espawnbio')
  } else if(type=='F'){
    x <- get_rep(fits, 'F')
  } else if(type=='recruit'){
    x <- get_rep(fits, 'recruit')
  } else { stop('something wrong with type=',type)}
  stopifnot(nrow(x)>0)
  x <- x %>%
    mutate(peel=as.numeric(gsub('peel','', version))) %>%
    group_by(year) %>%
    mutate(Pct_Diff=100*(value-value[peel==0])/value[peel==0]) %>%
    ungroup()
  if(!is.null(max_peels)) x <- filter(x, peel<=max_peels)
  ## Calculate Mohn's rho. For each terminal point in the peel SSB,
  ## compare to that same year in the full run and calculate
  ## relative differences. Already did this above so just grab the
  ## right years to average
  rho <- x %>% filter(max(year)-peel==year & year!=max(year)) %>%
    pull(Pct_Diff)
  message("Percentage range of annual differences:",round(min(rho),1)," to ", round(max(rho),1))
  message("Median absolute error of rho:",round(median(abs(rho))))
  rho <- mean(rho)/100
  rho.lab <- paste0("Mohn's rho= ", round(rho,3)," for ", type)
  message(rho.lab)
  return(round(rho,3))
}


#' Extract and melt a named object from a replist (beta)
#' @param replist Output from \code{get_rep}. Can be a single
#'   report object or a list of them.
#' @param slot The slot name to extract
#' @export
#' @return A data frame with named arguments sometimes
#'
mymelt <- function(replist, slot){

  ## multiple runs together or single?
  ## multi <- !is.pkfit(replist)
  ## if(multi){
  ##   if(!all(sapply(replist, is.pkfit)))
  ##     stop("Some elements of replist are not pkfit objects")
  ## }
  multi <- ifelse(length(replist)>100, FALSE,TRUE)
  if(!multi){
    ## Matrix already has dimnames for ages and years as appropriate
    y <- replist[[slot]]
    if(NROW(y)==0) stop("Slot '", slot, "' not found in report")
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
      temp <- data.frame(version=replist$version, value=y)
      if(nrow(temp)==length(replist$ages)) temp <- cbind(temp, age=replist$ages)
      if(nrow(temp)==length(replist$years)) temp <- cbind(temp, year=replist$years)
      temp
    }
    x <- cbind(temp, name=slot)
  }  else {
    ## use recursion to get individual ones
    x <- do.call(bind_rows, lapply(replist, function(z) mymelt(z, slot=slot)))
  }
  return(x)
}




