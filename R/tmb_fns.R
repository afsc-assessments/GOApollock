
#' Fit a GOA pollock model (BETA)
#' @param input Input list as returned by
#'   \code{prepare_pk_input}.
#' @param newtonsteps,control Optimizer settings
#' @param getsd Calculate and return sdreport?
#' @param do.fit Optimize or return obj? Used for testing.
#' @param use_bounds Whether to use bounds, slows optimization
#'   down but may be necessary. See \link{\code{get_bounds}}
#'   which is called internally.
#' @return A list object of class 'pkfit' which contains a
#'   "version" model name, rep, opt as returned by
#'   \code{TMBHelper::fit_tmb} but without the SD slot, sd if
#'   \code{getsd=TRUE}, and the obj.
#' @details This function is beta still.
#' @export
fit_pk <- function(input, getsd=TRUE, newtonsteps=1,
                   control=NULL, do.fit=TRUE,
                   use_bounds=FALSE){

  compile("../source/goa_pk_tmb.cpp")
  dyn.load('../source/goa_pk_tmb.dll')
  obj <- MakeADFun(data=input$dat, parameters=input$pars,
                   map=input$map, random=input$random,
                   DLL='goa_pk_tmb', silent=TRUE)
  obj <- MakeADFun(data=input$dat, parameters=input$pars,
                   map=input$map, random=input$random,
                   DLL='goa_pk_tmb', silent=TRUE)
  if(!do.fit) return(obj)
  if(use_bounds){
    lwr <- get_bounds(obj)$lwr
    upr <- get_bounds(obj)$upr
  } else {
    lwr <- -Inf; upr <- Inf
  }
  if(is.null(control))
    control <- list(eval.max=10000, iter.max=10000, trace=100)
  ## optimize and compare
  opt <- TMBhelper::fit_tmb(obj, control=control,
                            newtonsteps=newtonsteps, getsd=FALSE,
                            lower=lwr, upper=upr)
  sdrep <- sdreport(obj)
  rep <- obj$report()
  std <- with(sdrep, data.frame(name=names(value), est=value, se=sqrt(diag(cov)))) %>%
    group_by(name) %>% mutate(year=1969+1:n(), lwr=est-1.96*se, upr=est+1.96*se) %>% ungroup

  fit <- list(version=input$version, rep=rep, opt=opt, sd=std, obj=obj)
  class(fit) <- c('pkfit', 'list')
  return(fit)
}

#' Prepare inputs for the TMB GOA pollock model.
#' @param dat Data list as read in with
#'   \link{\code{read_pk_dat}}.
#' @param version A character string for the model name, used for
#'   plotting downstream
#' @param random A character vector declaring random effect
#'   vectors to be integrated. Defaults to NULL (fully penalized
#'   ML)
#' @return A list with the version, dat, pars, map and random
#'   which are used to build a TMB 'obj' in
#'   \link{\code{fit_pk}}.
#' @export
prepare_pk_input <- function(dat, version='none', random=NULL){
  ## Prepare the parameter list based on the data
  pars <- prepare_par(dat)
  map <- prepare_map(pars)
  out <- list(version=version, dat=dat, pars=pars, map=map, random=random)
  return(out)
}

#' Internal fucntion to prepare parameter list from dat
#' list. Uses initial values from the 2022 model for scalars, and
#' vectors are all 0.
prepare_par <- function(dat){
  ## build a default starting list, based off MLE in 2022 but
  ## truncated a bit
  nyrs <- length(dat$styr:dat$endyr)
  stopifnot(nyrs>30)
  pars <- list(dev_log_initN = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
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
