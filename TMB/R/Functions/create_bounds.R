#' Calculate bounds for parameters
#' @param obj Object
#' @return list of lower (lwr) and upper (upr) vectors matching
#'   the object
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
