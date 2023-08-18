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
  upr['log_slp2_srv1'] <- 5
  lwr['log_slp2_srv1'] <- -5
  upr['inf2_srv1'] <- 20#10
  lwr['inf2_srv1'] <- 3
  upr['log_slp1_srv2'] <- 5
  lwr['log_slp1_srv2'] <- -5
  upr['inf1_srv2'] <- 50
  lwr['inf1_srv2'] <- 1
  upr['log_q2_mean'] <- 10
  lwr['log_q2_mean'] <- -10

  # - Non-parametric fish pars
  upr['sel_rho_c'] <- 20
  lwr['sel_rho_c'] <- -20

  upr['sel_rho_a'] <- 20
  lwr['sel_rho_a'] <- -20

  upr['sel_rho_y'] <- 20
  lwr['sel_rho_y'] <- -20

  upr['mean_sel'] <- 50
  lwr['mean_sel'] <- -50

  upr['selpars_re'] <- 20
  lwr['selpars_re'] <- -20

  upr['ln_sel_sd'] <- 5
  lwr['ln_sel_sd'] <- -50
  return(list(lwr=lwr, upr=upr))
}
