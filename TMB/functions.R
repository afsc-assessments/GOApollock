#' Explore OSA residuals for multinomial composition data and
#' compare to Pearson
#' @param obs,exp,pearson the observed, expected and Pearson
#'   residual matrices with rows as years and columns as ages (or
#'   lengths)
#' @param ages,years vectors giving the ages and years
#' @param model Character for model version
#' @return returns nothing but creates a PDF file in the working
#' directory
#'
plot_osa_comps <- function(obs, exp, ages, years, Neff, model,
                           plot.type=c('default', 'ggplot')){
  require(compResidual) ## https://github.com/fishfollower/compResidual
  plot.type <- match.arg(plot.type)
  stopifnot(all.equal(nrow(obs), nrow(exp), length(years)))
  stopifnot(all.equal(ncol(obs), ncol(exp), length(ages)))
  Neff <- ceiling(Neff)
  o <- round(Neff*obs/rowSums(obs),0); p=exp/rowSums(exp)
  rowSums(o)
  rowSums(p)
  res <- resMulti(t(o), t(p))
  if(!all(is.finite(res))){
    warning("failed to calculate OSA residuals")
    browser()
    return(NULL)
  }
  if(plot.type=='default'){
    ## default output
    plot(res)
    mtext(text=paste('GOA pollock fishery OSA comps:', model),
          line=3)
  } else {
    require(ggplot2)
    mat <- t(matrix(res, nrow=nrow(res), ncol=ncol(res)))
    dimnames(mat) <- list(year=years, age=ages[-1])
    reslong <- reshape2::melt(mat, value.name='resid')
    g1 <- ggplot(reslong, aes(year, age, size=abs(resid),
                              color=resid>0)) + geom_point() +
      ggtitle(model) + ylim(range(ages))
    return(g1)
  }
  return(invisible(res))
}


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
  return(list(lwr=lwr, upr=upr))
}
