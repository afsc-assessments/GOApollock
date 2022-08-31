#' Plot spawning biomass
#'
#' @param x A list of data frames as read in from
#'   \link\code{read_pk_cor}
#' @param add_uncertainty Whether to add 95% confidence
#'   intervals. These are calculated in log space.
#' @param uselog Whether to calculate in log space or not
#' @param plotlog Whether to plot in log space or not
#' @param plot Whether to plot and return the ggplot object
#'   (TRUE) or return the data (FALSE).
#'
#' @return Either a ggplot object or the data used for it,
#' depending on \code{plot} argument.
#' @export
#'
plot_pk_ssb <- function(x, add_uncertainty=TRUE, plotlog=TRUE, uselog=FALSE, plot=TRUE){
  nmods <- length(x)
  x <- bind_rows(x)
  alpha1 <- .5 # for ribbon
  alpha2 <- .8 # for lines
  tmp <- 'Espawnbio'
  if(uselog)  tmp <- 'Espawnbio_log'
  if(uselog & nrow(filter(x, name==tmp))==0){
    tmp <- 'Espawnbio'
    uselog <- FALSE
    warning("Log of SSB not found in input, likely from old model run. Setting uselog=FALSE.")
  }
  if(uselog){
    ses <- x %>% filter(name==tmp) %>%
      mutate(est=exp(est), lwr=exp(lwr), upr=exp(upr))
  } else {
    ses <- x %>% filter(name==tmp) %>%
      mutate(lwr=est/exp(1.96*sqrt(log(1+(se/est)^2))),
             upr=est*exp(1.96*sqrt(log(1+(se/est)^2))))
  }
  g <- ggplot(ses, aes(year, est, ymin=lwr, ymax=upr, color=version, fill=version)) +
    geom_line(alpha=alpha2, lwd=1) +
    theme_bw()+ labs(x=NULL,y='Spawning biomass (M mt)')
  if(plotlog) g <- g+ scale_y_log10()
  if(add_uncertainty) g <-  g+ geom_ribbon(alpha=alpha1)
  if(plot)  {print(g); g}  else  return(ses)
}
