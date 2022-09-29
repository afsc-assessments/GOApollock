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

#' Plot Pearson residual matrix
#'
#' @param mat Matrix of Pearson residuals (years by ages) as
#'   returned by \code{read_pk_rep}
#' @param years Vector of years of observations
#'
#' @return A long data frame ready for ggplot
#' @export
plot_resids <- function(mat, years, minyr=NULL){
  melt_resids <- function(mat, years){
    resids <- mat %>% data.frame() %>% cbind(year=years)%>%
      pivot_longer(-year, names_prefix='X', names_to='age',
                   values_to='resid', names_transform=list(age=as.numeric))
    resids
  }
  x <- melt_resids(mat, years) %>% mutate(resid=ifelse(resid==0,NA,resid))
  if(!is.null(minyr)) x <- filter(x, year>=minyr)
  g <- ggplot(x, aes(year, age, size=sqrt(abs(resid)), fill=resid<0)) +
    geom_point(alpha=.8, pch=21, fg=gray(.5) ) +
    scale_y_continuous(breaks=1:10)
  rmin <- min(x$resid, na.rm=TRUE) %>% round(1)
  rmax <- max(x$resid, na.rm=TRUE) %>% round(1)
   g <- g+ggtitle(label=paste0("Pearson residual range: ", rmin, " to ", rmax))+
     labs(x=NULL, y="Age") + theme(legend.position='none')
  g + theme(plot.title=element_text(hjust = 1))
}


#' Plot observed vs expected compositions as lines
#'
#' @param mat Matrix of observed and expected compositions as
#'   returned by \code{read_pk_rep}
#' @param years Vector of years with data
#' @param ncol Number of columns for faceting
#' @param title Plot title
#' @param minyr Optional minimum year to plot
#' @param minage Optional minimum age to plot
#' @return Prints a ggplot
#' @export
plot_obs_exp <- function(mat, years, ncol=5, title, minyr=NULL, minage=NULL){
  get_pos <- function(x) x %>% group_by(year) %>% summarize(x=9, maxp=max(proportion))
  melt_obs_exp <- function(mat, years){
    x <- mat %>% data.frame() %>%
      setNames(c(paste0('obs', 1:10), paste0('exp', 1:10))) %>%
      cbind(year=years)
    x <- x %>% pivot_longer(-year, names_to=c('type'), values_to='proportion') %>%
      mutate(age=as.numeric(stringr::str_sub(type,4,5)),
             type=stringr::str_sub(type, 1,3))
    ## fake data point ot control ylimit height
    x0 <- group_by(x, year) %>% summarize(age=NA, proportion=max(proportion)*1.1) %>% ungroup
    bind_rows(x,x0)
  }
  x <- melt_obs_exp(mat, years)
  pos <- get_pos(x)                     # breaks if after filtering
  if(!is.null(minyr)) x <- filter(x, year>=minyr)
  if(!is.null(minage)) x$proportion[x$age < minage & !is.na(x$age)] <-  NA
  g <- ggplot(x, aes(age,proportion, color=type, lty=type, label=year)) +
    geom_line(lwd=.8)+
    geom_point(cex=1) +
    facet_wrap("year", dir='v',ncol=ncol, scales='free_y') +
    scale_x_continuous(breaks=1:10) +
    geom_text(data=pos, size=2.5, aes(x=x, lty=NULL,y=maxp), vjust = 1, color='black') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border =element_rect(color = gray(.8)),
          panel.spacing=unit(0,'cm'),
          strip.text.x = element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position='none')+
    labs(x="Age") + ggtitle(title)
  g
}
