
#' Plot overview of data sources
#'
#' @param datlist Datlist as returned by \code{read_pk_dat}.
#' @return Makes plot and invisibly returns the processed tidy
#'   data frame
#' @details This works well with png(..., width=7, height=5,
#'   units='in')
#' @export
plot_data_overview <- function(datlist){
  dd <- datlist
  x1 <- with(dd, data.frame(year=styr:endyr, size=cattot, survey='fishery', type='catch'))
  x2 <- with(dd, data.frame(year=fshyrs, size=multN_fsh, survey='fishery', type='ages'))
  x3 <- with(dd, data.frame(year=srvyrs1, size=indxsurv_log_sd1, survey='Shelikof', type='age 3+ index'))
  x4 <- with(dd, data.frame(year=srv_acyrs1, size=multN_srv1, survey='Shelikof', type='ages'))
  x5 <- with(dd, data.frame(year=srvyrs2, size=indxsurv_log_sd2, survey='NMFS BT', type='index'))
  x6 <- with(dd, data.frame(year=srv_acyrs2, size=multN_srv2, survey='NMFS BT', type='ages'))
  x7 <- with(dd, data.frame(year=srvyrs3, size=indxsurv_log_sd3, survey='ADF&G', type='index'))
  x8 <- with(dd, data.frame(year=srv_acyrs3, size=multN_srv3, survey='ADF&G', type='ages'))
  x9 <- with(dd, data.frame(year=srvyrs4, size=indxsurv_log_sd4, survey='Shelikof', type='age 1 index'))
  x10 <- with(dd, data.frame(year=srvyrs5, size=indxsurv_log_sd5, survey='Shelikof', type='age 2 index'))
  x11 <- with(dd, data.frame(year=srvyrs6, size=indxsurv_log_sd6, survey='Summer AT', type='index'))
  x12 <- with(dd, data.frame(year=srv_acyrs6, size=multN_srv6, survey='Summer AT', type='ages'))
  ## lengths special cases
  x13 <- with(dd, data.frame(year=srv_lenyrs2, size=multNlen_srv2, survey='NMFS BT', type='lengths'))
  x14 <- with(dd, data.frame(year=srv_lenyrs6, size=multNlen_srv6, survey='Summer AT', type='lengths'))
  dat <- rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14)
  dat <- group_by(dat, survey, type) %>% mutate(relsize=size/max(size)) %>% ungroup
  dat <- mutate(dat, id=paste(survey, type, sep='_'))
  ## png('data.png', width=7, height=5, units='in', res=400)
  size.cex <- 1
  maxsize <- 2
  ymax <- 19
  xmid <- mean(unique(dat$year))
  par(mar=c(1.75,.75,.75,.75), mgp=c(1.5,.35,0), tck=-.01)
  plot(0, xlim = c(min(dat$year), dd$endyr+14), ylim = c(-1, ymax+1), axes = FALSE,  yaxs = "i",
       type = "n", xlab = NA, ylab = "" )
  box()
  mycircles <- function(survey, type,y, color, lab=type){
    text(x=dd$endyr+2, y=ymax-y, label=lab, pos=4)
    xx <- dat[dat$survey==survey &  dat$type==type & dat$size>0,]
    if(nrow(xx)>0){
      symbols(x=xx$year, y=rep(ymax-y, length(xx$year)), circles=sqrt(xx$relsize),
              bg = adjustcolor(color, alpha.f=.6),
              add = TRUE, inches = .07)
    }
  }
  cols <- c(rgb(127,201,127,max=256),rgb(190,174,212,max=256),rgb(253,192,134,max=256),rgb(255,255,153,max=256),rgb(56,108,176, max=256))
  text(x=xmid, y=ymax-.15, labels='Fishery', font=2, cex=1.1)
  mycircles(survey='fishery', type='catch', y=1, color=cols[1], lab='Catch')
  mycircles(survey='fishery', type='ages', y=2, color=cols[1], lab='Age Comps')
  text(x=xmid, y=ymax-3.15, labels='Shelikof', font=2, cex=1.1)
  mycircles(survey='Shelikof', type='ages', y=4, color=cols[2], lab='Age comps')
  mycircles(survey='Shelikof', type='age 3+ index', y=5, color=cols[2], lab='Age 3+ index')
  mycircles(survey='Shelikof', type='age 1 index', y=6, color=cols[2], lab='Age 1 index')
  mycircles(survey='Shelikof', type='age 2 index', y=7, color=cols[2], lab='Age 2 index')
  text(x=xmid, y=ymax-8.15, labels='Summer AT', font=2, cex=1.1)
  mycircles(survey='Summer AT', type='index', y=9, color=cols[3], lab='Index')
  mycircles(survey='Summer AT', type='ages', y=10, color=cols[3], lab='Age comps')
  mycircles(survey='Summer AT', type='lengths', y=11, color=cols[3], lab='Length comps')
  text(x=xmid, y=ymax-12.15, labels='NMFS BT', font=2, cex=1.1)
  mycircles(survey='NMFS BT', type='index', y=13, color=cols[4], lab='Index')
  mycircles(survey='NMFS BT', type='ages', y=14, color=cols[4], lab='Age comps')
  mycircles(survey='NMFS BT', type='lengths', y=15, color=cols[4], lab='Length comps')
  text(x=xmid, y=ymax-16.15, labels='ADF&G BT', font=2, cex=1.1)
  mycircles(survey='ADF&G', type='index', y=17, color=cols[5], lab='Index')
  mycircles(survey='ADF&G', type='ages', y=18, color=cols[5], lab='Age comps')
                                       #abline(v=dd$endyr, type=3, col=gray(.5))
  axis(1, at=seq(1970,dd$endyr, by=5))
  return(invisible(dat))
}


#' Plot survey selectivities with +/- 1 SE
#' @param x A list of data frames as read in from
#'   \link{\code{read_pk_std}}
#'
#' @param add_uncertainty Whether to add +/- SE
#'   intervals. These are calculated in logit space.
#' @param plot_logit Whether to plot in logit space or not
#' @param plot Whether to plot and return the ggplot object
#'   (TRUE) or return the data (FALSE).
#'
#' @return Either a ggplot object or the data used for it,
#' depending on \code{plot} argument.
#' @export
#'
plot_pk_selex <- function(x, add_uncertainty=TRUE, add_fishery=TRUE,
                          plot_logit=FALSE, plot=TRUE){
  x <-  bind_rows(x)
  ## redefine these to clear up plot
  x <- mutate(x, lwr=est-se, upr=est+se)
  svys <- x %>%
    filter(name!='slctfsh_logit' & grepl('_logit',name)) %>%
    mutate(survey=gsub("slctsrv","Survey ", name),
           survey=gsub('_logit', '',survey)) %>%
    filter(is.finite(est))
  if(nrow(svys)==0)
    stop("No survey selex found, maybe from old model version?")
  if(add_fishery){
    fsh <-
      filter(x, name=='slctfsh_logit' & is.finite(est)) %>%
      mutate(survey='Fishery (2021)')
    if(nrow(fsh)==0)
      stop("No fishery selex found, maybe from old model version?")
    x <- bind_rows(svys,fsh)
  } else {
    x <- svs
  }
  ## plot in logit space?
  if(!plot_logit){
    f <- function(x) 1/(1+exp(-x))
    x <- mutate(x, est=f(est), lwr=f(lwr), upr=f(upr))
  }
  alpha1 <- .2 # for ribbon
  alpha2 <- 1 # for lines
  ylab <- ifelse(plot_logit, 'logit selectivity','Selectivity')
  if(add_uncertainty) ylab <- paste(ylab, '(+/- 1 SE)')
  g <- ggplot(x, aes(i,est, ymin=lwr, ymax=upr, color=survey, fill=survey))+
    geom_line(alpha=alpha2, lwd=1) + geom_point(alpha=alpha2)+
    facet_wrap('version')+   scale_x_continuous(breaks=1:10)+
    labs(fill=NULL, color=NULL, x='Age', y=ylab)
  if(add_uncertainty) g <-  g+ geom_ribbon(alpha=alpha1, color=NA)
  if(plot)  {print(g); g}  else  return(x)
}

#' Plot spawning biomass
#'
#' @param x A list of data frames as read in from
#'   \link{\code{read_pk_cor}}
#' @param add_uncertainty Whether to add 95% confidence
#'   intervals. These are calculated in log space.
#' @param uselog Whether to calculate in log space or not
#' @param plotlog Whether to plot in log space or not
#' @param plot Whether to plot and return the ggplot object
#'   (TRUE) or return the data (FALSE).
#' @param alpha1 Transparency for ribbons
#' @param alpha2 Transparency for lines
#'
#' @return Either a ggplot object or the data used for it,
#' depending on \code{plot} argument.
#' @export
#'
plot_pk_ssb <- function(x, add_uncertainty=TRUE, plotlog=TRUE,
  uselog=FALSE, plot=TRUE, alpha1=.5, alpha2=.8){
  nmods <- length(x)
  x <- bind_rows(x)
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
  if(!plotlog) g <- g+ylim(0,NA)
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
  if(!is.null(minyr)) x <- filter(x, year>=minyr)
  pos <- get_pos(x)                     # breaks if after filtering
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
