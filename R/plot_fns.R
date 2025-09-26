
#' Plot overview of data sources
#'
#' @param datlist Datlist as returned by \code{read_dat}.
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
#' @param fits A list of data frames as read in from
#'   \code{read_std}
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
plot_selex <- function(fits, add_uncertainty=TRUE, add_fishery=TRUE,
                          plot_logit=FALSE, plot=TRUE){
  if(is.pkfit(fits)){
    nmods <- 1
    labs <- fits$version
  } else {
    nmods <- length(fits)
    labs <- sapply(fits, function(x) x$version[1])
  }
  sel <- get_std(fits) %>% filter(grepl('_logit',name)) %>%
  ## redefine these to clear up plot
    mutate(lwr=est-se, upr=est+se) %>%
    group_by(version, name) %>%
    mutate(age=1:n()) %>% ungroup
  svys <- sel %>%
    filter(name!='slctfsh_logit') %>%
    mutate(survey=gsub("slctsrv","", name),
           survey=gsub('_logit', '',survey),
           surveyf=surveyf(as.numeric(survey))) %>%
    filter(is.finite(est))
  if(nrow(svys)==0)
    stop("No survey selex found, maybe from old model version?")
  if(add_fishery){
    thisyear <- max(get_rep(fits, 'endyr')$value)
    fsh <- filter(sel, name=='slctfsh_logit' & is.finite(est)) %>%
      mutate(surveyf=paste0('Fishery (',thisyear-1,')'))
    if(nrow(fsh)==0)
      stop("No fishery selex found, maybe from old model version?")
    x <- bind_rows(svys,fsh)
  } else {
    x <- svys
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
  g <- ggplot(x, aes(age,est, ymin=lwr, ymax=upr, color=surveyf, fill=surveyf))+
    geom_line(alpha=alpha2, lwd=1) + geom_point(alpha=alpha2)+
    facet_wrap('version')+   #scale_x_continuous(breaks=1:10)+
    labs(fill=NULL, color=NULL, x='Age', y=ylab)
  if(add_uncertainty) g <-  g+ geom_ribbon(alpha=alpha1, color=NA)
  if(plot)  {print(g); g}  else  return(x)
}

#' Plot spawning biomass
#'
#' @param fits A model fit or list of model fits as produced by
#'   \code{fit_tmb}.
#' @param add_uncertainty Whether to add 95% confidence
#'   intervals. These are calculated in log space depending on
#'   \code{uselog}.
#' @param uselog Whether to calculate in log space or not
#' @param plotlog Whether to plot SSB in log space or not
#' @param plot Whether to plot and return the ggplot object
#'   (TRUE) or return the data used in the plot(FALSE).
#' @param alpha1 Transparency for ribbons
#' @param alpha2 Transparency for lines
#'
#' @return Either a ggplot object or the data used for it,
#' depending on \code{plot} argument.
#' @export
#'
plot_ssb <- function(fits, add_uncertainty=TRUE, plotlog=FALSE,
                        uselog=FALSE, addproj=FALSE,
                        plot=TRUE, alpha1=.5, alpha2=.8){
  if(is.pkfit(fits)){
    nmods <- 1
    labs <- fits$version
  } else {
    nmods <- length(fits)
    labs <- sapply(fits, function(x) x$version[1])
  }
  ses <- get_std(fits)
  ## get the order to match what the user puts in the fits list
  ses$version <- factor(ses$version, levels=labs)
  tmp <- 'Espawnbio'
  if(addproj & uselog) {
    warning("Can't have uselog and addproj both, setting uselog=FALSE")
    uselog <- FALSE
  }
  if(uselog)  tmp <- 'Espawnbio_log'
  if(uselog & nrow(filter(ses, name==tmp))==0){
    tmp <- 'Espawnbio'
    uselog <- FALSE
    warning("Log of SSB not found in input, likely from old model run. Setting uselog=FALSE.")
  }
  if(uselog){
    ssb <- ses %>% filter(name==tmp) %>%
      mutate(est=exp(est), lwr=exp(lwr), upr=exp(upr))
  } else {
    ## assume they are lognormal
    ssb <- ses %>% filter(name==tmp) %>%
      mutate(lwr=est/exp(1.96*sqrt(log(1+(se/est)^2))),
             upr=est*exp(1.96*sqrt(log(1+(se/est)^2))))
  }
  ## tack on the proj years, have to get the proj years manually
  ## and by model since may differ later
  if(addproj){
    ssb <- ses %>% filter(name=='Espawnbio_proj') %>%
      mutate(year=1+year-min(year))%>%
      bind_rows(ssb) %>% group_by(version) %>%
      mutate(maxyear=max(year)) %>%
      mutate(year=ifelse(name=='Espawnbio', year, year+maxyear)) %>% ungroup
    maxyear <- filter(ssb, name=='Espawnbio') %>% summarize(max(year)) %>% pull
  }
  g <- ggplot(ssb, aes(year, est, ymin=lwr, ymax=upr, color=version, fill=version)) +
    geom_line(alpha=alpha2, lwd=1) +
    theme_bw()+ labs(x=NULL,y='Spawning biomass (Mt)')
  if(plotlog) g <- g+ scale_y_log10()
  if(!plotlog) g <- g+ylim(0,NA)
  if(addproj) g <-  g+geom_vline(xintercept=maxyear)
  if(add_uncertainty) g <-  g+ geom_ribbon(alpha=alpha1)
  if(plot)  {print(g); g}  else  return(ssb)
}

#' Plot Pearson residual matrix
#'
#' @param mat Matrix of Pearson residuals (years by ages) as
#'   returned by \code{read_rep}
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
#'   returned by \code{read_rep}
#' @param years Vector of years with data
#' @param ncol Number of columns for faceting
#' @param title Plot title
#' @param minyr Optional minimum year to plot
#' @param minage Optional minimum age to plot
#' @return Prints a ggplot
#' @export
plot_obs_exp <- function(mat, years, ncol=5, title, minyr=NULL,
                         minage=NULL){
  maxage <- ncol(mat)/2
  get_pos <- function(x) x %>% group_by(year) %>% summarize(x=9, maxp=max(proportion))
  melt_obs_exp <- function(mat, years){
    x <- mat %>% data.frame() %>%
      setNames(c(paste0('obs', 1:maxage), paste0('exp', 1:maxage))) %>%
      cbind(year=years)
    x <- x %>% pivot_longer(-year, names_to=c('type'), values_to='proportion') %>%
      mutate(age=as.numeric(stringr::str_sub(type,4,5)),
             type=stringr::str_sub(type, 1,3))
    ## fake data point ot control ylimit height
    x0 <- group_by(x, year) %>% summarize(age=NA, type='obs', proportion=max(proportion)*1.1) %>% ungroup
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
    scale_x_continuous(breaks=1:maxage) +
    geom_text(data=pos, size=2.5, aes(x=x, lty=NULL,y=maxp, type=NULL), vjust = 1, color='black') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border =element_rect(color = gray(.8)),
          panel.spacing=unit(0,'cm'),
          strip.text.x = element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position='right')+
    labs(x="Age") + ggtitle(title)
  g
}

#' Plot cohort contribution to SSB as a stacked ribbon plot
#'
#' @param fit A model fit as produced by \code{fit_tmb}
#' @param type The type of biomass to plot, options are 'ssb',
#'   'total', or 'summary' which are spawning, total, and age 2+,
#'   respectively. Default is SSB.
#' @param plot Whether to create a plot and invisibly return the
#'   ggplot object,  or return a data.frame of the data.
#'
#' @details The plus group breaks the idea of a cohort a little
#'   bit. Here it assumes that the cohort dies when hitting the
#'   plus group.
#' @export
#' @return A ggplot or data frame depending on \code{plot}.
plot_cohort_biomass <- function(fit, type=c('ssb', 'total', 'summary'),
                            plot=TRUE, print=FALSE){
  stopifnot(is.pkfit(fit))
  type <- match.arg(type)
  naa <- fit$rep$N
  swaa <- fit$input$dat$wt_srv1
  mat <- fit$input$dat$mat
  ssb <- t(apply(naa,1, function(x) x*mat))*swaa/2
  dimnames(ssb) <- dimnames(naa) <- list(year=fit$rep$years, age=1:ncol(naa))
  totb <- naa*fit$input$dat$wt_srv2
  sumb <- totb; sumb[,1] <- 0
  mat2longdf <- function(x){
    x <- as.data.frame(x)
    x$year <- row.names(x)
    y <- pivot_longer(x, -year) %>%
      mutate(age=as.numeric(name), year=as.numeric(year),
             cohort=year-age) %>% select(-name)
    ## need to merge on empty values for cohorts in each year
    z <- expand.grid(year=unique(y$year), cohort=unique(y$cohort), value2=0)
    tmp <- merge(y,z, all.y=TRUE) %>%
      mutate(value=ifelse(is.na(value),0,value), ssb=value+value2) %>%
      arrange(year,cohort) %>% group_by(year) %>%
      mutate(pct.ssb=ssb/sum(ssb)) %>% ungroup() %>%
      select(year, ssb, pct.ssb, cohort)
    tmp %>% pivot_longer(-c(year,cohort))
  }

  if(type=='ssb'){
    x <- mat2longdf(ssb)
    labs <- c('Spawning biomass (M t)', '% Spawning biomass')
  } else if(type=='total') {
    x <- mat2longdf(totb)
    labs <- c('Total biomass (M t)', '% Total biomass')
  } else {
    x <- mat2longdf(sumb)
    labs <- c('Summar biomass (M t)', '% Summary biomass')
  }
  x <- mutate(x, name=factor(name, levels=c('ssb', 'pct.ssb'),
                         labels=labs), version=fit$version)
  g <- ggplot(x, aes(x=year, y=value, group=factor(cohort))) +
    geom_area(color=1, linewidth=.25, fill='white') +
    facet_wrap('name', ncol=1, scales='free_y') +
    labs(x=NULL,y=NULL)
  if(plot) return(g)
  return(x)
}

#' Plot retrospective peels against the base model with CI
#' @param retros A list of pkfits
#' @param type The metric to use: 'ssb', 'F', 'recruit'
#' @export
plot_retros <- function(retros, type=c('ssb', 'F', 'recruit'), title=NULL){
  type <- match.arg(type)
  version <- retros[[1]]$version
  thisyear <- retros[[1]]$rep$endyr
  if(type=='ssb'){
    x <- get_rep(retros, 'Espawnbio')
    y <- get_std(retros, 'Espawnbio')
    lab <- 'Spawning biomass (Mt)'
  } else if(type=='F'){
    x <- get_rep(retros, 'F')
    y <- get_std(retros, 'F')
    lab <- 'Fishing Effort (F)'
  } else if(type=='recruit'){
    x <- get_rep(retros, 'recruit')
    y <- get_std(retros, 'recruit')
    lab <- 'Recruitment (billions)'
  } else { stop('something wrong with type=',type)}
  y <- mutate(y, lwr=pmax(lwr,0.001)) %>% filter(version=='peel0')
  stopifnot(nrow(x)>0)
  stopifnot(nrow(y)>0)
  rho <- calculate_rho(retros,type=type)
  rho.lab <- paste0("Mohn's rho= ", round(rho,3))
  ## Do some heavy processing to plot quickly
  x <- x %>%
    mutate(version=as.numeric(gsub('peel','', version))) %>%
    rename(peel=version) %>% group_by(year) %>%
    mutate(Pct_Diff=100*(value-value[peel==0])/value[peel==0]) %>%
    ungroup()
  ## Plot it
  g1 <- ggplot(y, aes(year, est, group=NULL, fill=NULL, color=NULL,
                      ymin=lwr, ymax=upr)) + geom_ribbon(alpha=.25)
  g1 <- g1+geom_line(data=x, aes(year, value, group=peel, ymax=NULL,
                                 ymin=NULL, fill=NULL, color=factor(peel)))
  x2 <- filter(x, thisyear-peel==year)
  g1 <- g1 +
    geom_point(data=x2, aes(y=value, color=factor(peel), ymin=NULL,ymax=NULL), size=2) +
    theme(legend.position='none') +
    annotate('label', x=2010,y=max(y$upr), label=rho.lab) +
    labs(x=NULL, y=lab, title=title)
  if(type=='recruit') g1 <- g1+scale_y_log10()
  g1
}

#' Plot effective sample sizes to compare Francis tuning to Dirichlet-multinomial (DM)
#' @param fitDM A model fit with the DM turned on
#' @param fitfrancis A model fit with Francis tuning
#' @param addISS Whether to add the input sample size (ISS) taken from the fitDM model.
#' @param plot Whether to print.
#' @return A ggplot object if plot is TRUE, otherwise the data
#' @export
plot_ess <- function(fitDM, fitfrancis=NULL, addISS=FALSE, plot=TRUE){
  d <- fitDM$input$dat
  r <- fitDM$rep
  f <- fitfrancis$input$dat
  x0 <- data.frame(survey=0, year=d$fshyrs, ISS=d$multN_fsh, DM=r$ESS_fsh)
  x1 <- data.frame(survey=1, year=d$srv_acyrs1, ISS=d$multN_srv1, DM=r$ESS_srv1)
  x2 <- data.frame(survey=2, year=d$srv_acyrs2, ISS=d$multN_srv2, DM=r$ESS_srv2)
  x3 <- data.frame(survey=3, year=d$srv_acyrs3, ISS=d$multN_srv3, DM=r$ESS_srv3)
  x6 <- data.frame(survey=6, year=d$srv_acyrs6, ISS=d$multN_srv6, DM=r$ESS_srv6)
  if(!is.null(fitfrancis)){
    x0 <- cbind(x0, Francis=f$multN_fsh)
    x1 <- cbind(x1, Francis=f$multN_srv1)
    x2 <- cbind(x2, Francis=f$multN_srv2)
    x3 <- cbind(x3, Francis=f$multN_srv3)
    x6 <- cbind(x6, Francis=f$multN_srv6)
  }
  x <-bind_rows(x0,x1,x2,x3,x6)
  if(!addISS) x <- select(x, -ISS)
  g <- x %>% pivot_longer(-c(survey,year)) %>%
    filter(value>0) %>% mutate(survey=surveyf(survey)) %>%
    ggplot(aes(year,value, color=name)) + geom_line() +
    facet_wrap('survey') + ylim(0,NA) +
    labs(x=NULL, y='Effecive sample size (ESS)')
  if(plot) return(g)
  return(x)
}

#' Plot fits to biomass indices for a single fit
#' @param fit A fitted model object
#' @export
plot_index_fits <- function(fit){
  survey.cex <- .9
  dd <- fit$input$dat
  rr <- fit$rep
  years <- rr$years
  thisyear <- max(years)
  ## index fits
  addsegs <- function(yrs, obs, CV){
    getlwr <- function(obs, CV) qlnorm(p=.025, meanlog=log(obs), sdlog=sqrt(log(1+CV^2)))
    getupr <- function(obs, CV) qlnorm(p=.975, meanlog=log(obs), sdlog=sqrt(log(1+CV^2)))
    segments(yrs, y0=getlwr(obs,CV), y1=getupr(obs,CV))
    points(yrs, obs, pch=22, bg='white')
  }
  par(mfrow=c(4,1), mar=c(.2,3,.2,.2), oma=c(2,0,0,0),
      mgp=c(1.5,.25,0), tck=-.01)
  plot(years, rr$Eindxsurv1, type='l', xlim=c(1987, thisyear),
       ylim=c(0,2.1),  xlab=NA, xaxt='n',
       ylab='')
  addsegs(yrs=dd$srvyrs1, obs=dd$indxsurv1, CV=dd$indxsurv_log_sd1)
  mtext('Shelikof acoustic survey',  line=-1.5, cex=survey.cex)
  plot(years, rr$Eindxsurv6, type='l',xaxt='n',
       xlim=c(1987, thisyear),
       ylim=c(0,2.7),  xlab=NA, ylab='')
  addsegs(yrs=dd$srvyrs6, obs=dd$indxsurv6, CV=dd$indxsurv_log_sd6)
  mtext('Summer acoustic survey',  line=-1.5, cex=survey.cex)
  plot(years, rr$Eindxsurv2, type='l', xlim=c(1987, thisyear),
       ylim=c(0,1.5),  xlab=NA, xaxt='n',
       ylab='')
  addsegs(yrs=dd$srvyrs2, obs=dd$indxsurv2, CV=dd$indxsurv_log_sd2)
  mtext('NMFS bottom trawl survey',  line=-1.5, cex=survey.cex)
  plot(years, rr$Eindxsurv3, type='l', xlim=c(1987, thisyear),
       ylim=c(0,.5), xlab=NA, ylab='')
  addsegs(yrs=dd$srvyrs3, obs=dd$indxsurv3, CV=dd$indxsurv_log_sd3)
  mtext('ADF&G bottom trawl survey',  line=-1.5, cex=survey.cex)
  legend('topright', legend=c('Expected', 'Observed'), lty=c(1,NA),
         pch=c(NA,0), bty='n')
  mtext('Biomass (Mt)', outer=TRUE, side=2, line=-1.5, cex=.8)
}

#' Plot OSA fits via afscOSA for fishery and all surveys
#' @param fit Fitted object
#' @export
#' @details This is just a wrapper around run_osa
plot_osa_fits <- function(fit){
  library(afscOSA)
  theta <- exp(fit$parList$log_DM_pars)
  fs <- run_osa(obs = fit$rep$res_fish[,2:10],
                exp = fit$rep$res_fish[,12:20],
                N=fit$input$dat$multN_fsh,
                index=2:10, index_label = 'Age',
                fleet='Fishery',
                theta = theta[1],
                year=fit$input$dat$fshyrs)

  s1 <- run_osa(obs = fit$rep$res_srv1[,3:10],
                exp = fit$rep$res_srv1[,13:20],
                N=fit$input$dat$multN_srv1,
                index=3:10, index_label = 'Age',
                fleet='Shelikof',
                theta = theta[2],
                year=fit$input$dat$srv_acyrs1)
  s2 <- run_osa(obs = fit$rep$res_srv2[,1:10],
                exp = fit$rep$res_srv2[,11:20],
                N=fit$input$dat$multN_srv2,
                index=1:10, index_label = 'Age',
                fleet='NMFS BT',
                theta = theta[3],
                year=fit$input$dat$srv_acyrs2)
  ind <- which(rowSums(fit$rep$res_srv3[,1:10])>0)
  s3 <- run_osa(obs = fit$rep$res_srv3[ind,1:10],
                exp = fit$rep$res_srv3[ind,11:20],
                N=fit$input$dat$multN_srv3[ind],
                index=1:10, index_label = 'Age',
                fleet='ADF&G BT',
                theta = theta[4],
                year=fit$input$dat$srv_acyrs3[ind])
  s6 <- run_osa(obs = fit$rep$res_srv6[,1:10],
                exp = fit$rep$res_srv6[,11:20],
                N=fit$input$dat$multN_srv6,
                index=1:10, index_label = 'Age',
                fleet='Summer AT',
                theta = theta[5],
                year=fit$input$dat$srv_acyrs6)
  afscOSA::plot_osa(list(fs, s1,s2,s3,s6))
}
