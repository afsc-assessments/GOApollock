gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Plot selectivity (persp, terminal, and time-varying plots)
#'
#' @param sdrep list or single sdreport object
#' @param model_names names of models in list
#'
#' @return
#' @export
#'
#' @examples
plot_selectivity <- function(sdrep, model_names = NULL){

  # Convert single one into a list
  if(class(sdrep) == "sdreport"){
    sdrep <- list(sdrep)
  }

  # Objects
  years <- 1970:(2022+5)
  ages <- 1:10

  sel_list <- list()
  sel_lwr_list <- list()
  sel_upr_list <- list()


  # Matrix to hold objects
  for(mod in 1:length(sdrep)){
    sel_lwr <- matrix(NA, length(years), length(ages))
    colnames(sel_lwr) <- ages
    rownames(sel_lwr) <- years
    logsel <- logsel_sd <- sel <- sel_sd <- sel_upr <- sel_lwr

    # - Get sel value and SD
    sel_ind <- which(names(sdrep[[mod]]$value) == "slctfsh_logit")
    logsel <- replace(logsel, values = sdrep[[mod]]$value[sel_ind])
    logsel_sd <- replace(logsel, values = sdrep[[mod]]$sd[sel_ind])

    # - Convert to natural scale
    invlogit <- function(x) 1/(1+exp(-x))

    sel = invlogit(logsel)
    sel_upr = invlogit(logsel + 1.96 * logsel_sd)
    sel_lwr = invlogit(logsel - 1.96 * logsel_sd)

    # FIXME: bounding issue
    sel_lwr[is.nan(sel_lwr)] <- 1
    sel_upr[is.nan(sel_upr)] <- 1
    sel[is.nan(sel)] <- 1

    # - Asign to list
    sel_list[[mod]] <- sel
    sel_upr_list[[mod]] <- sel_upr
    sel_lwr_list[[mod]] <- sel_lwr

    # - Sel on natural scale
    sel_ind <- which(names(sdrep[[mod]]$value) == "slctfsh")
    sel <- replace(sel, values = sdrep[[mod]]$value[sel_ind])
    sel_sd <- replace(sel_sd, values = sdrep[[mod]]$sd[sel_ind])
    sel_upr <- sel + 1.96 * sel_sd
    sel_lwr <- sel - 1.96 * sel_sd

    # PLOT 1: Perspective plot ----
    par(
      mar = c(3.2, 3.2 , 1 , 0.5) ,
      oma = c(0 , 0 , 0 , 0),
      tcl = -0.8,
      mgp = c(10, 0.6, 0)
    )

    persp(
      y = years,
      x =  ages,
      z = t(sel),
      col = "white",
      xlab = "Age",
      ylab = "\n\nYear",
      zlab = "\n\nSelectivity",
      expand = 0.5,
      box = TRUE,
      ticktype = "detailed",
      phi = 35,
      theta = -19,
      main = model_names[mod]
    )
  }


  # PLOT 2: Terminal selectivity ----

  # - Get range of terminal sel
  ylim_sel <- 0
  for (mod in 1:length(sdrep)) {
    ylim_sel = range(c(ylim_sel,
                       tail(sel_upr_list[[mod]], n = 1),
                       tail(sel_lwr_list[[mod]], n = 1),
                       tail(sel_list[[mod]], n = 1)))
  }


  # - Plot setup
  par(
    mar = c(3.2, 3.2 , 0.5 , 0.5) ,
    oma = c(0 , 0 , 0 , 0),
    tcl = -0.35,
    mgp = c(1.75, 0.5, 0)
  )

  plot(
    y = NA,
    x = NA,
    ylim = ylim_sel,
    xlim = c(min(ages),  max(ages)),
    xlab = "Age",
    ylab = "Terminal selectivity"
  )

  mod_colors <- gg_color_hue(length(sdrep))

  # - Selectivity
  for (mod in 1:length(sdrep)) {

    # - 95% CI
    polygon(
      x = c(ages, rev(ages)),
      y = c(tail(sel_upr_list[[mod]], n = 1), rev(tail(sel_lwr_list[[mod]], n = 1))),
      col = adjustcolor( mod_colors[mod], alpha.f = 0.8/2),
      border = NA
    )

    # - Mean
    lines(
      x = ages,
      y = tail(sel_list[[mod]], n = 1),
      lty = 1,
      lwd = 1.5,
      col = mod_colors[mod]
    )
  }

  # PLOT 3: Time-varying by age
  #TODO
  # g <- ggplot(tmp, aes(year, est, ymin=lwr, ymax=upr, color=version,
  #                      fill=version))+
  #   geom_ribbon(alpha=.5) + geom_line(alpha=.8) +
  #   facet_wrap('name',nrow=2) + labs(x=NULL,y='Selectivity',
  #                                    color=NULL, fill=NULL)
  # ggsave('plots/3_fshsel_sel.png', g, width=5, height=3.5)
  # sel <- mymelt(list(r1, r21_9), 'Fishery_selectivity') %>%
  #   filter(year %in% c(1980, 1993, 2006, 2020)) %>%
  #   mutate(yearf=as.factor(year))
  # g <- ggplot(sel, aes(age, value, color=model)) + geom_line() +
  #   facet_wrap('year', nrow=3) +
  #   geom_point()+
  #   geom_vline(xintercept=4, color=gray(.5), lty=3) +
  #   scale_x_continuous(breaks=1:10) +
  #   labs(y='Selectivity', color=NULL)+
  #   theme(legend.position='top')
  # ggsave('plots/3_fshsel_sel_examples.png', g, width=5, height=3.5)
}
