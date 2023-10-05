## Script to run the projections for various models
library(tidyverse)
source("~/GitHub/GOApollock/R/proj_fns.R", echo=TRUE)

## Final 2022 model

abc_calc_tmb <- function(replist, datlist){
  ## New automated way using the new 'spm' package. Confirmed the
  ## same as manually doing it with 'proj' below.
  write_spm_inputs_tmb(replist, datlist, path = 'TMB/R/proj_calcs/proj')
  prev_dir <- getwd()
  setwd('TMB/R/proj_calcs/proj')
  system("spm")
  setwd(prev_dir)
  bf <- read.csv('TMB/R/proj_calcs/proj/spm_detail.csv')
  exec_table <- get_exec_table_tmb(replist,datlist,bf)
  exec_tableF <- format_exec_table(exec_table)

  proj_scens <- bf %>% select(-Spp) %>%
    group_by(Alternative,Yr)%>%
    summarize_all('mean') %>% arrange(Alternative,Yr) %>% ungroup %>%
    as.data.frame()

  return(list(exec_table = exec_table, exec_tableF = exec_tableF, proj_scens = proj_scens))
}


write_spm_inputs_tmb <- function(replist, datlist, path=getwd()){
  ## delete old files in case it fails it'll stop there
  ff <- file.path(path, c("goa_wp.txt", "spm.dat", "spp_catch.dat"))
  trash <- lapply(ff, function(x) if(file.exists(x)) file.remove(x))
  ayr <- datlist$endyr
  write_spm_setup(ayr, path=path, catch=tail(datlist$cattot,1))
  write_spm_dynamics_tmb(replist, datlist, ayr, path)
  ## old unused file, just need dummy values so it runs
  ## tacpar <- readLines(tacpar.dat'); dput(tacpar)
  write.table(file=file.path(path,'tacpar.dat'),
              x=c("7", "6", " 2.60197 0.417 0.372 0.361 0.296 0.2733 0.125",
                  " -0.300856 -0.741664 -1.26797 -1.78614 -2.15198 -2.39595 -2.59939",
                  " -0.00703968 -0.956592 -1.48909 -2.09122 -2.26005 -2.31812 -2.45321",
                  " 0.372276 -1.35066 -1.76076 -2.22163 -2.31726 -2.28107 -2.34758",
                  " 0.505535 -1.19926 -1.68812 -2.39543 -2.40752 -2.5904 -2.58726",
                  " 0.368722 -1.59025 -1.79548 -2.67422 -2.61358 -2.41092 -2.70204",
                  " 0.308248 -1.78052 -2.23264 -2.85605 -3.39375 -3.05209 -2.94219",
                  " 0.180676 -1.7586 -1.80222 -3.09402 -2.2618 -3.43215 -3.06417"))
}



write_spm_dynamics_tmb <- function(replist, datlist, ayr, path){
  x <- c()
  x <-
    c(x,paste("# proj input file written by write_spm_inputs on", Sys.time()))
  x <- c(x, "goa_wp")
  x <- c(x, "1 # Flag to tell if this is a SSL forage species")
  x <- c(x, "0 # Flag to Dorn's version of a constant buffer")
  x <- c(x, "1 # number of fisheries")
  x <- c(x, "1 # number of sexes")
  x <- c(x, paste(mean(tail(replist$F,5)), " # average F over last 5 years"))
  x <- c(x, "1 # Author's F multiplier (to MaxPermissible)")
  x <- c(x, "0.4 # ABC SPR")
  x <- c(x, "0.35 # OFL SPR")
  x <- c(x, "3.52 # Month of spawning")
  x <- c(x, "10 # Number of ages")
  x <- c(x, "1 # Fratio")
  x <- c(x, "# Natural mortality")
  x <- c(x, paste(replist$M, collapse=' '))
  x <- c(x, "# Maturity")
  x <- c(x, paste(datlist$mat, collapse=' '))
  x <- c(x, "# Spawning WAA")
  x <- c(x, paste(colMeans(tail(datlist$wt_srv1,5))*1000, collapse=' '))
  x <- c(x, "# Fishery WAA")
  x <- c(x, paste(datlist$wt_fsh[datlist$nyrs,]*1000, collapse=' '))
  x <- c(x, "# Fishery selex, averaged over last 5 years")
  ## This is wrong b/c ignores most recent year
  ## x <- c(x, paste(colMeans(tail(replist$Fishery_selectivity,
  ## 5)), collapse= ' '))
  ##  x <- c(x, paste(replist$Mean_fishery_selectivity, collapse='
  ##  '))
  x <- c(x, paste(replist$slctfsh[datlist$nyrs+1,], collapse=' ')) # 2023 (nyrs+1)
  x <- c(x, "# current year starting NAA")
  x <- c(x, paste(tail(replist$N,1)*1000, collapse=' '))
  ind <- which(datlist$yrs %in% 1978:(ayr-1))
  recs <- replist$recruit[ind]*1000
  ssb <- replist$Espawnbio[ind-1]*1000
  x <- c(x, paste(length(recs), " # num of recruits"))
  x <- c(x, "# Recruits from 1978 to this year -1 due to uncertain last year")
  x <- c(x, paste(recs, collapse=' '))
  x <- c(x, '# SSB from 1977 to this year-2 to match recruits')
  x <- c(x, paste(ssb, collapse=' '))
  writeLines(x, con=file.path(path, 'goa_wp.txt'))
}


#' Get table from projection output
#' @param replist A list as read in by \link{read_pk_rep}
#' @param bigfile A data frame of full projection scenario
#'   outputs
#' @export
get_exec_table_tmb <- function(replist, datlist, bigfile, maxABCratio=1){
  ayr <- tail(datlist$yrs,1)
  M <- c(.3,.3)
  means <- bigfile %>% filter(Yr > ayr & Yr <= ayr+2) %>%
    group_by(Alternative, Yr, Spp, SpNo) %>%
    summarize_all(mean) %>% ungroup
  sumbio <- replist[[113]][1:2]*1e6
  ssb <- filter(means,  Alternative==1) %>% pull(SSB) %>%
    round(0)
  if('B0' %in% names(means))
    bfrac <- filter(means, Alternative==1) %>% select(B0,B40,B35) %>%
    round(-3) %>% t
  else {warning("no B0 col"); bfrac=matrix(NA,3,2)}
  fofl <- filter(means,  Alternative==2) %>% pull(FOFL) %>% round(3)
  fabc <-filter(means, Alternative==1) %>% pull(F) %>% round(3)
  maxfabc <- fabc
  ofl <- filter(means,  Alternative==2) %>% pull(OFL) %>% round(0)
  maxabc <- filter(means,  Alternative==1) %>% pull(Catch) %>% round(0)
  abc <- maxabc*maxABCratio
  tab <- rbind(M, sumbio, ssb, bfrac, fofl, maxfabc, fabc,ofl,maxabc,abc)
  tab <- as.data.frame(tab) %>% setNames(c(ayr+1,ayr+2)) %>%
    cbind(name=row.names(tab),.)
  row.names(tab) <- NULL
  return(tab)
}
