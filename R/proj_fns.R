
#' Calculate Fabc for a fitted object.
#'
#' @param fit Fitted TMB model
#' @param ratio The desired target biomass ratio
#' @param F An option fishing mortality rate which can be
#'   specified to get SPR(F), used for status figure
#' @return A list containing SSB100, fabc, abc, btarget
#'   (optimzed) and SPR at F=0, the optimized target, and what is specified by F.
#' @export
#'
get_fabc <- function(fit, ratio=.4, F=NULL){
  ## data inputs
  if(is.null(fit$input))
    stop("Function requires input slot -- add it manually or use updated package version")
  M <- fit$rep$M
  popwt <- fit$rep$wt_pop_proj # weights in kg
  fshwt <- fit$rep$wt_fsh_proj
  spawnwt <- fit$rep$wt_spawn_proj
  mat <- fit$input$dat$mat
  sel <- fit$rep$slctfsh_proj
  maxage <- max(fit$rep$ages)
  stopifnot(all.equal(length(popwt), length(fshwt), length(spawnwt), length(mat), length(sel), maxage))
  if(maxage != 10) warning("SPR function not tested for maxage !=10")
  R0 <- fit$rep$recruit_proj[1]# mean age-1 recruits (billion) used in proj
  stopifnot(!is.null(R0))
  spr <- function(F){
    Z <- M+ F*sel
    N <- rep(NA,maxage)
    N[1] <- R0/2
    for(i in 2:maxage) N[i] <- N[i-1]*exp(-Z[i-1])
    N[maxage] <- N[maxage]/(1-exp(-Z[maxage]))
    ssb <- sum(N*mat*spawnwt*exp(-.21*Z))
    catch <- 2*sum(N*fshwt*(sel*F/Z)*(1-exp(-Z)))
    return(list(ssb=ssb, catch=catch))
  }
  if(!is.null(F)){
    sprF <- spr(F)
  } else {
    sprF <- NULL
  }
  ssb100 <- spr(0)$ssb
  ssq <- function(ratio){
    f <- function(F) 100*(spr(F)$ssb-ssb100*ratio)^2
    opt <- optimize(f=f, interval=c(0,1))
    opt$minimum
  }
  fabc <- ssq(ratio)
  x <- spr(fabc)
  return(list(ssb100=ssb100,
              fabc=fabc,
              abc=x$catch,
              btarget=x$ssb/ssb100,
              spr_at_F0=ssb100/R0,
              spr_at_target=x$ssb/R0,
              spr_at_F=sprF$ssb/R0))
}

#' Get table from projection output
#' @param replist A list as read in by \link{read_rep}
#' @param bigfile A data frame of full projection scenario
#'   outputs
#' @export
get_exec_table <- function(replist, bigfile, maxABCratio=1){
  ayr <- tail(replist$years,1)
  M <- c(.3,.3)
  means <- bigfile %>% filter(Yr > ayr & Yr <= ayr+2) %>%
    group_by(Alternative, Yr, Spp, SpNo) %>%
    summarize_all(mean) %>% ungroup
  sumbio <- tail(replist$Esumbio,2)*1e6
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

#' Format the executive table
#' @param tab Table as read in by \link{get_exec_table}
#' @export
format_exec_table <- function(tab){
  ## kludgy way but its finicky
  tab0 <- tab
  for(i in c(1:5,9:11)){
    for(j in 2:3){
      tab[i,j] <- formatC(tab0[i,j], format='d', big.mark=',')
    }
  }
  tab
}

#' Write input files for the 'spm' projection module. Differs
#' slightly from the old projection module.
#'
#' @param replist A list as read in by \link{read_rep}
#' @param datlist A list as read in by \link{read_dat}
#' @param path Directory to write in
#'
#' @return Nothing but writes files spm.dat, goa_wp.txt, and
#'   tacpar.dat (dummy file) to 'path' based on values in the
#'   replist and datlist objections
#' @export
#'
write_spm_inputs <- function(replist, datlist, path=getwd()){
  ## delete old files in case it fails it'll stop there
  ff <- file.path(path, c("goa_wp.txt", "spm.dat", "spp_catch.dat"))
  trash <- lapply(ff, function(x) if(file.exists(x)) file.remove(x))
  ayr <- tail(replist$years,1)
  write_spm_setup(ayr, path=path, catch=tail(datlist$cattot,1))
  write_spm_dynamics(replist, datlist, ayr, path)
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


write_spm_dynamics <- function(replist, datlist, ayr, path){
  if(is.null(replist$Fishing_mortalities)){
    ## new TMB names don't match so fix here
    replist$Fishing_mortalities <- replist$F
    replist$Natural_mortality <- replist$M
    replist$Projection_spawning_weight_at_age <- replist$wt_spawn_proj
    replist$Projection_fishery_weight_at_age <- replist$wt_fsh_proj
    replist$Projection_fishery_selectivity <- colMeans(tail(replist$slctfsh,5)[-5,]) #replist$slctfsh_proj
    replist$Expected_spawning_biomass <- replist$Espawnbio
    replist$Recruits <- replist$recruit
    replist$Numbers_at_age <- replist$N
  }
  x <- c()
  x <-
  c(x,paste("# proj input file written by write_spm_inputs on", Sys.time()))
  x <- c(x, "goa_wp")
  x <- c(x, "1 # Flag to tell if this is a SSL forage species")
  x <- c(x, "0 # Flag to Dorn's version of a constant buffer")
  x <- c(x, "1 # number of fisheries")
  x <- c(x, "1 # number of sexes")
  x <- c(x, paste(mean(tail(replist$Fishing_mortalities,5)), " # average F over last 5 years"))
  x <- c(x, "1 # Author's F multiplier (to MaxPermissible)")
  x <- c(x, "0.4 # ABC SPR")
  x <- c(x, "0.35 # OFL SPR")
  x <- c(x, "3.52 # Month of spawning")
  x <- c(x, "10 # Number of ages")
  x <- c(x, "1 # Fratio")
  x <- c(x, "# Natural mortality")
  x <- c(x, paste(replist$Natural_mortality, collapse=' '))
  x <- c(x, "# Maturity")
  x <- c(x, paste(datlist$mat, collapse=' '))
  x <- c(x, "# Spawning WAA")
  x <- c(x, paste(replist$Projection_spawning_weight_at_age*1000, collapse=' '))
  x <- c(x, "# Fishery WAA")
  x <- c(x, paste(replist$Projection_fishery_weight_at_age*1000, collapse=' '))
  x <- c(x, "# Fishery selex, averaged over last 5 years")
  ## This is wrong b/c ignores most recent year
  ## x <- c(x, paste(colMeans(tail(replist$Fishery_selectivity,
  ## 5)), collapse= ' '))
 ##  x <- c(x, paste(replist$Mean_fishery_selectivity, collapse='
  ##  '))
  x <- c(x, paste(replist$Projection_fishery_selectivity, collapse=' '))
  x <- c(x, "# current year starting NAA")
  x <- c(x, paste(tail(replist$Numbers_at_age,1)*1000, collapse=' '))
  ind <- which(replist$years %in% 1978:(ayr-1))
  recs <- replist$Recruits[ind]*1000
  ssb <- replist$Expected_spawning_biomass[ind-1]*1000
  x <- c(x, paste(length(recs), " # num of recruits"))
  x <- c(x, "# Recruits from 1978 to this year -1 due to uncertain last year")
  x <- c(x, paste(recs, collapse=' '))
  x <- c(x, '# SSB from 1977 to this year-2 to match recruits')
  x <- c(x, paste(ssb, collapse=' '))
  writeLines(x, con=file.path(path, 'goa_wp.txt'))
}

write_spm_setup <- function(ayr, path, catch){
  x <- c()
  x <- c(x,paste("# proj input file written by write_proj_inputs on",
                 Sys.time()))
  x <- c(x, " std # run name")
  x <- c(x, " 3 # tier")
  x <- c(x, "7 # number of alt scenarios")
  x <- c(x, paste(1:7, collapse=' '))
  x <- c(x, "1 # flag to set TAC=ABC (1 is true)")
  x <- c(x, "2 # SR type (1 ricker, 2 bholt")
  x <- c(x, "1 # proj rec form (default: 1= use obs mean/SD")
  x <- c(x, "0 # SR conditioning (0 means no)")
  x <- c(x, "0.0 # turn off prior thing")
  x <- c(x, "1 # flag to write bigfile")
  x <- c(x, "14 # number of proj years")
  x <- c(x, "1000 # number of simulations")
  x <- c(x, paste(ayr, " # begin year"))
  x <- c(x, "1 # number of years with specified catch")
  x <- c(x, "1 # numberof species")
  x <- c(x, "0 # OY min")
  x <- c(x, "2e6 # OY max")
  x <- c(x, "goa_wp.txt # input file name for biology")
  x <- c(x, "1 # ABC multipliers")
  x <- c(x, "1 # pop scalars")
  x <- c(x, "0.75 # new alt 4 Fabc SPRs (unused?)")
  x <- c(x, "1 # num of TAC model categories")
  x <- c(x, "1 # TAC model indices")
  x <- c(x, paste(ayr,catch,"# Catch in each year starting w/ begining NAA"))
  writeLines(x, con=file.path(path, 'spm.dat'))
}


#' Write input files for the Projection module
#'
#' @param replist A list as read in by \link{read_rep}
#' @param datlist A list as read in by \link{read_dat}
#' @param path Directory to write in
#'
#' @return Nothing but writes three files
#' @export
write_proj_inputs <- function(replist, datlist, path=getwd()){
  ## delete old files in case it fails it'll stop there
  ff <- file.path(path, c("goa_wp.txt", "setup.dat",
                          "spp_catch.dat"))
  lapply(ff, function(x) if(file.exists(x)) file.remove(x))
  ayr <- tail(replist$years,1)
  write_proj_setup(ayr, path=path)
  write_proj_catch(replist, ayr, path=path)
  write_proj_dynamics(replist, datlist, ayr, path)
}
write_proj_setup <- function(ayr, path){
  x <- c()
  x <- c(x,paste("# proj input file written by write_proj_inputs on",
                 Sys.time()))
  x <- c(x, " std # run name")
  x <- c(x, "7 # number of alt scenarios")
  x <- c(x, paste(1:7, collapse=' '))
  x <- c(x, "1 # flag to set TAC=ABC (1 is true)")
  x <- c(x, "2 # SR type (1 ricker, 2 bholt")
  x <- c(x, "1 # proj rec form (default: 1= use obs mean/SD")
  x <- c(x, "0 # SR conditioning (0 means no)")
  x <- c(x, "0.0 # turn off prior thing")
  x <- c(x, "1 # flag to write bigfile")
  x <- c(x, "14 # number of proj years")
  x <- c(x, "1000 # number of simulations")
  x <- c(x, paste(ayr, " # begin year"))
  writeLines(x, con=file.path(path, 'setup.dat'))
}

write_proj_catch <- function(replist, ayr, path){
  x <- c()
  x <- c(x,paste("# proj input file written by write_proj_inputs on",
                 Sys.time()))
  x <- c(x, "1 # num years of catch")
  x <- c(x, "1 # num species")
  x <- c(x, "1343.248 # unused")
  x <- c(x, "1943.248 # unused")
  x <- c(x, "goa_wp.txt # input file name")
  x <- c(x, "1 # ABC multipliers")
  x <- c(x, "1 # pop scalars")
  x <- c(x, "0.75 # new alt 4 Fabc SPRs (unused?)")
  x <- c(x, "1 # num of TAC model categories")
  x <- c(x, "1 # TAC model indices")
  x <- c(x, "# Catch in each year starting w/ begining NAA")
  ## Not the last year b/c this is data and if using retro option
  ## it's full length
  ind <- which(replist$years == ayr)
  x <- c(x, paste(ayr, replist$Total_catch[ind]))
  writeLines(x, con=file.path(path, 'spp_catch.dat'))
}


write_proj_dynamics <- function(replist, datlist, ayr, path){
  x <- c()
  x <-
  c(x,paste("# proj input file written by write_proj_inputs on", Sys.time()))
  x <- c(x, "goa_wp")
  x <- c(x, "1 # Flag to tell if this is a SSL forage species")
  x <- c(x, "0 # Flag to Dorn's version of a constant buffer")
  x <- c(x, "1 # number of fisheries")
  x <- c(x, "1 # number of sexes")
  x <- c(x, paste(mean(tail(replist$Fishing_mortalities,5)), " # average F over last 5 years"))
  x <- c(x, "1 # Author's F multiplier (to MaxPermissible)")
  x <- c(x, "0.4 # ABC SPR")
  x <- c(x, "0.35 # OFL SPR")
  x <- c(x, "3.52 # Month of spawning")
  x <- c(x, "10 # Number of ages")
  x <- c(x, "1 # Fratio")
  x <- c(x, "# Natural mortality")
  x <- c(x, paste(replist$Natural_mortality, collapse=' '))
  x <- c(x, "# Maturity")
  x <- c(x, paste(datlist$mat, collapse=' '))
  x <- c(x, "# Spawning WAA")
  ## x <- c(x, paste(datlist$wt_spawn_proj*1000, collapse=' '))
  x <- c(x, paste(replist$Projection_spawning_weight_at_age*1000, collapse=' '))
  x <- c(x, "# Fishery WAA")
  ##  x <- c(x, paste(datlist$wt_fsh_proj*1000, collapse=' '))
  x <- c(x, paste(replist$Projection_fishery_weight_at_age*1000, collapse=' '))
  x <- c(x, "# Fishery selex, averaged over last 5 years")
  ## This is wrong b/c ignores most recent year
  ## x <- c(x, paste(colMeans(tail(replist$Fishery_selectivity,
  ## 5)), collapse= ' '))
 ##  x <- c(x, paste(replist$Mean_fishery_selectivity, collapse='
  ##  '))
  x <- c(x, paste(replist[[111]], collapse=' '))
  x <- c(x, "# current year starting NAA")
  x <- c(x, paste(tail(replist$Numbers_at_age,1)*1000, collapse=' '))
  ind <- which(replist$years %in% 1978:(ayr-1))
  recs <- replist$Recruits[ind]*1000
  ssb <- replist$Expected_spawning_biomass[ind-1]*1000
  x <- c(x, paste(length(recs), " # num of recruits"))
  x <- c(x, "# Recruits from 1978 to this year -1 due to uncertain last year")
  x <- c(x, paste(recs, collapse=' '))
  x <- c(x, '# SSB from 1977 to this year-2 to match recruits')
  x <- c(x, paste(ssb, collapse=' '))
  writeLines(x, con=file.path(path, 'goa_wp.txt'))
}
