
#' Write input files for the Projection module
#'
#' @param replist A list as read in by \link{read_pk_rep}
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
  x <- c(x, paste(datlist$wt_spawn_proj*1000, collapse=' '))
  x <- c(x, "# Fishery WAA")
  x <- c(x, paste(datlist$wt_fsh_proj*1000, collapse=' '))
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
