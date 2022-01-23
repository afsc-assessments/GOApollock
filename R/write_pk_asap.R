# function to write goa pollock dat and rep file to asap3 file for reading into wham

write_pk_asap <- function(dat, 
                          rep, 
                          surveys=c(1,2,3,6), 
                          path=NULL,
                          filename=NULL) {
  
  if(!is.null(path)){
    oldwd <- getwd()
    on.exit(setwd(oldwd))
    setwd(path)
  }
  
  if(is.null(filename)){ filename <- 'goa_pk_asap3.txt' }
  
  # Model dimensions ----
  n_years_model <- dat$endyr - dat$styr + 1
  n_ages <- dat$trmage - dat$rcrage + 1
  n_fleets <- 1
  n_indices <- length(surveys)
  
  cat("# Alaska GOA pollock ", dat$endyr, " assessment", "\n",
      "# Contact: cole.monnahan@noaa.gov", "\n",
      "# Number of Years", "\n",
      n_years_model, "\n",
      "# First year", "\n",
      dat$styr, "\n",
      "# Number of ages", "\n",
      n_ages, "\n",
      "# Number of fleets", "\n",
      n_fleets, "\n",
      "# Number of selectivity blocks", "\n",
      1, "\n",
      "# Number of available indices", "\n",
      n_indices, "\n",
      "# M matrix", "\n", 
      sep = "", file = filename)
  
  # Natural mortality ----
  write.table(matrix(data = rep$Natural_mortality, ncol = n_ages, nrow = n_years_model, byrow = TRUE), 
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("# Fecundity option", "\n",
      0, "\n", 
      "# Fraction of year elapsed before SSB calculation", "\n",
      0.21, "\n",  # from goa_pk.tpl line 760
      "# MATURITY matrix", "\n",
      sep = "", file = filename, append = TRUE)
  
  # Maturity ----
  write.table(matrix(data = dat$mat, ncol = n_ages, nrow = n_years_model, byrow = TRUE), 
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Weight-at-age ----
  
  cat("# Number of WAA matrices", "\n",
      4, "\n", 
      "# WAA matrix-1 fishery", "\n",
      sep = "", file = filename, append = TRUE)
  
  write.table(dat$wt_fsh, file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # round(dat$wt_srv1,3)==dat$wt_spawn
  cat("# WAA matrix-2 survey 1 shelikof also use for ssb", "\n",
      sep = " ", file = filename, append = TRUE)
  
  write.table(dat$wt_srv1, file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("# WAA matrix-3 survey 2 NMFS BTS also use for total pop biom", "\n",
      sep = " ", file = filename, append = TRUE)
  
  write.table(dat$wt_srv2, file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("# WAA matrix-4 survey 6 Summer Accoustic Trawl", "\n",
      sep = " ", file = filename, append = TRUE)
  
  write.table(dat$wt_srv6, file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("# WEIGHT AT AGE POINTERS", "\n",
      1, " # fleet 1 catch", "\n",
      1, " # fleet 1 discards (not used)", "\n", 
      1, " # total catch", "\n", 
      1, " # total discards (not used)", "\n", 
      2, " # SSB - use Shelikof", "\n", 
      3, " # Jan 1 - total pop use NMFS BTS", "\n", 
      "# Selectivity blocks (blocks within years)", "\n", 
      "# Fleet 1 Selectivity Block Assignment", "\n", 
      sep = "", file = filename, append = TRUE)
  
  write.table(matrix(data = 1, ncol = 1, nrow = n_years_model), file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("# Selectivity options for each block", "\n",
      3, " # fishery uses double logistic", "\n",
      "# Selectivity Block #1 Data", "\n", 
      sep = "", file = filename, append = TRUE)
  
  # Selectivity in ASAP dat file is structured as follows: four columns (Initial guess,
  # Phase, Lambda, and Coefficient of Variation). Number of rows = number of ages
  # + 2 parameters for logistic + 4 parameters for double logistic
  
  # Fishery selectivity averaged by age over all years ----
  
  # names(rep)[which(grepl('Sel|sel', names(rep)))]
  fsh_sel <- as.vector(colMeans(rep$Fishery_selectivity))
  fsh_sel[fsh_sel == 0] <- 0.0001 # can't be 1 or 0 otherwise it doesn't play nice with the logit trans
  fsh_sel[fsh_sel == 1] <- 0.9999
  
  fsh_sel_inits <- c(
    # age-specific
    fsh_sel,
    # double logistic
    c(rep$Selectivity_means[2], 1/exp(rep$Selectivity_means[1]),
      rep$Selectivity_means[4], 1/exp(rep$Selectivity_means[3])),
    # logistic
    c(rep$Selectivity_means[2], 1/exp(rep$Selectivity_means[1]))
  )
  
  # age-7 fixed at 1 in ADMB
  fsh_sel_phase <- c(1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  fsh_sel_lambda <- rep(0, length(fsh_sel_inits))
  fsh_sel_cv <- rep(1, length(fsh_sel_inits))
  
  write.table(matrix(data = c(fsh_sel_inits, fsh_sel_phase, 
                              fsh_sel_lambda, fsh_sel_cv), 
                     ncol = 4), 
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("# Selectivity start age by fleet", "\n",
      1,  "\n",
      "# Selectivity end age by by fleet", "\n", 
      10, "\n",
      "# Age range for average F", "\n",
      1, " ", 10, "\n",
      "# Average F report option", "\n", 
      1,  "\n",
      "# Use likelihood constants", "\n", 
      0,  "\n",
      "# Release mortality by fleet", "\n", 
      0,  "\n",
      "# Catch data", "\n", 
      "# Fleet-1 Catch data", "\n", 
      sep = "", file = filename, append = TRUE)
  
  # Fishery age comps + catch ----
  
  # names(dat)[which(grepl('cat|fsh', names(dat)))]
  # dat$cattot
  # dat$catp
  # length(dat$fshyrs)
  
  catp <- cbind(dat$fshyrs,dat$catp)
  colnames(catp) <- c('year', 1:n_ages)
  catmat <- merge(catp,
                  data.frame(year = dat$styr:dat$endyr,
                             catch = dat$cattot),
                  all = TRUE)
  catmat[is.na(catmat)] <- -999
  
  # Catch data, first comps then catch (dim: nrow = nyrs, ncols = nage+1); no year column
  write.table(unname(catmat[,-1]), 
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("# Discards at age by fleet", "\n", 
      "# Fleet-1 Discards data", "\n", 
      sep = "", file = filename, append = TRUE)
  
  # all 0's for discard data 
  write.table(matrix(data = 0, nrow = n_years_model, ncol = n_ages + 1), 
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("# Release proportion at age by fleet", "\n", 
      "# Fleet-1 Release data", "\n", 
      sep = "", file = filename, append = TRUE)
  
  # all 0's for release data 
  write.table(matrix(data = 0, nrow = n_years_model, ncol = n_ages + 1), 
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Survey units/WAA pointers ----
  cat("# Survey Index Data", "\n", 
      "# Index units", "\n", 
      sep = "", file = filename, append = TRUE)
  
  # units 1=biomass, 2=abundance
  cat(rep(1, n_indices), "\n", sep = " ", file = filename, append = TRUE)
  
  cat("# Index Age comp units", "\n", 
      sep = "", file = filename, append = TRUE)
  
  # units 1=biomass, 2=abundance (all numbers-at-age)
  cat(rep(2, n_indices), "\n", sep = " ", file = filename, append = TRUE)
  
  cat("# Index WAA matrix", "\n", 
      sep = "", file = filename, append = TRUE)
  
  waa_pointer_indices <- c(2,3,3,2,2,4)
  waa_pointer_indices <- waa_pointer_indices[surveys]
  cat(waa_pointer_indices, "\n", sep = " ", file = filename, append = TRUE)
  
  cat("# Index month", "\n", 
      sep = "", file = filename, append = TRUE)
  
  names(dat)[which(grepl('frct', names(dat)))]
  fracyr_indices <- 12 * c(dat$yrfrct_srv1[1], dat$yrfrct_srv2[1], dat$yrfrct_srv3[1], 
                           dat$yrfrct_srv1[1], dat$yrfrct_srv1[1], dat$yrfrct_srv6[1]) 
  fracyr_indices <- fracyr_indices[surveys]
  cat(fracyr_indices, "\n", sep = " ", file = filename, append = TRUE)
  
  # What is this? Just use a vector of -1 with length = n_indices
  cat("# Index link to fleet", "\n", 
      sep = "", file = filename, append = TRUE)
  cat(rep(-1, n_indices), "\n", sep = " ", file = filename, append = TRUE)
  
  # Survey selectivities ----
  cat("# Index selectivity option", "\n", 
      sep = "", file = filename, append = TRUE)
  index_sel_option <- c(
    # Survey 1 Shelikof accoustic: descending logistic with age-1 and age-2 fixed at
    # 0, age-3 fixed at 1; option 4 in WHAM. set to option 1 (age-specific) for now, 
    1,
    # Survey 2 NMFS BTS: logistic/asymptotic 
    2,
    # Survey 3 ADFG coastal survey: logistic/asymptotic
    2, 
    # Survey 4 Shelikof accoustic age-1: fixed at 1 for age-1, 0 for all other ages
    # - use age-specific option in wham, fix sel pars
    1,
    # Survey 5 Shelikof accoustic age-2: fixed at 1 for age-2, 0 for all other ages
    # - use age-specific option in wham, fix sel pars
    1,
    # Survey 6 Summer accoustic trawl (AT): fixed at 1 for all ages - use
    # age-specific option in wham, fix sel pars
    1)
  
  index_sel_option <- index_sel_option[surveys]
  cat(index_sel_option, "\n", sep = " ", file = filename, append = TRUE)
  
  cat("# Index start age", "\n", 
      sep = "", file = filename, append = TRUE)
  cat(rep(1, n_indices), "\n", sep = " ", file = filename, append = TRUE)
  
  cat("# Index end age", "\n", 
      sep = "", file = filename, append = TRUE)
  cat(rep(10, n_indices), "\n", sep = " ", file = filename, append = TRUE)
  
  # FLAG - is this how to estimate an index without any composition data?
  cat("# Index Estimate Proportion (YES=1)", "\n", 
      sep = "", file = filename, append = TRUE)
  estimate_index_proportion <- c(1, 1, 1, 0, 0, 1) # age-1 and age-2 shelikof accoustic surveys have no age comps
  estimate_index_proportion <- estimate_index_proportion[surveys]
  cat(estimate_index_proportion, "\n", sep = " ", file = filename, append = TRUE)
  
  # FLAG to turn off index fitting(?)
  cat("# Use Index", "\n", 
      sep = "", file = filename, append = TRUE)
  cat(rep(1, n_indices), "\n", sep = " ", file = filename, append = TRUE)
  
  # Survey 1 selectivity ----
  
  if(1 %in% surveys){
    
    cat("# Index-1 Selectivity Data Shelikof age-3 plus", "\n", 
        sep = "", file = filename, append = TRUE)
    
    names(rep)[which(grepl('Sel|sel', names(rep)))]
    
    # Survey 1 Shelikof accoustic: descending logistic with age-1 and age-2 fixed at
    # 0, age-3 fixed at 1; option 4 in WHAM
    srv1_sel_init <- rep$Survey_1_selectivity
    srv1_sel_init[srv1_sel_init == 0] <- 0.0001
    srv1_sel_init[srv1_sel_init == 1] <- 0.9999
    
    # placeholder until read_pk_rep() returns survey-specific parameter estimates
    srv1_k <- 0.61416 # estimate from admb
    srv1_a50 <- 9.60821
    
    # survey 1 selectivity matrix - desc logistic
    write.table(matrix(data = 
                         c(# initial guess/fixed values 
                           srv1_sel_init, # age-specific
                           srv1_a50, 1/exp(srv1_k), # logistic
                           3.5, 1.5, srv1_a50, 1/exp(srv1_k), #double logistic
                           # phases
                           -1, -1, -1, 1, 1, 1, 1 ,1, 1, 1, # age-specific
                           1, 1, # logistic
                           1, 1, 1, 1, #double logistic
                           rep(0, n_ages+6), # lambdas
                           rep(1, n_ages+6) # cv
                         ),
                       nrow = n_ages+6, ncol = 4),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Survey 2 selectivity ----
  
  if(2 %in% surveys){
    
    cat("# Index-2 Selectivity Data", "\n", 
        sep = "", file = filename, append = TRUE)
    # Survey 2 NMFS BTS: logistic/asymptotic
    srv2_sel_init <- rep$Survey_2_selectivity
    srv2_sel_init[srv2_sel_init == 0] <- 0.0001
    srv2_sel_init[srv2_sel_init == 1] <- 0.9999
    
    # placeholder until read_pk_rep() returns survey-specific parameter estimates
    srv2_k <- -0.304719 # estimate from admb
    srv2_a50 <- 3.58329
    
    write.table(matrix(data = 
                         c(# initial guess/fixed values (from par/rep file)
                           srv2_sel_init, # age-specific
                           srv2_a50, 1/exp(srv2_k), # logistic
                           srv2_a50, 1/exp(srv2_k), n_ages-0.0001, 1/exp(1), #double logistic
                           # phases
                           1, 1, 1, 1, 1, 1, 1, 1, 1, -1, # age-specific
                           1, 1, # logistic
                           1, 1, 1, 1, #double logistic
                           rep(0, n_ages+6), # lambdas
                           rep(1, n_ages+6) # cv
                         ),
                       nrow = n_ages+6, ncol = 4),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Survey 3 selectivity ----
  
  if(3 %in% surveys){
    
    cat("# Index-3 Selectivity Data", "\n", 
        sep = "", file = filename, append = TRUE)
    
    # Survey 3 ADFG coastal survey: logistic/asymptotic
    srv3_sel_init <- rep$Survey_3_selectivity
    srv3_sel_init[srv3_sel_init == 0] <- 0.0001
    srv3_sel_init[srv3_sel_init == 1] <- 0.9999
    
    # placeholder until read_pk_rep() returns survey-specific parameter estimates
    srv3_k <- 0.298298 # estimate from admb
    srv3_a50 <- 4.8761
    
    write.table(matrix(data = 
                         c(# initial guess/fixed values (from par/rep file)
                           srv3_sel_init, # age-specific
                           srv3_a50, 1/exp(srv3_k), # logistic
                           srv3_a50, 1/exp(srv3_k), n_ages-0.0001, 1/exp(1), #double logistic
                           # phases
                           1, 1, 1, 1, 1, 1, 1, 1, 1, -1, # age-specific
                           1, 1, # logistic
                           1, 1, 1, 1, #double logistic
                           rep(0, n_ages+6), # lambdas
                           rep(1, n_ages+6) # cv
                         ),
                       nrow = n_ages+6, ncol = 4),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  # Survey 4 selectivity ----
  if(4 %in% surveys){
    
    cat("# Index-4 Selectivity Data", "\n", 
        sep = "", file = filename, append = TRUE)
    
    # Survey 4 Shelikof accoustic age-1: fixed at 1 for age-1, 0 for all other ages
    # - use age-specific option in wham, fix sel pars
    srv4_sel_init <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    srv4_sel_init[srv4_sel_init == 0] <- 0.0001
    srv4_sel_init[srv4_sel_init == 1] <- 0.9999
    
    write.table(matrix(data =
                         c(# initial guess/fixed values
                           srv4_sel_init, # age-specific
                           4.8761, 0.742080167, # logistic
                           4.8761, 0.742080167, n_ages-0.0001, 1/exp(1), #double logistic
                           # phases
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, # age-specific
                           1, 1, # logistic
                           1, 1, 1, 1, #double logistic
                           rep(0, n_ages+6), # lambdas
                           rep(1, n_ages+6) # cv
                         ),
                       nrow = n_ages+6, ncol = 4),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Survey 5 selectivity -----
  
  if(5 %in% surveys){
    
    cat("# Index-5 Selectivity Data", "\n",
        sep = "", file = filename, append = TRUE)
    
    # Survey 5 Shelikof accoustic age-2: fixed at 1 for age-2, 0 for all other ages
    # - use age-specific option in wham
    srv5_sel_init <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
    srv5_sel_init[srv5_sel_init == 0] <- 0.0001
    srv5_sel_init[srv5_sel_init == 1] <- 0.9999
    
    write.table(matrix(data =
                         c(# initial guess/fixed values
                           srv5_sel_init, # age-specific
                           4.8761, 0.742080167, # logistic
                           4.8761, 0.742080167, n_ages-0.0001, 1/exp(1), #double logistic
                           # phases
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, # age-specific
                           1, 1, # logistic
                           1, 1, 1, 1, #double logistic
                           rep(0, n_ages+6), # lambdas
                           rep(1, n_ages+6) # cv
                         ),
                       nrow = n_ages+6, ncol = 4),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Survey 6 selectivity -----
  if(6 %in% surveys){
    
    cat("# Index-6 Selectivity Data", "\n", 
        sep = "", file = filename, append = TRUE)
    
    # Survey 6 Summer accoustic trawl (AT): fixed at 1 for all ages - use
    # age-specific option in wham
    srv6_sel_init <- rep$Survey_6_selectivity
    srv6_sel_init[srv6_sel_init == 0] <- 0.0001
    srv6_sel_init[srv6_sel_init == 1] <- 0.9999
    
    write.table(matrix(data = 
                         c(# initial guess/fixed values 
                           srv6_sel_init, # age-specific
                           4.8761, 0.742080167, # logistic
                           4.8761, 0.742080167, n_ages-0.0001, 1/exp(1), #double logistic
                           # phases
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, # age-specific
                           1, 1, # logistic
                           1, 1, 1, 1, #double logistic
                           rep(0, n_ages+6), # lambdas
                           rep(1, n_ages+6) # cv
                         ),
                       nrow = n_ages+6, ncol = 4),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Survey data ----
  
  # nrow = n_years_model, ncol = n_ages+4
  cat("# Index data matrices (year, index, index_sigma, age comps, neff)", "\n", 
      sep = "", file = filename, append = TRUE)
  
  # Survey 1 data ----
  names(dat)[which(grepl('srv1|surv|yr', names(dat)))]
  names(dat)[which(grepl('multN', names(dat)))]
  
  if(1 %in% surveys){
    
    cat("# Index 1 Shelikof age-3 plus 1992-pres.", "\n", 
        sep = "", file = filename, append = TRUE)
    
    srv1_indx <- merge(expand.grid(data.frame(year = dat$styr:dat$endyr)),
                       data.frame(year = dat$srvyrs1,
                                  index = dat$indxsurv1 * 1e6,
                                  cv = dat$indxsurv_log_sd1),
                       all.x = TRUE)
    
    srv1_ac <- cbind(cbind(data.frame(year = dat$srv_acyrs1),
                           dat$srvp1),
                     data.frame(Neff = dat$multN_srv1))
    
    srv1_dat <- merge(srv1_indx, srv1_ac, all.x = TRUE)
    srv1_dat[is.na(srv1_dat)] <- -999
    
    write.table(unname(as.matrix(srv1_dat)),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Survey 2 data ----
  
  if(2 %in% surveys){
    
    cat("# Index 2 NMFS BTS", "\n", 
        sep = "", file = filename, append = TRUE)
    
    srv2_indx <- merge(expand.grid(data.frame(year = dat$styr:dat$endyr)),
                       data.frame(year = dat$srvyrs2,
                                  index = dat$indxsurv2 * 1e6,
                                  cv = dat$indxsurv_log_sd2),
                       all.x = TRUE)
    
    srv2_ac <- cbind(cbind(data.frame(year = dat$srv_acyrs2),
                           dat$srvp2),
                     data.frame(Neff = dat$multN_srv2))
    
    srv2_dat <- merge(srv2_indx, srv2_ac, all.x = TRUE)
    srv2_dat[is.na(srv2_dat)] <- -999
    # cbind(asap3$dat$IAA_mats[[2]][,c(1,2)] , srv2_dat[, c(1,2)])
    
    write.table(unname(as.matrix(srv2_dat)),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Survey 3 data ----
  
  if(3 %in% surveys){
    cat("# Index 3 ADFG coastal survey", "\n", 
        sep = "", file = filename, append = TRUE)
    
    srv3_indx <- merge(expand.grid(data.frame(year = dat$styr:dat$endyr)),
                       data.frame(year = dat$srvyrs3,
                                  index = dat$indxsurv3 * 1e6,
                                  cv = dat$indxsurv_log_sd3),
                       all.x = TRUE)
    
    srv3_ac <- cbind(cbind(data.frame(year = dat$srv_acyrs3),
                           dat$srvp3),
                     data.frame(Neff = dat$multN_srv3))
    
    srv3_dat <- merge(srv3_indx, srv3_ac, all.x = TRUE)
    srv3_dat[is.na(srv3_dat)] <- -999
    # cbind(asap3$dat$IAA_mats[[3]][,c(1,2)] , srv3_dat[, c(1,2)])
    
    write.table(unname(as.matrix(srv3_dat)),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Survey 4 data ----
  
  if(4 %in% surveys){
    cat("# Index 4 Shelikof age-1", "\n",
        sep = "", file = filename, append = TRUE)
    
    srv4_indx <- merge(expand.grid(data.frame(year = dat$styr:dat$endyr)),
                       data.frame(year = dat$srvyrs4,
                                  index = dat$indxsurv4 * 1e6,
                                  cv = dat$indxsurv_log_sd4),
                       all.x = TRUE)
    
    srv4_ac <- matrix(data = -999, ncol = n_ages + 1, nrow = n_years_model)
    srv4_dat <- cbind(srv4_indx, srv4_ac)
    srv4_dat[is.na(srv4_dat)] <- -999
    # dim(srv4_dat) == c(n_years_model, n_ages + 4)
    
    write.table(unname(as.matrix(srv4_dat)),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Survey 5 data ----
  
  if(5 %in% surveys){
    cat("# Index 5 Shelikof age-2", "\n",
        sep = "", file = filename, append = TRUE)
    
    srv5_indx <- merge(expand.grid(data.frame(year = dat$styr:dat$endyr)),
                       data.frame(year = dat$srvyrs5,
                                  index = dat$indxsurv5 * 1e6,
                                  cv = dat$indxsurv_log_sd5),
                       all.x = TRUE)
    
    srv5_ac <- matrix(data = -999, ncol = n_ages + 1, nrow = n_years_model)
    srv5_dat <- cbind(srv5_indx, srv5_ac)
    srv5_dat[is.na(srv5_dat)] <- -999
    # dim(srv5_dat) == c(n_years_model, n_ages + 4)
    
    write.table(unname(as.matrix(srv5_dat)),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Survey 6 data ----
  
  if(6 %in% surveys){
    cat("# Index 6 Summer accoustic trawl (AT)", "\n", 
        sep = "", file = filename, append = TRUE)
    
    srv6_indx <- merge(expand.grid(data.frame(year = dat$styr:dat$endyr)),
                       data.frame(year = dat$srvyrs6,
                                  index = dat$indxsurv6 * 1e6,
                                  cv = dat$indxsurv_log_sd6),
                       all.x = TRUE)
    
    srv6_ac <- cbind(cbind(data.frame(year = dat$srv_acyrs6),
                           dat$srvp6),
                     data.frame(Neff = dat$multN_srv6))
    
    srv6_dat <- merge(srv6_indx, srv6_ac, all.x = TRUE)
    srv6_dat[is.na(srv6_dat)] <- -999
    # cbind(asap3$dat$IAA_mats[[4]][,c(1,2)] , srv6_dat[, c(1,2)])
    
    write.table(unname(as.matrix(srv6_dat)),
                file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Phases ----
  cat("#########################################", "\n",
      "# Phase data", "\n",
      "# Phase for Fmult in 1st year", "\n",
      1, "\n",
      "# Phase for Fmult deviations", "\n",
      3, "\n",
      "# Phase for recruitment deviations", "\n", 
      3, "\n",
      "# Phase for N in 1st year ", "\n",
      2, "\n",
      "# Phase for catchability in 1st year", "\n",
      1, "\n",
      "# Phase for catchability deviations", "\n",
      -1, "\n",
      "# Phase for stock recruit relationship", "\n",
      1, "\n",
      "# Phase for steepness", "\n",
      -2, "\n",
      "#########################################", "\n",
      "# Lambdas and CVs","\n",
      "# Recruitment CV by year", "\n",
      sep = "", file = filename, append = TRUE)
  
  write.table(matrix(data = 1, nrow = n_years_model, ncol = 1),
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Lambdas and CVs ----
  cat("# Lambda for each index", "\n", sep = "", file = filename, append = TRUE)
  cat(rep(1, n_indices), "\n", sep = " ", file = filename, append = TRUE)
  
  cat("# Lambda for Total catch in weight by fleet", "\n", sep = "", file = filename, append = TRUE)
  cat(rep(1, n_fleets), "\n", sep = " ", file = filename, append = TRUE)
  
  cat("# Lambda for total discards at age by fleet", "\n", sep = "", file = filename, append = TRUE)
  cat(rep(0, n_fleets), "\n", sep = " ", file = filename, append = TRUE)
  
  cat("# Catch Total CV by year and fleet", "\n", sep = "", file = filename, append = TRUE)
  write.table(dat$cattot_log_sd, 
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("# Discard Total CV by year and fleet", "\n", sep = "", file = filename, append = TRUE)
  write.table(matrix(data = 0, ncol = n_fleets, nrow = n_years_model),
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Neff for catch comps ----
  cat_Neff <- cbind(dat$fshyrs, dat$multN_fsh)
  colnames(cat_Neff) <- c('year', 'Neff')
  cat_Neff <- merge(cat_Neff,
                    data.frame(year = dat$styr:dat$endyr),
                    all.y = TRUE)
  cat_Neff[is.na(cat_Neff)] <- -999
  
  cat("# Input effective sample size for catch at age by year and fleet", "\n", sep = "", file = filename, append = TRUE)
  write.table(unname(cat_Neff[-1]),
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("# Input effective sample size for discards at age by year and fleet", "\n", sep = "", file = filename, append = TRUE)
  write.table(matrix(data = 0, ncol = n_fleets, nrow = n_years_model),
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # More lambdas and CVs ----
  cat("# Lambda for Fmult in first year by fleet", "\n",
      0, "\n",
      "# CV for Fmult in first year by fleet", "\n",
      1, "\n",
      "# Lambda for Fmult deviations", "\n",
      0, "\n",
      "# CV for Fmult deviations", "\n",
      1, "\n",
      "# Lambda for N in 1st year deviations", "\n", 
      0, "\n",
      "# CV for N in 1st year deviations", "\n",
      1, "\n",
      "# Lambda for recruitment deviations", "\n",
      1, "\n",
      "# Lambda for catchability in first year by index", "\n",
      sep = "", 
      file = filename, append = TRUE)
  
  # Lambda for catchability in first year by index - could this be used to fix q?
  cat(rep(0, n_indices), "\n",file = filename, sep = " ", append = TRUE)
  
  # FLAG CV for catchability - could this be used to fix the odd ball survey qs?
  cat("# CV for catchability in first year by index", "\n",
      file = filename, sep = "", append = TRUE)
  cat(rep(1, n_indices), "\n", file = filename, sep = " ", append = TRUE)
  
  cat("# Lambda for catchability deviations by index", "\n",
      file = filename, sep = "", append = TRUE)
  cat(rep(1, n_indices), "\n", file = filename, sep = " ", append = TRUE)
  
  cat("# CV for catchability deviations by index", "\n",
      file = filename, sep = "", append = TRUE)
  cat(rep(1, n_indices), "\n", file = filename, sep = " ", append = TRUE)
  
  cat("# Lambda for deviation from initial steepness", "\n",
      0, "\n",
      "# CV for deviation from initial steepness", "\n",
      1, "\n",
      "# Lambda for deviation from initial SSB0", "\n", 
      0, "\n",
      "# CV for deviation from initial SSB0", "\n", 
      1, "\n",
      "# NAA Deviations flag (1=   , 0=  )", "\n", 
      1, "\n",
      "###########################################", "\n",
      "###  Initial Guesses", "\n",
      "# NAA for year1", "\n",
      sep = "", file = filename, append = TRUE)
  
  # Inits for scaling pars ----
  # names(rep)[which(grepl('Numbers', names(rep)))]
  init_naa <- unname(rep$Numbers_at_age[1,] * 1e6)
  cat(init_naa, "\n", file = filename, sep = " ", append = TRUE)
  
  # names(rep)[which(grepl('Fish', names(rep)))]
  init_Fmort <- rep$Fishing_mortalities[1]
  cat("# Fmult in 1st year by fleet", "\n", 
      init_Fmort, "\n", 
      "# Catchability in 1st year by index", "\n",
      file = filename, sep = "", append = TRUE)
  
  # names(rep)[which(grepl('q', names(rep)))]
  q_inits <- c(mean(rep$Survey_1_q[-1]),     # FLAG - what is q_bs?
               mean(rep$Survey_2_q), 
               mean(rep$Survey_3_q), 
               rep$Survey_4_q, 
               rep$Survey_5_q, 
               rep$Survey_6_q)
  q_inits <- q_inits[surveys]
  cat(q_inits,
      "\n", file = filename, sep = " ", append = TRUE)
  
  # S-R, projections, MCMC ----
  cat("# S-R Unexploited specification (1=   0=)", "\n",
      1, "\n",
      "# Unexploited initial guess", "\n",
      1e+07, "\n",
      "# Steepness initial guess", "\n",
      1, "\n",
      "# Maximum F (upper bound on Fmult)", "\n",
      5, "\n",
      "# Ignore guesses", "\n",
      0, "\n",
      "###########################################", "\n",
      "###  Projection Control data", "\n",
      "# Do projections", "\n",
      0, "\n",
      "# Fleet directed flag", "\n",
      1, "\n",
      "# Final year of projections", "\n",
      dat$endyr+2, "\n",
      "# Year, projected recruits, what projected, target, non-directed Fmult", "\n",
      sep = "", file = filename, append = TRUE)
  
  # I don't know what these inputs are exactly or they're even used in wham. not in asap documentation.
  write.table(matrix(data = c(dat$endyr+1, dat$endyr+2, -1, -1, 1, 3, 150, -99, 0, 0),
                     ncol = 5, nrow = 2),
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("###########################################", "\n",
      "###  MCMC Control data", "\n",
      "# do mcmc", "\n",
      0, "\n",
      "# MCMC nyear option", "\n",
      0, "\n",
      "# MCMC number of saved iterations desired", "\n",
      1000, "\n",
      "# MCMC thinning rate", "\n",
      200, "\n",
      "# MCMC random number seed", "\n",
      5230547, "\n",
      "###########################################", "\n",
      "###  A few AGEPRO specs", "\n",
      "# R in agepro.bsn file", "\n",
      0, "\n",
      "# Starting year for calculation of R", "\n",
      2004, "\n",
      "# Ending year for calculation of R", "\n",
      2014, "\n",
      "# Export to R flag (1=  0=)", "\n",
      1, "\n",
      "# test value", "\n",
      -23456, "\n",
      "###########################################", "\n",
      "###### FINIS ######", "\n",
      "# Fleet Names", "\n",
      "Trawl_fleet", "\n",
      "# Survey Names", "\n",
      sep = "", file = filename, append = TRUE)
  
  survey_names <- c("Shelikof_age3plus", 
                    "NMFS_BTS", 
                    "ADFG_Survey", 
                    "Shelikof_age1", 
                    "Shelikof_age2", 
                    "Summer_AT")
  survey_names <- survey_names[surveys]
  write.table(matrix(data = survey_names,
                     ncol = 1, nrow = length(surveys)),
              file = filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
}


