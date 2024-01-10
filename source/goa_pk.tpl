// GOApollock age-structured model originally developed by
// Martin Dorn (NFMS AFSC), taken over by Cole Monnahan (NMFS,
// AFSC) in 2021 (pk20_8 renamed to goa_pk).

// One fishery, three surveys: acoustic, bottom trawl, ADFG crab/groundfish  
// Double logistic selectivity for fishery
// Random walks in selectivity
// Logistic or double logistic for surveys
// Five year projection
// Shelikof Strait EIT split into three catchability periods
// Biosonics MillerFreeman and OscarDyson 
// Implements 2 corrections to coding errors noticed by Teresa in Spring 2014
// Namely: Uses 20 iterations rather than 10 to get the correct spawning biomass in the HCR, AND implements the bias corrected log likelihood for survey biomass correctly.

// Command line option '-retro N' added to allow for easy
// retrospective analyses. In this mode RMSE for indices are not
// calculated (fixed at zero).

// Added log biomass sdreport vectors. Will use those for
// uncertainties moving forward.

// Cleaned out old inputs  9/22.

// Model 19.1a: add sigmaR=1.3 for all devs, turn on descending
// logistic selectivity for survey 6, including some broad priors
// to stabilize estimation. Also add more sdreport variables to
// output.

DATA_SECTION
  // Command line argument to do a retrospective peel
  // To use this flag, run the model using: "-retro n" to peel n years.
  //
  // Modified in 2021 by code Steve Martell wrote for ATF. The
  // general idea is to read in the data with dimensions as is,
  // then redefine endyr in the parameter section to be shorter
  // by n years. To avoid over-running arrays it breaks out of
  // for loops early and has conditionals based on isretro value
  int retro_yrs 
  int isretro
  !! int on,opt;
  !! retro_yrs=0;
  !! isretro=0;
  !! if((on=option_match(ad_comm::argc,ad_comm::argv,"-retro",opt))>-1){
  !!   retro_yrs=atoi(ad_comm::argv[on+1]);
  !!   isretro=1;
  !!   cout << "\n!!  Implementing retro w/ peels="<< retro_yrs << " !!\n";
  !! }
  
  !!CLASS ofstream report1("mceval.dat")

  init_int styr                                  // Starting year for population model
  init_int endyr                                 // Ending year for population model
  init_int rcrage                                // Recruitment age
  init_int trmage                                // Last modeled age
  init_int nbins1                                // Number of length bins in transitiom matrix 1
  init_int nbins2                                // Number of length bins in transitiom matrix 2
  init_int nbins3                                // Number of length bins in transitiom matrix 3

  // Fishery
  init_vector cattot(styr,endyr)                 // Total catch in tons
  init_vector cattot_log_sd(styr,endyr)          // Total catch (cv) = sdev of log(cattot)
  init_int nyrs_fsh                              // Number of fishery age comps
  init_ivector fshyrs(1,nyrs_fsh)                // Years for the fishery age comps
  init_vector multN_fsh(1,nyrs_fsh)              // Multinomial sample size by year
  init_ivector ac_yng_fsh(1,nyrs_fsh)            // Accumulation of lower ages
  init_ivector ac_old_fsh(1,nyrs_fsh)            // Accumulation of upper ages
  init_int nyrslen_fsh                           // Number of fishery length comps
  init_ivector fshlenyrs(1,nyrslen_fsh)          // Years for the fishery length comps
  init_vector multNlen_fsh(1,nyrslen_fsh)        // Multinomial sample size by year
  init_vector rwlk_sd(styr,endyr-1)              // Random walk stdevs
  init_matrix catp(1,nyrs_fsh,rcrage,trmage)     // Catch proportions at age
  init_matrix lenp(1,nyrslen_fsh,1,nbins1)       // Catch proportions at age
  init_matrix wt_fsh(styr,endyr,rcrage,trmage)   // Weight at age by year
                                                 // For matrices indices for rows, then col
  // Survey 1 (Acoustic) EK500 (biosonic deleted in 2022)
  init_int nyrs_srv1                             // Number of survey biomass estimates
  init_ivector srvyrs1(1,nyrs_srv1)           // Years in which surveys occured
  init_vector indxsurv1(1,nyrs_srv1)           // Survey index
  init_vector indxsurv_log_sd1(1,nyrs_srv1)    // Survey index (cv) = sdev of log(indxsurv)
  init_vector q1_rwlk_sd(styr,endyr-1)              // Random walk stdevs
  init_vector yrfrct_srv1(styr,endyr)            // Fraction of year to midpoint of survey
  init_int nyrsac_srv1                           // Number of survey age comps
  init_ivector srv_acyrs1(1,nyrsac_srv1)         // Years for the survey age comp
  init_vector multN_srv1(1,nyrsac_srv1)          // Multinomial sample size by year
  init_ivector ac_yng_srv1(1,nyrsac_srv1)        // Accumulation of lower ages
  init_ivector ac_old_srv1(1,nyrsac_srv1)        // Accumulation of upper ages
  init_int nyrslen_srv1                          // Number of survey length comps
  init_ivector srv_lenyrs1(1,nyrslen_srv1)       // Years for the survey length comps
  init_vector multNlen_srv1(1,nyrslen_srv1)      // Multinomial sample size by year
  init_matrix srvp1(1,nyrsac_srv1,rcrage,trmage) // Survey proportions at age
  init_matrix srvlenp1(1,nyrslen_srv1,1,nbins3)  // Survey proportions at length
  init_matrix wt_srv1(styr,endyr,rcrage,trmage)  // Survey weights at age

  // Survey 2 (Bottom trawl)
  init_int nyrs_srv2                             // Number of surveys
  init_ivector srvyrs2(1,nyrs_srv2)              // Years in which surveys occured
  init_vector indxsurv2(1,nyrs_srv2)             // Survey index
  init_vector indxsurv_log_sd2(1,nyrs_srv2)      // Survey index (cv) = sdev of log(indxsurv)
  init_vector q2_rwlk_sd(styr,endyr-1)              // Random walk stdevs
  init_vector yrfrct_srv2(styr,endyr)            // Fraction of year to midpoint of survey
  init_int nyrsac_srv2                           // Number of survey age comps
  init_ivector srv_acyrs2(1,nyrsac_srv2)         // Years for the survey age comp
  init_vector multN_srv2(1,nyrsac_srv2)          // Multinomial sample size by year
  init_ivector ac_yng_srv2(1,nyrsac_srv2)        // Accumulation of lower ages
  init_ivector ac_old_srv2(1,nyrsac_srv2)        // Accumulation of upper ages
  init_int nyrslen_srv2                          // Number of survey length comps
  init_ivector srv_lenyrs2(1,nyrslen_srv2)       // Years for the survey length comps
  init_vector multNlen_srv2(1,nyrslen_srv2)      // Multinomial sample size by year
  init_matrix srvp2(1,nyrsac_srv2,rcrage,trmage) // Survey proportions at age
  init_matrix srvlenp2(1,nyrslen_srv2,1,nbins2)  // Survey proportions at length
  init_matrix wt_srv2(styr,endyr,rcrage,trmage)  // Survey weights at age

  // Survey 3 (ADFG coastal survey)
  init_int nyrs_srv3                             // Number of survey biomass estimates
  init_ivector srvyrs3(1,nyrs_srv3)              // Years in which surveys occured
  init_vector indxsurv3(1,nyrs_srv3)             // Survey index
  init_vector indxsurv_log_sd3(1,nyrs_srv3)      // Survey index (cv) = sdev of log(indxsurv)
  init_vector q3_rwlk_sd(styr,endyr-1)              // Random walk stdevs
  init_vector yrfrct_srv3(styr,endyr)            // Fraction of year to midpoint of survey
  init_int nyrsac_srv3                           // Number of survey age comps
  init_ivector srv_acyrs3(1,nyrsac_srv3)         // Years for the survey age comps
  init_vector multN_srv3(1,nyrsac_srv3)          // Multinomial sample size by year
  init_int nyrslen_srv3                          // Number of survey length comps
  init_ivector srv_lenyrs3(1,nyrslen_srv3)       // Years for the survey length comps
  init_vector multNlen_srv3(1,nyrslen_srv3)      // Multinomial sample size by year
  init_matrix srvp3(1,nyrsac_srv3,rcrage,trmage) // Survey proportions at age
  init_matrix srvlenp3(1,nyrslen_srv3,1,nbins2)  // Survey proportions at length
  init_matrix wt_srv3(styr,endyr,rcrage,trmage)  // Survey weights at age

  // Survey 4 (Age 1 acoustic)
  init_int nyrs_srv4                             // Number of surveys
  init_ivector srvyrs4(1,nyrs_srv4)              // Years in which surveys occured
  init_vector indxsurv4(1,nyrs_srv4)             // Survey index
  init_vector indxsurv_log_sd4(1,nyrs_srv4)      // Survey index (cv) = sdev of log(indxsurv)
 
  // Survey 5 (Age 2 acoustic)
  init_int nyrs_srv5                             // Number of surveys
  init_ivector srvyrs5(1,nyrs_srv5)              // Years in which surveys occured
  init_vector indxsurv5(1,nyrs_srv5)             // Survey index
  init_vector indxsurv_log_sd5(1,nyrs_srv5)      // Survey index (cv) = sdev of log(indxsurv)
 
  // Survey 6 (Summer acoustic)
  init_int nyrs_srv6                             // Number of surveys
  init_ivector srvyrs6(1,nyrs_srv6)              // Years in which surveys occured
  init_vector indxsurv6(1,nyrs_srv6)             // Survey index
  init_vector indxsurv_log_sd6(1,nyrs_srv6)      // Survey index (cv) = sdev of log(indxsurv)
  init_vector yrfrct_srv6(styr,endyr)            // Fraction of year to midpoint of survey
  init_int nyrsac_srv6                           // Number of survey age comps
  init_ivector srv_acyrs6(1,nyrsac_srv6)         // Years for the survey age comp
  init_vector multN_srv6(1,nyrsac_srv6)          // Multinomial sample size by year
  init_ivector ac_yng_srv6(1,nyrsac_srv6)        // Accumulation of lower ages
  init_ivector ac_old_srv6(1,nyrsac_srv6)        // Accumulation of upper ages
  init_int nyrslen_srv6                          // Number of survey length comps
  init_ivector srv_lenyrs6(1,nyrslen_srv6)       // Years for the survey length comps
  init_vector multNlen_srv6(1,nyrslen_srv6)      // Multinomial sample size by year
  init_matrix srvp6(1,nyrsac_srv6,rcrage,trmage) // Survey proportions at age
  init_matrix srvlenp6(1,nyrslen_srv6,1,nbins2)  // Survey proportions at length
  init_matrix wt_srv6(styr,endyr,rcrage,trmage)  // Survey weights at age

  // Ageing error transition matrix
  init_matrix age_trans(rcrage,trmage,rcrage,trmage)
  // Age to length transition matrix, to calculate expected length comps
  init_matrix len_trans1(rcrage,trmage,1,nbins1)
  init_matrix len_trans2(rcrage,trmage,1,nbins2)
  init_matrix len_trans3(rcrage,trmage,1,nbins3)

  // Population vectors
  matrix wt_pop(styr,endyr,rcrage,trmage)   // Population weight at age
  matrix wt_spawn(styr,endyr,rcrage,trmage) // Population weight at age at spawning (April 15)
  !! wt_pop=wt_srv2;
  !! wt_spawn=wt_srv1;
  init_vector mat(rcrage,trmage)                 // Proportion mature

  // Projection parameters
  vector wt_pop_proj(rcrage,trmage);
  vector wt_spawn_proj(rcrage,trmage);
  vector wt_fsh_proj(rcrage,trmage);
  vector wt_srv_proj(rcrage,trmage);
  LOCAL_CALCS
   // for projections take averages of WAA only from recent survey years with data
    wt_pop_proj.initialize();
    wt_spawn_proj.initialize();
    for(int a=rcrage;a<=trmage;a++){
      for(int i=1;i<=3;i++)
          wt_pop_proj(a)+=wt_srv2(srv_acyrs2(nyrsac_srv2-i+1),a)/3;
      for(int i=1;i<=5;i++){
          wt_spawn_proj(a)+=wt_srv1(srv_acyrs1(nyrsac_srv1-i+1),a)/5;
      // for predicting what the survey will see next year, Shelikof for now
          wt_srv_proj(a)+=wt_srv1(srv_acyrs1(nyrsac_srv1-i+1),a)/5;
      }	  
     }
     wt_fsh_proj=wt_fsh(endyr);
  END_CALCS
  
  init_vector Ftarget(endyr+1,endyr+5)
  init_number B40		                 // mean log recruitment
  init_number log_mean_recr_proj 
  init_number sigmasq_recr                       // Variance for log recr, recruitment indices

  init_number check
  !! if(check != -999){ cerr << "Failed to parse .dat file, final value=" << check <<endl; ad_exit(1);}
  int styr_avg_slct
  int endyr_avg_slct
  int i                                          // Index for year
  int j                                          // Index for age
  int loop

  number o                                       // A small number

 // After all the data is read in, redefine the endyr to be the
 // retroyr, then **carefully* redefine calculations throughout
 // to not overindex the data inputs
 int endyr0;			// original endyr
 LOC_CALCS
  if(retro_yrs<0){cerr << "bad peels in -retro option" << endl; ad_exit(1);};
  endyr0=endyr;
  endyr=endyr-retro_yrs;
 END_CALCS


INITIALIZATION_SECTION

  mean_log_initN      0.0   // Mean log initial age composition
  mean_log_recruit    0.0   // Mean log recruitment 
  mean_log_F        -1.6   // Mean log fishing mortality

  M                  0.30  // Natural mortality
  log_q1_mean        0.0   // Survey 1 catchability
  log_q1_dev         0.0   // Survey 1 catchability
  log_q2_mean        0.0   // Survey 2 catchability
  log_q2_dev         0.0  
  log_q3_mean       -1.6   // Survey 3 catchability
  log_q3_dev         0.0
  log_q6             0.0   // Survey 6 catchability
  
  slp1_fsh_dev       0.0   // Selectivity deviance terms
  inf1_fsh_dev       0.0
  slp2_fsh_dev       0.0
  inf2_fsh_dev       0.0

//Starting values for selectivity curves  
  log_slp1_fsh_mean  1.0 
  inf1_fsh_mean      4.0
  log_slp2_fsh_mean  1.0
  inf2_fsh_mean      8.0

//  log_slp1_srv1     -0.8
//  inf1_srv1          4.0
  log_slp2_srv1      0.5
  inf2_srv1          9.0

//  log_slp1_srv2      0.0
//  inf1_srv2          4.0
//  log_slp2_srv2      0.0
//  inf2_srv2          8.0

  log_slp1_srv2     -0.8
  inf1_srv2          4.0
  log_slp2_srv2     1
  inf2_srv2          20


  log_slp1_srv3      0.0
  inf1_srv3          5.0

  log_slp1_srv6     4.9
  inf1_srv6          0.5
//  log_slp2_srv6     0
//  inf2_srv6          7
  log_slp2_srv6     1
  inf2_srv6          20
  sigmaR 1.3 
  natMscalar   1
  
PARAMETER_SECTION

 // Population parameters
  //  init_bounded_number M(0.1,0.5,-1)
  init_bounded_vector M(rcrage,trmage,0.1,5.0,-1)
  init_bounded_number mean_log_initN(-15,15,-1)
  init_bounded_dev_vector dev_log_initN(rcrage+1,trmage,-15,15,-2)
  vector initN(rcrage+1,trmage)
  init_bounded_number mean_log_recruit(-15,15,1)
  init_bounded_dev_vector dev_log_recruit(styr,endyr,-15,15,3)
  init_bounded_number sigmaR(0,5,-1);
  sdreport_vector recruit(styr,endyr)
  sdreport_vector log_recruit(styr,endyr);

 // Forward projections 
  init_bounded_vector log_recr_proj(endyr+1,endyr+5,-5,5,10)
  sdreport_vector recruit_proj(endyr+1,endyr+5)  
  matrix N_proj(endyr+1,endyr+5,rcrage,trmage)
  vector F_proj(endyr+1,endyr+5)
  matrix Z_proj(endyr+1,endyr+5,rcrage,trmage)
  matrix C_proj(endyr+1,endyr+5,rcrage,trmage)
  matrix Nsrv_proj(endyr+1,endyr+5,rcrage,trmage)
  vector slctfsh_proj(rcrage,trmage)
  vector Ecattot_proj(endyr+1,endyr+5)  
  sdreport_vector Esumbio_proj(endyr+1,endyr+5)  
  sdreport_vector Espawnbio_proj(endyr+1,endyr+5) 
  sdreport_vector Esrv_proj(endyr+1,endyr+5) 
  sdreport_vector Exrate_proj(endyr+1,endyr+5)
  number sbio

 // Selectivity parameters

 // Fishery selectivity
  init_bounded_number log_slp1_fsh_mean(-5,5,4)
  init_bounded_number inf1_fsh_mean(1,5,4)
  init_bounded_number log_slp2_fsh_mean(-5,5,4)
  init_bounded_number inf2_fsh_mean(7,20,4)
  //    init_bounded_number log_slp2_fsh_mean(-5,5,-1)	  
  //    init_bounded_number inf2_fsh_mean(7,20,-1)	  
  init_bounded_dev_vector slp1_fsh_dev(styr,endyr,-5,5,7)
  init_bounded_dev_vector inf1_fsh_dev(styr,endyr,-5,5,7)
  //    init_bounded_dev_vector slp2_fsh_dev(styr,endyr,-5,5,6)
  //    init_bounded_dev_vector inf2_fsh_dev(styr,endyr,-5,5,6)
  init_bounded_dev_vector slp2_fsh_dev(styr,endyr,-5,5,-1)
  init_bounded_dev_vector inf2_fsh_dev(styr,endyr,-5,5,-1)

  vector slp1_fsh(styr,endyr)
  vector inf1_fsh(styr,endyr)
  vector slp2_fsh(styr,endyr)
  vector inf2_fsh(styr,endyr)

 // Acoustic survey selectivity
  //  init_bounded_number log_slp1_srv1(-10,5,7)
  //  init_bounded_number inf1_srv1(1,20,7)
  init_bounded_number log_slp2_srv1(-5,5,7)
  init_bounded_number inf2_srv1(3,10,7)
  // Trawl selectivity
  init_bounded_number log_slp1_srv2(-5,5,6)
  init_bounded_number inf1_srv2(1,50,6)
  //    init_bounded_number inf1_srv2(1,50,-1)
  init_bounded_number log_slp2_srv2(-5,5,-1)
  init_bounded_number inf2_srv2(5,25,-1)
  // ADFG selectivity
  init_bounded_number log_slp1_srv3(-5,5,9)
  init_bounded_number inf1_srv3(1,20,9)
  //    init_bounded_number log_slp2_srv3(-5,5,3)
  //    init_bounded_number inf2_srv3(3,20,3)
  // Summer acoustic selectivity
  init_bounded_number log_slp1_srv6(-5,5,-1)
  init_bounded_number inf1_srv6(0,50,-1)
  //    init_bounded_number inf1_srv6(1,50,-1)
  //    init_bounded_number inf1_srv6(1,50,-1)
  init_bounded_number log_slp2_srv6(-5,5,7)
  init_bounded_number inf2_srv6(5,25,7)

 // Fishing mortality and survey catchablility
  init_bounded_number mean_log_F(-10,10,1)
  init_bounded_dev_vector dev_log_F(styr,endyr,-10,10,2)
  vector F(styr,endyr)
  init_bounded_number log_q1_mean(-10,10,5)
  init_bounded_dev_vector log_q1_dev(styr,endyr,-5,5,5) 
   //  init_bounded_number log_q2(-10,10,-1)
  init_bounded_number log_q2_mean(-10,10,5)
  init_bounded_dev_vector log_q2_dev(styr,endyr,-5,5,-1)  
  init_bounded_number log_q3_mean(-10,10,6)
  init_bounded_dev_vector log_q3_dev(styr,endyr,-5,5,5) 
  init_bounded_number log_q4(-10,10,6)
  //  init_bounded_number q4_pow(-10,10,6)
  init_bounded_number q4_pow(-10,10,-6)
  init_bounded_number log_q5(-10,10,6)
  //  init_bounded_number q5_pow(-10,10,6)  
  init_bounded_number q5_pow(-10,10,-6)
  init_bounded_number log_q6(-10,10,5)
  // This scales M vector below so that M={M}*natMscalar. If 1 does nothing. 
  init_bounded_number natMscalar(0,5,-5)

 // Dependent parameters
  vector q1(styr,endyr)
  vector q2(styr,endyr)
  vector q3(styr,endyr)
  sdreport_vector log_q1(styr,endyr)
  sdreport_vector log_q2(styr,endyr)
  sdreport_vector log_q3(styr,endyr)
  number q4
  number q5
  number q6
  matrix N(styr,endyr,rcrage,trmage)
  sdreport_vector endN(rcrage,trmage)  
  matrix Z(styr,endyr,rcrage,trmage)
  matrix C(styr,endyr,rcrage,trmage)
  matrix Nsrv1(styr,endyr,rcrage,trmage)
  sdreport_vector slctsrv1(rcrage,trmage)
  matrix Nsrv2(styr,endyr,rcrage,trmage)
  sdreport_vector slctsrv2(rcrage,trmage)
  matrix Nsrv3(styr,endyr,rcrage,trmage)
  sdreport_vector slctsrv3(rcrage,trmage)
  matrix Nsrv6(styr,endyr,rcrage,trmage)
  vector slctsrv6(rcrage,trmage)
  matrix slctfsh(styr,endyr,rcrage,trmage)
  sdreport_vector slctsrv1_logit(rcrage,trmage)
  sdreport_vector slctsrv2_logit(rcrage,trmage)
  sdreport_vector slctsrv3_logit(rcrage,trmage)
  sdreport_vector slctsrv6_logit(rcrage,trmage)
  sdreport_vector slctfsh_logit(rcrage,trmage);


// Expected values
  vector Eecocon(styr,endyr)
  matrix Eec(styr,endyr,rcrage,trmage)
  vector Ecattot(styr,endyr)
  matrix Ecatp(styr,endyr,rcrage,trmage)
  matrix Elenp(styr,endyr,1,nbins1)
  vector Eindxsurv1(styr,endyr)
  matrix Esrvp1(styr,endyr,rcrage,trmage)
  matrix Esrvlenp1(styr,endyr,1,nbins3)
  vector Eindxsurv2(styr,endyr)
  matrix Esrvp2(styr,endyr,rcrage,trmage)
  matrix Esrvlenp2(styr,endyr,1,nbins2)
  vector Eindxsurv3(styr,endyr)
  matrix Esrvp3(styr,endyr,rcrage,trmage)
  matrix Esrvlenp3(styr,endyr,1,nbins2)
  vector Eindxsurv4(styr,endyr) 
  vector Eindxsurv5(styr,endyr) 
  vector Eindxsurv6(styr,endyr)
  matrix Esrvp6(styr,endyr,rcrage,trmage)
  matrix Esrvlenp6(styr,endyr,1,nbins2)
  sdreport_vector Espawnbio(styr,endyr)
  sdreport_vector Esumbio(styr,endyr)
  // log versions make more sense for asymptotics
  sdreport_vector Espawnbio_log(styr,endyr)
  sdreport_vector Esumbio_log(styr,endyr)
  vector Espawnbio_2plus(styr,endyr)
  vector Etotalbio(styr,endyr)
  // log likelihood containers
  vector loglik(1,24)
  vector llcatp(1,nyrs_fsh)
  vector lllenp(1,nyrslen_fsh)
  vector llsrvp1(1,nyrsac_srv1)
  vector llsrvlenp1(1,nyrslen_srv1)
  vector llsrvp2(1,nyrsac_srv2)
  vector llsrvlenp2(1,nyrslen_srv2)
  vector llsrvp3(1,nyrsac_srv3)
  vector llsrvlenp3(1,nyrslen_srv3)
  vector llsrvp6(1,nyrsac_srv6)
  vector llsrvlenp6(1,nyrslen_srv6)
 

  //residual output matrices
  matrix res_fish(1,nyrs_fsh,rcrage,2*trmage-rcrage+1)
  matrix res_srv1(1,nyrsac_srv1,rcrage,2*trmage-rcrage+1)
  matrix res_srv2(1,nyrsac_srv2,rcrage,2*trmage-rcrage+1)
  matrix res_srv3(1,nyrsac_srv3,rcrage,2*trmage-rcrage+1)
  matrix res_srv3len(1,nyrslen_srv3,1,2*nbins2)
  matrix res_srv6(1,nyrsac_srv6,rcrage,2*trmage-rcrage+1)
  matrix pearson_fish(1,nyrs_fsh,rcrage,trmage)
  matrix pearson_srv1(1,nyrsac_srv1,rcrage,trmage)
  matrix pearson_srv2(1,nyrsac_srv2,rcrage,trmage)
  matrix pearson_srv3(1,nyrsac_srv3,rcrage,trmage)
  matrix pearson_srv3len(1,nyrslen_srv3,1,nbins2)
  matrix pearson_srv6(1,nyrsac_srv6,rcrage,trmage)

  vector effN_fsh(1,nyrs_fsh)
  vector effN_srv1(1,nyrsac_srv1)
  vector effN_srv2(1,nyrsac_srv2)
  vector effN_srv3(1,nyrsac_srv3)
  vector effN_srv6(1,nyrsac_srv6)

  number RMSE_srv1
  number RMSE_srv2
  number RMSE_srv3
  number RMSE_srv4
  number RMSE_srv5
  number RMSE_srv6

 // Objective function
  likeprof_number var_prof
  objective_function_value objfun

PRELIMINARY_CALCS_SECTION
 // Do upper and lower accumulations in age composition data
 // Fishery
  for (i=1;i<=nyrs_fsh;i++) {
   for (j=rcrage;j<=trmage;j++) {
     if(j<ac_yng_fsh(i)) {
       catp(i,ac_yng_fsh(i)) += catp(i,j);
       catp(i,j) = 0;
     }
     if(j>ac_old_fsh(i)) {
       catp(i,ac_old_fsh(i)) += catp(i,j);
       catp(i,j) = 0;
     }
    }
  }

// Survey 1
 for (i=1;i<=nyrsac_srv1;i++) {
   for (j=rcrage;j<=trmage;j++) {
     if(j<ac_yng_srv1(i)) {
       srvp1(i,ac_yng_srv1(i)) += srvp1(i,j);
       srvp1(i,j) = 0;
     }
     if(j>ac_old_srv1(i)) {
       srvp1(i,ac_old_srv1(i)) += srvp1(i,j);
       srvp1(i,j) = 0;
     }
   }
 }

// Survey 2
  for (i=1;i<=nyrsac_srv2;i++) {
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_srv2(i)) {
	srvp2(i,ac_yng_srv2(i)) += srvp2(i,j);
	srvp2(i,j) = 0;
      }
      if(j>ac_old_srv2(i)) {
	srvp2(i,ac_old_srv2(i)) += srvp2(i,j);
	srvp2(i,j) = 0;
      }
    }
  }
	
  // Survey 6
  for (i=1;i<=nyrsac_srv6;i++) {
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_srv6(i)) {
	srvp6(i,ac_yng_srv6(i)) += srvp6(i,j);
	srvp6(i,j) = 0;
      }
      if(j>ac_old_srv6(i)) {
	srvp6(i,ac_old_srv6(i)) += srvp6(i,j);
	srvp6(i,j) = 0;
      }
    }
  }

  o = 0.00001;
  var_prof.set_stepnumber(30);
  //  var_prof.set_stepnumber(10);
  var_prof.set_stepsize(0.1);


PROCEDURE_SECTION

  Convert_log_parameters();
  Selectivity();
  Mortality();
  Numbers_at_age();
  Catch_at_age();
  Expected_values();
  if(last_phase())
  {
    Projections();
  } 
  Objective_function();
  MCMC_output();

FUNCTION Convert_log_parameters
 // Assume close to equilibrium at f=0 at the start   
 //Use the initial recruitment dev to fill out the initial age composition
  initN(rcrage+1) = mfexp(mean_log_recruit +  dev_log_recruit(styr) - M(rcrage));
  for (j=rcrage+2;j<=trmage;j++) {
    initN(j) = initN(j-1)*mfexp(-M(j));
  }
  initN(trmage) /= (1.0 - mfexp(-M(trmage)));
 //devs for initN are turned off 
 for (j=rcrage+1;j<=trmage;j++) {
   initN(j) = initN(j)*mfexp(dev_log_initN(j));
 }
 log_recruit=mean_log_recruit +  dev_log_recruit;
 recruit = mfexp(log_recruit);

 // Acoustic survey random walk setup
 for (i=styr;i<=endyr;i++) {
   log_q1(i)=log_q1_mean+log_q1_dev(i);
   log_q2(i)=log_q2_mean+log_q2_dev(i);
   log_q3(i)=log_q3_mean+log_q3_dev(i); 
   q1(i)=mfexp(log_q1(i));
   q2(i)=mfexp(log_q2(i));
   q3(i)=mfexp(log_q3(i));
 }
 F = mfexp(mean_log_F + dev_log_F);
 q4 = mfexp(log_q4);
 q5 = mfexp(log_q5);
 q6 = mfexp(log_q6);
 
FUNCTION Selectivity

 // Fishery selectivity
 for (i=styr;i<=endyr;i++) {
   slp1_fsh(i)=mfexp(log_slp1_fsh_mean+slp1_fsh_dev(i));
   inf1_fsh(i)=inf1_fsh_mean+inf1_fsh_dev(i);
   slp2_fsh(i)=mfexp(log_slp2_fsh_mean+slp2_fsh_dev(i));
   inf2_fsh(i)=inf2_fsh_mean+inf2_fsh_dev(i);
   for (j=rcrage;j<=trmage;j++) {
     slctfsh(i,j) = (1/(1+mfexp(-(slp1_fsh(i))*(double(j)-(inf1_fsh(i))))))*
       (1-1/(1+mfexp(-(slp2_fsh(i))*(double(j)-(inf2_fsh(i))))));
   }
 // The plan would be to check and adjust the max selected age as needed
   slctfsh(i)=slctfsh(i)/slctfsh(i,7);
 }

 M(1)=1.39; M(2)=0.69; M(3)=0.48; M(4)=0.37; M(5)=0.34;
 M(6)=0.30; M(7)=0.30; M(8)=0.29; M(9)=0.28; M(10)=0.29;
 for(int i=1; i<=10; i++) M(i)*=natMscalar;

//Survey 1 selectivity
 for (j=rcrage;j<=trmage;j++) {
   slctsrv1(j) = (1-1/(1+mfexp(-mfexp(log_slp2_srv1)*(double(j)-inf2_srv1))));
 }
 slctsrv1=slctsrv1/slctsrv1(3);
 slctsrv1(rcrage)=0;
 slctsrv1(rcrage+1)=0;	
 
//Survey 2 selectivity
 for (j=rcrage;j<=trmage;j++) {
   slctsrv2(j) = (1/(1+mfexp(-mfexp(log_slp1_srv2)*(double(j)-inf1_srv2))))
     *(1-1/(1+mfexp(-mfexp(log_slp2_srv2)*(double(j)-inf2_srv2))));
 }
 slctsrv2=slctsrv2/slctsrv2(10);

//Survey 3 selectivity
  for (j=rcrage;j<=trmage;j++) {
    slctsrv3(j) = (1/(1+mfexp(-mfexp(log_slp1_srv3)*(double(j)-inf1_srv3))));
  }
  slctsrv3=slctsrv3/slctsrv3(10);

 // Survey 6 selectivity
  for (j=rcrage;j<=trmage;j++) {
    slctsrv6(j) = (1/(1+mfexp(-mfexp(log_slp1_srv6)*(double(j)-inf1_srv6))))
      *(1-1/(1+mfexp(-mfexp(log_slp2_srv6)*(double(j)-inf2_srv6))));
  }
  slctsrv6=slctsrv6/slctsrv6(1);

  // Calculate uncertainty in logit space where it makes more
  // sense. But it does break when selex=0 or 1 identically so
  // need to be careful there.
  slctsrv1_logit=-log(1/(slctsrv1-1e-10)-1);
  slctsrv2_logit=-log(1/(slctsrv2-1e-10)-1);
  slctsrv3_logit=-log(1/(slctsrv3-1e-10)-1);
  slctsrv6_logit=-log(1/(slctsrv6-1e-10)-1);
  slctfsh_logit=-log(1/(slctfsh(endyr-1)-1e-10)-1);


FUNCTION Mortality
  for (i=styr;i<=endyr;i++) {
    for (j=rcrage;j<=trmage;j++) {
      Z(i,j)=(F(i)*slctfsh(i,j))+M(j);
    }
  }  


FUNCTION Numbers_at_age
  N(styr)(rcrage+1,trmage)=initN;
  for (i=styr;i<=endyr;i++) {
    N(i,rcrage)=recruit(i);
  }
  for (i=styr;i<endyr;i++) {
    for (j=rcrage;j<trmage;j++) {
      N(i+1,j+1)=N(i,j)*mfexp(-Z(i,j));
    }  
    N(i+1,trmage)+=N(i,trmage)*mfexp(-Z(i,trmage));
  }
  endN=N(endyr);


FUNCTION Catch_at_age
 // Catch at age and survey numbers at age
  for (i=styr;i<=endyr;i++) {
    for (j=rcrage;j<=trmage;j++) {
      C(i,j)=N(i,j)*((F(i)*slctfsh(i,j))/Z(i,j))*(1-mfexp(-Z(i,j)));
      Eec(i,j)=N(i,j)*(M(j)/Z(i,j))*(1-mfexp(-Z(i,j)));
      Nsrv1(i,j)=slctsrv1(j)*N(i,j)*mfexp(-yrfrct_srv1(i)*Z(i,j));
      Nsrv2(i,j)=slctsrv2(j)*N(i,j)*mfexp(-yrfrct_srv2(i)*Z(i,j));
      Nsrv3(i,j)=slctsrv3(j)*N(i,j)*mfexp(-yrfrct_srv3(i)*Z(i,j));
      Nsrv6(i,j)=slctsrv6(j)*N(i,j)*mfexp(-yrfrct_srv6(i)*Z(i,j));
     }
  }

FUNCTION Expected_values
 for (i=styr;i<=endyr;i++){
   Ecattot(i) = 1000000*sum(elem_prod(C(i),wt_fsh(i)));
   Eecocon(i) = 1000000*sum(elem_prod(Eec(i),wt_pop(i)));
   Ecatp(i) = (C(i)/sum(C(i)))*age_trans;
   Elenp(i) = Ecatp(i) * len_trans1;
   Eindxsurv1(i)= q1(i)*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv1(i)*Z(i))),slctsrv1),wt_srv1(i)));
   Esrvp1(i) = (Nsrv1(i)/sum(Nsrv1(i)))*age_trans;
   Esrvlenp1(i) = Esrvp1(i) * len_trans3;
   Eindxsurv2(i)= q2(i)*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv2(i)*Z(i))),slctsrv2),wt_srv2(i)));
   Esrvp2(i) = (Nsrv2(i)/sum(Nsrv2(i)))*age_trans;
   Esrvlenp2(i) = Esrvp2(i) * len_trans2;
   Eindxsurv3(i)= q3(i)*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv3(i)*Z(i))),slctsrv3),wt_srv3(i)));
   Esrvp3(i) = (Nsrv3(i)/sum(Nsrv3(i)))*age_trans;
   Esrvlenp3(i) = Esrvp3(i) * len_trans2;
   Eindxsurv4(i)= q4*pow(N(i,1),(q4_pow+1));
   Eindxsurv5(i)= q5*pow(N(i,2),(q5_pow+1));	
   Eindxsurv6(i)= q6*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv6(i)*Z(i))),slctsrv6),wt_srv6(i)));
   Esrvp6(i) = (Nsrv6(i)/sum(Nsrv6(i)))*age_trans;
   Esrvlenp6(i) = Esrvp6(i) * len_trans2;
   // 3+ biomass
   Esumbio(i)= N(i)(rcrage+2,trmage)*wt_pop(i)(rcrage+2,trmage);
   // Total biomass
   Etotalbio(i)= N(i)(rcrage,trmage)*wt_pop(i)(rcrage,trmage);
   // 2+ biomass at spawning (use for apportioning to management area)
   Espawnbio_2plus(i)= sum(elem_prod(elem_prod(N(i)(rcrage+1,trmage),mfexp(-yrfrct_srv1(i)*Z(i)(rcrage+1,trmage))),wt_srv1(i)(rcrage+1,trmage)));
   // //1+ biomass at spawning
   //     Esumbio(i)= sum(elem_prod(elem_prod(N(i)(rcrage,trmage),mfexp(-yrfrct_srv1(i)*Z(i)(rcrage,trmage))),wt_srv1(i)(rcrage,trmage)));
   Espawnbio(i)= sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-0.21*Z(i))),wt_spawn(i)),0.5*mat));
 }
 Espawnbio_log=log(Espawnbio);
 Esumbio_log=log(Esumbio);
 

 // Do upper and lower accumulation in expected age composition

 // Fishery
  for (i=1;i<=nyrs_fsh;i++) {
    if(fshyrs(i)>endyr) break;
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_fsh(i)) {
	Ecatp(fshyrs(i),ac_yng_fsh(i)) += Ecatp(fshyrs(i),j);
	Ecatp(fshyrs(i),j) = 0;
      }
      if(j>ac_old_fsh(i)) {
	Ecatp(fshyrs(i),ac_old_fsh(i)) += Ecatp(fshyrs(i),j);
	Ecatp(fshyrs(i),j) = 0;
      }
    }
  } 
 // Survey 1
  for (i=1;i<=nyrsac_srv1;i++) {
    if(srv_acyrs1(i)>endyr) break;
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_srv1(i)) {
	Esrvp1(srv_acyrs1(i),ac_yng_srv1(i)) += Esrvp1(srv_acyrs1(i),j);
	Esrvp1(srv_acyrs1(i),j) = 0;
      }
      if(j>ac_old_srv1(i)) {
	Esrvp1(srv_acyrs1(i),ac_old_srv1(i)) += Esrvp1(srv_acyrs1(i),j);
	Esrvp1(srv_acyrs1(i),j) = 0;
      }
    }
  }

  // Survey 2
  for (i=1;i<=nyrsac_srv2;i++) {
    if(srv_acyrs2(i)>endyr) break;
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_srv2(i)) {
	Esrvp2(srv_acyrs2(i),ac_yng_srv2(i)) += Esrvp2(srv_acyrs2(i),j);
	Esrvp2(srv_acyrs2(i),j) = 0.;
      }
      if(j>ac_old_srv2(i)) {
	Esrvp2(srv_acyrs2(i),ac_old_srv2(i)) += Esrvp2(srv_acyrs2(i),j);
	Esrvp2(srv_acyrs2(i),j) = 0;
      }
    }
  }

  // Survey 6
  for (i=1;i<=nyrsac_srv6;i++) {
    if(srv_acyrs6(i)>endyr) break;
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_srv6(i)) {
	Esrvp2(srv_acyrs6(i),ac_yng_srv6(i)) += Esrvp6(srv_acyrs6(i),j);
	Esrvp6(srv_acyrs6(i),j) = 0.;
      }
      if(j>ac_old_srv6(i)) {
	Esrvp6(srv_acyrs6(i),ac_old_srv6(i)) += Esrvp6(srv_acyrs6(i),j);
	Esrvp6(srv_acyrs6(i),j) = 0;
      }
    }
  }
	
FUNCTION Projections

//Recruitments 
 for (i=endyr+1;i<=endyr+5;i++) {
   //  recruitment with bias correction
   //    recruit_proj(i)=mfexp(log_recr_proj(i)+(sigmasq_recr/2));
   // for MCMC projections to get the prob < B20
   //   recruit_proj(i)=mfexp(log_recr_proj(i));
   // or just use average recruitment after 1977
   // note that endyr-1 is the last year for mean
   // standard method
   recruit_proj(i)=mean(recruit(1978,endyr-1));
 }
 for (i=endyr+1;i<=endyr+5;i++) {
   N_proj(i,rcrage)=recruit_proj(i);
 }
//Initialize the age composition

//  Standard projection
 for (j=rcrage;j<trmage;j++) {
   N_proj(endyr+1,j+1)=N(endyr,j)*mfexp(-Z(endyr,j));
 }  
 N_proj(endyr+1,trmage)+=N(endyr,trmage)*mfexp(-Z(endyr,trmage));

 // Set 2007 year class to mean
 // note that endyr-1 is the last year for mean
 //   for (j=rcrage;j<trmage;j++)
 //   {
 //   N_proj(endyr+1,j+1)=N(endyr,j)*mfexp(-Z(endyr,j));
 //   }  
 //   N_proj(endyr+1,rcrage+1)=mean(recruit(1979,endyr-1))*mfexp(-Z(endyr,rcrage));
 //   N_proj(endyr+1,trmage)+=N(endyr,trmage)*mfexp(-Z(endyr,trmage));

 // Averaging window for selectivity 
 endyr_avg_slct=endyr-1;
 styr_avg_slct=endyr_avg_slct-4;

 for (j=rcrage;j<=trmage;j++) {
   slctfsh_proj(j) = 0;
   for (i=styr_avg_slct;i<=endyr_avg_slct;i++) {
     slctfsh_proj(j) += slctfsh(i,j);
   }
 }
 slctfsh_proj=slctfsh_proj/max(slctfsh_proj);
 


 //Forward projections
  for (i=endyr+1;i<=endyr+5;i++) {
    // have to get Z twice
    //  Tuning loop to get the spawning biomass adjustment right
    F_proj(i)=Ftarget(i+retro_yrs);
    for (loop=1;loop<=20;loop++) {
      for (j=rcrage;j<=trmage;j++) {
	Z_proj(i,j)=(F_proj(i)*slctfsh_proj(j))+M(j);
      }  
      sbio = sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-0.21*Z_proj(i))),wt_spawn_proj),0.5*mat));
      //  Set the fishing mortality rate
      F_proj(i)=Ftarget(i+retro_yrs);
      if (sbio < B40) {
	F_proj(i)=Ftarget(i+retro_yrs)*(((sbio/B40)-0.05)/(1-0.05));
	// SSL control rule
	//   F_proj(i)=Ftarget(i+retro_yrs)*(((sbio/B40)-0.2)/(1-0.2));
      }
    }
  // Total mortality
    for (j=rcrage;j<=trmage;j++) {
      Z_proj(i,j)=(F_proj(i)*slctfsh_proj(j))+M(j);
    }  
  //  Numbers at age
   if(i<endyr+5) {
     for (j=rcrage;j<trmage;j++) {
       N_proj(i+1,j+1)=N_proj(i,j)*mfexp(-Z_proj(i,j));
     }  
     N_proj(i+1,trmage)+=N_proj(i,trmage)*mfexp(-Z_proj(i,trmage)); 
   }
  // Catches
  for (j=rcrage;j<=trmage;j++) {
    C_proj(i,j)=N_proj(i,j)*((F_proj(i)*slctfsh_proj(j))/Z_proj(i,j))*(1-mfexp(-Z_proj(i,j)));
    //    Nsrv_proj(i,j)=q2*slctsrv1(j)*N_proj(i,j)*mfexp(-yrfrct_srv1(endyr)*Z_proj(i,j));  
    Nsrv_proj(i,j)=N_proj(i,j)*mfexp(-yrfrct_srv6(endyr)*Z_proj(i,j));  
  }
  //  Total catches and biomass
  Ecattot_proj(i) = 1000000*sum(elem_prod(C_proj(i),wt_fsh_proj));
  // 3+ biomass
  Esumbio_proj(i)= N_proj(i)(rcrage+2,trmage)*wt_pop_proj(rcrage+2,trmage);
  // Alternative: 2+ biomass
  //    Esumbio_proj(i)= N_proj(i)(rcrage+1,trmage)*wt_pop_proj(rcrage+1,trmage);
  Exrate_proj(i)=Ecattot_proj(i)/(1000000*Esumbio_proj(i));
  Espawnbio_proj(i)= sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-0.21*Z_proj(i))),wt_spawn_proj),0.5*mat));
  //Summer acoustic
  //    Esrv_proj(i)= q6*sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-yrfrct_srv6(endyr)*Z_proj(i))),slctsrv6),wt_srv_proj));
  //Winter acoutic 
  Esrv_proj(i)= q1(endyr)*sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-yrfrct_srv1(endyr)*Z_proj(i))),slctsrv1),wt_srv_proj));
  }

FUNCTION Objective_function
 // Fishery likelihoods
 //Total catch
  //loglik(1) = -.5*norm2(elem_div((log(cattot)-log(Ecattot)),cattot_log_sd));
 loglik(1)=0;
 for(i=styr; i<=endyr;i++){
   if(i>endyr) break;
   loglik(1) += -.5*square((log(cattot(i))-log(Ecattot(i)))/cattot_log_sd(i));
 }	  

//Age composition
 loglik(2)=0;
 for (i=1;i<=nyrs_fsh;i++) {
    if(fshyrs(i)>endyr) break; 	// ignore data after retroyear
    llcatp(i) = 0;
    for (j=ac_yng_fsh(i);j<=ac_old_fsh(i);j++) {
      llcatp(i) += multN_fsh(i)*(catp(i,j)+o)*log((Ecatp(fshyrs(i),j)+o)/(catp(i,j)+o));
      res_fish(i,j)=catp(i,j);
      res_fish(i,trmage-rcrage+j+1)=Ecatp(fshyrs(i),j);
      if(multN_fsh(i)>0) {
	pearson_fish(i,j)=(catp(i,j)-Ecatp(fshyrs(i),j))/sqrt((Ecatp(fshyrs(i),j)*(1.-Ecatp(fshyrs(i),j)))/multN_fsh(i));	
      }
    }
    if(multN_fsh(i)>0) {	
      effN_fsh(i) = sum(elem_prod(Ecatp(fshyrs(i)),(1-Ecatp(fshyrs(i)))))/sum(square(catp(i)-Ecatp(fshyrs(i))));
    }
    loglik(2) += llcatp(i);
 }
 //Length composition
 loglik(3)=0;
 for (i=1;i<=nyrslen_fsh;i++) {
   if(fshlenyrs(i)>endyr) break;
   lllenp(i) = 0;
   for (j=1;j<=nbins1;j++) {
     lllenp(i) += multNlen_fsh(i)*(lenp(i,j)+o)*log((Elenp(fshlenyrs(i),j)+o)/(lenp(i,j)+o));
   }
   loglik(3) += lllenp(i);
 }

 // Survey 1 likelihoods
 // Total biomass  
  loglik(4)=0;
  for(i=1; i<=nyrs_srv1;i++){
    if(srvyrs1(i)>endyr) break;
    loglik(4)+=-.5*square((log(indxsurv1(i))-log(Eindxsurv1(srvyrs1(i)))+square(indxsurv_log_sd1(i))/2.)/indxsurv_log_sd1(i));
  }

  RMSE_srv1=0;
  if(!isretro)
    RMSE_srv1= sqrt(norm2(log(indxsurv1)-log(Eindxsurv1(srvyrs1))+square(indxsurv_log_sd1)/2.)/nyrs_srv1);
  
 //Age composition
  loglik(5)=0;
  for (i=1;i<=nyrsac_srv1;i++) {
    if(srv_acyrs1(i)>endyr) break;
    llsrvp1(i) = 0;
    for (j=ac_yng_srv1(i);j<=ac_old_srv1(i);j++) {
      llsrvp1(i) += multN_srv1(i)*(srvp1(i,j)+o)*log((Esrvp1(srv_acyrs1(i),j)+o)/(srvp1(i,j)+o));
      res_srv1(i,j)=srvp1(i,j);
      res_srv1(i,trmage-rcrage+j+1)=Esrvp1(srv_acyrs1(i),j);
      if(multN_srv1(i)>0) {
	pearson_srv1(i,j)=(srvp1(i,j)-Esrvp1(srv_acyrs1(i),j))/sqrt((Esrvp1(srv_acyrs1(i),j)*(1.-Esrvp1(srv_acyrs1(i),j)))/multN_srv1(i));	
      }
    }
    if(multN_srv1(i)>0) {
      effN_srv1(i) = sum(elem_prod(Esrvp1(srv_acyrs1(i)),(1-Esrvp1(srv_acyrs1(i)))))/sum(square(srvp1(i)-Esrvp1(srv_acyrs1(i))));
    }
    loglik(5) += llsrvp1(i);
  }
 //length composition
 loglik(6)=0;
 for (i=1;i<=nyrslen_srv1;i++) {
   if(srv_lenyrs1(i)>endyr) break;
   llsrvlenp1(i) = 0;
   for (j=1;j<=nbins3;j++) {
     llsrvlenp1(i) += multNlen_srv1(i)*(srvlenp1(i,j)+o)*log((Esrvlenp1(srv_lenyrs1(i),j)+o)/(srvlenp1(i,j)+o));
   }
   loglik(6) += llsrvlenp1(i);
 }

 // Survey 2 likelihoods
 //Total biomass    
  loglik(7) =0;
  for (i=1;i<=nyrs_srv2;i++){
    if(srvyrs2(i)>endyr) break;
    loglik(7)+=-.5*square((log(indxsurv2(i))-log(Eindxsurv2(srvyrs2(i)))+square(indxsurv_log_sd2(i))/2.)/indxsurv_log_sd2(i));
  }
  RMSE_srv2=0;
  if(!isretro)
    RMSE_srv2= sqrt(norm2(log(indxsurv2)-log(Eindxsurv2(srvyrs2))+square(indxsurv_log_sd2)/2.)/nyrs_srv2);
   
 // Age composition
  loglik(8)=0;
  for (i=1;i<=nyrsac_srv2;i++) {
    llsrvp2(i) = 0;
   if(srv_acyrs2(i)>endyr) break;
   for (j=ac_yng_srv2(i);j<=ac_old_srv2(i);j++) {
     llsrvp2(i) += multN_srv2(i)*(srvp2(i,j)+o)*log((Esrvp2(srv_acyrs2(i),j)+o)/(srvp2(i,j)+o));
     res_srv2(i,j)=srvp2(i,j);
     res_srv2(i,trmage-rcrage+j+1)=Esrvp2(srv_acyrs2(i),j);
     if(multN_srv2(i)>0) {
       pearson_srv2(i,j)=(srvp2(i,j)-Esrvp2(srv_acyrs2(i),j))/sqrt((Esrvp2(srv_acyrs2(i),j)*(1.-Esrvp2(srv_acyrs2(i),j)))/multN_srv2(i));	
     }  
   }
   if(multN_srv2(i)>0) {	
     effN_srv2(i) = sum(elem_prod(Esrvp2(srv_acyrs2(i)),(1-Esrvp2(srv_acyrs2(i)))))/sum(square(srvp2(i)-Esrvp2(srv_acyrs2(i))));
   }
   loglik(8) += llsrvp2(i);
  }
  // length composition
  loglik(9)=0;
  for (i=1;i<=nyrslen_srv2;i++) {
    if(srv_lenyrs2(i) > endyr) break;
    llsrvlenp2(i) = 0;
    for (j=1;j<=nbins2;j++) {
      llsrvlenp2(i) += multNlen_srv2(i)*(srvlenp2(i,j)+o)*log((Esrvlenp2(srv_lenyrs2(i),j)+o)/(srvlenp2(i,j)+o));
    }
    loglik(9) += llsrvlenp2(i);
  }
  loglik(10) = 0;

// Survey 3 likelihoods
//Total biomass
    loglik(11)=0;
    for(i=1; i<=nyrs_srv3;i++){
      if(srvyrs3(i)>endyr) break;
      loglik(11) += -.5*square((log(indxsurv3(i))-log(Eindxsurv3(srvyrs3(i)))+square(indxsurv_log_sd3(i))/2.)/indxsurv_log_sd3(i));
    }
    RMSE_srv3=0;
    if(!isretro)
      RMSE_srv3= sqrt(norm2(log(indxsurv3)-log(Eindxsurv3(srvyrs3))+square(indxsurv_log_sd3)/2.)/nyrs_srv3);

 // age composition
    loglik(12)=0;
  for (i=1;i<=nyrsac_srv3;i++) {
    if(srv_acyrs3(i)>endyr) break;
    llsrvp3(i) = 0;
    for (j=rcrage;j<=trmage;j++) {
      llsrvp3(i) += multN_srv3(i)*(srvp3(i,j)+o)*log((Esrvp3(srv_acyrs3(i),j)+o)/(srvp3(i,j)+o));
      res_srv3(i,j)=srvp3(i,j);
      res_srv3(i,trmage-rcrage+j+1)=Esrvp3(srv_acyrs3(i),j);
      if(multN_srv3(i)>0) {
	pearson_srv3(i,j)=(srvp3(i,j)-Esrvp3(srv_acyrs3(i),j))/sqrt((Esrvp3(srv_acyrs3(i),j)*(1.-Esrvp3(srv_acyrs3(i),j)))/multN_srv3(i));	
      }
    }
    if(multN_srv3(i)>0) {		
      effN_srv3(i) = sum(elem_prod(Esrvp3(srv_acyrs3(i)),(1-Esrvp3(srv_acyrs3(i)))))/sum(square(srvp3(i)-Esrvp3(srv_acyrs3(i))));	
    }
    loglik(12) += llsrvp3(i);
  }
 // length composition
  loglik(13)=0;
  for (i=1;i<=nyrslen_srv3;i++) {
    if(srv_lenyrs3(i)>endyr) break;
    llsrvlenp3(i) = 0;
    for (j=1;j<=nbins2;j++) {
      llsrvlenp3(i) += multNlen_srv3(i)*(srvlenp3(i,j)+o)*log((Esrvlenp3(srv_lenyrs3(i),j)+o)/(srvlenp3(i,j)+o));
      res_srv3len(i,j)=srvlenp3(i,j);
      res_srv3len(i,nbins2+j)=Esrvlenp3(srv_lenyrs3(i),j);
      if(multNlen_srv3(i)>0) {
	pearson_srv3len(i,j)=(srvlenp3(i,j)-Esrvlenp3(srv_lenyrs3(i),j))/sqrt((Esrvlenp3(srv_lenyrs3(i),j)*(1.-Esrvlenp3(srv_lenyrs3(i),j)))/multNlen_srv3(i));	
      }
    }
    loglik(13) += llsrvlenp3(i);
  }

// Survey 4 and 5 likelihoods
  loglik(14)=0; loglik(15)=0;
  for(i=1; i<=nyrs_srv4;i++){ 	// assuming srv4 and srv5 have identical structure
    if(srvyrs4(i) >endyr) break;
    loglik(14) += -.5*square((log(indxsurv4(i))-log(Eindxsurv4(srvyrs4(i)))+square(indxsurv_log_sd4(i))/2.)/indxsurv_log_sd4(i));
   //   loglik(14) = 0;
    loglik(15) += -.5*square((log(indxsurv5(i))-log(Eindxsurv5(srvyrs5(i)))+square(indxsurv_log_sd5(i))/2.)/indxsurv_log_sd5(i));
  }
   RMSE_srv5=0;
   if(!isretro)
     RMSE_srv5= sqrt(norm2(log(indxsurv5)-log(Eindxsurv5(srvyrs5))+square(indxsurv_log_sd5)/2.)/nyrs_srv5);
   RMSE_srv4=0;
   if(!isretro)
     RMSE_srv4= sqrt(norm2(log(indxsurv4)-log(Eindxsurv4(srvyrs4))+square(indxsurv_log_sd4)/2.)/nyrs_srv4);

// Survey 6 likelihoods
  //Total biomass    
 loglik(16)=0;
 for(i=1;i<=nyrs_srv6;i++){
   if(srvyrs6(i)>endyr) break;
   loglik(16)+=-.5*square((log(indxsurv6(i))-log(Eindxsurv6(srvyrs6(i)))+square(indxsurv_log_sd6(i))/2.)/indxsurv_log_sd6(i));
 }
 RMSE_srv6=0;
 if(!isretro)
   RMSE_srv6= sqrt(norm2(log(indxsurv6)-log(Eindxsurv6(srvyrs6))+square(indxsurv_log_sd6)/2.)/nyrs_srv6);
   
 // Age composition
  loglik(17)=0;
 for (i=1;i<=nyrsac_srv6;i++) {
   if(srv_acyrs6(i)>endyr) break;
   llsrvp6(i) = 0;
   for (j=ac_yng_srv6(i);j<=ac_old_srv6(i);j++) {
     llsrvp6(i) += multN_srv6(i)*(srvp6(i,j)+o)*log((Esrvp6(srv_acyrs6(i),j)+o)/(srvp6(i,j)+o));
     res_srv6(i,j)=srvp6(i,j);
     res_srv6(i,trmage-rcrage+j+1)=Esrvp6(srv_acyrs6(i),j);
     if(multN_srv6(i)>0) {
       pearson_srv6(i,j)=(srvp6(i,j)-Esrvp6(srv_acyrs6(i),j))/sqrt((Esrvp6(srv_acyrs6(i),j)*(1.-Esrvp6(srv_acyrs6(i),j)))/multN_srv6(i));	
     }  
   }
   if(multN_srv6(i)>0) {	
     effN_srv6(i) = sum(elem_prod(Esrvp6(srv_acyrs6(i)),(1-Esrvp6(srv_acyrs6(i)))))/sum(square(srvp6(i)-Esrvp6(srv_acyrs6(i))));
   }
   loglik(17) +=llsrvp6(i);
 }
// length composition
  for (i=1;i<=nyrslen_srv6;i++) {
    if(srv_lenyrs6(i)>endyr) break;
    llsrvlenp6(i) = 0;
  for (j=1;j<=nbins2;j++) {
    llsrvlenp6(i) += multNlen_srv6(i)*(srvlenp6(i,j)+o)*log((Esrvlenp6(srv_lenyrs6(i),j)+o)/(srvlenp6(i,j)+o));
  }
  loglik(17) += llsrvlenp6(i);
  }

  // Constraints on recruitment. Assumed sigmaR=1.3 for all devs
  loglik(18)= 0;
  loglik(18) += -0.5*norm2(dev_log_recruit/sigmaR);

 // Normal process error on selectivity deviations. Note
 // rwlk_sd(styr,endyr-1) b/c if using retro they will be too
 // long since read in as data with the original endyr value
  loglik(19)  = -0.5*norm2(elem_div(first_difference(slp1_fsh_dev),rwlk_sd(styr,endyr-1))); 
  loglik(19) += -0.5*norm2(elem_div(first_difference(inf1_fsh_dev),4.0*rwlk_sd(styr,endyr-1)));
  loglik(19) += -0.5*norm2(elem_div(first_difference(slp2_fsh_dev),rwlk_sd(styr,endyr-1)));
  loglik(19) += -0.5*norm2(elem_div(first_difference(inf2_fsh_dev),rwlk_sd(styr,endyr-1)));
 
 // Recruitment in projection mode, but skip it if doing retro
 // peels b/c it makes no sense and breaks this code
  if(last_phase()) {
    loglik(20) =  -(1/(2.0*sigmasq_recr))*norm2(log_recr_proj - log_mean_recr_proj);
  } else {
    loglik(20)=0;
  }
 // Normal process error on catchability deviations. Note
 // rwlk_sd(styr,endyr-1) b/c if using retro they will be too
 // long since read in as data with the original endyr value
  loglik(21)  = -0.5*norm2(elem_div(first_difference(log_q1_dev),q1_rwlk_sd(styr,endyr-1))); 
  loglik(21)  += -0.5*norm2(elem_div(first_difference(log_q2_dev),q2_rwlk_sd(styr,endyr-1))); 
  loglik(21)  += -0.5*norm2(elem_div(first_difference(log_q3_dev),q3_rwlk_sd(styr,endyr-1))); 
  loglik(22)= 0;

 // Prior on trawl catchability       
 loglik(23) = -.5*square((log_q2_mean-log(0.85))/0.1);
 loglik(24) = 0;
 // these broad priors stabilize estimation, particularly for
 // retros, and imply uniform 1 selex for this survey
 loglik(24) -= dnorm(log_slp2_srv6, 0,2, true);
 loglik(24) -= dnorm(inf2_srv6, 10,3, true);
 objfun = -sum(loglik);

 // Variable to do a likelihood profile over
  var_prof=Esumbio(endyr);


FUNCTION MCMC_output
 if(mceval_phase()) {
   // Dumping data to mceval.dat
   report1<<mean(recruit(1978,endyr-1));
   report1<<" ";
   report1<<endyr-2;
   report1<<" ";
   report1<<F(endyr-2);
   report1<<" ";
   report1<<Espawnbio(endyr-2);
   report1<<" ";
   report1<<recruit(endyr-2);
   report1<<" ";
   
   report1<<endyr-1;
   report1<<" ";
   report1<<F(endyr-1);
   report1<<" ";
   report1<<Espawnbio(endyr-1);
   report1<<" ";
   report1<<recruit(endyr-1);
   report1<<" ";
   
   report1<<endyr;
   report1<<" ";
   report1<<F(endyr);
   report1<<" ";
   report1<<Espawnbio(endyr);
   report1<<" ";
   report1<<recruit(endyr);
   report1<<" ";
   
   report1<<endyr+1;
   report1<<" ";
   report1<<F_proj(endyr+1);
   report1<<" ";
   report1<<Espawnbio_proj(endyr+1);
   report1<<" ";
   report1<<recruit_proj(endyr+1);
   report1<<" ";
   
   report1<<endyr+2;
   report1<<" ";
   report1<<F_proj(endyr+2);
   report1<<" ";
   report1<<Espawnbio_proj(endyr+2);
   report1<<" ";
   report1<<recruit_proj(endyr+2);
   report1<<" ";
   
   report1<<endyr+3;
   report1<<" ";
   report1<<F_proj(endyr+3);
   report1<<" ";
   report1<<Espawnbio_proj(endyr+3);
   report1<<" ";
   report1<<recruit_proj(endyr+3);
   report1<<" ";
   
   report1<<endyr+4;
   report1<<" ";
   report1<<F_proj(endyr+4);
   report1<<" ";
   report1<<Espawnbio_proj(endyr+4);
   report1<<" ";
   report1<<recruit_proj(endyr+4);
   report1<<" ";
   
   report1<<endyr+5;
   report1<<" ";
   report1<<F_proj(endyr+5);
   report1<<" ";
   report1<<Espawnbio_proj(endyr+5);
   report1<<" ";
   report1<<recruit_proj(endyr+5);
   report1 << endl;
  }


RUNTIME_SECTION

  convergence_criteria 1.e0, 1.e-1, 1.e-4, 1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7
  maximum_function_evaluations 1000, 1000, 1000, 1000

TOP_OF_MAIN_SECTION
 arrmblsize = 3000000;
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);


REPORT_SECTION
  report << "Objective function" << endl << objfun << endl;
  report << "Likelihood components" << endl << loglik << endl;
  report << "Recruits" << endl <<  recruit << endl; 
  report << "Expected spawning biomass" << endl << Espawnbio << endl;
  report << "Expected total biomass" << endl << Etotalbio << endl;
  report << "Expected summary (age 3+) biomass" << endl << Esumbio << endl;
  report << "Expected spawning biomass age 2+" << endl << Espawnbio_2plus << endl;
  report << "Fishing mortalities" << endl <<  F << endl ; 
  report << "Ecosystem comsumption" << endl <<  Eecocon << endl; 
  report << "Initial age comp" << endl <<  initN << endl; 
  report << "Total catch" << endl <<  cattot << endl;  
  report << "Expected total catch" << endl <<  Ecattot << endl;  
  report << "Natural mortality" << endl <<  M << endl; 
  report << "Numbers at age" << endl << N << endl << endl;


  report << "Fishery selectivity means" << endl;  
  report << log_slp1_fsh_mean << " " <<  inf1_fsh_mean << " " << log_slp2_fsh_mean  << " " << inf2_fsh_mean << endl; 
  report << "Fishery selectivity deviances" << endl;  
  report << slp1_fsh_dev << endl; 
  report << inf1_fsh_dev << endl; 
  report << slp2_fsh_dev << endl; 
  report << inf2_fsh_dev << endl; 
  report << "Fishery selectivity vectors" << endl;  
  report <<  slp1_fsh << endl; 
  report <<  inf1_fsh << endl; 
  report <<  slp2_fsh << endl; 
  report <<  inf2_fsh << endl; 
  report << "Fishery selectivity" << endl << slctfsh << endl;
  report << "Fishery age composition likelihoods" << endl << llcatp << endl;
  report << "Fishery age composition" << endl << catp << endl;
  report << "Fishery expected age composition" << endl << Ecatp << endl;
  report << "Fishery observed and expected age comp" << endl << res_fish << endl;
  report << "Fishery Pearson residuals age comp" << endl << pearson_fish << endl; 
  report << "Fishery input N" << endl << multN_fsh << endl;  
  report << "Fishery effective N age comp" << endl << effN_fsh << endl;   
  report << "Fishery length composition likelihoods" << endl << lllenp << endl;
  report << "Fishery length composition" << endl << lenp << endl;
  report << "Fishery expected length composition" << endl << Elenp << endl << endl;

  report << "Survey 1 q" << endl << q1 << endl;
  report << "Survey 1 selectivity parameters" << endl;  
  report << log_slp2_srv1 << " " << inf2_srv1 << endl; 
  report << "Survey 1 selectivity" << endl << slctsrv1 << endl;
  report << "Survey 1 expected index" << endl << Eindxsurv1 << endl;
  report << "Survey 1 RMSE" << endl << RMSE_srv1 << endl;
  report << "Survey 1 age composition likelihoods" << endl << llsrvp1 << endl;
  report << "Survey 1 age composition" << endl << srvp1 << endl;
  report << "Survey 1 expected age composition" << endl << Esrvp1 << endl;
  report << "Survey 1 observed and expected age comp" << endl << res_srv1 << endl;
  report << "Survey 1 Pearson residuals age comp" << endl << pearson_srv1 << endl; 
  report << "Survey 1 input N" << endl << multN_srv1 << endl;  
  report << "Survey 1 effective N age comp" << endl << effN_srv1 << endl;   
  report << "Survey 1 length composition likelihoods" << endl << llsrvlenp1 << endl;
  report << "Survey 1 length composition" << endl << srvlenp1 << endl;
  report << "Survey 1 expected length composition" << endl << Esrvlenp1 << endl << endl;

  report << "Survey 2 q" << endl << q2 << endl;
  report << "Survey 2 selectivity parameters" << endl;  
  report << log_slp1_srv2 << " " << inf1_srv2 << " " << log_slp2_srv2 << " " << inf2_srv2 << endl; 
  report << "Survey 2 selectivity" << endl << slctsrv2 << endl;
  report << "Survey 2 expected index" << endl << Eindxsurv2 << endl;
  report << "Survey 2 RMSE" << endl << RMSE_srv2 << endl;
  report << "Survey 2 age composition likelihoods" << endl << llsrvp2 << endl;
  report << "Survey 2 age composition" << endl << srvp2 << endl;
  report << "Survey 2 expected age composition" << endl << Esrvp2 << endl;
  report << "Survey 2 observed and expected age comp" << endl << res_srv2 << endl;
  report << "Survey 2 Pearson residuals age comp" << endl << pearson_srv2 << endl;  
  report << "Survey 2 input N" << endl << multN_srv2 << endl;  
  report << "Survey 2 effective N age comp" << endl << effN_srv2 << endl;   
  report << "Survey 2 length composition likelihoods" << endl << llsrvlenp2 << endl;
  report << "Survey 2 length composition" << endl << srvlenp2 << endl;
  report << "Survey 2 expected length composition" << endl << Esrvlenp2 << endl << endl;

  report << "Survey 3 q" << endl << q3 << endl;
  report << "Survey 3 selectivity parameters" << endl;  
  report << log_slp1_srv3 << " " <<  inf1_srv3 << endl; 
  report << "Survey 3 selectivity" << endl << slctsrv3 << endl;
  report << "Survey 3 expected index" << endl << Eindxsurv3 << endl;
  report << "Survey 3 RMSE" << endl << RMSE_srv3 << endl;
  report << "Survey 3 age composition likelihoods" << endl << llsrvp3 << endl;
  report << "Survey 3 age composition" << endl << srvp3 << endl;
  report << "Survey 3 expected age composition" << endl << Esrvp3 << endl;
  report << "Survey 3 observed and expected age comp" << endl << res_srv3 << endl;
  report << "Survey 3 Pearson residuals age comp" << endl << pearson_srv3 << endl; 
  report << "Survey 3 input N" << endl << multN_srv3 << endl;
  report << "Survey 3 effective N age comp" << endl << effN_srv3 << endl;   
  report << "Survey 3 length composition likelihoods" << endl << llsrvlenp3 << endl;
  report << "Survey 3 length composition" << endl << srvlenp3 << endl;
  report << "Survey 3 expected length composition" << endl << Esrvlenp3 << endl;
  report << "Survey 3 observed and expected length comp" << endl << res_srv3len << endl;
  report << "Survey 3 Pearson residuals length comp" << endl << pearson_srv3len << endl << endl;
  
  report << "Survey 4 q" << endl << q4 << endl;
  report << "Survey 4 power" << endl << q4_pow << endl; 
  report << "Survey 4 expected index" << endl << Eindxsurv4 << endl;
  report << "Survey 4 RMSE" << endl << RMSE_srv4 << endl << endl;
  report << "Survey 5 q" << endl << q5 << endl;
  report << "Survey 5 power" << endl << q5_pow << endl; 
  report << "Survey 5 expected index" << endl << Eindxsurv5 << endl;
  report << "Survey 5 RMSE" << endl << RMSE_srv5 << endl << endl;
   
  report << "Survey 6 q" << endl << q6 << endl;
  report << "Survey 6 selectivity parameters" << endl;  
  report << log_slp1_srv6 << " " << inf1_srv6 << " " << log_slp2_srv6 << " " << inf2_srv6 << endl; 
  report << "Survey 6 selectivity" << endl << slctsrv6 << endl;
  report << "Survey 6 expected index" << endl << Eindxsurv6 << endl;
  report << "Survey 6 RMSE" << endl << RMSE_srv6 << endl;
  report << "Survey 6 age composition likelihoods" << endl << llsrvp6 << endl;
  report << "Survey 6 age composition" << endl << srvp6 << endl;
  report << "Survey 6 expected age composition" << endl << Esrvp6 << endl;
  report << "Survey 6 observed and expected age comp" << endl << res_srv6 << endl;
  report << "Survey 6 Pearson residuals age comp" << endl << pearson_srv6 << endl;  
  report << "Survey 6 input N" << endl << multN_srv6 << endl;  
  report << "Survey 6 effective N age comp" << endl << effN_srv6 << endl;   
  report << "Survey 6 length composition likelihoods" << endl << llsrvlenp6 << endl;
  report << "Survey 6 length composition" << endl << srvlenp6 << endl;
  report << "Survey 6 expected length composition" << endl << Esrvlenp6 << endl << endl;

  report << "Projection selectivity" << endl << slctfsh_proj << endl;
  report << "Projection recruits" << endl << recruit_proj << endl;
  report << "Projection log mean recruitment" << endl << log_recr_proj << endl;
  report << "Projection recruitment variability" << endl <<  sigmasq_recr<< endl;
  report << "Projection numbers at age" << endl << N_proj << endl;
  report << "Projection total mortality" << endl << Z_proj << endl;	  
  report << "Projection catch at age" << endl << C_proj << endl; 
  report << "Projection survey numbers at age" << endl << Nsrv_proj << endl;
  report << "Projection population weight at age" << endl <<  wt_pop_proj << endl;  
  report << "Projection spawning weight at age" << endl << wt_spawn_proj << endl;  
  report << "Projection fishery weight at age" << endl << wt_fsh_proj << endl; 
  report << "Projection fishery selectivity" << endl << slctfsh_proj << endl; 
  report << "Projection total catches" << endl << Ecattot_proj << endl;  
  report << "Projection summary biomass" << endl << Esumbio_proj << endl;  
  report << "Projection spawning biomass" << endl << Espawnbio_proj << endl;  
  report << "Projection survey biomass" << endl << Esrv_proj << endl;  
  report << "Projection Ftarget" << endl << Ftarget << endl;
  report << "Projection B40" << endl <<  B40 << endl;  
  report << "Projection fishing mortality" << endl <<  F_proj << endl;  
  

