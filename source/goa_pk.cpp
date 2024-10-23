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

// Added log biomass sdreport vectors. Will use those for
// uncertainties moving forward.

// Cleaned out old inputs  9/22.

// Model 19.1a: add sigmaR=1.3 for all devs, turn on descending
// logistic selectivity for survey 6, including some broad priors
// to stabilize estimation. Also add more sdreport variables to
// output.

// July 2023 started conversion to TMB
// Oct 2023 took out retro option (done in R)

// Jan 2024 added internal projection module back in. Now has
// dynamic length depending on length of Ftarget input

// Feb 2024 switched to dnorm and dmultinom and added SIMULATE

// Oct 2024 add Ecov link to q1 (Rogers et al. 2024)



#include <TMB.hpp>
#include <iostream>
//#include "helper_functions.hpp"

template <class Type>
Type square(Type x){return x*x;};
template <class Type>
Type norm2(vector<Type> x){return (x*x).sum();};
// from TMB examples, modified name, constrains between -1 and 1
template <class Type>
Type rho_trans(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}
#define see(object) std::cout << #object ":\n" << object << "\n";
// rmultinom taken from WHAM on 2024-03-22
// https://github.com/timjmiller/wham/blob/master/src/age_comp_sim.hpp#L1C1-L19C2
template<class Type>
vector<Type> rmultinom(Type N, vector<Type> p)
{
  //multinomial
  int dim = p.size();
  vector<Type> x(dim);
  int Nint = CppAD::Integer(N);
  x.setZero();
  for(int i = 0; i < Nint; i++)
  {
    Type y = runif(0.0,1.0);
    for(int a = 0; a < dim; a++) if(y < p.head(a+1).sum())
    {
      x(a) += 1.0;
      break;
    }
  }
  return x;
}

template <class Type>
vector<Type> rdirichlet(vector<Type> alpha){
  vector<Type> x=rgamma(alpha,Type(1));
  return x/sum(x);
}

template<class Type>
vector<Type> rdirmultinom(Type N, vector<Type> alpha, vector<int> ages) //dirichlet generated from iid gammas
{
  vector<Type> dp = rdirichlet(alpha);
  vector<Type> obs = rmultinom(N,dp);
  return(obs);
}


template<class Type>
Type ddirmultinom(vector<Type> obs, vector<Type> alpha, int do_log)
{
  int dim = obs.size();
  Type N = obs.sum();
  Type phi=sum(alpha);
  Type ll = lgamma(N + 1.0) + lgamma(phi) - lgamma(N + phi);
  for(int a = 0; a < dim; a++) ll += -lgamma(obs(a) + 1.0) + lgamma(obs(a) + alpha(a)) - lgamma(alpha(a));
  if(do_log == 1) return ll;
  else return exp(ll);
}
template<class Type>
Type objective_function<Type>::operator() ()
{


  // !!CLASS ofstream report1("mceval.dat")
  DATA_INTEGER(styr);	    // Starting year for population model
  DATA_INTEGER(endyr);	    // Ending year for population model
  DATA_INTEGER(rcrage);	    // Recruitment age
  DATA_INTEGER(trmage);	    // Last modeled age
  DATA_INTEGER(nbins1);	// Number of length bins in transitiom matrix 1
  DATA_INTEGER(nbins2);	// Number of length bins in transitiom matrix 2
  DATA_INTEGER(nbins3);	// Number of length bins in transitiom matrix 3

  // Fishery
  DATA_VECTOR(cattot);			// Total catch in tons
  DATA_VECTOR(cattot_log_sd); // Total catch (cv) = sdev of log(cattot)
  DATA_INTEGER(nyrs_fsh);     // Number of fishery age comps
  DATA_IVECTOR(fshyrs);       // Years for the fishery age comps
  DATA_VECTOR(multN_fsh);       // Multinomial sample size by year
  DATA_IVECTOR(ac_yng_fsh);	// Accumulation of lower ages
  DATA_IVECTOR(ac_old_fsh);	// Accumulation of upper ages
  DATA_INTEGER(nyrslen_fsh);	// Number of fishery length comps
  DATA_IVECTOR(fshlenyrs);          // Years for the fishery length comps
  DATA_VECTOR(multNlen_fsh); // Multinomial sample size by year
  DATA_VECTOR(rwlk_sd); // Random walk stdevs
  DATA_MATRIX(catp); // Catch proportions at age
  DATA_MATRIX(lenp); // Catch proportions at age
  DATA_MATRIX(wt_fsh); // Weight at age by year
  // For matrices indices for rows, then col
  // Survey 1 (Acoustic) EK500 (biosonic deleted in 2022)
  DATA_INTEGER(nyrs_srv1); // Number of survey biomass estimates
  DATA_IVECTOR(srvyrs1);   // Years in which surveys occured
  DATA_VECTOR(indxsurv1);  // Survey index
  DATA_VECTOR(indxsurv_log_sd1); // Survey index (cv) = sdev of log(indxsurv)
  DATA_VECTOR(q1_rwlk_sd);	 // Random walk stdevs
  DATA_VECTOR(yrfrct_srv1); // Fraction of year to midpoint of survey
  DATA_INTEGER(nyrsac_srv1);	// Number of survey age comps
  DATA_IVECTOR(srv_acyrs1);	// Years for the survey age comp
  DATA_VECTOR(multN_srv1);     // Multinomial sample size by year
  DATA_IVECTOR(ac_yng_srv1);   // Accumulation of lower ages
  DATA_IVECTOR(ac_old_srv1);   // Accumulation of upper ages
  DATA_INTEGER(nyrslen_srv1);	// Number of survey length comps
  DATA_IVECTOR(srv_lenyrs1); // Years for the survey length comps
  DATA_VECTOR(multNlen_srv1);	// Multinomial sample size by year
  DATA_MATRIX(srvp1);		// Survey proportions at age
  DATA_MATRIX(srvlenp1);	// Survey proportions at length
  DATA_MATRIX(wt_srv1);	// Survey weights at age

  // Survey 2 (Bottom trawl)
  DATA_INTEGER(nyrs_srv2);	    // Number of surveys
  DATA_IVECTOR(srvyrs2);	// Years in which surveys occured
  DATA_VECTOR(indxsurv2);	// Survey index
  DATA_VECTOR(indxsurv_log_sd2); // Survey index (cv) = sdev of log(indxsurv)
  DATA_VECTOR(q2_rwlk_sd);	 // Random walk stdevs
  DATA_VECTOR(yrfrct_srv2); // Fraction of year to midpoint of survey
  DATA_INTEGER(nyrsac_srv2);	// Number of survey age comps
  DATA_IVECTOR(srv_acyrs2);	// Years for the survey age comp
  DATA_VECTOR(multN_srv2);     // Multinomial sample size by year
  DATA_IVECTOR(ac_yng_srv2);   // Accumulation of lower ages
  DATA_IVECTOR(ac_old_srv2);   // Accumulation of upper ages
  DATA_INTEGER(nyrslen_srv2);	// Number of survey length comps
  DATA_IVECTOR(srv_lenyrs2); // Years for the survey length comps
  DATA_VECTOR(multNlen_srv2);  // Multinomial sample size by year
  DATA_MATRIX(srvp2);	       // Survey proportions at age
  DATA_MATRIX(srvlenp2);       // Survey proportions at length
  DATA_MATRIX(wt_srv2);	       // Survey weights at age

  // Survey 3 (ADFG coastal survey)
  DATA_INTEGER(nyrs_srv3); // Number of survey biomass estimates
  DATA_IVECTOR(srvyrs3);   // Years in which surveys occured
  DATA_VECTOR(indxsurv3);  // Survey index
  DATA_VECTOR(indxsurv_log_sd3); // Survey index (cv) = sdev of log(indxsurv)
  DATA_VECTOR(q3_rwlk_sd);	 // Random walk stdevs
  DATA_VECTOR(yrfrct_srv3); // Fraction of year to midpoint of survey
  DATA_INTEGER(nyrsac_srv3);	// Number of survey age comps
  DATA_IVECTOR(srv_acyrs3);	// Years for the survey age comps
  DATA_VECTOR(multN_srv3);     // Multinomial sample size by year
  DATA_INTEGER(nyrslen_srv3);	// Number of survey length comps
  DATA_IVECTOR(srv_lenyrs3); // Years for the survey length comps
  DATA_VECTOR(multNlen_srv3);  // Multinomial sample size by year
  DATA_MATRIX(srvp3);	       // Survey proportions at age
  DATA_MATRIX(srvlenp3);       // Survey proportions at length
  DATA_MATRIX(wt_srv3);	       // Survey weights at age

  // Survey 4 (Age 1 acoustic)
  DATA_INTEGER(nyrs_srv4);	    // Number of surveys
  DATA_IVECTOR(srvyrs4);	// Years in which surveys occured
  DATA_VECTOR(indxsurv4);	// Survey index
  DATA_VECTOR(indxsurv_log_sd4); // Survey index (cv) = sdev of log(indxsurv)

  // Survey 5 (Age 2 acoustic)
  // DATA_INTEGER(nyrs_srv5);	    // Number of surveys
  DATA_IVECTOR(srvyrs5);	// Years in which surveys occured
  DATA_VECTOR(indxsurv5);	// Survey index
  DATA_VECTOR(indxsurv_log_sd5); // Survey index (cv) = sdev of log(indxsurv)

  // Survey 6 (Summer acoustic)
   DATA_INTEGER(nyrs_srv6);	    // Number of surveys
  DATA_IVECTOR(srvyrs6);	// Years in which surveys occured
  DATA_VECTOR(indxsurv6);	// Survey index
  DATA_VECTOR(indxsurv_log_sd6); // Survey index (cv) = sdev of log(indxsurv)
  DATA_VECTOR(yrfrct_srv6); // Fraction of year to midpoint of survey
  DATA_INTEGER(nyrsac_srv6);	// Number of survey age comps
  DATA_IVECTOR(srv_acyrs6);	// Years for the survey age comp
  DATA_VECTOR(multN_srv6);   // Multinomial sample size by year
  DATA_IVECTOR(ac_yng_srv6); // Accumulation of lower ages
  DATA_IVECTOR(ac_old_srv6); // Accumulation of upper ages
  DATA_INTEGER(nyrslen_srv6);	// Number of survey length comps
  DATA_IVECTOR(srv_lenyrs6); // Years for the survey length comps
  DATA_VECTOR(multNlen_srv6);	// Multinomial sample size by year
  DATA_MATRIX(srvp6);		// Survey proportions at age
  DATA_MATRIX(srvlenp6);	// Survey proportions at length
  DATA_MATRIX(wt_srv6);	// Survey weights at age

  // Ageing error transition matrix
  DATA_MATRIX(age_trans);

  // Age to length transition matrix, to calculate expected length comps
  DATA_MATRIX(len_trans1);
  DATA_MATRIX(len_trans2);
  DATA_MATRIX(len_trans3);

  int nyrs=endyr-styr+1;
  int nages=trmage-rcrage+1;
  vector<int> ages(nages);
  vector<int> years(nyrs);
  for(int i=0;i<nages;i++) ages(i)=i+1;
  for(int i=0;i<nyrs;i++) years(i)=1970+i;

  // Population vectors
  matrix<Type> wt_pop(nyrs,nages);   // Population weight at age
  matrix<Type> wt_spawn(nyrs,nages); // Population weight at age at spawning (April 15)
  wt_pop=wt_srv2;
  wt_spawn=wt_srv1;
  DATA_VECTOR(mat);		// Proportion mature

  DATA_VECTOR(Ftarget);		// determines length of projections
  int nyears_proj = Ftarget.size();


  DATA_SCALAR(B40); // mean log recruitment
  // DATA_SCALAR(log_mean_recr_proj);
  DATA_SCALAR(sigmasq_recr);                       // Variance for log recr, recruitment indices
  DATA_IVECTOR(Ecov_obs_year);
  DATA_VECTOR(Ecov_obs);
  DATA_INTEGER(complike); // 1 is multinomial, 2 is D-M
  bool isDM=false;
  if(complike==2) isDM=true;
  REPORT(isDM);
  // int styr_avg_slct;
  // int endyr_avg_slct;
  int i;                                          // Index for year
  int j;                                          // Index for age
  //int loop=0;

  Type o=0.00001;                                       // A small number

   // Projection parameters
  vector<Type> wt_pop_proj(nages);
  vector<Type> wt_spawn_proj(nages);
  vector<Type> wt_fsh_proj(nages);
  vector<Type> wt_srv_proj(nages);

  // for projections take averages of WAA only from recent survey years with data
  wt_pop_proj.setZero();
  wt_spawn_proj.setZero();
  wt_fsh_proj.setZero();
  // // have to be careful with retro option or it will grab data
  // // that shouldn't exist, walk backward through the survey
  // // data only counting years included in the likelihood and
  // // accumulate a mean
  // int counter=0;
  // for(int i=nyrsac_srv1; i>=1; i--){
  //   if(srv_acyrs1(i) <= endyr){
  //     counter++;
  //     for(int a=rcrage;a<=trmage;a++) wt_spawn_proj(a)+=wt_srv1(srv_acyrs1(i),a)/5;
  //     if(counter==5) break;
  //   }
  // }
  // counter=0;
  // for(int i=nyrsac_srv2; i>=1; i--){
  //   if(srv_acyrs2(i) <= endyr){
  //     counter++;
  //     for(int a=rcrage;a<=trmage;a++) wt_pop_proj(a)+=wt_srv2(srv_acyrs2(i),a)/3;
  //     if(counter==3) break;
  //   }
  // }
  // for predicting what the survey will see next year, Shelikof for now
  //  wt_srv_proj=wt_spawn_proj;
  // wt_fsh_proj=wt_fsh(nyrs);

  // Population parameters
  //  init_bounded_number M(0.1,0.5,-1)
  // PARAMETER_VECTOR(M);
  // PARAMETER(mean_log_initN);
  PARAMETER_VECTOR(dev_log_initN);
  vector<Type> initN(nages-1); // goes from 2 to 10
  PARAMETER(mean_log_recruit);
  PARAMETER_VECTOR(dev_log_recruit);
  PARAMETER(sigmaR);
  vector<Type> recruit(nyrs);
  vector<Type> log_recruit(nyrs);

  // Forward projections
  PARAMETER_VECTOR(log_recr_proj);
  vector<Type> recruit_proj(nyears_proj);
  matrix<Type> N_proj(nyears_proj,nages);
  vector<Type> F_proj(nyears_proj);
  matrix<Type> Z_proj(nyears_proj,nages);
  matrix<Type> C_proj(nyears_proj,nages);
  matrix<Type> Nsrv_proj(nyears_proj,nages);
  vector<Type> slctfsh_proj(nages);
  vector<Type> Ecattot_proj(nyears_proj);
  vector<Type> Esumbio_proj(nyears_proj);
  vector<Type> Espawnbio_proj(nyears_proj);
  vector<Type> Esrv_proj(nyears_proj);
  vector<Type> Exrate_proj(nyears_proj);

  // Selectivity parameters

  // Fishery selectivity
  PARAMETER(log_slp1_fsh_mean);
  PARAMETER(inf1_fsh_mean);
  PARAMETER(log_slp2_fsh_mean);
  PARAMETER(inf2_fsh_mean);
  //    PARAMETER(log_slp2_fsh_mean);
  //    PARAMETER(inf2_fsh_mean);
  PARAMETER_VECTOR(slp1_fsh_dev);
  PARAMETER_VECTOR(inf1_fsh_dev);
  //    PARAMETER_VECTOR(slp2_fsh_dev);
  //    PARAMETER_VECTOR(inf2_fsh_dev);
  PARAMETER_VECTOR(slp2_fsh_dev);
  PARAMETER_VECTOR(inf2_fsh_dev);

  vector<Type> slp1_fsh(nyrs);
  vector<Type> inf1_fsh(nyrs);
  vector<Type> slp2_fsh(nyrs);
  vector<Type> inf2_fsh(nyrs);

  // Acoustic survey selectivity
  //  PARAMETER(log_slp1_srv1);
  //  PARAMETER(inf1_srv1);
  PARAMETER(log_slp2_srv1);
  PARAMETER(inf2_srv1);
  // Trawl selectivity
  PARAMETER(log_slp1_srv2);
  PARAMETER(inf1_srv2);
  //    PARAMETER inf1_srv2(1,50,-1)
  PARAMETER(log_slp2_srv2);
  PARAMETER(inf2_srv2);
  // ADFG selectivity
  PARAMETER(log_slp1_srv3);
  PARAMETER(inf1_srv3);
  //    PARAMETER(log_slp2_srv3);
  //    PARAMETER(inf2_srv3);
  // Summer acoustic selectivity
  PARAMETER(log_slp1_srv6);
  PARAMETER(inf1_srv6);
  //    PARAMETER(inf1_srv6);
  //    PARAMETER(inf1_srv6);
  PARAMETER(log_slp2_srv6);
  PARAMETER(inf2_srv6);

  // Fishing mortality and survey catchablility
  PARAMETER(mean_log_F);
  PARAMETER_VECTOR(dev_log_F);
  vector<Type> F(nyrs);
  PARAMETER(log_q1_mean);
  PARAMETER_VECTOR(log_q1_dev);;
  //  PARAMETER(log_q2);
  PARAMETER(log_q2_mean);
  PARAMETER_VECTOR(log_q2_dev);;
  PARAMETER(log_q3_mean);
  PARAMETER_VECTOR(log_q3_dev);;
  PARAMETER(log_q4);
  //  PARAMETER(q4_pow);
  PARAMETER(q4_pow);
  PARAMETER(log_q5);
  //  PARAMETER(q5_pow);;
  PARAMETER(q5_pow);
  PARAMETER(log_q6);
  // This scales M vector below so that M={M}*natMscalar. If 1 does nothing.
  PARAMETER(natMscalar);
  // Ecov effects on q1
  PARAMETER(transf_rho);
  PARAMETER_VECTOR(Ecov_exp);
  PARAMETER(log_Ecov_obs_sd);
  PARAMETER(log_Ecov_sd);
  PARAMETER(Ecov_beta);
  // linear Dirichlet-multinomial parameters for age comps (fsh,
  // survs 1, 2, 3, 6), NOT length comps. Simulation built in but
  // not yet tested thoroughly
  PARAMETER_VECTOR(log_DM_pars);
  vector<Type> DM_pars = exp(log_DM_pars);
  vector<Type> ESS_fsh=multN_fsh*invlogit(log_DM_pars(0));
  REPORT(ESS_fsh);
  vector<Type> ESS_srv1=multN_srv1*invlogit(log_DM_pars(1));
  REPORT(ESS_srv1);
  vector<Type> ESS_srv2=multN_srv2*invlogit(log_DM_pars(2));
  REPORT(ESS_srv2);
  vector<Type> ESS_srv3=multN_srv3*invlogit(log_DM_pars(3));
  REPORT(ESS_srv3);
  vector<Type> ESS_srv6=multN_srv6*invlogit(log_DM_pars(4));
  REPORT(ESS_srv6);

  // Dependent parameters
  vector<Type> q1(nyrs);
  vector<Type> q2(nyrs);
  vector<Type> q3(nyrs);
  vector<Type> log_q1(nyrs);
  vector<Type> log_q2(nyrs);
  vector<Type> log_q3(nyrs);
  Type q4;
  Type q5;
  Type q6;
  matrix<Type> N(nyrs,nages);
  vector<Type> endN(nages);
  matrix<Type> Z(nyrs,nages);
  matrix<Type> C(nyrs,nages);
  matrix<Type> Nsrv1(nyrs,nages);
  vector<Type> slctsrv1(nages);
  matrix<Type> Nsrv2(nyrs,nages);
  vector<Type> slctsrv2(nages);
  matrix<Type> Nsrv3(nyrs,nages);
  vector<Type> slctsrv3(nages);
  matrix<Type> Nsrv6(nyrs,nages);
  vector<Type> slctsrv6(nages);
  matrix<Type> slctfsh(nyrs,nages);
  vector<Type> slctsrv1_logit(nages);
  vector<Type> slctsrv2_logit(nages);
  vector<Type> slctsrv3_logit(nages);
  vector<Type> slctsrv6_logit(nages);
  vector<Type> slctfsh_logit(nages);


  // Expected values
  vector<Type> Eecocon(nyrs);
  matrix<Type> Eec(nyrs,nages);
  vector<Type> Ecattot(nyrs);
  matrix<Type> Ecatp(nyrs,nages);
  matrix<Type> Elenp(nyrs,nbins1);
  vector<Type> Eindxsurv1(nyrs);
  matrix<Type> Esrvp1(nyrs,nages);
  matrix<Type> Esrvlenp1(nyrs,nbins3);
  vector<Type> Eindxsurv2(nyrs);
  matrix<Type> Esrvp2(nyrs,nages);
  matrix<Type> Esrvlenp2(nyrs,nbins2);
  vector<Type> Eindxsurv3(nyrs);
  matrix<Type> Esrvp3(nyrs,nages);
  matrix<Type> Esrvlenp3(nyrs,nbins2);
  vector<Type> Eindxsurv4(nyrs);
  vector<Type> Eindxsurv5(nyrs);
  vector<Type> Eindxsurv6(nyrs);
  matrix<Type> Esrvp6(nyrs,nages);
  matrix<Type> Esrvlenp6(nyrs,nbins2);
  vector<Type> Espawnbio(nyrs);
  vector<Type> Esumbio(nyrs);
  // log versions make more sense for asymptotics
  vector<Type> Espawnbio_log(nyrs);
  vector<Type> Esumbio_log(nyrs);
  vector<Type> Espawnbio_2plus(nyrs);
  vector<Type> Etotalbio(nyrs);
  // log likelihood containers
  vector<Type> loglik(25);
  vector<Type> llcatp(nyrs_fsh);
  vector<Type> lllenp(nyrslen_fsh);
  vector<Type> llsrvp1(nyrsac_srv1);
  vector<Type> llsrvlenp1(nyrslen_srv1);
  vector<Type> llsrvp2(nyrsac_srv2);
  vector<Type> llsrvlenp2(nyrslen_srv2);
  vector<Type> llsrvp3(nyrsac_srv3);
  vector<Type> llsrvlenp3(nyrslen_srv3);
  vector<Type> llsrvp6(nyrsac_srv6);
  vector<Type> llsrvlenp6(nyrslen_srv6);


  //residual output matrices
  matrix<Type> res_fish(nyrs_fsh,2*nages);
  matrix<Type> res_srv1(nyrsac_srv1,2*nages);
  matrix<Type> res_srv2(nyrsac_srv2,2*nages);
  matrix<Type> res_srv3(nyrsac_srv3,2*nages);
  matrix<Type> res_srv3len(nyrslen_srv3,2*nbins2);
  matrix<Type> res_srv6(nyrsac_srv6,2*nages);
  matrix<Type> pearson_fish(nyrs_fsh,nages);
  matrix<Type> pearson_srv1(nyrsac_srv1,nages);
  matrix<Type> pearson_srv2(nyrsac_srv2,nages);
  matrix<Type> pearson_srv3(nyrsac_srv3,nages);
  matrix<Type> pearson_srv3len(nyrslen_srv3,nbins2);
  matrix<Type> pearson_srv6(nyrsac_srv6,nages);
  pearson_fish.setZero();
  pearson_srv1.setZero();
  pearson_srv2.setZero();
  pearson_srv3.setZero();
  pearson_srv3len.setZero();
  pearson_srv6.setZero();
  res_fish.setZero();
  res_srv1.setZero();
  res_srv2.setZero();
  res_srv3.setZero();
  res_srv3len.setZero();
  res_srv6.setZero();

  // These are old, probably used for Francis tuning but not
  // anymore. Left blank for now but consider taking out later
  vector<Type> effN_fsh(nyrs_fsh);
  vector<Type> effN_srv1(nyrsac_srv1);
  vector<Type> effN_srv2(nyrsac_srv2);
  vector<Type> effN_srv3(nyrsac_srv3);
  vector<Type> effN_srv6(nyrsac_srv6);
  effN_fsh.setZero();
  effN_srv1.setZero();
  effN_srv2.setZero();
  effN_srv3.setZero();
  effN_srv6.setZero();

  Type RMSE_srv1;
  Type RMSE_srv2;
  Type RMSE_srv3;
  Type RMSE_srv4;
  Type RMSE_srv5;
  Type RMSE_srv6;

  // Objective function
  //likeprof_number var_prof;
  Type objfun=0.0;


  // PROCEDURE_SECTION

  //   Convert_log_parameters();
  //   Selectivity();
  //   Mortality();
  //   Numbers_at_age();
  //   Catch_at_age();
  //   Expected_values();
  //   if(last_phase())
  //   {
  //     Projections();
  //   }
  //   Objective_function();
  //   MCMC_output();

  // Make some new vectors to help simplify indexing from 0
  int y0=0;
  int y1=nyrs-1;
  int a0=0;
  int a1=nages-1;
  vector<int> ifshyrs=fshyrs-styr;
  vector<int> isrvyrs1=srvyrs1-styr;
  vector<int> isrvyrs2=srvyrs2-styr;
  vector<int> isrvyrs3=srvyrs3-styr;
  vector<int> isrvyrs4=srvyrs4-styr;
  vector<int> isrvyrs5=srvyrs5-styr;
  vector<int> isrvyrs6=srvyrs6-styr;
  vector<int> isrv_acyrs1=srv_acyrs1-styr;
  vector<int> isrv_acyrs2=srv_acyrs2-styr;
  vector<int> isrv_acyrs3=srv_acyrs3-styr;
  // vector<int> isrv_acyrs4=srv_acyrs4-styr;
  // vector<int> isrv_acyrs5=srv_acyrs5-styr;
  vector<int> isrv_acyrs6=srv_acyrs6-styr;
  vector<int> isrv_lenyrs1=srv_lenyrs1-styr;
  vector<int> isrv_lenyrs2=srv_lenyrs2-styr;
  vector<int> isrv_lenyrs3=srv_lenyrs3-styr;
  vector<int> isrv_lenyrs6=srv_lenyrs6-styr;
  vector<int> iac_yng_fsh=ac_yng_fsh-1;
  vector<int> iac_old_fsh=ac_old_fsh-1;
  vector<int> iac_yng_srv1=ac_yng_srv1-1;
  vector<int> iac_old_srv1=ac_old_srv1-1;

  // for projections take averages of WAA only from recent survey years with data
  wt_pop_proj.setZero();
  wt_spawn_proj.setZero();
  wt_srv_proj.setZero();
  for(int a=a0;a<=a1;a++){
    for(int i=0;i<=2;i++)
     wt_pop_proj(a)+=wt_srv2(isrv_acyrs2(nyrsac_srv2-i-1),a)/3;
    for(int i=0;i<=4;i++){
      wt_spawn_proj(a)+=wt_srv1(isrv_acyrs1(nyrsac_srv1-i-1),a)/5;
      // for predicting what the survey will see next year, Shelikof for now
      wt_srv_proj(a)+=wt_srv1(isrv_acyrs1(nyrsac_srv1-i-1),a)/5;
    }
  }
  wt_fsh_proj=wt_fsh.row(y1);

  vector<Type> M(nages);
  M(0)=1.39; M(1)=0.69; M(2)=0.48; M(3)=0.37; M(4)=0.34;
  M(5)=0.30; M(6)=0.30; M(7)=0.29; M(8)=0.28; M(9)=0.29;
  for(j=a0;j<=a1;j++) M(j)*=natMscalar;


  // Assume close to equilibrium at f=0 at the start Use the
  // initial recruitment dev to fill out the initial age
  // composition. Need to be really careful with indexing here
  // since initN is of length 9 and represents ages 2 to 10 but
  // has limits (0,8) unlike the other vectors
  initN.setZero();
    initN(0) = exp(mean_log_recruit +  dev_log_recruit(y0) - M(a0)); // age 2
  for (j=1;j<=8;j++) {
    // j is a different age between initN and M
    initN(j) = initN(j-1)*exp(-M(j+1));
  }
  initN(8) /= (1.0 - exp(-M(a1)));
  //devs for initN are turned off
  for (j=0;j<=8;j++) {
    initN(j) = initN(j)*exp(dev_log_initN(j));
  }
  log_recruit=mean_log_recruit+dev_log_recruit;
  recruit = exp(log_recruit);

  // Acoustic survey random walk setup
  for (i=y0;i<=y1;i++) {
    log_q1(i)=log_q1_mean+log_q1_dev(i)+Ecov_beta*Ecov_exp(i);
    log_q2(i)=log_q2_mean+log_q2_dev(i);
    log_q3(i)=log_q3_mean+log_q3_dev(i);
    q1(i)=exp(log_q1(i));
    q2(i)=exp(log_q2(i));
    q3(i)=exp(log_q3(i));
  }
  F = exp(mean_log_F + dev_log_F);
  q4 = exp(log_q4);
  q5 = exp(log_q5);
  q6 = exp(log_q6);


  // Fishery selectivity
  for (i=y0;i<=y1;i++) {
    slp1_fsh(i)=exp(log_slp1_fsh_mean+slp1_fsh_dev(i));
    inf1_fsh(i)=inf1_fsh_mean+inf1_fsh_dev(i);
    slp2_fsh(i)=exp(log_slp2_fsh_mean+slp2_fsh_dev(i));
    inf2_fsh(i)=inf2_fsh_mean+inf2_fsh_dev(i);
    for (j=a0;j<=a1;j++) {
      slctfsh(i,j) = (1/(1+exp(-(slp1_fsh(i))*(double(j+1)-(inf1_fsh(i))))))*
	(1-1/(1+exp(-(slp2_fsh(i))*(double(j+1)-(inf2_fsh(i))))));
    }
    // The plan would be to check and adjust the max selected age as needed
    slctfsh.row(i)=slctfsh.row(i)/slctfsh(i,6);
  }


  //Survey 1 selectivity
  for (j=a0;j<=a1;j++) {
    slctsrv1(j) = (1-1/(1+exp(-exp(log_slp2_srv1)*(double(j+1)-inf2_srv1))));
  }
  slctsrv1=slctsrv1/slctsrv1(2);
  slctsrv1(a0)=0;
  slctsrv1(a0+1)=0;

  //Survey 2 selectivity
  for (j=a0;j<=a1;j++) {
    slctsrv2(j) = (1/(1+exp(-exp(log_slp1_srv2)*(double(j+1)-inf1_srv2))))
      *(1-1/(1+exp(-exp(log_slp2_srv2)*(double(j+1)-inf2_srv2))));
  }
  slctsrv2=slctsrv2/slctsrv2(a1);

  //Survey 3 selectivity
  for (j=a0;j<=a1;j++) {
    slctsrv3(j) = (1/(1+exp(-exp(log_slp1_srv3)*(double(j+1)-inf1_srv3))));
  }
  slctsrv3=slctsrv3/slctsrv3(a1);

  // Survey 6 selectivity
  for (j=a0;j<=a1;j++) {
    slctsrv6(j) = (1/(1+exp(-exp(log_slp1_srv6)*(double(j+1)-inf1_srv6))))
      *(1-1/(1+exp(-exp(log_slp2_srv6)*(double(j+1)-inf2_srv6))));
  }
  slctsrv6=slctsrv6/slctsrv6(0);

  // Calculate uncertainty in logit space where it makes more
  // sense. But it does break when selex=0 or 1 identically so
  // need to be careful there.
  slctsrv1_logit=-log(1/(slctsrv1-1e-10)-1);
  slctsrv2_logit=-log(1/(slctsrv2-1e-10)-1);
  slctsrv3_logit=-log(1/(slctsrv3-1e-10)-1);
  slctsrv6_logit=-log(1/(slctsrv6-1e-10)-1);
  // last year of fishery selex that has data
  vector<Type> tmp=slctfsh.row(y1-1);
  slctfsh_logit=-log(1/(tmp -1e-10)-1);

  for (i=y0;i<=y1;i++) {
    for (j=a0;j<=a1;j++) {
      Z(i,j)=(F(i)*slctfsh(i,j))+M(j);
    }
  }

  N.setZero();
  for(j=a0+1;j<=a1;j++) N(y0,j)=initN(j-1); // ages 2 to 10 only
  for (i=y0;i<=y1;i++) {
    N(i,a0)=recruit(i);
  }
  for (i=y0;i<y1;i++) {
    for (j=a0;j<a1;j++) {
      N(i+1,j+1)=N(i,j)*exp(-Z(i,j));
    }
    N(i+1,a1)+=N(i,a1)*exp(-Z(i,a1));
  }
  endN=N.row(y1);

  // Catch at age and survey numbers at age
  for (i=y0;i<=y1;i++) {
    for (j=a0;j<=a1;j++) {
      C(i,j)=N(i,j)*((F(i)*slctfsh(i,j))/Z(i,j))*(1-exp(-Z(i,j)));
      Eec(i,j)=N(i,j)*(M(j)/Z(i,j))*(1-exp(-Z(i,j)));
      Nsrv1(i,j)=slctsrv1(j)*N(i,j)*exp(-yrfrct_srv1(i)*Z(i,j));
      Nsrv2(i,j)=slctsrv2(j)*N(i,j)*exp(-yrfrct_srv2(i)*Z(i,j));
      Nsrv3(i,j)=slctsrv3(j)*N(i,j)*exp(-yrfrct_srv3(i)*Z(i,j));
      Nsrv6(i,j)=slctsrv6(j)*N(i,j)*exp(-yrfrct_srv6(i)*Z(i,j));
    }
  }

  Ecattot.setZero(); Eecocon.setZero();
  Esumbio.setZero(); Espawnbio.setZero();
  Etotalbio.setZero(); Espawnbio_2plus.setZero();
  Eindxsurv1.setZero(); Eindxsurv2.setZero();
  Eindxsurv3.setZero(); Eindxsurv6.setZero();
  for (i=y0;i<=y1;i++){
    for(j=a0;j<=a1;j++){
      Ecattot(i) += 1000000*(C(i,j)*wt_fsh(i,j));
      Eecocon(i) += 1000000*(Eec(i,j)*wt_pop(i,j));
      Eindxsurv1(i) += q1(i)*N(i,j)*exp(-yrfrct_srv1(i)*Z(i,j))*slctsrv1(j)*wt_srv1(i,j);
      Eindxsurv2(i) += q2(i)*N(i,j)*exp(-yrfrct_srv2(i)*Z(i,j))*slctsrv2(j)*wt_srv2(i,j);
      Eindxsurv3(i) += q3(i)*N(i,j)*exp(-yrfrct_srv3(i)*Z(i,j))*slctsrv3(j)*wt_srv3(i,j);
      Eindxsurv6(i) += q6*N(i,j)*exp(-yrfrct_srv6(i)*Z(i,j))*slctsrv6(j)*wt_srv6(i,j);
      // 3+ biomass
      if(j>=2) Esumbio(i)+= N(i,j)*wt_pop(i,j);
    // Total biomass
      Etotalbio(i)+= N(i,j)*wt_pop(i,j);
    // 2+ biomass at spawning (use for apportioning to management area)
      if(j>=1)
	Espawnbio_2plus(i)+= N(i,j)*exp(-yrfrct_srv1(i)*Z(i,j))*wt_srv1(i,j);
    // //1+ biomass at spawning
    //     Esumbio(i)= sum(elem_prod(elem_prod(N(i)(a0,a1),exp(-yrfrct_srv1(i)*Z(i)(a0,a1))),wt_srv1(i)(a0,a1)));
      Espawnbio(i)+= N(i,j)*exp(-0.21*Z(i,j))*wt_spawn(i,j)*0.5*mat(j);
    }
    Ecatp.row(i) = (C.row(i)/C.row(i).sum())*age_trans;
    Elenp.row(i) = Ecatp.row(i) * len_trans1;
    Esrvp1.row(i) = (Nsrv1.row(i)/Nsrv1.row(i).sum())*age_trans;
    Esrvlenp1.row(i) = Esrvp1.row(i) * len_trans3;
    Esrvp2.row(i) = (Nsrv2.row(i)/Nsrv2.row(i).sum())*age_trans;
    Esrvlenp2.row(i) = Esrvp2.row(i) * len_trans2;
    Esrvp3.row(i) = (Nsrv3.row(i)/Nsrv3.row(i).sum())*age_trans;
    Esrvlenp3.row(i) = Esrvp3.row(i) * len_trans2;
    Eindxsurv4(i)= q4*pow(N(i,0),(q4_pow+1));
    Eindxsurv5(i)= q5*pow(N(i,1),(q5_pow+1));
    Esrvp6.row(i) = (Nsrv6.row(i)/Nsrv6.row(i).sum())*age_trans;
    Esrvlenp6.row(i) = Esrvp6.row(i) * len_trans2;
  }
  Espawnbio_log=log(Espawnbio);
  Esumbio_log=log(Esumbio);


   // ------------------------------------------------------------
   // Projections

   // Averaging window for selectivity
  slctfsh_proj.setZero();
  int endyr_avg_slct=y1-1;
  int y0_avg_slct=endyr_avg_slct-4;
   for (j=a0;j<=a1;j++) {
     slctfsh_proj(j) = 0;
     for (i=y0_avg_slct;i<=endyr_avg_slct;i++) {
       slctfsh_proj(j) += slctfsh(i,j);
     }
   }
   slctfsh_proj=slctfsh_proj/max(slctfsh_proj);

  N_proj.setZero();
  recruit_proj.setZero();
  Z_proj.setZero();
  C_proj.setZero();
  F_proj.setZero();
  Ecattot_proj.setZero();
  Esumbio_proj.setZero();
  Exrate_proj.setZero();
  Espawnbio_proj.setZero();
  Esrv_proj.setZero();
  Type sbio;

  // recruits from 1978 to previous year, index starting at 0,
  // with the first arg being the index start, and the second the
  // number of years
  vector<Type> rectmp=recruit.segment(1978-styr,endyr-1978);
   for (i=0;i<nyears_proj;i++)  N_proj(i,a0)=rectmp.mean();
   recruit_proj=N_proj.col(0);
  //  Standard projection to get NAA in Jan-1 of first proj year
  //  which is first row of projections, comes from last row of
  //  standard quantities (this year)
   for (j=a0;j<a1;j++)  N_proj(0,j+1)=N(y1,j)*exp(-Z(y1,j));
   N_proj(0,a1)+=N(y1,a1)*exp(-Z(y1,a1)); // plus group

  // Forward projections
  for (i=0; i<nyears_proj; i++) {
    // have to get Z twice
    //  Tuning loop to get the spawning biomass adjustment right
    F_proj(i)=Ftarget(i);
    for (int loop=1;loop<=20;loop++) {
      sbio=0.0;
      for (j=a0;j<=a1;j++) {
	Z_proj(i,j)=F_proj(i)*slctfsh_proj(j)+M(j);
	sbio += N_proj(i,j)*exp(-0.21*Z_proj(i,j))*wt_spawn_proj(j)*0.5*mat(j);
      }
      //  Set the fishing mortality rate, adjusting by the harvest control rule
      F_proj(i)=Ftarget(i);
      if (sbio < B40) 	F_proj(i)=Ftarget(i)*(((sbio/B40)-0.05)/(1-0.05));
    } // end tuning loop so F_proj is right

    // Finalize total mortality
    for (j=a0;j<=a1;j++) Z_proj(i,j)=F_proj(i)*slctfsh_proj(j) + M(j);
    //  Calculate numbers at age given F
    if(i<nyears_proj-1) {
      for (j=a0;j<a1;j++) N_proj(i+1,j+1)=N_proj(i,j)*exp(-Z_proj(i,j));
               N_proj(i+1,a1)+=N_proj(i,a1)*exp(-Z_proj(i,a1));
    }
    // Calculate all derived quantities given F and N
    for (j=a0;j<=a1;j++) {
      C_proj(i,j)=N_proj(i,j)*((F_proj(i)*slctfsh_proj(j))/Z_proj(i,j))*(1-exp(-Z_proj(i,j)));
      Nsrv_proj(i,j)=N_proj(i,j)*exp(-yrfrct_srv6(y1)*Z_proj(i,j));
      Esrv_proj(i)  += q1(y1)*N_proj(i,j)*exp(-yrfrct_srv1(y1)*Z_proj(i,j))*slctsrv1(j)*wt_srv_proj(j);
      Ecattot_proj(i) += 1000000*C_proj(i,j)*wt_fsh_proj(j);
      if(j>=2) Esumbio_proj(i) += N_proj(i,j)*wt_pop_proj(j);
      Espawnbio_proj(i) += N_proj(i,j)*exp(-0.21*Z_proj(i,j))*wt_spawn_proj(j)*0.5*mat(j);
    }
    Exrate_proj(i)=Ecattot_proj(i)/(1000000*Esumbio_proj(i));
  } // end proj years

    // age accumulation (add 1st and 2nd ages) for fishery, turned off for the surveys
    for (i=0;i<nyrs_fsh;i++) {
    for (j=a0;j<=a1;j++) {
      if(j<iac_yng_fsh(i)) {
	Ecatp(ifshyrs(i),iac_yng_fsh(i)) += Ecatp(ifshyrs(i),j);
	Ecatp(ifshyrs(i),j) = 0;
	catp(i,iac_yng_fsh(i)) += catp(i,j);
	catp(i,j) = 0;
      }
    }
  }
    for (i=0;i<nyrs_srv1;i++) {
     for (j=a0;j<=a1;j++) {
      if(j<iac_yng_srv1(i)) {
	Esrvp1(isrvyrs1(i),iac_yng_srv1(i)) += Esrvp1(isrvyrs1(i),j);
	Esrvp1(isrvyrs1(i),j) = 0;
	srvp1(i,iac_yng_srv1(i)) += srvp1(i,j);
	srvp1(i,j) = 0;
      }
    }
  }

    // ------------------------------------------------------------
    // Calculate the likelihood components
    loglik.setZero();

    // Fishery likelihoods
    // Total catch
    loglik(0) = dnorm(log(cattot), log(Ecattot), cattot_log_sd, true).sum();
    SIMULATE {
      cattot=exp(rnorm(log(Ecattot), cattot_log_sd));
      REPORT(cattot);
    }

    //Age composition
    res_fish.setZero();
    pearson_fish.setZero();
    llcatp.setZero();
    for (i=0;i<nyrs_fsh;i++) {
      int nagestmp=1+iac_old_fsh(i)-iac_yng_fsh(i); // .segment takes initial place and vector length
      if(multN_fsh(i)>0) { // this allows turning off contribution easily from dat file
        // Not sure why I have to do this in two steps
        vector<Type> otmp= catp.row(i).segment(iac_yng_fsh(i), nagestmp);
        vector<Type> etmp=Ecatp.row(ifshyrs(i)).segment(iac_yng_fsh(i), nagestmp);
        otmp+=o; otmp*=multN_fsh(i); etmp+=o;
        vector<Type> alphas = sum(otmp) * etmp * DM_pars(0);
        if(isDM){
          llcatp(i) = ddirmultinom(otmp, alphas,  true);
        } else {
          llcatp(i) = dmultinom(otmp, etmp,  true);
        }
        SIMULATE {
          if(isDM){
            otmp=rdirmultinom(multN_fsh(i), alphas, ages);
          } else {
            otmp=rmultinom(multN_fsh(i),etmp);
          }
          catp.row(i).segment(iac_yng_fsh(i), nagestmp)=otmp/otmp.sum();
        }
        // Residuals and outputs which get used by R functions
        for (j=iac_yng_fsh(i);j<=iac_old_fsh(i);j++) {
          res_fish(i,j)=catp(i,j); res_fish(i,a1-a0+j+1)=Ecatp(ifshyrs(i),j);
          pearson_fish(i,j)=(catp(i,j)-Ecatp(ifshyrs(i),j))/sqrt((Ecatp(ifshyrs(i),j)*(1.-Ecatp(ifshyrs(i),j)))/multN_fsh(i));
        }
      }
    }
    REPORT(catp);
    loglik(1) = llcatp.sum();
    // //Length compositions no longer implemented
    lllenp.setZero();
    loglik(2) =0; // unused for now, placeholder for fishery length comps

    // Survey 1 likelihoods
    // Total biomass
    RMSE_srv1=0;
    for(i=0; i<nyrs_srv1;i++){
      // Add bias correction on
      Type Eindxsurv1_tmp=log(Eindxsurv1(isrvyrs1(i)));
      Eindxsurv1_tmp -= square(indxsurv_log_sd1(i))/2.;
      if(indxsurv_log_sd1(i)>0){
        loglik(3)+= dnorm(log(indxsurv1(i)), Eindxsurv1_tmp, indxsurv_log_sd1(i), true);
        SIMULATE  indxsurv1(i)=exp(rnorm(Eindxsurv1_tmp, indxsurv_log_sd1(i)));
      }
    }
    // Age compositions
    res_srv1.setZero(); pearson_srv1.setZero(); llsrvp1.setZero();
    for (i=0;i<nyrsac_srv1;i++) {
      int nagestmp=1+iac_old_srv1(i)-iac_yng_srv1(i); // .segment takes initial place and vector length
      if(multN_srv1(i)>0) { // this allows turning off contribution easily from dat file
        // Not sure why I have to do this in two steps.
        vector<Type> otmp=srvp1.row(i).segment(iac_yng_srv1(i), nagestmp);
        vector<Type> etmp=Esrvp1.row(isrv_acyrs1(i)).segment(iac_yng_srv1(i), nagestmp);
        otmp+=o; otmp*=multN_srv1(i); etmp+=o;
        vector<Type> alphas = sum(otmp) * etmp * DM_pars(1);
        if(isDM){
          llsrvp1(i) = ddirmultinom(otmp, alphas,  true);
        } else {
          llsrvp1(i) = dmultinom(otmp,etmp,true);
        }
        SIMULATE {
          if(isDM){
            otmp=rdirmultinom(multN_srv1(i), alphas, ages);
          } else {
            otmp=rmultinom(multN_srv1(i), etmp);
          }
          srvp1.row(i).segment(iac_yng_srv1(i), nagestmp)=otmp/otmp.sum();
        }
        //Residuals and outputs which get used by R functions
        for (j=iac_yng_srv1(i);j<=iac_old_srv1(i);j++) {
          res_srv1(i,j)=srvp1(i,j); res_srv1(i,a1-a0+j+1)=Esrvp1(isrv_acyrs1(i),j);
          pearson_srv1(i,j)=(srvp1(i,j)-Esrvp1(isrv_acyrs1(i),j))/sqrt((Esrvp1(isrv_acyrs1(i),j)*(1.-Esrvp1(isrv_acyrs1(i),j)))/multN_srv1(i));
        }
      }
    }
    REPORT(srvp1);
    REPORT(indxsurv1);
    loglik(4)=llsrvp1.sum();
    loglik(5)=0; // unused for now, placeholder for survey 1 length comps

    // Survey 2 likelihoods
    // Total biomass
    RMSE_srv2=0;
    for(i=0; i<nyrs_srv2;i++){
      // Add bias correction on
      Type Eindxsurv2_tmp=log(Eindxsurv2(isrvyrs2(i)));
      Eindxsurv2_tmp -= square(indxsurv_log_sd2(i))/2.;
      if(indxsurv_log_sd2(i)>0){
        loglik(6)+= dnorm(log(indxsurv2(i)), Eindxsurv2_tmp, indxsurv_log_sd2(i), true);
        SIMULATE indxsurv2(i)=exp(rnorm(Eindxsurv2_tmp, indxsurv_log_sd2(i)));
      }
    }
    REPORT(indxsurv2);
    // Age compositions
    res_srv2.setZero(); pearson_srv2.setZero(); llsrvp2.setZero();
    for (i=0;i<nyrsac_srv2;i++) {
      if(multN_srv2(i)>0) { // this allows turning off contribution easily from dat file
        // Not sure why I have to do this in two steps.
        vector<Type> otmp=srvp2.row(i);
        vector<Type> etmp=Esrvp2.row(isrv_acyrs2(i));
        otmp+=o; otmp*=multN_srv2(i); etmp+=o;
        vector<Type> alphas = sum(otmp) * etmp * DM_pars(2);
        if(isDM){
          llsrvp2(i) = ddirmultinom(otmp, alphas,  true);
        } else {
          llsrvp2(i) = dmultinom(otmp, etmp,  true);
        }
        SIMULATE {
          if(isDM){
            otmp=rdirmultinom(multN_srv2(i), alphas, ages);
          } else {
            otmp=rmultinom(multN_srv2(i), etmp);
          }
          srvp2.row(i)=otmp/otmp.sum();
        }
        //Residuals and outputs which get used by R functions
        for (j=a0;j<=a1;j++) {
          res_srv2(i,j)=srvp2(i,j); res_srv2(i,a1-a0+j+1)=Esrvp2(isrv_acyrs2(i),j);
          pearson_srv2(i,j)=(srvp2(i,j)-Esrvp2(isrv_acyrs2(i),j))/sqrt((Esrvp2(isrv_acyrs2(i),j)*(1.-Esrvp2(isrv_acyrs2(i),j)))/multN_srv2(i));
        }
      }
    }
    REPORT(srvp2);
    loglik(7)=llsrvp2.sum();
    llsrvlenp2.setZero();
    for (i=0;i<nyrslen_srv2;i++) {
      if(multNlen_srv2(i)>0) { // this allows turning off contribution easily from dat file
        // Not sure why I have to do this in two steps.
        vector<Type> otmp=srvlenp2.row(i);
        vector<Type> etmp=Esrvlenp2.row(isrv_lenyrs2(i));
        otmp+=o; otmp*=multNlen_srv2(i); etmp+=o;
        llsrvlenp2(i) = dmultinom(otmp, etmp, true);
        SIMULATE {
          otmp=rmultinom(multNlen_srv2(i), etmp);
          srvlenp2.row(i)=otmp/otmp.sum();
        }
      }
    }
    REPORT(srvlenp2);
    loglik(8)=llsrvlenp2.sum();

    // Survey 3 likelihoods
    // Total biomass
    RMSE_srv3=0;
    for(i=0; i<nyrs_srv3;i++){
      // Add bias correction on
      Type Eindxsurv3_tmp=log(Eindxsurv3(isrvyrs3(i)));
      Eindxsurv3_tmp -= square(indxsurv_log_sd3(i))/2.;
      if(indxsurv_log_sd3(i)>0){
        loglik(9)+= dnorm(log(indxsurv3(i)), Eindxsurv3_tmp, indxsurv_log_sd3(i), true);
        SIMULATE indxsurv3(i)=exp(rnorm(Eindxsurv3_tmp, indxsurv_log_sd3(i)));
      }
    }
    REPORT(indxsurv3);
    // Age compositions
    llsrvp3.setZero(); res_srv3.setZero(); pearson_srv3.setZero();
    for (i=0;i<nyrsac_srv3;i++) {
      if(multN_srv3(i)>0) { // this allows turning off contribution easily from dat file
        // Not sure why I have to do this in two steps.
        vector<Type> otmp=srvp3.row(i);
        vector<Type> etmp=Esrvp3.row(isrv_acyrs3(i));
        otmp+=o; otmp*=multN_srv3(i); etmp+=o;
        vector<Type> alphas = sum(otmp) * etmp * DM_pars(3);
        if(isDM){
          llsrvp3(i) = ddirmultinom(otmp, alphas,  true);
        } else {
          llsrvp3(i) = dmultinom(otmp, etmp,  true);
        }
        SIMULATE {
          if(isDM){
            otmp=rdirmultinom(multN_srv3(i), alphas, ages);
          } else {
            otmp=rmultinom(multN_srv3(i), etmp);
          }
          srvp3.row(i)=otmp/otmp.sum();
        }
        //Residuals and outputs which get used by R functions
        for (j=a0;j<=a1;j++) {
          res_srv3(i,j)=srvp3(i,j); res_srv3(i,a1-a0+j+1)=Esrvp3(isrv_acyrs3(i),j);
          pearson_srv3(i,j)=(srvp3(i,j)-Esrvp3(isrv_acyrs3(i),j))/sqrt((Esrvp3(isrv_acyrs3(i),j)*(1.-Esrvp3(isrv_acyrs3(i),j)))/multN_srv3(i));
        }
      }
    }
    REPORT(srvp3);
    loglik(10)=llsrvp3.sum();
    loglik(11)=0; // unused for now, placeholder for survey 3 length comps

    // Survey 4 and 5 likelihoods
    RMSE_srv5=0;  RMSE_srv4=0;
    for(i=0; i<nyrs_srv4;i++){
      // Add bias correction on
      Type Eindxsurv4_tmp=log(Eindxsurv4(isrvyrs4(i)));
      Eindxsurv4_tmp -= square(indxsurv_log_sd4(i))/2.;
      Type Eindxsurv5_tmp=log(Eindxsurv5(isrvyrs5(i)));
      Eindxsurv5_tmp -= square(indxsurv_log_sd5(i))/2.;
      if(indxsurv_log_sd4(i)>0) {
        loglik(12)+= dnorm(log(indxsurv4(i)), Eindxsurv4_tmp, indxsurv_log_sd4(i), true);
        SIMULATE indxsurv4(i)=exp(rnorm(Eindxsurv4_tmp, indxsurv_log_sd4(i)));
      }
      if(indxsurv_log_sd5(i)>0){
        loglik(13)+= dnorm(log(indxsurv5(i)), Eindxsurv5_tmp, indxsurv_log_sd5(i), true);
        SIMULATE indxsurv4(i)=exp(rnorm(Eindxsurv4_tmp, indxsurv_log_sd4(i)));
      }
    }
    REPORT(indxsurv4); REPORT(indxsurv5);

    // Survey 6 likelihoods
    // Total biomass
    RMSE_srv6=0;
    for(i=0; i<nyrs_srv6;i++){
      // Add bias correction on
      Type Eindxsurv6_tmp=log(Eindxsurv6(isrvyrs6(i)));
      Eindxsurv6_tmp -= square(indxsurv_log_sd6(i))/2.;
      if(indxsurv_log_sd6(i)>0){
        loglik(14)+= dnorm(log(indxsurv6(i)), Eindxsurv6_tmp, indxsurv_log_sd6(i), true);
        SIMULATE indxsurv6(i)=exp(rnorm(Eindxsurv6_tmp, indxsurv_log_sd6(i)));
      }
    }
    REPORT(indxsurv6);
    // Age compositions
    llsrvp6.setZero(); res_srv6.setZero(); pearson_srv6.setZero();
    for (i=0;i<nyrsac_srv6;i++) {
      if(multN_srv6(i)>0) { // this allows turning off contribution easily from dat file
        // Not sure why I have to do this in two steps.
        vector<Type> otmp=srvp6.row(i);
        vector<Type> etmp=Esrvp6.row(isrv_acyrs6(i));
        otmp+=o; otmp*=multN_srv6(i);	etmp+=o;
        vector<Type> alphas = sum(otmp) * etmp * DM_pars(4);
        if(isDM){
          llsrvp6(i) = ddirmultinom(otmp, alphas,  true);
        } else {
          llsrvp6(i) = dmultinom(otmp, etmp,  true);
        }
        SIMULATE {
          if(isDM){
            otmp=rdirmultinom(multN_srv6(i), alphas, ages);
          } else {
            otmp=rmultinom(multN_srv6(i), etmp);
          }
          srvp6.row(i)=otmp/otmp.sum();
        }
        //Residuals and outputs which get used by R functions
        for (j=a0;j<=a1;j++) {
          res_srv6(i,j)=srvp6(i,j); res_srv6(i,a1-a0+j+1)=Esrvp6(isrv_acyrs6(i),j);
          pearson_srv6(i,j)=(srvp6(i,j)-Esrvp6(isrv_acyrs6(i),j))/sqrt((Esrvp6(isrv_acyrs6(i),j)*(1.-Esrvp6(isrv_acyrs6(i),j)))/multN_srv6(i));
        }
      }
    }
    REPORT(srvp6);
    loglik(15)=llsrvp6.sum();
    // length comps which are only used in the current year and otherwise have N=0
    llsrvlenp6.setZero();
    for (i=0;i<nyrslen_srv6;i++) {
      if(multNlen_srv6(i)>0) { // this allows turning off contribution easily from dat file
        // Not sure why I have to do this in two steps.
        vector<Type> otmp=srvlenp6.row(i);
        vector<Type> etmp=Esrvlenp6.row(isrv_lenyrs6(i));
        otmp+=o; otmp*=multNlen_srv6(i); etmp+=o;
        llsrvlenp6(i) = dmultinom(otmp, etmp, true);
        SIMULATE {
          otmp=rmultinom(multNlen_srv6(i), etmp);
          srvlenp6.row(i)=otmp/otmp.sum();
        }
      }
    }
    REPORT(srvlenp6);
    loglik(16)+=llsrvlenp6.sum();
    // End of data likelihoods
    // ------------------------------------------------------------

    // ------------------------------------------------------------
    // Priors and random effect penalties
    // TODO: move this above data and build in RE simulation

    // Constraints on recruitment. Assumed sigmaR=1.3 for all devs
    loglik(17) += dnorm(dev_log_recruit, Type(0.0), sigmaR, true).sum();
    // random walk penalties
    for(i=y0+1;i<=y1;i++){
      loglik(18) += dnorm(slp1_fsh_dev(i)-slp1_fsh_dev(i-1), Type(0.0), rwlk_sd(i-1), true);
      loglik(18) += dnorm(inf1_fsh_dev(i)-inf1_fsh_dev(i-1), Type(0.0), 4*rwlk_sd(i-1), true);
      loglik(18) += dnorm(slp2_fsh_dev(i)-slp2_fsh_dev(i-1), Type(0.0), rwlk_sd(i-1), true);
      loglik(18) += dnorm(inf2_fsh_dev(i)-inf2_fsh_dev(i-1), Type(0.0), 4*rwlk_sd(i-1), true);
      loglik(19) += dnorm(log_q1_dev(i)-log_q1_dev(i-1), Type(0.0), q1_rwlk_sd(i-1), true);
      loglik(19) += dnorm(log_q2_dev(i)-log_q2_dev(i-1), Type(0.0), q2_rwlk_sd(i-1), true);
      loglik(19) += dnorm(log_q3_dev(i)-log_q3_dev(i-1), Type(0.0), q3_rwlk_sd(i-1), true);
    }
    loglik(20)=0; // placeholder for projection recruitment
    // // Prior on trawl catchability
    loglik(21) = dnorm(log_q2_mean, log(Type(0.85)), Type(0.1), true);

    // these broad priors on selex stabilize estimation and were
    // determined using prior pushforward checks. Any log_slp >10
    // or < .01 is too flat to be differentiated with data so
    // used as thresholds here.
    loglik(22) += dnorm(log_slp1_fsh_mean, Type(-1.0),Type(1.5), true);
    loglik(22) += dnorm(log_slp2_fsh_mean, Type(-1.0),Type(1.5), true);
    loglik(22) += dnorm(inf1_fsh_mean, Type(0.0),Type(3.0), true);
    loglik(22) += dnorm(inf2_fsh_mean, Type(10.0),Type(3.0), true);
    // loglik(22) += dnorm(log_slp1_srv1, Type(-1.0),Type(1.5), true);
    loglik(22) += dnorm(log_slp2_srv1, Type(-1.0),Type(1.5), true);
    // loglik(22) += dnorm(inf1_srv1, Type(0.0),Type(3.0), true);
    loglik(22) += dnorm(inf2_srv1, Type(10.0),Type(3.0), true);
    loglik(22) += dnorm(log_slp1_srv2, Type(-1.0),Type(1.5), true);
    loglik(22) += dnorm(log_slp2_srv2, Type(-1.0),Type(1.5), true);
    loglik(22) += dnorm(inf1_srv2, Type(0.0),Type(3.0), true);
    loglik(22) += dnorm(inf2_srv2, Type(10.0),Type(3.0), true);
    loglik(22) += dnorm(log_slp1_srv3, Type(-1.0),Type(1.5), true);
    // loglik(22) += dnorm(log_slp2_srv3, Type(-1.0),Type(1.5), true);
    loglik(22) += dnorm(inf1_srv3, Type(0.0),Type(3.0), true);
    // loglik(22) += dnorm(inf2_srv3, Type(10.0),Type(3.0), true);
    loglik(22) += dnorm(log_slp1_srv6, Type(-1.0),Type(1.5), true);
    loglik(22) += dnorm(log_slp2_srv6, Type(-1.0),Type(1.5), true);
    loglik(22) += dnorm(inf1_srv6, Type(0.0),Type(3.0), true);
    loglik(22) += dnorm(inf2_srv6, Type(10.0),Type(3.0), true);

    Type rho=rho_trans(transf_rho);
    Type Ecov_obs_sd=exp(log_Ecov_obs_sd);
    Type Ecov_sd=exp(log_Ecov_sd);//pow(pow(exp(log_Ecov_sd),2)/(1-pow(rho,2)),0.5);
    using namespace density;
    // process error
    loglik(23)= -SCALE(AR1(rho),Ecov_sd)(Ecov_exp);
    // observation error
    for(int i=0; i<Ecov_obs.size(); i++){
      loglik(24)+=dnorm(Ecov_obs(i),Ecov_exp(Ecov_obs_year(i)-styr), Ecov_obs_sd,true);
      SIMULATE Ecov_obs(i)=rnorm(Ecov_exp(i), Ecov_obs_sd);
     }

  objfun = -sum(loglik);

  REPORT(objfun);
  REPORT(styr);
  REPORT(endyr);
  REPORT(ages);
  REPORT(years);
  REPORT(loglik);
  REPORT(recruit);
  REPORT(Espawnbio);
  REPORT(Etotalbio);
  REPORT(Esumbio);
  REPORT(Espawnbio_2plus);
  REPORT(F);
  ADREPORT(F);
  REPORT(Eecocon);
  REPORT(initN);
  REPORT(cattot);
  REPORT(Ecattot);
  REPORT(M);
  REPORT(N);
  REPORT(inf2_fsh_mean);
  REPORT(slp1_fsh_dev);
  REPORT(inf1_fsh_dev);
  REPORT(slp2_fsh_dev);
  REPORT(inf2_fsh_dev);
  REPORT(slp1_fsh);
  REPORT(inf1_fsh);
  REPORT(slp2_fsh);
  REPORT(inf2_fsh);
  REPORT(slctfsh);
  REPORT(llcatp);
  REPORT(catp);
  REPORT(Ecatp);
  REPORT(res_fish);
  REPORT(pearson_fish);
  REPORT(multN_fsh);
  REPORT(effN_fsh);
  REPORT(lllenp);
  REPORT(lenp);
  REPORT(Elenp);
  REPORT(q1);
  REPORT(log_slp2_srv1);
  REPORT(inf2_srv1);
  REPORT(slctsrv1);
  REPORT(Eindxsurv1);
  REPORT(RMSE_srv1);
  REPORT(llsrvp1);
  REPORT(srvp1);
  REPORT(Esrvp1);
  REPORT(res_srv1);
  REPORT(pearson_srv1);
  REPORT(multN_srv1);
  REPORT(effN_srv1);
  REPORT(llsrvlenp1);
  REPORT(srvlenp1);
  REPORT(Esrvlenp1);
  REPORT(q2);
  REPORT(log_slp1_srv2);
  REPORT(inf1_srv2);
  REPORT(log_slp2_srv2);
  REPORT(inf2_srv2);
  REPORT(slctsrv2);
  REPORT(Eindxsurv2);
  REPORT(RMSE_srv2);
  REPORT(llsrvp2);
  REPORT(srvp2);
  REPORT(Esrvp2);
  REPORT(res_srv2);
  REPORT(pearson_srv2);
  REPORT(multN_srv2);
  REPORT(effN_srv2);
  REPORT(llsrvlenp2);
  REPORT(srvlenp2);
  REPORT(Esrvlenp2);
  REPORT(q3);
  REPORT(log_slp1_srv3);
  REPORT(inf1_srv3);
  REPORT(slctsrv3);
  REPORT(Eindxsurv3);
  REPORT(RMSE_srv3);
  REPORT(llsrvp3);
  REPORT(srvp3);
  REPORT(Esrvp3);
  REPORT(res_srv3);
  REPORT(pearson_srv3);
  REPORT(multN_srv3);
  REPORT(effN_srv3);
  REPORT(llsrvlenp3);
  REPORT(srvlenp3);
  REPORT(Esrvlenp3);
  REPORT(res_srv3len);
  REPORT(pearson_srv3len);
  REPORT(q4);
  REPORT(q4_pow);
  REPORT(Eindxsurv4);
  REPORT(RMSE_srv4);
  REPORT(q5);
  REPORT(q5_pow);
  REPORT(Eindxsurv5);
  REPORT(RMSE_srv5);
  REPORT(q6);
  REPORT(log_slp1_srv6);
  REPORT(inf1_srv6);
  REPORT(log_slp2_srv6);
  REPORT(inf2_srv6);
  REPORT(slctsrv6);
  REPORT(Eindxsurv6);
  REPORT(RMSE_srv6);
  REPORT(llsrvp6);
  REPORT(srvp6);
  REPORT(Esrvp6);
  REPORT(res_srv6);
  REPORT(pearson_srv6);
  REPORT(multN_srv6);
  REPORT(effN_srv6);
  REPORT(llsrvlenp6);
  REPORT(srvlenp6);
  REPORT(Esrvlenp6);
  REPORT(slctfsh_proj);
  REPORT(recruit_proj);
  REPORT(log_recr_proj);
  REPORT(sigmasq_recr);
  REPORT(N_proj);
  REPORT(C_proj);
  REPORT(Nsrv_proj);
  REPORT(wt_pop_proj);
  REPORT(wt_spawn_proj);
  REPORT(wt_fsh_proj);
  REPORT(slctfsh_proj);
  REPORT(Ecattot_proj);
  REPORT(Esumbio_proj);
  REPORT(Espawnbio_proj);
  REPORT(Esrv_proj);
  REPORT(Ftarget);
  REPORT(B40);
  REPORT(F_proj);
  REPORT(Z_proj);
  ADREPORT(recruit);
  ADREPORT(log_recruit);
  ADREPORT(recruit_proj);
  ADREPORT(Esumbio_proj);
  ADREPORT(Espawnbio_proj);
  ADREPORT(Esrv_proj);
  ADREPORT(Exrate_proj);
  ADREPORT(log_q1);
  ADREPORT(log_q2);
  ADREPORT(log_q3);
  ADREPORT(log_q4);
  ADREPORT(log_q5);
  ADREPORT(log_q6);
  ADREPORT(endN);
  ADREPORT(slctsrv1);
  ADREPORT(slctsrv2);
  ADREPORT(slctsrv3);
  ADREPORT(slctsrv1_logit);
  ADREPORT(slctsrv2_logit);
  ADREPORT(slctsrv3_logit);
  ADREPORT(slctsrv6_logit);
  ADREPORT(slctfsh_logit);
  ADREPORT(Espawnbio);
  ADREPORT(Esumbio);
  ADREPORT(Espawnbio_log);
  ADREPORT(Esumbio_log);
  REPORT(Ecov_obs_sd);
  ADREPORT(rho);
  REPORT(rho);
  ADREPORT(Ecov_sd);
  REPORT(Ecov_sd);
  REPORT(Ecov_exp);
  REPORT(Ecov_beta);
  return(objfun);
}
