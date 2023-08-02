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

// July 2023 started conversion to TMB

#include <TMB.hpp>
#include <iostream>
//#include "helper_functions.hpp"


// -----------------------------------------------
// FUNCTIONS ----
// -----------------------------------------------
template <class Type>
Type square(Type x){return x*x;};
template <class Type>
Type norm2(vector<Type> x){return (x*x).sum();};
#define see(object) std::cout << #object ":\n" << object << "\n";

// transformation to ensure correlation parameters are between -1 and 1
// 2/(1+exp(-2*x)) - 1
template <class Type>
Type rho_trans(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}


// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Function to assemble sparse precision matrix
template<class Type>
// @description: Function that constructs a precision matrix, separable along the
// year, age, and cohort axis. Var_Param allows users to switch between conditional
// variance, and marginal variance.
Eigen::SparseMatrix<Type> construct_Q(int n_years, // Integer of years
                                      int n_ages, // Integer of ages
                                      matrix<Type> ay_Index, // Index matrix to construct
                                      Type rho_y, // Partial correlation by years
                                      Type rho_a, // Partial correlation by ages
                                      Type rho_c, // Partial correlation by cohort
                                      Type log_sigma2, // Variance parameter governing GMRF
                                      int Var_Param // Parameterization of Variance ==0 (Conditional), == 1(Marginal)
) {

  // Dimension to construct matrices
  int total_n = n_years * n_ages;

  // Construct matrices for precision matrix
  Eigen::SparseMatrix<Type> B(total_n,total_n); // B matrix
  Eigen::SparseMatrix<Type> I(total_n,total_n); // Identity matrix
  I.setIdentity(); // Set I to identity matrix
  Eigen::SparseMatrix<Type> Omega(total_n,total_n); // Omega matrix (variances)
  Eigen::SparseMatrix<Type> Q_sparse(total_n, total_n); // Precision matrix

  for(int n = 0; n < total_n; n++) {

    // Define year and age objects
    Type age = ay_Index(n,0);
    Type year = ay_Index(n,1);

    // Constructing B matrix to determine where the correlation pars should go
    if(age > 1) {

      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year)
          B.coeffRef(n, n1) = rho_y;
      } // n1 loop

    } // end age > 1

    if(year > 1) {

      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age && ay_Index(n1, 1) == year - 1)
          B.coeffRef(n, n1) = rho_a;
      } // n1 loop

    } // if year > 1

    if(year > 1 && age > 1) {

      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year - 1)
          B.coeffRef(n,n1) = rho_c; // correlation by cohort
      } // n1 loop

    } // if both year and age > 1

  } // end n loop

  // Fill in Omega matrix here (variances)
  if(Var_Param == 0) { // Conditional variance

    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/exp(log_sigma2);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop

  } // end if conditional variance

  if(Var_Param == 1) { // Marginal Variance

    // Construct container objects
    matrix<Type> L(total_n, total_n); // L Matrix
    matrix<Type> tmp_I_B = I-B; // Temporary Matrix to store I-B
    L =  tmp_I_B.inverse(); // Invert to get L
    vector<Type> d(total_n); // Store variance calculations

    for(int n = 0; n < total_n; n++) {
      if(n == 0) {
        d(n) = exp(log_sigma2); // marginal variance parameter
      } else{

        Type cumvar = 0; // Cumulative Variance Container

        for(int n1 = 0; n1 < n; n1++) {
          cumvar += L(n, n1) * d(n1) * L(n, n1);
        } // n1 loop

        // Calculate diagonal values for omega
        d(n) = (exp(log_sigma2) - cumvar) / pow(L(n, n), 2);

      } // else loop
    } // n loop

    // Now fill in our diagonals for Omega
    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/d(i);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop

  } // end if marginal variance

  // Now, do calculations to construct (Q = (I - t(B)) %*% Omega %*% (I-B))
  Eigen::SparseMatrix<Type> B_transpose = B.transpose(); // transpose B matrix

  // Calculate Precision Matrix
  Q_sparse = (I - B_transpose) * Omega * (I-B);

  return(Q_sparse);

} // end construct_Q function


// -----------------------------------------------
// POLLOCK MODEL ----
// -----------------------------------------------
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  // -----------------------------------------------
  // DATA INPUTS ----
  // -----------------------------------------------
  DATA_INTEGER(styr);	    // Starting year for population model
  DATA_INTEGER(endyr);	    // Ending year for population model
  DATA_INTEGER(rcrage);	    // Recruitment age
  DATA_INTEGER(trmage);	    // Last modeled age
  DATA_INTEGER(nbins1);	// Number of length bins in transitiom matrix 1
  DATA_INTEGER(nbins2);	// Number of length bins in transitiom matrix 2
  DATA_INTEGER(nbins3);	// Number of length bins in transitiom matrix 3

  // Fishery
  DATA_MATRIX(ay_Index); // (n_years * n_ages), 2
  DATA_INTEGER(sel_vartype);  // Parameterization of 3D AR Precision Matrix: == 0 (Conditional), == 1(Marginal)
  DATA_INTEGER(seltype);    // Fishery selectivity form (1 = double logistic, 2 = age-specific AR1)
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

  // Population vectors
  matrix<Type> wt_pop(nyrs,nages);   // Population weight at age
  matrix<Type> wt_spawn(nyrs,nages); // Population weight at age at spawning (April 15)
  wt_pop=wt_srv2;
  wt_spawn=wt_srv1;
  DATA_VECTOR(mat);		// Proportion mature

  DATA_VECTOR(Ftarget);
  DATA_SCALAR(B40); // mean log recruitment
  // DATA_SCALAR(log_mean_recr_proj);
  DATA_SCALAR(sigmasq_recr);                       // Variance for log recr, recruitment indices

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
  wt_srv_proj=wt_spawn_proj;
  wt_fsh_proj=wt_fsh(nyrs);

  // -----------------------------------------------
  // POPULATION PARAMETERS
  // -----------------------------------------------
  // RECRUITMENT PARAMETERS
  PARAMETER_VECTOR(dev_log_initN);
  vector<Type> initN(nages-1); // goes from 2 to 10
  PARAMETER(mean_log_recruit);
  PARAMETER_VECTOR(dev_log_recruit);
  PARAMETER(sigmaR);
  vector<Type> recruit(nyrs);
  vector<Type> log_recruit(nyrs);

  // Forward projections
  PARAMETER_VECTOR(log_recr_proj);
  vector<Type> recruit_proj(5);
  matrix<Type> N_proj(5,nages);
  vector<Type> F_proj(5);
  matrix<Type> Z_proj(5,nages);
  matrix<Type> C_proj(5,nages);
  matrix<Type> Nsrv_proj(5,nages);
  vector<Type> slctfsh_proj(nages);
  vector<Type> Ecattot_proj(5);
  vector<Type> Esumbio_proj(5);
  vector<Type> Espawnbio_proj(5);
  vector<Type> Esrv_proj(5);
  vector<Type> Exrate_proj(5);

  // SELECTIVITY PARAMETERS
  // - Fishery selectivity
  PARAMETER(mean_sel); // Mean selectivity
  PARAMETER_ARRAY(selpars_re); // AR selectivity parameters
  PARAMETER(sel_rho_a); // AR1 age correlation
  PARAMETER(sel_rho_y); // AR1 year correlation
  PARAMETER(sel_rho_c); // AR1 cohort correlation
  PARAMETER(log_slp1_fsh_mean);
  PARAMETER(inf1_fsh_mean);
  PARAMETER(log_slp2_fsh_mean);
  PARAMETER(inf2_fsh_mean);
  PARAMETER_VECTOR(slp1_fsh_dev);
  PARAMETER_VECTOR(inf1_fsh_dev);
  PARAMETER_VECTOR(slp2_fsh_dev);
  PARAMETER_VECTOR(inf2_fsh_dev);
  PARAMETER(ln_sel_sd); // SD of sel parms

  Type Sigma_sig_sel = 0;
  Type sel_sd = exp(ln_sel_sd);
  Type rho_a = rho_trans(sel_rho_a); // Scale from -1 to 1
  Type rho_y = rho_trans(sel_rho_y);
  Type rho_c = rho_trans(sel_rho_c);
  vector<Type> slp1_fsh(nyrs);
  vector<Type> inf1_fsh(nyrs);
  vector<Type> slp2_fsh(nyrs);
  vector<Type> inf2_fsh(nyrs);

  // -- Define precision matrix for GMRF
  Eigen::SparseMatrix<Type> Q_sparse(nages * nyrs, nages * nyrs); // Precision matrix

  // -- Construct precision matrix here
  Q_sparse = construct_Q(nyrs, nages, ay_Index,
                         rho_y, rho_a, rho_c,
                         log(square(sel_sd)), sel_vartype);

  // - Acoustic survey selectivity
  PARAMETER(log_slp2_srv1);
  PARAMETER(inf2_srv1);

  // - Trawl selectivity
  PARAMETER(log_slp1_srv2);
  PARAMETER(inf1_srv2);
  PARAMETER(log_slp2_srv2);
  PARAMETER(inf2_srv2);

  // - ADFG selectivity
  PARAMETER(log_slp1_srv3);
  PARAMETER(inf1_srv3);

  // - Summer acoustic selectivity
  PARAMETER(log_slp1_srv6);
  PARAMETER(inf1_srv6);
  PARAMETER(log_slp2_srv6);
  PARAMETER(inf2_srv6);

  // - Fishing mortality
  PARAMETER(mean_log_F);
  PARAMETER_VECTOR(dev_log_F);
  vector<Type> F(nyrs);

  // - Survey catchablility
  PARAMETER(log_q1_mean);
  PARAMETER_VECTOR(log_q1_dev);
  PARAMETER(log_q2_mean);
  PARAMETER_VECTOR(log_q2_dev);
  PARAMETER(log_q3_mean);
  PARAMETER_VECTOR(log_q3_dev);
  PARAMETER(log_q4);
  PARAMETER(q4_pow);
  PARAMETER(log_q5);
  PARAMETER(q5_pow);
  PARAMETER(log_q6);

  // This scales M vector below so that M={M}*natMscalar. If 1 does nothing.
  PARAMETER(natMscalar);

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


  // -----------------------------------------------
  // Expected values
  // -----------------------------------------------
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
  vector<Type> loglik(24);
  loglik.setZero();
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

  vector<Type> effN_fsh(nyrs_fsh);
  vector<Type> effN_srv1(nyrsac_srv1);
  vector<Type> effN_srv2(nyrsac_srv2);
  vector<Type> effN_srv3(nyrsac_srv3);
  vector<Type> effN_srv6(nyrsac_srv6);

  // -----------------------------------------------
  // Objective function
  // -----------------------------------------------
  Type objfun=0.0;

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
  vector<int> iac_yng_fsh=ac_yng_fsh-1;
  vector<int> iac_old_fsh=ac_old_fsh-1;
  vector<int> iac_yng_srv1=ac_yng_srv1-1;
  vector<int> iac_old_srv1=ac_old_srv1-1;


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
    log_q1(i)=log_q1_mean+log_q1_dev(i);
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

  switch(seltype){
  // - Double logistic with random effects on parameters
  case 1:
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
    break;

    // - Double logistic with age AR1 random effects on selectivity function
  case 2:
    for (i=y0;i<=y1;i++) {
      slp1_fsh(i)=exp(log_slp1_fsh_mean);
      inf1_fsh(i)=inf1_fsh_mean;
      slp2_fsh(i)=exp(log_slp2_fsh_mean);
      inf2_fsh(i)=inf2_fsh_mean;
      for (j=a0;j<=a1;j++) {
        slctfsh(i,j) = (1/(1+exp(-(slp1_fsh(i))*(double(j+1)-(inf1_fsh(i))))))*
          (1-1/(1+exp(-(slp2_fsh(i))*(double(j+1)-(inf2_fsh(i)))))) +  selpars_re(j,0);
      }
      // The plan would be to check and adjust the max selected age as needed
      slctfsh.row(i)=slctfsh.row(i)/slctfsh(i,6);
    }
    break;

    // - Double logistic with age,year 2D AR1 random effects on selectivity function
  case 3:
    for (i=y0;i<=y1;i++) {
      slp1_fsh(i)=exp(log_slp1_fsh_mean);
      inf1_fsh(i)=inf1_fsh_mean;
      slp2_fsh(i)=exp(log_slp2_fsh_mean);
      inf2_fsh(i)=inf2_fsh_mean;
      for (j=a0;j<=a1;j++) {
        slctfsh(i,j) = (1/(1+exp(-(slp1_fsh(i))*(double(j+1)-(inf1_fsh(i))))))*
          (1-1/(1+exp(-(slp2_fsh(i))*(double(j+1)-(inf2_fsh(i)))))) + selpars_re(j,i);
      }
      // The plan would be to check and adjust the max selected age as needed
      slctfsh.row(i)=slctfsh.row(i)/slctfsh(i,6);
    }
    break;

    // Non-parametric AR1 on age
  case 4:
    for (i=y0;i<=y1;i++) {
      for (j=a0;j<=a1;j++) {
        slctfsh(i,j) = 1 / (1 + exp(-(mean_sel + selpars_re(j,0)))); // Random effects are constant across years and cohorts
      }
      slctfsh.row(i)=slctfsh.row(i)/slctfsh(i,6);
    }
    break;

    // Non-parametric 2D-AR1 on age and year
  case 5:
    for (i=y0;i<=y1;i++) {
      for (j=a0;j<=a1;j++) {
        slctfsh(i,j) = 1 / (1 + exp(-(mean_sel + selpars_re(j,i)))); // Random effects are constant across years and cohorts
      }
      //slctfsh.row(i)=slctfsh.row(i)/slctfsh(i,6);
    }
    break;

    // Non-parametric 3D-AR1 on age, year, and cohort
  case 6:
    for (i=y0;i<=y1;i++) {
      for (j=a0;j<=a1;j++) {
        slctfsh(i,j) = 1 / (1 + exp(-(mean_sel + selpars_re(j,i)))); // Random effects are constant across years and cohorts
      }
    }
    break;
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
  slctfsh_logit=-log(1/(slctfsh(y1-1)-1e-10)-1);

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


  // -----------------------------------------------
  // LIKELIHOODS
  // -----------------------------------------------
  loglik.setZero();

  // age accumulation (add 1st and 2nd ages) for fishery, turned off for the surveys
  for (i=0;i<nyrs_fsh;i++) {
    if(ifshyrs(i)>endyr) break;
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
    if(isrvyrs1(i)>endyr) break;
    for (j=a0;j<=a1;j++) {
      if(j<iac_yng_srv1(i)) {
        Esrvp1(isrvyrs1(i),iac_yng_srv1(i)) += Esrvp1(isrvyrs1(i),j);
        Esrvp1(isrvyrs1(i),j) = 0;
        srvp1(i,iac_yng_srv1(i)) += srvp1(i,j);
        srvp1(i,j) = 0;
      }
    }
  }

  // Fishery likelihoods
  //Total catch
  for(i=y0; i<=y1;i++){
    if(i>y1) break;
    loglik(0) += -.5*square((log(cattot(i))-log(Ecattot(i)))/cattot_log_sd(i));
  }


  // Age composition
  for (i=0;i<nyrs_fsh;i++) {
    // if(ifshyrs(i)>y1) break; 	// ignore data after retroyear
    llcatp(i) = 0;
    for (j=iac_yng_fsh(i);j<=iac_old_fsh(i);j++) {
      llcatp(i) += multN_fsh(i)*(catp(i,j)+o)*log((Ecatp(ifshyrs(i),j)+o)/(catp(i,j)+o));
      res_fish(i,j)=catp(i,j);
      res_fish(i,a1-a0+j+1)=Ecatp(ifshyrs(i),j);
      if(multN_fsh(i)>0) {
        pearson_fish(i,j)=(catp(i,j)-Ecatp(ifshyrs(i),j))/sqrt((Ecatp(ifshyrs(i),j)*(1.-Ecatp(ifshyrs(i),j)))/multN_fsh(i));
      }
    }
    if(multN_fsh(i)>0) {
      //   effN_fsh(i) = sum(Ecatp(ifshyrs(i))*(1-Ecatp(ifshyrs(i))))/sum(square(catp(i)-Ecatp(ifshyrs(i))));
    }
    loglik(1) += llcatp(i);
  }


  // Survey 1 (shelikof acoustic) likelihoods
  // - Total biomass
  loglik(3)=0;
  for(i=0; i<nyrs_srv1;i++){
    if(isrvyrs1(i)>y1) break;
    loglik(3)+=-.5*square((log(indxsurv1(i))-log(Eindxsurv1(isrvyrs1(i)))+square(indxsurv_log_sd1(i))/2.)/indxsurv_log_sd1(i));
  }

  // - Age composition
  loglik(4)=0;
  for (i=0;i<nyrsac_srv1;i++) {
    if(isrv_acyrs1(i)>y1) break;
    llsrvp1(i) = 0;
    for (j=a0;j<=a1;j++) {
      llsrvp1(i) += multN_srv1(i)*(srvp1(i,j)+o)*log((Esrvp1(isrv_acyrs1(i),j)+o)/(srvp1(i,j)+o));
      res_srv1(i,j)=srvp1(i,j);
      res_srv1(i,a1-a0+j+1)=Esrvp1(isrv_acyrs1(i),j);
      if(multN_srv1(i)>0) {
        pearson_srv1(i,j)=(srvp1(i,j)-Esrvp1(isrv_acyrs1(i),j))/sqrt((Esrvp1(isrv_acyrs1(i),j)*(1.-Esrvp1(isrv_acyrs1(i),j)))/multN_srv1(i));
      }
    }
    if(multN_srv1(i)>0) {
      //effN_srv1(i) = sum(Esrvp1(isrv_acyrs1(i))*(1-Esrvp1(isrv_acyrs1(i))))/sum(square(srvp1(i)-Esrvp1(isrv_acyrs1(i))));
    }
    loglik(4) += llsrvp1(i);
  }

  // Survey 2 (bottom trawl) likelihoods
  // - Total biomass
  loglik(6) =0;
  for (i=0;i<nyrs_srv2;i++){
    if(isrvyrs2(i)>y1) break;
    loglik(6)+=-.5*square((log(indxsurv2(i))-log(Eindxsurv2(isrvyrs2(i)))+square(indxsurv_log_sd2(i))/2.)/indxsurv_log_sd2(i));
  }

  // - Age composition
  loglik(7)=0;
  for (i=0;i<nyrsac_srv2;i++) {
    llsrvp2(i) = 0;
    if(isrv_acyrs2(i)>y1) break;
    for (j=a0;j<=a1;j++) {
      llsrvp2(i) += multN_srv2(i)*(srvp2(i,j)+o)*log((Esrvp2(isrv_acyrs2(i),j)+o)/(srvp2(i,j)+o));
      res_srv2(i,j)=srvp2(i,j);
      res_srv2(i,a1-a0+j+1)=Esrvp2(isrv_acyrs2(i),j);
      if(multN_srv2(i)>0) {
        pearson_srv2(i,j)=(srvp2(i,j)-Esrvp2(isrv_acyrs2(i),j))/sqrt((Esrvp2(isrv_acyrs2(i),j)*(1.-Esrvp2(isrv_acyrs2(i),j)))/multN_srv2(i));
      }
    }
    if(multN_srv2(i)>0) {
      //effN_srv2(i) = sum(Esrvp2(isrv_acyrs2(i))*(1-Esrvp2(isrv_acyrs2(i))))/sum(square(srvp2(i)-Esrvp2(isrv_acyrs2(i))));
    }
    loglik(7) += llsrvp2(i);
  }
  loglik(9) = 0;

  // Survey 3 (ADFG coastal survey) likelihoods
  // - Total biomass
  loglik(10)=0;
  for(i=0; i<nyrs_srv3;i++){
    if(isrvyrs3(i)>y1) break;
    loglik(10) += -.5*square((log(indxsurv3(i))-log(Eindxsurv3(isrvyrs3(i)))+square(indxsurv_log_sd3(i))/2.)/indxsurv_log_sd3(i));
  }

  // - Age composition
  loglik(11)=0;
  for (i=0;i<nyrsac_srv3;i++) {
    if(isrv_acyrs3(i)>y1) break;
    llsrvp3(i) = 0;
    for (j=a0;j<=a1;j++) {
      llsrvp3(i) += multN_srv3(i)*(srvp3(i,j)+o)*log((Esrvp3(isrv_acyrs3(i),j)+o)/(srvp3(i,j)+o));
      res_srv3(i,j)=srvp3(i,j);
      res_srv3(i,a1-a0+j+1)=Esrvp3(isrv_acyrs3(i),j);
      if(multN_srv3(i)>0) {
        pearson_srv3(i,j)=(srvp3(i,j)-Esrvp3(isrv_acyrs3(i),j))/sqrt((Esrvp3(isrv_acyrs3(i),j)*(1.-Esrvp3(isrv_acyrs3(i),j)))/multN_srv3(i));
      }
    }
    if(multN_srv3(i)>0) {
      //effN_srv3(i) = sum(Esrvp3(isrv_acyrs3(i))*(1-Esrvp3(isrv_acyrs3(i))))/sum(square(srvp3(i)-Esrvp3(isrv_acyrs3(i))));
    }
    loglik(11) += llsrvp3(i);
  }

  // Survey 4 and 5 likelihoods (age-1 and -2 acoustic)
  for(i=0; i<nyrs_srv4;i++){ 	// assuming srv4 and srv5 have identical structure
    if(isrvyrs4(i) >y1) break;
    loglik(13) += -.5*square((log(indxsurv4(i))-log(Eindxsurv4(isrvyrs4(i)))+square(indxsurv_log_sd4(i))/2.)/indxsurv_log_sd4(i));
    loglik(14) += -.5*square((log(indxsurv5(i))-log(Eindxsurv5(isrvyrs5(i)))+square(indxsurv_log_sd5(i))/2.)/indxsurv_log_sd5(i));
  }

  // Survey 6 likelihoods (summer acoustic)
  // - Total biomass
  loglik(15)=0;
  for(i=0;i<nyrs_srv6;i++){
    if(isrvyrs6(i)>y1) break;
    loglik(15)+=-.5*square((log(indxsurv6(i))-log(Eindxsurv6(isrvyrs6(i)))+square(indxsurv_log_sd6(i))/2.)/indxsurv_log_sd6(i));
  }

  // - Age composition
  loglik(16)=0;
  for (i=0;i<nyrsac_srv6;i++) {
    if(isrv_acyrs6(i)>y1) break;
    llsrvp6(i) = 0;
    for (j=a0;j<=a1;j++) {
      llsrvp6(i) += multN_srv6(i)*(srvp6(i,j)+o)*log((Esrvp6(isrv_acyrs6(i),j)+o)/(srvp6(i,j)+o));
      res_srv6(i,j)=srvp6(i,j);
      res_srv6(i,a1-a0+j+1)=Esrvp6(isrv_acyrs6(i),j);
      if(multN_srv6(i)>0) {
        pearson_srv6(i,j)=(srvp6(i,j)-Esrvp6(isrv_acyrs6(i),j))/sqrt((Esrvp6(isrv_acyrs6(i),j)*(1.-Esrvp6(isrv_acyrs6(i),j)))/multN_srv6(i));
      }
    }
    loglik(16) +=llsrvp6(i);
  }

  // Constraints on recruitment. Assumed sigmaR=1.3 for all devs
  loglik(17) += -0.5*norm2(dev_log_recruit)/square(sigmaR);

  // Fishery selectivity
  // - Double logistic with deviates
  if(seltype == 1){
    for(i=y0+1;i<=y1;i++){
      loglik(18) += -0.5*square( (slp1_fsh_dev(i)-slp1_fsh_dev(i-1))/sel_sd);
      loglik(18) += -0.5*square( (inf1_fsh_dev(i)-inf1_fsh_dev(i-1))/(4*sel_sd));
      loglik(18) += -0.5*square( (slp2_fsh_dev(i)-slp2_fsh_dev(i-1))/sel_sd);
      loglik(18) += -0.5*square( (inf2_fsh_dev(i)-inf2_fsh_dev(i-1))/sel_sd);
    }
  }

  // - AR1 on age
  if((seltype == 2) | (seltype == 4)){
    vector<Type> tmp_AR1 = selpars_re.col(0); // Random effects are constant across years and cohorts
    Sigma_sig_sel = pow(pow(sel_sd,2) / (1-pow(rho_a,2)),0.5);
    loglik(18) -= SCALE(AR1(rho_a), Sigma_sig_sel)(tmp_AR1);
  }

  //- 2D AR1 on age/year
  if((seltype == 3 )| (seltype == 5)){
    Sigma_sig_sel = pow(pow(sel_sd,2) / ((1-pow(rho_y,2))*(1-pow(rho_a,2))),0.5);
    loglik(18) -= SCALE(SEPARABLE(AR1(rho_a),AR1(rho_y)), Sigma_sig_sel)(selpars_re);
  }

  //- 3D AR1 on age/year
  if((seltype == 6 )){
      loglik(18) -= GMRF(Q_sparse)(selpars_re);
  }


  // Random walk for q devs
  for(i=y0+1;i<=y1;i++){
    loglik(20) += -0.5*square( (log_q1_dev(i)-log_q1_dev(i-1))/q1_rwlk_sd(i-1));
    loglik(20) += -0.5*square( (log_q3_dev(i)-log_q3_dev(i-1))/q3_rwlk_sd(i-1));
  }

  // // Prior on trawl catchability
  loglik(22) = -.5*square((log_q2_mean-log(0.85))/0.1);
  // these broad priors stabilize estimation, particularly for
  // retros, and imply uniform 1 selex for this survey
  loglik(23) += dnorm(log_slp2_srv6, Type(0.0),Type(2.0), true);
  loglik(23) += dnorm(inf2_srv6, Type(10.0),Type(3.0), true);
  objfun = -sum(loglik);



  // -----------------------------------------------
  // REPORT
  // -----------------------------------------------
  REPORT(objfun);
  REPORT(loglik);
  REPORT(recruit);
  REPORT(Espawnbio);
  REPORT(Etotalbio);
  REPORT(Esumbio);
  REPORT(Espawnbio_2plus);
  REPORT(F);
  REPORT(Eecocon);
  REPORT(initN);
  REPORT(cattot);
  REPORT(Ecattot);
  REPORT(M);
  REPORT(N);

  // Fish sel
  REPORT(Q_sparse);
  REPORT(sel_sd);
  REPORT(inf2_fsh_mean);
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
  REPORT(q5);
  REPORT(q5_pow);
  REPORT(Eindxsurv5);
  REPORT(q6);
  REPORT(log_slp1_srv6);
  REPORT(inf1_srv6);
  REPORT(log_slp2_srv6);
  REPORT(inf2_srv6);
  REPORT(slctsrv6);
  REPORT(Eindxsurv6);
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
  ADREPORT(endN);
  ADREPORT(slctsrv1);
  ADREPORT(slctsrv2);
  ADREPORT(slctsrv3);
  ADREPORT(slctsrv1_logit);
  ADREPORT(slctsrv2_logit);
  ADREPORT(slctsrv3_logit);
  ADREPORT(slctsrv6_logit);
  ADREPORT(slctfsh_logit);
  ADREPORT(slctfsh);
  ADREPORT(Espawnbio);
  ADREPORT(Esumbio);
  ADREPORT(Espawnbio_log);
  ADREPORT(Esumbio_log);

  return(objfun);
}
