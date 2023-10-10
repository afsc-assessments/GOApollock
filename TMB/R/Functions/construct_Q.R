#
# png('cors.png', width=7, height=3, units='in', res=500)
# par(mfrow=c(1,3))
# plot_cor(3,6, .9,0,.0)
# plot_cor(3,6, 0,.9,0)
# plot_cor(3,6, 0,0,-1.9)
# dev.off()

plot_cor <- function(n_years, n_ages, rho_y, rho_a, rho_c, log_sigma){
  ay_Index <- expand.grid(age=1:n_ages, year=1:n_years)
  total_n <- n_years * n_ages;
  B <- Omega <- Q_sparse <- matrix(0, total_n,total_n); ## B matrix
  I <- diag(total_n)
  for(n in 1:total_n){
    ## Define year and age objects
    age = ay_Index[n,1];
    year = ay_Index[n,2];
    ## Constructing B matrix to determine where the correlation pars should go
    if(age > 1) {
      ## Get column index for years
      for(n1 in 1:total_n) {
        if(ay_Index[n1, 1] == age - 1 && ay_Index[n1, 2] == year)
          ##B.coeffRef[n, n1] = rho_y;
          B[n, n1] <- rho_y;
      } ## n1 loop
    } ## end age > 1
    if(year > 1) {
      ## Get column index for years
      for(n1 in 1:total_n) {
        if(ay_Index[n1,1] == age && ay_Index[n1, 2] == year - 1)
          B[n, n1] = rho_a;
      } ## n1 loop
    } ## if year > 1
    if(year > 1 && age > 1) {
      ## Get column index for years
      for(n1 in 1:total_n) {
        if(ay_Index[n1,1] == age - 1 && ay_Index[n1, 2] == year - 1)
          B[n,n1] = rho_c; ## correlation by cohort
      } ## n1 loop
    } ## if both year and age > 1
  } ## end n loop

  ## Fill in Omega matrix here (variances)
  Var_Param <- 0
  if(Var_Param == 0) { ## Conditional variance
    for(i in 1:total_n) {
      for(j in 1:total_n) {
        if(i == j){
          Omega[i,j] = 1/exp(log_sigma)^2;
        }  else {
          Omega[i,j] = 0
        }
      } ## j loop
    } ## i loop
  } ## end if conditional variance

  ## if(Var_Param == 1) { ## Marginal Variance
  ##   ## Construct container objects
  ##   matrix<Type> L(total_n, total_n); ## L Matrix
  ##   matrix<Type> tmp_I_B = I-B; ## Temporary Matrix to store I-B
  ##   L =  tmp_I_B.inverse(); ## Invert to get L
  ##   vector<Type> d(total_n); ## Store variance calculations
  ##   for(int n = 0; n < total_n; n++) {
  ##     if(n == 0) {
  ##       d(n) = exp(log_sigma2); ## marginal variance parameter
  ##     } else{
  ##       Type cumvar = 0; ## Cumulative Variance Container
  ##       for(int n1 = 0; n1 < n; n1++) {
  ##         cumvar += L(n, n1) * d(n1) * L(n, n1);
  ##       } ## n1 loop
  ##       ## Calculate diagonal values for omega
  ##       d(n) = (exp(log_sigma2) - cumvar) / pow(L(n, n), 2);
  ##     } ## else loop
  ##   } ## n loop
  ##   ## Now fill in our diagonals for Omega
  ##   for(int i = 0; i < total_n; i++) {
  ##     for(int j = 0; j < total_n; j++) {
  ##       if(i == j) Omega.coeffRef(i,j) = 1/d(i);
  ##       else Omega.coeffRef(i,j) = Type(0.0);
  ##     } ## j loop
  ##   } ## i loop
  ## } ## end if marginal variance

  ## Now, do calculations to construct (Q = (I - t(B)) %*% Omega %*% (I-B))
  B_transpose <- t(B); ## transpose B matrix

  ## Calculate Precision Matrix
  Q_sparse = (I - B_transpose) %*% Omega %*% (I-B);

  covar <- solve(Q_sparse)
  library(corrplot)
  x <- paste0('a=', ay_Index[,1], ', y=',ay_Index[,2])
  corr <- cov2cor(covar)
  dimnames(corr) <- list(x,x)
  corrplot(corr, type='lower')
}
