library(ggplot2)
library(mvtnorm)
library(dplyr)

logistic <- function(age, par1 = 0.5, par2 = 5){
  return(1/(1+exp(-par1 * (age - par2))))
}

sim_dat <- data.frame(model = "Parametric", Year = 1, Age = 1:15, Selectivity = logistic(1:15))

# - Semi-parametric
par1_vec <- seq(0.2,0.8, length.out = 6)
for(i in 1:7){
  sim_dat <- rbind(sim_dat, data.frame(model = "Semi-parametric", Year = i,  Age = 1:15,
                                       Selectivity = logistic(1:15, par1 = par1_vec[i]))
  )
}

# - Non-parametric
age_var_cov <- matrix(0, 15, 15)
year_var_cov <- matrix(0, 7, 7)
rho_age = 0.5
for(a1 in 1:15){
  for(a2 in 1:15){
    age_var_cov[a1,a2] = rho_age^abs(a1-a2)
  }
}
rho_yr = 0.2
for(y1 in 1:7){
  for(y2 in 1:7){
    year_var_cov[y1, y2] = rho_yr^abs(y1-y2)
  }
}
diag(age_var_cov) <- 1
diag(year_var_cov) <- 1
sigma_yr <- 0.6
sigma_age = 0.7
age_var_cov <- diag(sigma_age,15) %*% age_var_cov %*% diag(sigma_age,15) # Get covariance
year_var_cov <- diag(sigma_yr,7) %*% year_var_cov %*% diag(sigma_yr,7) # Get covariance

colnames(age_var_cov) <- rownames(age_var_cov) <- paste0("Age",1:15)
colnames(year_var_cov) <- rownames(year_var_cov) <- paste0("Yr",1:7)

npdat <- expand.grid( Year = 1:7, Age = 1:15)
npdat$model = "Non-parametric"
var_cov_full <- kronecker(year_var_cov, age_var_cov)
npdat$Selectivity <- as.numeric(rmvnorm(1,  mean = rep(c(-0.6, - 0.5, -0.25, 0, 0.25, 0, -0.25) * 4, each = 15), var_cov_full))
npdat$Selectivity <- 1/(1+exp(-(npdat$Selectivity )))

npdat <- npdat %>%
  select(model, Year, Age, Selectivity)




# plot
sim_dat <- rbind(sim_dat, npdat)
sim_dat <- sim_dat %>%
  mutate(model = factor(model, levels = c("Parametric", "Semi-parametric", "Non-parametric")))
ggplot(sim_dat, aes(x = Age, y = Selectivity, colour = Year, group = Year)) +
  geom_line() +
  facet_wrap(~model, nrow = 3) +
  theme_classic()
