library(GOApollock)

test_that("retros work",{
  skip_if(TRUE)
  setwd('../model_20_8/')
  stds <- reps <- list()
  f1 <- 'goa_pk.std'
  trash <- if(file.exists(f1)) file.remove(f1)
  ## run it once without retros to make sure -0 option is
  ## identical
  message("Running base model")
  system(paste("goa_pk -display 0"))
  std0 <- read_pk_std(version=0, endyr=2021)
  k <- 0
  for(i in 0:5){
    trash <- if(file.exists(f1)) file.remove(f1)
    system(paste("goa_pk -display 0 -retro",i))
    if(file.exists(f1)){
      k <- k+1
      stds[[i+1]] <- read_pk_std(version=i, endyr=2021-i)
      reps[[i+1]] <- read_pk_rep(version=i, endyr=2021-i)
    } else {
      stop("Failed to produce .std file for peel ", i)
    }
  if(i==0) expect_equal(std0,stds[[i+1]])
  }
  expect_true(k==6)
  rho <- calculate_rho(reps)
  expect_equal(-0.104,rho)
})


