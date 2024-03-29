
test_that("simulator works",{
  setwd('../model_20_8')
  datlist <- read_dat('goa_pk.dat')
  replist <- read_pk_rep(endyr=2021, version='simulated')
  set.seed(231512)
  file.remove('goa_pk.rep')
  sim_dat(datlist, replist, fileout='goa_pk_sim.dat')
  system("goa_pk -nohess -goa_pk_sim.dat", ignore.stdout = TRUE)
  x <- read_pk_rep(endyr=2021, version='simulated')
  file.remove('goa_pk_sim.dat')
  setwd('../testthat')
  ## saveRDS(x,  "_expect_sim_rep.rds")
  oldrep <- readRDS('_expect_sim_rep.rds')
  expect_equal(x, oldrep)
})
