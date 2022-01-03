
test_that("write_dat works",{
  setwd('../model_20_8')
  system("goa_pk -ind goa_pk.dat -maxfn 0", ignore.stdout = TRUE)
  rep1 <- read_pk_rep(endyr=2020, version='test')
  dat <- read_dat('goa_pk.dat')
  write_dat(datlist=dat, fileout = 'goa_pk2.dat')
  system("goa_pk -ind goa_pk2.dat -maxfn 0", ignore.stdout = TRUE)
  rep2 <- read_pk_rep(endyr=2020, version='test')
  clean_pk_dir()
  setwd('../testthat')
  expect_equal(rep1$Objective_function,rep2$Objective_function)
  expect_equal(rep1$Total_catches,rep2$Total_catches)
})

