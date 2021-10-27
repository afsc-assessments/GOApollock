
test_that("read rep and cor and mymelt functions work",{
  r1 <- read_pk_rep(path='../model_20_8', endyr=2020)
  c1 <- read_pk_cor(path='../model_20_8', endyr=2020)
  r2 <- r1; r2$version='none2'
  r <- list(r1,r2)
  x <- mymelt(r, 'Numbers_at_age')
  ## saveRDS(x,  "_expect_read_rep.rds")
  oldrep <- readRDS('_expect_read_rep.rds')
  expect_equal(x, oldrep)
  ## saveRDS(r,  "_expect_read_cor.rds")
  oldcor <- readRDS('_expect_read_cor.rds')
  expect_equal(r, oldcor)
})
