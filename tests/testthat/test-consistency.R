
message("Recompiling and estimating model 20_8...")
setwd('..')
file.copy('../source/goa_pk.tpl', to='model_20_8/goa_pk.tpl')
setwd('model_20_8')
system("admb goa_pk")
system("goa_pk -nox -iprint 500")
clean_pk_dir()
myrep <- read_pk_rep('goa_pk.rep')
setwd('../testthat')

## saveRDS(myrep, '_expect_rep_20_8.rds')
oldrep <- readRDS('_expect_rep_20_8.rds')
expect_equal(myrep, oldrep)
