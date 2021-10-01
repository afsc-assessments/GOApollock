
test_that("model 20_8 is unchanged",{
message("\nRecompiling and estimating model 20_8...\n")
setwd('..')
file.copy('../source/goa_pk.tpl', to='model_20_8/goa_pk.tpl', overwrite = TRUE)
setwd('model_20_8')
system("admb goa_pk", ignore.stdout = TRUE)
system("goa_pk -nox -nohess -iprint 500", ignore.stdout = TRUE)
clean_pk_dir()
setwd('../testthat')
})
