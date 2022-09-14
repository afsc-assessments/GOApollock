library(testthat)
library(GOApollock)

message("\nRecompiling  model 20_8...\n")
setwd('..')
file.copy('../source/goa_pk.tpl', to='model_20_8/goa_pk.tpl', overwrite = TRUE)
setwd('model_20_8')
clean_pk_dir(full=TRUE)
system("admb goa_pk", ignore.stdout = TRUE)

print(getwd())
setwd('..')
test_check("GOApollock")
