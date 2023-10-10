
## Get MLE parameters into list form for TMB
## mle <- adnuts:::.read_mle_fit('goa_pk', path='C:/Users/cole.monnahan/Work/assessments/GOA_pollock/2022/model_runs/m22_13_2022_final')
## parvec <- mle$est[1:mle$nopar]
## parnames <- mle$par.names[1:mle$nopar]
## parnames <- lapply(strsplit(parnames, split='\\['), function(x) x[1]) |> unlist()
## names(parvec) <- parnames
## active_pars <- split(unname(parvec),names(parvec))
## saveRDS(active_pars, 'active_pars.RDS')
active_pars <- readRDS("active_pars.RDS")

## Have to initialize inactive pars and then map off
nyrs <- 53
## initialize some parameters
some_pars <- list(dev_log_initN=rep(0,10),
                 q4_pow=0, q5_pow=0, mean_log_recruit=0.0,
                 mean_log_F=-1.6, log_q1_mean=0.0,
                 log_q1_dev=rep(0,nyrs), log_q2_mean=0.0,
                 log_q2_dev=rep(0,nyrs), log_q3_mean=-1.6,
                 log_q3_dev=rep(0,nyrs), log_q6=0.0,
                 slp1_fsh_dev=rep(0,nyrs),
                 inf1_fsh_dev=rep(0,nyrs),
                 slp2_fsh_dev=rep(0,nyrs),
                 inf2_fsh_dev=rep(0,nyrs),
                 log_slp1_fsh_mean=1.0,
                 inf1_fsh_mean=4.0, log_slp2_fsh_mean=1.0,
                 inf2_fsh_mean=8.0, log_slp2_srv1=0.5,
                 inf2_srv1=9.0, log_slp1_srv2=-0.8,
                 inf1_srv2=4.0, log_slp2_srv2=1, inf2_srv2=20,
                 log_slp1_srv3=0.0, inf1_srv3=5.0,
                 log_slp1_srv6=4.9, inf1_srv6=0.5,
                 log_slp2_srv6=1, inf2_srv6=20, sigmaR=1.3,
                 natMscalar=1)
inactive_pars <- some_pars[!names(some_pars) %in% names(active_pars)]
map <- lapply(inactive_pars, function(x) as.factor(x*NA))
## turn off projection stuff for now (broken in cpp)
map$log_recr_proj <- as.factor(rep(NA,5))
## turn off the "means" which are really the initial value for
## the RW. Needed b/c ADMB had dev_vectors and I dropped those
## and added a df so take it away. Fixing them at arbitrary value
## (previous MLE).
map$mean_log_F <- map$log_q3_mean <- map$log_q1_mean <- factor(NA)
map$log_slp1_fsh_mean <- map$inf1_fsh_mean <- factor(NA)
pars <- c(active_pars, inactive_pars)

## for some strnage reason I had to manually reorder this to
## match the template, otherwise it was totally wrong (??!!)
ind <- c("dev_log_initN", "mean_log_recruit", "dev_log_recruit",
         "sigmaR", "log_recr_proj", "log_slp1_fsh_mean",
         "inf1_fsh_mean", "log_slp2_fsh_mean", "inf2_fsh_mean",
         "slp1_fsh_dev", "inf1_fsh_dev", "slp2_fsh_dev",
         "inf2_fsh_dev", "log_slp2_srv1", "inf2_srv1",
         "log_slp1_srv2", "inf1_srv2", "log_slp2_srv2",
         "inf2_srv2", "log_slp1_srv3", "inf1_srv3",
         "log_slp1_srv6", "inf1_srv6", "log_slp2_srv6",
         "inf2_srv6", "mean_log_F", "dev_log_F",
         "log_q1_mean", "log_q1_dev", "log_q2_mean",
         "log_q2_dev", "log_q3_mean", "log_q3_dev",
         "log_q4", "q4_pow", "log_q5", "q5_pow",
         "log_q6", "natMscalar" )
pars <- pars[ind]
saveRDS(pars, 'pars.RDS')
saveRDS(map, 'map.RDS')
