library(wham)
# https://timjmiller.github.io/wham/articles/ex4_selectivity.html

wham.dir <- find.package("wham")
file.path(wham.dir, "example_scripts")
source(file.path(wham.dir, "example_scripts", "ex4_selectivity.R"))

write.dir <- "TMB/Tests/WHAM" # otherwise will be saved in working directory
dir.create(write.dir)
setwd(write.dir)

# Get data
file.copy(from=file.path(wham.dir,"extdata","ex1_SNEMAYT.dat"), to=write.dir, overwrite=FALSE)
asap3 <- read_asap3_dat("ex1_SNEMAYT.dat")

# Specify models
# m1-m5 logistic, m6-m9 age-specific
sel_model <- c(rep("age-specific",5))

# time-varying options for each of 3 blocks (b1 = fleet, b2-3 = indices)
sel_re <- list(c("none","none","none"), # m1-m5 age-specific
               c("iid","none","none"),
               c("ar1","none","none"),
               c("ar1_y","none","none"),
               c("2dar1","none","none"))
n.mods <- length(sel_re)

# summary data frame
df.mods <- data.frame(Model=paste0("m",1:n.mods),
                      Selectivity=sel_model, # Selectivity model (same for all blocks)
                      Block1_re=sapply(sel_re, function(x) x[[1]])) # Block 1 random effects
rownames(df.mods) <- NULL

df.mods


# Fit models
mod_list <- list()
for(m in 1:n.mods){

  # Often need to fix selectivity = 1 for at least one age per age-specific block: ages 4-6 / 4 / 3-6
  input <- prepare_wham_input(asap3, model_name=paste(paste0("Model ",m), sel_model[m], paste(sel_re[[m]], collapse="-"), sep=": "), recruit_model=2,
                              selectivity=list(model=rep("age-specific",3), re=sel_re[[m]],
                                               initial_pars=list(c(0.1,0.5,0.5,1,1,1),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,0.5,1,1,1,1)),
                                               fix_pars=list(4:6,4,3:6)),
                              NAA_re = list(sigma='rec+1',cor='iid'),
                              age_comp = "logistic-normal-miss0") # logistic normal, treat 0 obs as missing


  # fit model
  mod_list[[m]] <- fit_wham(input, do.check=T, do.osa=F, do.proj=F, do.retro=F, MakeADFun.silent = TRUE)
  # saveRDS(mod_list[[m]], file=paste0("TMB/Tests/WHAM/m",m,".rds"))
}

mod <- 3 # Look at AR1 on age (4 is AR1 on year, 5 is 2DAR1)
selpars_re <- mod_list[[mod]]$parList$selpars_re # Random effects
selpars_re
length(unique(selpars_re)) # 3 random effects
mod_list[[mod]]$parList$logit_selpars[1,]  # 3 fixed effects
mod_list[[mod]]$parList$sel_repars[1,] # var and rho of AR1

logitfun <- function(x) {1/(1+exp(- x))} # Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1)

rep <- mod_list[[mod]]$report()
head(rep$selpars[[mod]]) # Sel pars (same as sel)
logitfun(mod_list[[3]]$parList$logit_selpars[1,1:6] + c(unique(mod_list[[3]]$parList$selpars_re), rep(0,3)))
head(rep$selAA[[mod]]) # Selectivity
