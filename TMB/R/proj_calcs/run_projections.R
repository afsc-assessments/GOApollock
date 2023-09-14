## Script to run the projections for various models
library(tidyverse)
library(GOApollock)

## Final 2022 model
replist <- readRDS('../results/repfile.RDS')
datlist <- readRDS("../results/datfile.RDS")

## New automated way using the new 'spm' package. Confirmed the
## same as manually doing it with 'proj' below.
write_spm_inputs(replist, datlist, 'spm_test2')
setwd('spm_test2')
system("spm")
setwd("..")
bf <- read.csv('spm_test2/spm_detail.csv')
(exec_table <- get_exec_table(replist,bf))
(exec_tableF <- format_exec_table(exec_table))
## write.csv(exec_tableF, file='../results/exec_tableF.csv', row.names=FALSE)
## Need last year's table too to make full exec table
replist0 <- read_pk_rep(path='../model_runs/m21_9_final/', endyr=2021)
bf0 <- read.table('data_2021_final/bigfile.out', header=TRUE)
exec_table0 <- get_exec_table(replist0,bf0)
## only need to do this in 2022
exec_table0[4,-1] <- 1000*c(430,430)
exec_table0[5,-1] <- 1000*c(172,172)
exec_table0[6,-1] <- 1000*c(150,150)
exec_table <- cbind(exec_table0[,-1], exec_table[,-1])
write.csv(exec_table, file='../results/exec_table.csv', row.names=FALSE)

proj_scens <- bf %>% select(-Spp) %>%
  group_by(Alternative,Yr)%>%
  summarize_all('mean') %>% arrange(Alternative,Yr) %>% ungroup
write.csv(proj_scens, file='../results/proj_scens.csv', row.names=FALSE)

## the old "proj" way done by hand from model outputs
bfproj <- read.table('data_2022_final/bigfile.out', header=TRUE)
tabspm <- get_proj_table(replist,bf)
tabproj <- get_proj_table(replist,bfproj)
cbind(tabspm, tabproj)
cbind(tabspm- tabproj) # super close, must be rounding thing?



## ### Old experiments for proj for the data additions
## dir.create('data_22_1')
## setwd("data_22_1/")
## file.copy('../main.exe', 'main.exe')
## datlist <- read_dat(filename='pk21_9.txt', path='../../model_runs/m22_1_add_sigmaR/')
## replist <- read_pk_rep(path='../../model_runs/m22_1_add_sigmaR/', endyr=2021)
## file.copy('../data_2021_final/tacpar.dat', 'tacpar.dat')
## write_proj_inputs(replist, datlist)
## system("main")
## setwd('..')

## dir.create('data_22_2')
## setwd("data_22_2/")
## file.copy('../main.exe', 'main.exe')
## datlist <- read_dat(filename='pk21_9.txt', path='../../model_runs/m22_2_no_q_prior/')
## replist <- read_pk_rep(path='../../model_runs/m22_2_no_q_prior/', endyr=2021)
## file.copy('../data_2021_final/tacpar.dat', 'tacpar.dat')
## write_proj_inputs(replist, datlist)
## system("main")
## setwd('..')

## dir.create('data_22_3')
## setwd("data_22_3/")
## file.copy('../main.exe', 'main.exe')
## datlist <- read_dat(filename='pk21_9.txt', path='../../model_runs/m22_3_sel4/')
## replist <- read_pk_rep(path='../../model_runs/m22_3_sel4/', endyr=2021)
## file.copy('../data_2021_final/tacpar.dat', 'tacpar.dat')
## write_proj_inputs(replist, datlist)
## system("main")
## setwd('..')

## dir.create('data_22_6')
## setwd("data_22_6/")
## file.copy('../main.exe', 'main.exe')
## datlist <- read_dat(filename='pk21_9.txt', path='../../model_runs/m22_6_est_AT_selex/')
## replist <- read_pk_rep(path='../../model_runs/m22_6_est_AT_selex/', endyr=2021)
## file.copy('../data_2021_final/tacpar.dat', 'tacpar.dat')
## write_proj_inputs(replist, datlist)
## system("main")
## setwd('..')


## ### test if can match manual way
## ## system("admb main")
## ## ## dir.create('test')
## ## ## file.copy('main.exe', 'test/main.exe')
## ## setwd("test")
## ## dat0 <- read_dat(filename='pk21_9.txt', path='../../model_runs/m21_9_final/')
## ## rep0 <- read_pk_rep(path='../../model_runs/m21_9_final/', endyr=2021)
## ## file.copy('../data_2021_final/tacpar.dat', 'tacpar.dat')
## ## write_proj_inputs(rep0, dat0)
## ## system("main")
## ## setwd('..')
## ## ## now check identical
## ## bf0 <- read_table('data_2021_final/bigfile.out') %>% as.data.frame()
## ## bf1 <- read_table('test/bigfile.out') %>% as.data.frame()
## ## all.equal(bf0,bf1) ## pretty close it must be a rounding thing


## ## ## Ingrids way of formatting
## ## ABC_tyF=formatC(ABC_ty,format="d",big.mark=",")
## ## ABC_nyF=formatC(ABC_ny,format="d",big.mark=",")
## ## OFL_tyF=formatC(OFL_ty,format="d",big.mark=",")
## ## OFL_nyF=formatC(OFL_ny,format="d",big.mark=",")
## ## TotBio_tyF=formatC(TotBio_ty,format="d",big.mark=",")
## ## TotBio_nyF=formatC(TotBio_ny,format="d",big.mark=",")
## ## FSB_tyF=formatC(FSB_ty,format="d",big.mark=",")
## ## FSB_nyF=formatC(FSB_ny,format="d",big.mark=",")
## ## B100_tyF=formatC(percentiles$SB0,format="d",big.mark=",")
## ## B40_tyF=formatC(percentiles$SB40,format="d",big.mark=",")
## ## B35_tyF=formatC(percentiles$SB35,format="d",big.mark=",")
## ## rbind(
## ##   c(FSB_tyF, FSB_nyF),
## ##   c(B100_tyF, B100_tyF),
## ##   c(B40_tyF, B40_tyF),
## ##   c(B35_tyF, B35_tyF),
## ##   c(FOFL_ty, FOFL_ny),
## ##   c(FABC_ty, FABC_ny),
## ##   c(FABC_ty, FABC_ny),
## ##   c(OFL_tyF, OFL_nyF),
## ##   c(ABC_tyF, ABC_nyF))

