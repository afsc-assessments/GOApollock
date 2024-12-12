

# globals
addfig <- function(x){
  ff <- file.path('C:/Users/cole.monnahan/Work/assessments/GOA_pollock/2024/writeup/figures',x)
  if(!file.exists(ff)){
    stop(paste("File",ff, "does not exists, subbing in dummy one for now"))
    addfig('historical.png')
  } else {
    knitr::include_graphics(ff)
  }
}
table.path <-   'C:/Users/cole.monnahan/Work/assessments/GOA_pollock/2024/writeup/2024_GOApollock tables.xlsx'
gettab <- function(sheet, skip=1, maxcols=NULL, ...){
  p <- table.path
  stopifnot(file.exists(p))
  cap <- suppressMessages(readxl::read_xlsx(path=p, sheet=sheet, col_names=FALSE, range='A1', col_types = 'text'))
  cap <- as.character(cap[1,1]) # ditch merged cell stuff
  tab <- suppressMessages(readxl::read_xlsx(path=p, sheet=sheet,
                                            skip=skip, ...))
  if(!is.null(maxcols)) tab <- tab[,1:maxcols]
  return(list(cap=cap, tab=tab))
}

#species = params$species
year = 2024
ayears <- 1970:year
date = 'October 22, 2024'
model = "base"
end_proj = year + 15
best_f = 999 # from the best_f function in groundfishr package
spr0 <- 0.0761106
myres <- 'C:/Users/cole.monnahan/Work/assessments/GOA_pollock/2024/results'

## Read in the key output files from this year
finalfit <- readRDS(file.path(myres, 'fitfinal.RDS'))
oldfit <- readRDS(file.path(myres, '2023_tmbfit.RDS'))
repfile <- finalfit$rep
datfile <- finalfit$input$dat
stdfile <- finalfit$sd
## also need the old repfile for figures
repfile.last <- oldfit$rep

nov_abc_table <- read.csv(file.path(myres, 'nov_abc_table.csv'))

spm_detail <- readRDS(file.path(myres, 'spm_detail.RDS'))
proj_scens <- spm_detail %>% dplyr::select(-Stock) %>%
  dplyr::group_by(Alt,Year)%>%
  dplyr::summarize_all('mean') %>% dplyr::arrange(Alt,Year) %>%
  dplyr::ungroup()
cc <- dplyr::filter(proj_scens, Year %in% (year+0:2) & Alt==1) %>% dplyr::pull(Catch) %>% round(0)
c1 = tail(datfile$cattot,1) # estimated catch in current year
c2 = cc[1] # proj catch year + 1
c3 = cc[2] # proj catch year + 2

## This exec_table is constructed externally to have the old and
## new estimates
exec_table <- read.csv(file.path(myres, 'exec_table.csv'))
F <- function(x, digits=0) formatC(round(x,digits), format='d', big.mark=',')
Finv <- function(x) as.numeric(gsub(',','',x))
sumbio <- F(exec_table[2,3])
ssb <- F(exec_table[3,3])
b100 <- F(exec_table[4,3])
b40 <- F(exec_table[5,3])
b35 <- F(exec_table[6,3])
fofl <- exec_table[7,3]
fabc <- exec_table[8,3]
ofl <- F(exec_table[10,3])
maxabc <- F(exec_table[11,3])
abc <- F(exec_table[12,3])
pct.status <- paste0(format(round(100*Finv(ssb)/Finv(b100),1),nsmall=1), "%")

sumbio.last <- F(exec_table[2,1])
ssb.last <- F(exec_table[3,1])
b100.last <- F(exec_table[4,1])
b40.last <- F(exec_table[5,1])
b35.last <- F(exec_table[6,1])
ofl.last <- F(exec_table[10,1])
maxabc.last <- F(exec_table[11,1])
abc.last <- F(exec_table[12,1])
pct.b40.change <- (Finv(b40)-Finv(b40.last))/Finv(b40.last)
pct.b40.changeF <- (paste0(format(round(100*pct.b40.change,1),nsmall=1),ifelse(pct.b40.change>0, '% increase ', '% decrease ')))
pct.abc.change <- (Finv(abc)-Finv(abc.last))/Finv(abc.last)
pct.abc.changeF <- (paste0(format(round(100*pct.abc.change,1),nsmall=1),ifelse(pct.abc.change>0, '% increase ', '% decrease ')))

abcN <- Finv(abc)
abcN.last <- Finv(abc.last)

avgrec <- with(repfile, round(mean(recruit[years %in% 1978:(tail(years,1)-1)]),3))
CVrec <- with(repfile, round(sd(recruit[years %in% 1978:(tail(years,1)-1)])/avgrec,1))
avgrec.last <- with(repfile.last, round(mean(recruit[years %in% 1978:(tail(years,1)-1)]),3))
pct.avgrec.change <- round(100*(avgrec-avgrec.last)/avgrec.last,1)

# todo: build a function/switch to auto run the tier(a/b) - also for overfishing/overfished
