#' Calculate ADSB (avg diff in SSB) for use in determining model
#' names. If ADSB>0.1 then it's considered a major version.
#'
#' @param base The base model std df or rep list
#' @param new The propsed model std df or rep list. Can have more
#'   years than base model
#' @return ADSB value and message
#' @export
calc_adsb <- function(base, new){
  if(is_tibble(base)){
    ## tibble from read_pk_cor
    yr0 <- filter(base, year >= 1977 & name=='Espawnbio') %>% pull('year')
    sb0 <- filter(base, year >= 1977 & name=='Espawnbio') %>% pull('est')
  } else {
    ## replist from read_pk_rep
    yr0 <- base$years
    sb0 <- base$Expected_spawning_biomass
    stopifnot(length(yr0)==length(sb0))
    ind <- which(yr0>=1977)
    yr0 <- yr0[ind]; sb0 <- sb0[ind]
  }
  if(is_tibble(new)){
    sb1 <- filter(new, year >= 1977 & name=='Espawnbio') %>% pull('est')
  } else {
    yr1 <- new$years
    sb1 <- new$Expected_spawning_biomass
    stopifnot(length(yr1)==length(sb1))
    ind <- which(yr1>=1977)
    yr1 <- yr1[ind]; sb1 <- sb1[ind]
  }
  stopifnot(length(sb0)>=length(sb1))
  adsb <- sqrt(mean((sb1[1:length(sb0)]/sb0-1)^2))
  message('This is a ', ifelse(adsb>.1, '**major**',
                               '**minor**'), ' model version (ADSB=',round(adsb,2),')')
  return(adsb)
}
