
#' Clean temporary files from a directory containing a model run
#'
#' @param path Path to directory to clean
#' @param full Whether to delete everything except .dat and .tpl files
#' @return Nothing
#' @export
clean_pk_dir <- function(path=getwd(), full=FALSE){
  stopifnot(is.logical(full))
  if(full){
    r <- list.files(path)
    r <- r[-grep(pattern='.tpl|.dat', r)]
    if(length(r)>0) trash <- file.remove(r)
  } else {
    r <- list.files(path, pattern='r0|b0|p0')
    r <- c(r,list.files(path, pattern='\\.cpp|\\.obj|\\.log|\\.eva|\\.bar|\\.htp|\\.dep'))
    if(length(r)>0) trash <- file.remove(file.path(path, r))
    s <- file.size(r <- file.path(path, 'mceval.dat'))
    if(!is.na(s)){
      if(s<.0001){
        trash <- file.remove(r)
      } else {
        message("Not removing mceval.dat becuase it appears to have results")
      }
    }
    r <- list.files(path, '\\.psv')
    s <- file.size(file.path(path,r))
    if(length(s)>0){
      if(s<.0001){
        trash <- file.remove(r)
      } else {
        message("Not removing .psv file becuase it appears to have results")
      }
    }
  }
}


#' Copy and build ADMB source from package
#'
#' @param path Path to folder
#' @param name Executable name (default goa_pk)
#' @param compile Whether to compile (default is TRUE)
#' @export
setup_exe <- function(path=getwd(), name='goa_pk', compile=TRUE){
  tpl <- 'C:/Users/cole.monnahan/GOApollock/source/goa_pk.tpl'
  stopifnot( file.exists(tpl))
  dir.exists(path)
  test <- file.copy(from=tpl, to=file.path(path,paste0(name, '.tpl')), overwrite=TRUE)
  if(!test) warning("Failed to copy file")
  if(compile){
    message("Compiling model..")
    old.wd <- getwd(); on.exit(setwd(old.wd))
    setwd(path)
    system(paste('admb', name))
  }
}
