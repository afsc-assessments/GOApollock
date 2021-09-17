
#' Clean temporary files from a directory containing a model run
#'
#' @param path Path to directory to clean
#' @return Nothing
clean_pk_dir <- function(path=getwd()){
  r <- list.files(path, pattern='r0|b0|p0')
  r <- c(r,list.files(path, pattern='\\.cpp|\\.obj|\\.log|\\.eva|\\.bar|\\.std|\\.htp|\\.dep'))
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
      file.remove(r)
    } else {
      message("Not removing .psv file becuase it appears to have results")
    }
  }
}
