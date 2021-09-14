
#' Clean temporary files from a directory containing a model run
#'
#' @param path
#'
clean_pk_dir <- function(path){
  r <- list.files(path, pattern='r0|b0|p0')
  file.remove(r)
  r <- file.path(path, paste0('goa_pk', c('.obj', '.log', '.eva',
                                    '.bar', '.std', '.htp', '.cpp')))
  file.remove(r)
  file.remove(file.path(path, 'admodel.dep'))
  s <- file.size(r <- file.path(path, 'mceval.dat'))
  if(s<.001){
    file.remove(r)
  } else {
   message("Not removing mceval.dat becuase it appears to have results")
  }
}
