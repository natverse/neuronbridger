.onAttach <- function(libname, pkgname){
  db = file.path(getwd(),"neuronbridger")
  options(neuronbridger=db)
}
