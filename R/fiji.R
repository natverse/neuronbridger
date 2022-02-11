# hidden, from jimpipeline
fiji <- function (fijiPath = NULL)
{
  if (!is.null(fijiPath)) {
    if (!file.exists(fijiPath))
      stop("fiji is not at: ", fijiPath)
  }
  else {
    fijiPath = getOption("jimpipeline.fiji")
    if (!is.null(fijiPath)) {
      if (!file.exists(fijiPath))
        stop("fiji is not at: ", fijiPath, " as specified by options('jimpipeline.fiji')!")
    }
    else {
      if (!nzchar(fijiPath <- Sys.which("fiji"))) {
        macapp = "/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx"
        if (file.exists(macapp))
          fijiPath = macapp
        else stop("Unable to find fiji!", "Set options('jimpipeline.fiji') to point to the fiji command line executable!")
      }
    }
  }
  options(jimpipeline.fiji = fijiPath)
  normalizePath(fijiPath)
}

# hidden, from jimpipeline
runFijiMacro <- function (macro = "", macroArg = "", headless = FALSE, batch = TRUE,
                          MinMem = MaxMem, MaxMem = "2500m", IncrementalGC = TRUE,
                          Threads = NULL, fijiArgs = NULL, javaArgs = NULL, ijArgs = NULL,
                          fijiPath = fiji(), DryRun = FALSE) {
  os <- get_os()
  if (os == "windows"){
    fijiPath <- paste0('"',fijiPath,'"')
  }
  if(is.null(Threads)|macroArg==""){
    simple = TRUE
  }else{
    simple = FALSE
  }
  if (headless)
    fijiArgs = c(fijiArgs, "--headless")
  fijiArgs = paste(fijiArgs, collapse = " ")
  javaArgs = c(paste("-Xms", MinMem, sep = ""), paste("-Xmx",
                                                      MaxMem, sep = ""), javaArgs)
  if (IncrementalGC)
    javaArgs = c(javaArgs, "-Xincgc")
  javaArgs = paste(javaArgs, collapse = " ")
  threadAdjust = ifelse(is.null(Threads), "", paste("run(\"Memory & Threads...\", \"parallel=",
                                                    Threads, "\");", sep = ""))
  ijArgs = paste(c(ijArgs, ifelse(batch, "-batch", "")), collapse = " ")
  if(simple){
    macroCall = sprintf(' -macro "%s"', macro)
    cmd <- paste(fijiPath, javaArgs, fijiArgs, macroCall,
                 ijArgs)
  }else{macro
    macroCall = paste(" -eval '", threadAdjust, "runMacro(\"",
                      macro, "\",\"", macroArg, "\");' ", sep = "")
    cmd <- paste(fijiPath, javaArgs, fijiArgs, "--", macroCall,
                 ijArgs)
  }
  if (DryRun)
    return(cmd)
  return(0 == system(cmd))
}

# hidden
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin")
      os <- "osx"
  }
  else {
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
