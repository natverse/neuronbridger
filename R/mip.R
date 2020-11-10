#' @title Get and plot colour MIPs for neurons fron NeuronBridge
#'
#' @description These functions allow one to get (\code{neuronbridge_mip}) and plot (\code{plot_mip}) colour MIPs (maximum projection images).
#' MIPs are not 'stacks' of images, but a single image plane where colour encodes depth. The bluer the hue, the more anterior (front of the brain) the signal.
#' The redder the hue, the more posterior (back fo the brain). There may be differences/difficulties plotting images between different operating systems.
#' Only extensively tested on MacOS Catalina.
#'
#' @inherit neuronbridge_info params
#' @param mip either a full path to \code{.png} file or an array describing a image file, e.g. as comes from \code{png::readPNG}.
#' @param images.url the online database from which NeuronBridge gets its MIP files.
#' Higher res .png files come from https://s3.amazonaws.com/janelia-flylight-color-depth.
#' Lower res .jpg thumbnails come from https://s3.amazonaws.com/janelia-flylight-color-depth-thumbnails.
#' @param db the directory in which image files are saved (\code{neuronbridge_mip}), and from which they are read (\code{plot_mip}). This defaults to a temporary
#' directory: \code{file.path(tempdir(),"neuronbridger")}.
#' @param file.delete logical, whether or not to delete the saved \code{.png} file after plotting.
#' @param mips either a data frame of 'hits' (e.g. as porduced by \code{\link{neuronbridge_search}}), a vector of file paths to saved MIP .png files
#' or ids for data items (i.e. from \code{\link{neuronbridge_ids}}).
#' @param no.hits the number of hits to visualise if \code{type="hits"}.
#' @param type the type of argument given as \code{mips}.
#' @param sleep defaults to \code{NULL}. User needs to press a key to move to the nex MIP. If a numeric value, then MIPs progress automatically and
#' a pause of \code{sleep} seconds occurs before the next MIP is shown.
#' @param interactive logical. When using \code{neuronbridge_remove_mips}, whether or not to be accosted with an 'are you sure'?.
#'
#' @details By default MIPs are saved to your R session's temporary directory. However, you cna specify an alternative directory (argument: \code{db}).
#' MIPs are saved as \code{.png} files, where the file names are given as \code{id_number.png}.
#'
#' FlyLight file names contain metadata as follows:
#' \code{'[Publishing Name]-[Slide Code]-[Driver]-[Gender]-[Objective]-[Area/Tile]-[Alignment Space]-CDM_[Channel].png'}
#'
#' Find information on these meta data entries by examining \code{?neuronbridge_hits}.
#'
#' @return a named list of \code{.png} (\code{neuronbridge_mip}) or visualisation of a \code{.png} MIP in an \code{rgl} display (\code{plot_mip}).
#'
#' @examples
#'
#' # Let us now see the related MIP file
#' mip = neuronbridge_mip("542634818")
#' ## This gets every MIP file associated with id
#' ## In this case, just one
#'
#' # Plot the MIP image in an rgl viewer
#' plot_mip(mip)
#'
#' @seealso \code{\link{neuronbridge_ids}},
#'   \code{\link{neuronbridge_info}},
#'   \code{\link{neuronbridge_search}}
#' @export
neuronbridge_mip <- function(id,
                             dataset = c("detect","by_line","by_body"),
                             images.url = "https://s3.amazonaws.com/janelia-flylight-color-depth",
                             db = options("neuronbridger"),
                             version = "v2_1_1"){
  # check dataset
  dataset=match.arg(dataset)

  # Get path to ID info
  nb.info=neuronbridge_info(id=id, dataset = dataset)

  # Get images
  img.list=get_mip(nb.info, images.url=images.url,db=db)

  # return
  img.list
}

#' @rdname neuronbridge_mip
#' @export
plot_mip <- function(mip,
                     db = options("neuronbridger"),
                     file.delete = FALSE){
  if(!dir.exists(unlist(db)) ){
    dir.create(unlist(db))
  }
  if(is.list(mip)){
    if(length(mip)>1){
      warning("Multiple MIPs given, selecting only the first")
    }
    nam = names(mip)[1]
    if(is.null(nam)){
      nam = "temp"
      file.delete = TRUE
    }else{
      if(grepl("_",nam)){
        nam = gsub("\\_[0-9]+\\.png",".png",nam)
      }
    }
    mip = mip[[1]]
  }
  if(grepl("\\.png|\\_",mip[1])){
    temp = mip
  }else{
    temp = file.path(db,paste0(nam,".png"))
  }
  if(length(mip)){
    if(!file.exists(temp)){
      png::writePNG(image = mip, target = temp)
    }
    rgl::bg3d(texture = temp, col = "white")
    if(file.delete){
      rm = file.remove(temp)
    }
  }
}

#' @rdname neuronbridge_mip
#' @export
scan_mip <- function(mips,
                     no.hits = 10,
                     type = c("hits","files","ids","images"),
                     sleep = NULL,
                     images.url = "https://s3.amazonaws.com/janelia-flylight-color-depth",
                     db = options("neuronbridger")){
  type = match.arg(type)

  # Get images
  if(type=="hits"){
    mips = mips[1:no.hits,]
    img.list = get_mip(mips, images.url=images.url,db=db)
  }else if(type=="ids"){
    img.list = neuronbridge_mip(mips)
  }else {
    img.list = mips
  }

  # Scan through images
  for(i in 1:length(img.list)){
    img = img.list[i]
    if(type=="hits"){
      message("name: ", mips[i,"publishedName"])
      message("score: ",mips[i,"normalizedScore"])
    }
    rgl::clear3d()
    plot_mip(img)
    if(!is.null(sleep)){
      Sys.sleep(sleep)
    }else{
      p = readline("Hit any key to continue ")
    }
  }

}

# hidden
get_mip <- function(nb.info,
                    images.url = "https://s3.amazonaws.com/janelia-flylight-color-depth",
                    db = options("neuronbridger")){
  if(!dir.exists(unlist(db)) ){
    dir.create(unlist(db))
  }
  img.list = list()
  pb <- progress::progress_bar$new(format = "  downloading [:bar] :current/:total eta: :eta", total = nrow(nb.info), clear = FALSE, show_after = 1)
  for(i in 1:nrow(nb.info)){
    pb$tick()
    img.path = file.path(images.url,nb.info[i,"imageURL"])
    extension = ifelse(grepl("\\.png",img.path),".png",".jpg")
    nam = rownames(nb.info)[i]
    if(grepl("_",nam)){
      nam = gsub("\\_[0-9]+\\.png",".png",nam)
    }
    temp = file.path(db,paste0(nam,extension))
    if(!file.exists(temp)){
      tryCatch(utils::download.file(url = img.path, destfile = temp, quiet = TRUE),
               error = function(e){
                 warning(as.character(e))
                 next
                 })
    }
    if(extension==".png"){
      img = png::readPNG(temp)
    }else{
      img = jpeg::readJPEG(temp)
    }
    names(img) = rownames(nb.info)[i]
    img.list[[rownames(nb.info)[i]]] = img
  }
  img.list
}

# remove MIPs
#' @rdname neuronbridge_mip
#' @export
neuronbridge_remove_mips <- function(db = options("neuronbridger"),
                                     interactive = TRUE){
  files = list.files(as.character(db), path = "\\.png|\\.jpg", full.names = TRUE)
  if(interactive){
    p = "unknown"
    while(p%in%c("y","yes","n","no")){
      p = readline(prompt =   paste0("Do you with to remove ", length(files), " image files in ", db, " ?"))
    }
    if(p%in%c("yes","y")){
      file.remove(files)
    }
  }else{
    file.remove(files)
  }
}




