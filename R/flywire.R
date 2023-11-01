#' @title Make a .nrrd file from neurons in a given flywire scene
#'
#' @description Creates a common image file (a .nrrd file) readable in Fiji from a \href{https://flywire.ai/}{flywire} scene. User must have access to the
#' flywire.ai production data set. See R package \href{https://natverse.org/fafbseg/}{fafbseg} scene for more information on accessing flywire data from R.
#'
#' @param scene a neuroglancer link to a flywire scene. You can get this by hitting the 'share' button in the flywire web browser GUI.
#' @param root_ids a vector of root IDs that specify flywire meshes.
#' @param reference the templatebrain to transform the flywire neuron(s) into before creating the .nrrd file. Defaults to \code{"JRC2018U_HR"}, which is
#' the brain space into which a lot of data for neuronbridge is formatted.
#' @param savefolder a folder in which to save .nrrd files. defaults to a folder called 'nrrd' that it will create in your working directory.
#' @param file the name of your .nrrd file. Defaults to the root IDs for the given flywire neurons, concatenated with the name of the reference brain space.
#'
#' @return a .nrrd file.
#'
#' @examples
#' \donttest{
#' # root_id_to_nrrd(root_ids = c("720575940632668823", "720575940644278051"))
#' # ascending neurons,
#' # if this link does not work for you,
#' # you do not yet have access to the data set
#' }
#' @seealso \code{\link{nrrd_to_mip}}
#' @export
ngl_scene_to_nrrd <- function(scene,
                                  reference="JRC2018U_HR",
                                  savefolder = file.path(getwd(),"nrrd"),
                                  file = NULL){
  # Get flywire IDs
  fw.decoded = fafbseg::ngl_decode_scene(scene)

  # Get meshes
  fw.meshes = fafbseg::read_cloudvolume_meshes(fw.decoded)

  #transform hemibrain neuron to template space
  mesh_to_nrrd(fw.meshes=fw.meshes,
               reference=reference,
               savefolder=savefolder,
               file=file)
}

#' @rdname ngl_scene_to_nrrd
#' @export
root_id_to_nrrd <- function(root_ids,
                            reference="JRC2018U_HR",
                            savefolder = file.path(getwd(),"nrrd"),
                            file = NULL){
  # Get meshes
  fw.meshes = fafbseg::read_cloudvolume_meshes(root_ids)

  #transform hemibrain neuron to template space
  mesh_to_nrrd(fw.meshes=fw.meshes,
               reference=reference,
               savefolder=savefolder,
               file=file)
}

# hidden
mesh_to_nrrd <- function(fw.meshes,
                         reference="JRC2018U_HR",
                         savefolder = file.path(getwd(),"nrrd"),
                         file = NULL){
  # Download requried registrations
  # nat.jrcbrains::download_saalfeldlab_registrations()
  fafbseg:::check_package_available('nat.flybrains')
  fafbseg:::check_package_available('nat.templatebrains')
  if(grepl("JRC2018",reference)){
    nat.jrcbrains::register_saalfeldlab_registrations()
  }
  if(reference == "JRC2018U_HR"){
    brain = "JRC2018U"
    janelia.mip = TRUE
  }else{
    janelia.mip = FALSE
    brain = reference
  }

  #transform hemibrain neuron to template space
  fw.ids = names(fw.meshes)
  fw.reg = nat.templatebrains::xform_brain(fw.meshes, reference=brain, sample="FAFB14")

  # make im3d
  x <- get(brain, pos = c("package:nat.flybrains"))
  if(janelia.mip){
    JRC2018U_HR = nat.templatebrains::templatebrain("JRC2018U_HR", dims = c(1210,566,174), voxdims = c(0.519,0.519,1), units = "microns")
    points=nat::xyzmatrix(fw.reg)
    I=nat::as.im3d(points,JRC2018U_HR)
  }else{
    points=nat::xyzmatrix(fw.reg)
    I=nat::as.im3d(points,x)
  }

  # write .nrrd
  dir.create(savefolder, showWarnings = FALSE)
  if(is.null(file)){
    file = paste0(paste0(fw.ids, collapse = '_'),"_",reference,".nrrd")
  }
  savefile = file.path(savefolder, file)
  nat::write.nrrd(I, savefile)
  message(paste0("saved: ", savefile))
}

#' @title Open Fiji on your computer and use it to convert your .nrrd files into color MIPs
#'
#' @description Open Fiji on your computer and use it to convert your .nrrd files into color MIPs for searches. You must have downloaded Fiji onto your machine and also have followed the instructions
#' \href{https://github.com/JaneliaSciComp/ColorMIP_Mask_Search}{here (search tool)} and \href{https://github.com/JaneliaSciComp/ColorMIP_Mask_Search/tree/master/ColorDepthMIP_Generator}{here (MIP generator)} to install the Fiji plugins
#' from Scientific Computing at Janelia Research Campus, for the Color MIP Search. When running the MIP generator, select directory that contains 3D stacks to be color depth MIP
#' then select directory for saving MIPs. Pay attention to the parameters you can then adjust. Skeleton bolded MIP can be useful for data coming from skeleton or even mesh sources.
#'
#' @param fiji.path path to your system's Fiji. If NULL, attempts to find using \code{getOption("jimpipeline.fiji")}.
#' @param macro path to the Fiji macro used to create the MIP file. Defaults to the macro contained within this R package: \code{Color_Depth_MIP_batch_0308_2021.ijm}.
#' @param MinMem lower memory limit
#' @param MaxMem upper memory limit
#'
#' @return a MIP files, saved as .tif, converted from a folder of selected .nrrd files
#'
#' @examples
#' \donttest{
#' nrrd_to_mip()
#' }
#' @seealso \code{\link{nrrd_to_mip}}
#' @export
nrrd_to_mip <- function(fiji.path = NULL,
                        macro = system.file(file.path("exdata","macros"),"Color_Depth_MIP_batch_0308_2021.ijm", package="neuronbridger"),
                        MinMem = MaxMem,
                        MaxMem = "2500m"){
  fiji.path = fiji(fijiPath  = fiji.path) # e.g. fiji.path = 'C:\\Program Files\\Fiji.app\\ImageJ-win64.exe'
  message('fiji: ', fiji.path)
  if(is.null(macro)){
    os <- get_os()
    if (os == "windows"){
      macro = normalizePath(file.path(dirname(fiji.path), 'plugins\\Macros\\Color_Depth_MIP_batch_0308_2021.ijm'))
    }else{
      macro = normalizePath('/Applications/Fiji.app/plugins/Macros/Color_Depth_MIP_batch_0308_2021.ijm')
    }
  }else{
    macro = normalizePath(macro)
  }
  runFijiMacro(
    macro = macro,
    macroArg = "",
    headless = FALSE,
    batch = FALSE,
    MinMem = MaxMem,
    MaxMem = "2500m",
    IncrementalGC = TRUE,
    Threads = NULL,
    fijiArgs = NULL,
    javaArgs = NULL,
    ijArgs = NULL,
    fijiPath = fiji.path,
    DryRun = FALSE
  )
}


