#' @title Get identifiers for neurons/lines with NeuronBridge results
#'
#' @description This function gets all of the neuron identifiers or line codes (e.g.R24A02 is a generation 1 GMR GLA4 line, LH173 is a split-GAL4 line, 5813055865 is a hemibrain neuron)
#' for published neurons/lines for which NeuronBridge has precomputed matching results. If you query is not among these identifiers you will not be able to get results for it.
#'
#' @inherit neuronbridge_info params
#'
#' @return a \code{vector} of identifiers
#'
#' @examples
#' \donttest{
#' all.nb.ids = neuronbridge_ids()
#' length(all.nb.ids)
#' }
#' @seealso \code{\link{neuronbridge_hits}},
#'   \code{\link{neuronbridge_info}},
#'   \code{\link{neuronbridge_search}}
#' @export
neuronbridge_ids <- function(version = "v2_1_1"){
  path=sprintf("https://janelia-neuronbridge-data-prod.s3.amazonaws.com/%s/publishedNames.txt?x-id=GetObject",version)
  nb=neuronbridge_fetch(path, parse = "none")
  text=httr::content(nb, as = "text", encoding = "UTF-8")
  ids=unlist(strsplit(text, split = "\n"))
  ids
}

# hidden
guess_dataset <- function(id){
  all.ids = neuronbridge_ids()
  good.id = intersect(id,all.ids)
  if(length(good.id)!=length(id)){
    bad.id = setdiff(id,all.ids)
    warning("The given ID(s): ", paste(bad.id,collapse=", "),"  not in the neuronbridge database")
  }
  if(!length(good.id)||is.na(good.id)||good.id==0){
    stop("Given ID not in the neuronbridge database")
  }
  line = suppressWarnings(is.na(as.integer(good.id)))[1]
  ifelse(line,"by_line","by_body")
}



