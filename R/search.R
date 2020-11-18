#' @title Get information on neurons and lines used by NeuronBridge
#'
#' @description This function gets information stored on \href{https://neuronbridge.janelia.org/}{NeuronBridge}
#' for *D. melanogaster* neurons and genetic driver lines that have been used in EM-LM searches.
#'
#' @param id character vector. An identifier/identifiers for the neuron(s)/line(s) you wish to search. This can either be a line code
#' for a GAL4 line (e.g. \href{https://v2.virtualflybrain.org/org.geppetto.frontend/geppetto?id=VFBexp_FBtp0057908&i=VFB_00017894}{R16F12})
#' or a split GAL4 line (e.g. \href{https://splitgal4.janelia.org/cgi-bin/view_splitgal4_imagery.cgi?line=LH173}{LH173})
#' or a 'body ID' for a neuron from the \href{https://neuprint.janelia.org/help/videos?dataset=hemibrain}{hemibrain connectome},
#' (e.g. \href{https://neuprint.janelia.org/view?bodyid=1537331894}{1537331894}).
#' @param dataset whether the ID given is a body ID (\code{"by_body"}) from the hemibrain connectome
#' or a genetic driver line (\code{"by_line"}). If left at \code{"detect"} then \code{neuronbridger} tries to guess to which \code{id} belongs by using
#' \code{neuronbridger:::guess_dataset}.
#' @param version the precomputed scores to search. For example, \code{"v2_1_1"} refers to this \href{https://neuronbridge.janelia.org/releasenotes/DATA}{release}.
#' @param ... methods passed to \code{neuronbridger:::neuronbridge_fetch}.
#'
#' @return a \code{data.frame} where each row corresponds to a MIP file that has been used to compute matching scores.
#' Each relates to the \code{id} you gave in the function call, but each also has its own unique \code{nb.id}.
#' \itemize{
#'
#'   \item{"publishedName"} { - the \code{id} you gave in the function call.}
#'
#'   \item{"libraryName"}{ - the data set from which this data item came.}
#'
#'   \item{"imageURL"}{ - the path on https://s3.amazonaws.com/, at which one can find a 'high res' .png of the MIP file for this data item.}
#'
#'   \item{"thumbnailURL"}{ - the path on https://s3.amazonaws.com/, at which one can find a 'low res' .jpg thumbnail of the MIP file for this data item.}
#'
#'   \item{"gender"}{ - the sex of the fly brain which this data item derives}
#'
#'   \item{"nb.id"}{ - the 'NeuronBridge ID' for the MIP file.}
#'
#'   \item{"no.images"}{ - the total number of MIP files (each with its own nb.id) that related to the \code{id} for which you searched.}
#'
#' }
#' @examples
#' \donttest{
#' \dontrun{
#' # Get information on a mushroombody related GMR GAL4 driver line
#' nb.info.lm = neuronbridge_info("MB543B")
#'
#' # Get information on the hemibrain neuron,
#' nb.info.em = neuronbridge_info("542634818")
#' }}
#' @seealso \code{\link{neuronbridge_hits}},
#'   \code{\link{neuronbridge_line_contents}},
#'   \code{\link{neuronbridge_mip}}
#' @export
neuronbridge_info <- function(id,
                              dataset = c("detect","by_line","by_body"),
                              version = "v2_1_1",
                              ...){
  id = as.character(id)
  dataset = match.arg(dataset)
  nb.df.ids = data.frame()
  rnames = c()
  pb <- progress::progress_bar$new(format = "  retrieving MIP information [:bar] :current/:total eta: :eta", total = length(id), clear = FALSE, show_after = 1)
  for(item in unique(id)){

    # check dataset
    if(dataset=="detect"){
      dataset = guess_dataset(item)
    }
    pb$tick()

    # Get path to ID info
    path=sprintf("https://janelia-neuronbridge-data-prod.s3.amazonaws.com/%s/metadata/%s/%s.json?x-id=GetObject",
                 version,
                 dataset,
                 item)
    nb=tryCatch(neuronbridge_fetch(path, parse = "json", ...), error = function(e) NULL)
    if(is.null(nb)){
      dataset = c("by_line","by_body")[!c("by_line","by_body")%in%dataset]
      warning("dataset field might be problematic, trying ", dataset)
      path=sprintf("https://janelia-neuronbridge-data-prod.s3.amazonaws.com/%s/metadata/%s/%s.json?x-id=GetObject",
                   version,
                   dataset,
                   id)
      nb=neuronbridge_fetch(path)
    }
    no = length(nb$results)
    nb.df = data.frame()
    for(r in nb$results){
      nbr = as.data.frame(t(unnest_df(r)), stringsAsFactors = FALSE)
      nb.df = plyr::rbind.fill(nb.df,nbr)
    }
    nb.df$nb.id = nb.df$id
    nb.df$id = NULL
    nb.df$no.images = no
    rnames = c(rnames, rep_name(nb.df$publishedName))
    nb.df.ids = plyr::rbind.fill(nb.df.ids,nb.df)
  }
  rnames[duplicated(rnames)] = rep_name(x=rnames[duplicated(rnames)], sep = "_")
  rownames(nb.df.ids) = rnames
  nb.df.ids
}

#' @title Get EM-LM neuron matches for fly brain neurons from neuronbridge.janelia.org
#'
#' @description The function \code{neuronbridge_hits} retrieves a ranked list of hits for a MIP file search.
#' Specify the MIP file you want to search for by using its \code{nb.id}. This can be found by
#' using \code{\link{neuronbridge_info}}. \code{neuronbridge_search} will run \code{neuronbridge_info} and then \code{neuronbridge_hits}, to make things easier.
#' The function \code{neuronbridge_search} can also handle being given multiple IDs, i.e. searching for multiple data items at once.
#' MIP options can be fetched and visualised with \code{\link{neuronbridge_mip}}.
#'
#' @param nb.id an internal ID used by \href{https://neuronbridge.janelia.org/}{neuronbridge.janelia.org} to identify neurons/lines. This
#' differs from \code{id} used by other functions, and can be found using \code{\link{neuronbridge_info}}.
#' @param threshold LM-EM matches with a \code{normalizedScore} below this value are not returned. If set to \code{NULL}, the results are not filtered.
#' @inherit neuronbridge_info params
#' @return a \code{data.frame} of hits. Each row idnciates a separate MIP file with its own \code{nb.id}. The \code{data.frame} is already ranked by \code{normalizedScore}.
#' Top scores (better match) are at the top of the data frame. The columns mean:
#' \itemize{
#'
#'   \item{"publishedName"} { - the \code{id} for the potential hit neuron/line. I.e. specififes a genetic driver resource or a connectome neuron.
#'   these are the same ids that can be seen with \code{\link{neuronbridge_ids}}.}
#'
#'   \item{"libraryName"}{ - the data set from which this data item came.}
#'
#'   \item{"imageURL"}{ - the path on https://s3.amazonaws.com/, at which one can find a 'high res' .png of the MIP file for this data item.}
#'
#'   \item{"thumbnailURL"}{ - the path on https://s3.amazonaws.com/, at which one can find a 'low res' .jpg thumbnail of the MIP file for this data item.}
#'
#'   \item{"slideCode"}{ - the unique identifier for the sample from which the MIP came. The first number indicates the date the image was taken by FlyLight.}
#'
#'   \item{"objective"}{ - the magnification under which the imasge was taken.}
#'
#'   \item{"gender"}{ - the sex of the fly brain which this data item derives. f = female, m = male.}
#'
#'   \item{"anatomicalArea"}{ - the gross part of the nervous system images, e.g. brain or ventral nervous system.}
#'
#'   \item{"alignmentSpace"}{ - the template brain to which the image that formed this MIP, was aligned.
#'   Typically, this is the JR2018 standard template brain from \href{https://doi.org/10.1101/376384}{Bogovic et al. 2018}.}
#'
#'   \item{"channel"}{ - number of the channel from the aligned image stack that is represented by this MIP.}
#'
#'   \item{"mountingProtocol"}{ - the protocol used to prepare brain sample for imaging.}
#'
#'   \item{"matchingPixels"}{ - the number of overlapping pixels between query (\code{searched.id}) and taregt (\code{nb.id}).}
#'
#'   \item{"gradientAreaGap "}{ - unsure, seeking clarification from NeuronBridge}
#'
#'   \item{"normalizedGapScore"}{ - unsure, seeking clarification from NeuronBridge}
#'
#'   \item{"normalizedScore"}{ - the matching score, created by examining the overlapped pixel number and color depth.
#'   If the color and xy position of the pixel match between the mask and the searching data,
#'   then the approach here will count it as a positive matching score}
#'
#'   \item{"searched.id"}{ - the \code{nb.id} you searched with, i.e. given to the function call}
#'
#'   \item{"nb.id"}{ - the 'NeuronBridge ID' for the MIP file.}
#'
#'}
#'
#' @examples
#'\donttest{
#' \dontrun{
#' nb.info.em = neuronbridge_info("542634818")
#' nb.hits = neuronbridge_hits(nb.info.em$nb.id[1])
#' # View(nb.hits[1:10,]) # see the top 10 hits
#'
#' # One can also do this in one with:
#' nb.search = neuronbridge_search("542634818")
#' ## However, this differs from the above in that every MIP
#' ## associated with the given id, is searched. For a connectome
#' ## neuron this is just 1, but for a GAL4 line with MCFO data it
#' ## can be many quite different images.
#'
#' # Note that nb.info.em here is actually an attribute you can see
#' nb.info.em = attr(nb.search,"search")
#' }}
#' @seealso \code{\link{neuronbridge_info}},
#'   \code{\link{neuronbridge_mip}},
#'   \code{\link{neuronbridge_ids}}
#' @export
neuronbridge_hits <- function(nb.id,
                              version = "v2_1_1"){
  # Get information on MIP hits
  nb.id = as.character(nb.id)
  nb.hits = data.frame()
  rnames = c()
  pb <- progress::progress_bar$new(format = "  finding hits for MIPs [:bar] :current/:total eta: :eta", total = length(nb.id), clear = FALSE, show_after = 1)
  for(n in nb.id){
    pb$tick()
    path2=sprintf("https://janelia-neuronbridge-data-prod.s3.amazonaws.com/%s/metadata/cdsresults/%s.json",
                  version,
                  n)
    nb2=neuronbridge_fetch(path2)
    nb.list = list()
    for(r in nb2$results){
      nbr = as.data.frame(t(unnest_df(r)), stringsAsFactors = FALSE)
      nb.list[[nbr$id]] = nbr
    }
    nb.df = do.call(plyr::rbind.fill, nb.list)
    nb.df$searched.id = n
    nb.df$nb.id = nb.df$id
    nb.df$id = NULL
    for(column in c("matchingPixels", "matchingRatio", "gradientAreaGap", "normalizedGapScore", "normalizedScore")){
      nb.df[,column] = as.numeric(unlist(nb.df[,column]))
    }
    nb.df$normalizedScore = as.numeric(nb.df$normalizedScore)
    nb.df = nb.df[order(nb.df$normalizedGapScore, decreasing = TRUE),]
    nb.hits = plyr::rbind.fill(nb.hits,nb.df)
    rnames = c(rnames, rep_name(nb.df$publishedName))
  }
  rnames[duplicated(rnames)] = rep_name(x=rnames[duplicated(rnames)], sep = "_")
  rownames(nb.hits) = rnames
  nb.hits
}

#' @rdname neuronbridge_hits
#' @export
neuronbridge_search <- function(id,
                                version = "v2_1_1",
                                dataset = c("detect","by_line","by_body"),
                                threshold = 10000){
  # check dataset
  id = as.character(id)
  dataset = match.arg(dataset)
  if(dataset=="detect"){
    dataset = guess_dataset(id)
  }

  # Get path to ID info
  nb.info=neuronbridge_info(id=id, dataset = dataset, version = version)

  # Get information on MIP hits
  nb.hits.all = list()
  rnames = c()
  pb <- progress::progress_bar$new(format = "  finding hits for MIPs [:bar] :current/:total eta: :eta", total =  nrow(nb.info), clear = FALSE, show_after = 1)
  for(i in 1:nrow(nb.info)){
    pb$tick()
    nb.id = nb.info[i,"nb.id"]
    the.id = nb.info[i,"publishedName"]
    nb.hits = neuronbridge_hits(nb.id = nb.id, version = version)
    nb.hits$searched = the.id
    if(!is.null(threshold)){
      nb.hits = subset(nb.hits, nb.hits$normalizedScore >= threshold)
    }
    rnames = c(rnames, rownames(nb.hits))
    nb.hits.all[[i]] = nb.hits
  }
  nb.hits.all = do.call(plyr::rbind.fill, nb.hits.all)
  rnames[duplicated(rnames)] = rep_name(x=rnames[duplicated(rnames)], sep = "_")
  rownames(nb.hits.all) = rnames

  # Attach some meta data
  attr(nb.hits.all,"search") = nb.info
  class(nb.hits.all) = c("neuronbridge", class(nb.hits.all))

  # return
  nb.hits.all
}
