#' @title Get LM-EM matches for one set of neurons while avoiding another set
#'
#' @description This function uses \code{\link{neuronbridge_search}} to get potential matches for data items (i.e. neurons/lines) that you wish to search
#' (\code{search}). It also fetches the same information on data items you would rather not have in your 'hits', i.e.
#' hemibrain connectome neurons that you do not want to 'also'-be in the lines that come up when you search for you neurons in \code{search}. Depending
#' on your use-case, you may want to read a series of table with \code{\link{neuronbridge_search}} rather than use this function, to filter your matches.
#' See details below.
#'
#' @param search data item IDs (in \code{neuronbridge_ids}) for which you want to find hits.
#' @param avoid data item IDs (in \code{neuronbridge_ids}) that you would rather not have in your hits.
#' @param threshold LM-EM matches with a \code{normalizedScore} below this value are not returned.
#'
#' @inherit neuronbridge_info params
#'
#' @return a \code{data.frame} of hits. Each row indicates a separate MIP file with its own \code{nb.id}. The \code{data.frame} is already ranked by \code{normalizedScore}.
#' Top scores (better match) are at the top of the data frame. The columns mean:
#' \itemize{
#'
#'   \item{"publishedName"} { - the \code{id} for the potential hit neuron/line. I.e. specifies a genetic driver resource or a connectome neuron.
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
#'   \item{"objective"}{ - the magnification under which the image was taken.}
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
#'   \item{"matchingPixels"}{ - the number of overlapping pixels between query (\code{searched.id}) and target (\code{nb.id}).}
#'
#'   \item{"gradientAreaGap "}{ - unsure, seeking clarification from NeuronBridge}
#'
#'   \item{"normalizedGapScore"}{ - unsure, seeking clarification from NeuronBridge}
#'
#'   \item{"normalizedScore"}{ - the matching score, created by examining the overlapped pixel number and colourdepth.
#'   If the colourand xy position of the pixel match between the mask and the searching data,
#'   then the approach here will count it as a positive matching score}
#'
#'   \item{"searched.id"}{ - the \code{nb.id} you searched with, i.e. given to the function call}
#'
#'   \item{"nb.id"}{ - the 'NeuronBridge ID' for the MIP file.}
#'
#'}
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # Get helpful package with some classe hemibrain neuron IDs
#' if (!require("hemibrainr")) remotes::install_github("flyconnectome/hemibrainr")
#'
#' # Load package
#' library(hemibrainr)
#'
#' # So this is the 'olfactory PN' neuron we want
#' search = "542634818"
#'
#' # And we do not want these other PNs
#' avoid = setdiff(hemibrainr::pn.ids,search)
#'
#' # Let's see what we get for ot
#' hits = neuronbridge_avoid(search = search, avoid = avoid)
#' # View(hits[1:20,])
#'
#' }}
#' @seealso \code{\link{neuronbridge_info}},
#'   \code{\link{neuronbridge_mip}},
#'   \code{\link{neuronbridge_search}}
#' @export
neuronbridge_avoid <- function(search,
                         avoid,
                         version = "v2_1_1",
                         threshold = 23000){

  # First run a search for the favoured IDs
  search.search = neuronbridge_search(search, version = version, threshold = threshold)

  # And one for the unfavoured ones
  avoid.search = neuronbridge_search(avoid, version = version, threshold = threshold)

  # Aggregate
  search.filtered = stats::aggregate(list(normalizedScore=search.search$normalizedScore),
                                          list(publishedName=search.search$publishedName,
                                               searched=search.search$searched,
                                               libraryName=search.search$libraryName),
                                          FUN = function(x) max(x, na.rm = TRUE))
  avoid.filtered = stats::aggregate(list(normalizedScore=avoid.search$normalizedScore),
                                     list(
                                       avoid=avoid.search$searched,
                                       publishedName=avoid.search$publishedName,
                                       libraryName=avoid.search$libraryName),
                                     FUN = function(x) max(x, na.rm = TRUE))

  # Get the number of 'avoid' neurons
  avoid.count= stats::aggregate(list(avoid.count=avoid.filtered$avoid),
                                    list(publishedName=avoid.filtered$publishedName),
                                    FUN = function(x) length(unique(x)))
  avoid.contents = stats::aggregate(list(avoid.contents=avoid.filtered$avoid),
                                list(publishedName=avoid.filtered$publishedName),
                                FUN = function(x) paste(sort(unique(x)),collapse=", "))

  # Combine
  search.results = merge(search.filtered, avoid.count)
  search.final = merge(search.results, avoid.contents)
  search.final = search.final[order(search.final$avoid.count, decreasing = TRUE),]
  search.final = search.final[order(search.final$normalizedScore, decreasing = TRUE),]
  search.final

}


#' @title Use NeuronBridge to estimate which hemibrain connectome neurons are taregted by a genetic driver line
#'
#' @description This function uses \code{\link{neuronbridge_search}} to get potential connectime matches for the given genetic driver line (\code{line}).
#' Each line may have multiple images associated with it, likely because it has undergone stochastic labelling to better see what the line contains. Each
#' of these images is searched and potential connectome hits returned. Each is assigned its highest score. The final output is filtered to that the LM-EM matching score
#' is below the given \code{threshold}.
#'
#' @param line a single ID for a genetic driver line (in \code{neuronbridge_ids}) for which you want to find hits.
#' @param line1 a line you want to try combining with \code{line2}. The aim is to see which neurons are in both \code{line1} and \code{line2}
#' @param line2 a line you want to try combining with \code{line1}. The aim is to see which neurons are in both \code{line1} and \code{line2}
#' @param neuprintr logical, whether or not to use the package \href{https://github.com/natverse/neuprintr}{neuprintr}
#' to fetch meta data (e.g. cell body fiber tract, cell type) on the connectome neurons returned. If \code{TRUE} then
#' \code{neuprintr::neuprint_get_meta} is used, and additional columns from this call are added to the returned \code{data.frame}.
#' @param threshold LM-EM matches with a \code{normalizedScore} below this value are not returned.
#' @param threshold1 LM-EM matches with a \code{normalizedScore} below this value are not returned. Applied to 'hits' for \code{line1}.
#' @param threshold2 LM-EM matches with a \code{normalizedScore} below this value are not returned. Applied to 'hits' for \code{line2}.
#'
#' @inherit neuronbridge_info params
#'
#' @return a \code{data.frame} of hits. Each row idnciates a separate MIP file with its own \code{nb.id}. The \code{data.frame} is already ranked by \code{normalizedScore}.
#' Top scores (better match) are at the top of the data frame. The columns mean:
#' \itemize{
#'
#'   \item{"publishedName"} { - the \code{id} for the potential hit neuron/line. I.e. specififes a genetic driver resource or a connectome neuron.
#'   these are the same ids that can be seen with \code{\link{neuronbridge_ids}}.}
#'
#'   \item{"libraryName"}{ - the data set from which this data item came.}
#'#'
#'   \item{"normalizedScore"}{ - the matching score, created by examining the overlapped pixel number and colourdepth.
#'   If the colourand xy position of the pixel match between the mask and the searching data,
#'   then the approach here will count it as a positive matching score}
#'
#'   \item{"searched"}{ - the \code{line} you searched with, i.e. given to the function call}
#'
#'}
#' @examples
#' \donttest{
#' \dontrun{
#' # Interesting line that labels mushroom body neurons
#' ## But which ones?
#' line = "MB543B"
#'
#' # Let's find out
#' contents = neuronbridge_line_contents(line = line, threshold = 23000)
#' ## Note! Choosing the right threshold for the nomarlsied MIP-comparison score
#' ### Is critical. This may take some titring for lines you are really interested in
#' #### Though this value is a good first pass.
#'
#' # So this is what is likely in this line!
#' ## Note we have meta-data from neuprint! So we can see neuron types and names!!
#' # View(contents)
#'
#' # Now let's check what these neurons look like
#' if (!require("neuprintr")) remotes::install_github("natverse/neuprintr")
#'
#' # Load package
#' library(neuprintr)
#' ## Note you need to 'login' to neuPrint through R
#' ### Look at the package README and/or examine: ?neuprint_login
#'
#' # Let's see the EM neurons all together
#' neurons = neuprintr::neuprint_read_neurons(unique(contents$publishedName))
#' plot3d(neurons)
#'
#' # And compare with the LM data:
#' open3d()
#' mips = neuronbridge_mip(line)
#' scan_mip(mips,type="images", sleep = 5)
#'
#' }}
#' @seealso \code{\link{neuronbridge_info}},
#'   \code{\link{neuronbridge_avoid}},
#'   \code{\link{neuronbridge_search}}
#' @export
neuronbridge_line_contents <- function(line,
                                       neuprintr = TRUE,
                                       version = "v2_1_1",
                                       threshold = 23000){
  # Just one line!
  if(length(line)>1){
    stop("line must be of length 1")
  }

  # And one for the unfavoured ones
  line.search = neuronbridge_search(line, version = version, dataset = "by_line")

  # Find the
  line.search.filtered = stats::aggregate(list(normalizedScore=line.search$normalizedScore),
                                          list(publishedName=line.search$publishedName,
                                               libraryName=line.search$libraryName,
                                               searched=line.search$searched),
                                          function(x) max(x, na.rm = TRUE))
  line.search.filtered=line.search.filtered[order(line.search.filtered$normalizedScore,decreasing = TRUE),]

  # Just return IDs
  line.searched = subset(line.search.filtered, line.search.filtered$normalizedScore >= threshold)
  line.searched = line.searched[!is.na(line.searched$publishedName),]

  # Get metadata from neuprint
  if(neuprintr){
    if(!requireNamespace("neuprintr", quietly = TRUE)) {
      stop("Please install neuprintr using:\n", call. = FALSE,
           "remotes::install_github('natverse/neuprintr')")
    }
    hemi  = subset(line.searched, grepl("Hemibrain|hemibrain", line.searched$libraryName))
    hemi2 = neuprintr::neuprint_get_meta(hemi$publishedName)
    hemi2$publishedName = hemi2$bodyid
    line.searched = merge(line.searched,hemi2, all.x = TRUE, all.y = FALSE)
  }

  # Return
  line.searched
}


#' @rdname neuronbridge_line_contents
#' @export
neuronbridge_predict_split <-function(line1,
                                      line2,
                                      neuprintr = TRUE,
                                      version = "v2_1_1",
                                      threshold1 = 23000,
                                      threshold2 = 23000){
  # Just one line1!
  if(length(line1)>1){
    stop("line1 must be of length 1")
  }
  # Just one line2!
  if(length(line2)>1){
    stop("line2 must be of length 1")
  }

  # Contents of line 1
  l1 = neuronbridge_line_contents(line = line1, neuprintr = neuprintr, version = version, threshold = threshold1)

  # Contents of line 2
  l2 = neuronbridge_line_contents(line = line2, neuprintr = neuprintr, version = version, threshold = threshold2)

  # Intersect
  split = l1[l1$publishedName%in%l2$publishedName,]
  split$line1 = line1
  split$line2 = line2
  split$normalizedScore.line1 = split$normalizedScore
  split$normalizedScore.line2 = l2[match(split$publishedName,l2$publishedName),"normalizedScore"]
  split$searched = split$normalizedScore = NULL
  split

}


