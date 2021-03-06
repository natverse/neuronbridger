% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mip.R
\name{neuronbridge_mip}
\alias{neuronbridge_mip}
\alias{plot_mip}
\alias{scan_mip}
\alias{neuronbridge_remove_mips}
\title{Get and plot colour MIPs for neurons from NeuronBridge}
\usage{
neuronbridge_mip(
  id,
  dataset = c("detect", "by_line", "by_body"),
  images.url = "https://s3.amazonaws.com/janelia-flylight-color-depth",
  db = options("neuronbridger"),
  version = "v2_1_1"
)

plot_mip(mip, db = options("neuronbridger"), file.delete = FALSE)

scan_mip(
  mips,
  no.hits = 10,
  type = c("hits", "files", "ids", "images"),
  sleep = NULL,
  images.url = "https://s3.amazonaws.com/janelia-flylight-color-depth",
  db = options("neuronbridger")
)

neuronbridge_remove_mips(db = options("neuronbridger"), interactive = TRUE)
}
\arguments{
\item{id}{character vector. An identifier/identifiers for the neuron(s)/line(s) you wish to search. This can either be a line code
for a GAL4 line (e.g. \href{https://v2.virtualflybrain.org/org.geppetto.frontend/geppetto?id=VFBexp_FBtp0057908&i=VFB_00017894}{R16F12})
or a split GAL4 line (e.g. \href{https://splitgal4.janelia.org/cgi-bin/view_splitgal4_imagery.cgi?line=LH173}{LH173})
or a 'body ID' for a neuron from the \href{https://neuprint.janelia.org/help/videos?dataset=hemibrain}{hemibrain connectome},
(e.g. \href{https://neuprint.janelia.org/view?bodyid=1537331894}{1537331894}).}

\item{dataset}{whether the ID given is a body ID (\code{"by_body"}) from the hemibrain connectome
or a genetic driver line (\code{"by_line"}). If left at \code{"detect"} then \code{neuronbridger} tries to guess to which \code{id} belongs by using
\code{neuronbridger:::guess_dataset}.}

\item{images.url}{the online database from which NeuronBridge gets its MIP files.
Higher res .png files come from https://s3.amazonaws.com/janelia-flylight-color-depth.
Lower res .jpg thumbnails come from https://s3.amazonaws.com/janelia-flylight-color-depth-thumbnails.}

\item{db}{the directory in which image files are saved (\code{neuronbridge_mip}), and from which they are read (\code{plot_mip}). This defaults to a temporary
directory: \code{file.path(tempdir(),"neuronbridger")}.}

\item{version}{the precomputed scores to search. For example, \code{"v2_1_1"} refers to this \href{https://neuronbridge.janelia.org/releasenotes/DATA}{release}.}

\item{mip}{either a full path to \code{.png} file or an array describing a image file, e.g. as comes from \code{png::readPNG}.}

\item{file.delete}{logical, whether or not to delete the saved \code{.png} file after plotting.}

\item{mips}{either a data frame of 'hits' (e.g. as produced by \code{\link{neuronbridge_search}}), a vector of file paths to saved MIP .png files
or ids for data items (i.e. from \code{\link{neuronbridge_ids}}).}

\item{no.hits}{the number of hits to visualise if \code{type="hits"}.}

\item{type}{the type of argument given as \code{mips}.}

\item{sleep}{defaults to \code{NULL}. User needs to press a key to move to the next MIP. If a numeric value, then MIPs progress automatically and
a pause of \code{sleep} seconds occurs before the next MIP is shown.}

\item{interactive}{logical. When using \code{neuronbridge_remove_mips}, whether or not to be accosted with an 'are you sure'?.}
}
\value{
a named list of \code{.png} (\code{neuronbridge_mip}) or visualisation of a \code{.png} MIP in an \code{rgl} display (\code{plot_mip}).
}
\description{
These functions allow one to get (\code{neuronbridge_mip}) and plot (\code{plot_mip}) colour MIPs (maximum projection images).
MIPs are not 'stacks' of images, but a single image plane where colour encodes depth. The bluer the hue, the more anterior (front of the brain) the signal.
The redder the hue, the more posterior (back of the brain). There may be differences/difficulties plotting images between different operating systems.
Only extensively tested on MacOS Catalina.
}
\details{
By default MIPs are saved to your R session's temporary directory. However, you can specify an alternative directory (argument: \code{db}).
MIPs are saved as \code{.png} files, where the file names are given as \code{id_number.png}.

FlyLight file names contain metadata as follows:
\code{'[Publishing Name]-[Slide Code]-[Driver]-[Gender]-[Objective]-[Area/Tile]-[Alignment Space]-CDM_[Channel].png'}

Find information on these meta data entries by examining \code{?neuronbridge_hits}.
}
\examples{
\donttest{
\dontrun{
# Let us now see the related MIP file
mip = neuronbridge_mip("542634818")
## This gets every MIP file associated with id
## In this case, just one

# Plot the MIP image in an rgl viewer
plot_mip(mip)
}}
}
\seealso{
\code{\link{neuronbridge_ids}},
  \code{\link{neuronbridge_info}},
  \code{\link{neuronbridge_search}}
}
