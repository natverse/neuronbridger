#' @title Resolve split-GAL4 lines to their AD and DBD hemidriver halves
#'
#' @description
#' Split-GAL4 driver lines are built by combining two hemidriver transgenes:
#' an activation-domain (AD) half and a DNA-binding-domain (DBD) half. Each
#' half is itself an enhancer-driven transgene, named for the enhancer
#' fragment (a GMR \code{R*} or VT \code{VT*} identifier). Knowing the two
#' halves matters for cross-referencing: a split can \dQuote{share} one half
#' with a published Gen1 line that may itself have independent expression
#' data (e.g. in the KDRC KGutProject).
#'
#' \code{split_halves()} looks up each \code{line} in the FlyLight Split-GAL4
#' image-release metadata on S3
#' (\url{https://s3.amazonaws.com/janelia-flylight-imagery/}). For each
#' \code{SS*}, \code{IS*}, \code{MB*}, \code{OL*} or \code{SL*} line that is
#' present in any of the published split-GAL4 collections, the FlyLight
#' metadata JSON carries \code{ad} and \code{dbd} string fields (e.g.
#' \code{"72F01-p65ADZp in attP40"}) which we parse into canonical enhancer
#' codes (\code{R72F01}, \code{VT048352}, etc.).
#'
#' Gen1 lines (plain \code{R*} or \code{VT*}) are returned as-is with
#' \code{ad} and \code{dbd} both \code{NA}.
#'
#' @param lines character vector of line identifiers.
#' @param timeout seconds for each individual S3 request.
#' @param quiet logical; if \code{FALSE} (the default) emit a progress line
#'   every 25 inputs.
#'
#' @return a \code{data.frame} with one row per input line and columns
#'   \itemize{
#'     \item{\code{line}}{ - the query string.}
#'     \item{\code{is_split}}{ - logical; \code{TRUE} if the line looks like
#'       a split-GAL4 (prefix SS/IS/MB/OL/SL).}
#'     \item{\code{ad}}{ - resolved AD-half enhancer code (\code{R*} or
#'       \code{VT*}), or \code{NA}.}
#'     \item{\code{ad_full}}{ - the raw \code{ad:} field from the S3 metadata
#'       (includes landing site and balancer info), or \code{NA}.}
#'     \item{\code{dbd}}{ - resolved DBD-half enhancer code, or \code{NA}.}
#'     \item{\code{dbd_full}}{ - the raw \code{dbd:} field, or \code{NA}.}
#'     \item{\code{source_release}}{ - the FlyLight image release in which
#'       metadata was found (e.g. \dQuote{Split-GAL4 Omnibus Broad}), or
#'       \code{NA}.}
#'   }
#'
#' @details
#' The S3 release folders checked are (in order, first hit wins):
#' \itemize{
#'   \item \code{"Split-GAL4 Omnibus Broad"} (Meissner et al. 2024 Initial
#'         Split \code{IS*} lines and the broad-screen \code{SS*} lines)
#'   \item \code{"Split-GAL4 Omnibus Rescreen"} (2024 rescreen)
#'   \item \code{"Split-GAL4"} (legacy Aso 2014 and similar \code{SS*}/
#'         \code{MB*} releases)
#'   \item \code{"Split GAL4"} (alternate folder name used by some releases)
#' }
#'
#' If none of these contain a metadata JSON for a queried line, \code{ad}
#' and \code{dbd} are returned as \code{NA}. As a fallback for unlisted
#' lines, users can paste the line name into
#' \url{https://splitgal4.janelia.org/cgi-bin/view_splitgal4_imagery.cgi} or
#' the Bloomington Drosophila Stock Center stock page to find the two
#' components by hand.
#'
#' @examples
#' \donttest{
#' \dontrun{
#' split_halves(c("SS36097", "IS34128", "MB543B", "R20A02"))
#' }}
#' @seealso \code{\link{kdrc_lookup_lines}}
#' @export
split_halves <- function(lines, timeout = 10, quiet = FALSE) {
  lines <- as.character(lines)
  out <- lapply(seq_along(lines), function(i) {
    row <- split_resolve_one(lines[i], timeout = timeout)
    if (!quiet && (i %% 25 == 0 || i == length(lines)))
      message(sprintf("  [split_halves] %d/%d", i, length(lines)))
    row
  })
  do.call(rbind, out)
}

# hidden
# FlyLight split-GAL4 release folders on s3://janelia-flylight-imagery/,
# ordered by roughly how many lines each hosts. `split_resolve_one()` picks
# the first release that has a metadata JSON for the query line.
SPLIT_RELEASES <- c(
  "Split-GAL4 Omnibus Broad",
  "Split-GAL4 Omnibus Rescreen",
  "MB Paper 2014",
  "Aso 2021",
  "Descending Neurons 2018",
  "Descending Neurons 2025",
  "Dorsal VNC 2023",
  "Cheong, Eichler, Stuerner 2023",
  "Cheong, Boone, Bennett 2023",
  "Ascending Neurons 2023",
  "Nern et al 2024",
  "Wolff 2018",
  "Wolff et al 2024",
  "Lateral Horn 2019",
  "SEZ 2021",
  "oviDN 2020",
  "oviDN 2023",
  "Lillvis 2022",
  "Lillvis 2023",
  "Schretter et al 2024",
  "Schretter et al. 2020",
  "Hampel 2015",
  "LC Paper",
  "Aso&Rubin 2016"
)

# hidden
split_is_split_like <- function(line) {
  grepl("^(SS|IS|MB|OL|SL)[0-9]", line)
}

# hidden
split_resolve_one <- function(line, timeout = 10) {
  template <- function(ad = NA_character_, ad_full = NA_character_,
                       dbd = NA_character_, dbd_full = NA_character_,
                       release = NA_character_)
    data.frame(line = line, is_split = split_is_split_like(line),
               ad = ad, ad_full = ad_full,
               dbd = dbd, dbd_full = dbd_full,
               source_release = release, stringsAsFactors = FALSE)

  if (!split_is_split_like(line)) return(template())

  for (release in SPLIT_RELEASES) {
    key <- split_find_metadata_key(release, line, timeout = timeout)
    if (is.null(key)) next
    meta <- split_fetch_metadata(key, timeout = timeout)
    if (is.null(meta)) next
    ad_raw  <- meta$ad  %||% NA_character_
    dbd_raw <- meta$dbd %||% NA_character_
    return(template(
      ad      = split_parse_enhancer(ad_raw),
      ad_full = ad_raw,
      dbd     = split_parse_enhancer(dbd_raw),
      dbd_full = dbd_raw,
      release = release))
  }
  template()
}

# hidden
split_find_metadata_key <- function(release, line, timeout = 10,
                                    max_keys = 200) {
  url <- sprintf(
    "https://s3.amazonaws.com/janelia-flylight-imagery/?list-type=2&prefix=%s/%s/&max-keys=%d",
    utils::URLencode(release, reserved = TRUE), line, max_keys)
  txt <- split_fetch_text(url, timeout = timeout)
  if (is.null(txt)) return(NULL)
  keys <- regmatches(txt,
                     gregexpr("<Key>[^<]+metadata\\.json</Key>", txt))[[1]]
  if (!length(keys)) return(NULL)
  sub("^<Key>(.*)</Key>$", "\\1", keys[1])
}

# hidden
split_fetch_text <- function(url, timeout = 10) {
  tryCatch({
    resp <- httr::GET(url, httr::timeout(timeout))
    if (httr::status_code(resp) != 200) return(NULL)
    httr::content(resp, as = "text", encoding = "UTF-8")
  }, error = function(e) NULL)
}

# hidden
split_fetch_metadata <- function(key, timeout = 10) {
  url <- paste0("https://s3.amazonaws.com/janelia-flylight-imagery/",
                utils::URLencode(key))
  txt <- split_fetch_text(url, timeout = timeout)
  if (is.null(txt)) return(NULL)
  tryCatch(jsonlite::fromJSON(txt, simplifyVector = TRUE),
           error = function(e) NULL)
}

# hidden
# Parse a FlyLight AD/DBD raw string into a canonical enhancer / driver code.
# Handles the common forms we see in the S3 metadata:
#   "72F01-p65ADZp in attP40; MKRS/TM6B" -> "R72F01"
#   "VT048352-ZpGDBD in attP2"           -> "VT048352"
#   "26D04-p65ADZp in VK00027"           -> "R26D04"
#   "Ddc-p65ADZp in VK00027"             -> "Ddc"   (gene-specific enhancer)
#   "HAV5; Sco/CyO; 13E04-ZpGdbd in attP2" -> "R13E04"  (balancer prefix)
#   "publishing_name"                    -> NA  (FlyLight placeholder
#                                                 when components were not
#                                                 populated for this sample)
#   NA / ""                              -> NA
split_parse_enhancer <- function(raw) {
  if (is.null(raw) || is.na(raw) || !nzchar(raw)) return(NA_character_)
  if (identical(raw, "publishing_name")) return(NA_character_)
  # Enhancer codes first: GMR or VT numeric, optionally preceded by balancer
  # cruft separated by ';'.
  m <- regmatches(raw, regexec("(VT[0-9]+|[0-9]+[A-Z][0-9]+)-(p65AD|ZpG)",
                               raw))[[1]]
  if (length(m) >= 2) {
    code <- m[2]
    return(if (startsWith(code, "VT")) code else paste0("R", code))
  }
  # Gene-name splits (e.g. "Ddc-p65ADZp ..."). Return the bare gene token.
  m <- regmatches(raw, regexec("^([A-Za-z][A-Za-z0-9_-]*)-(p65AD|ZpG)",
                               raw))[[1]]
  if (length(m) >= 2) return(m[2])
  NA_character_
}
