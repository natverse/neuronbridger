#' @title Query the KDRC KGutProject database for line expression in the
#'   adult \emph{Drosophila} gut
#'
#' @description
#' The \href{https://flyinfo.kr/FlexForm_KGutProject_Load.do?tmpl=template/E}{KDRC
#' KGutProject} hosts a screen of 353 GAL4 / split-GAL4 driver lines for
#' expression across ten regions of the adult \emph{Drosophila} gut
#' (\code{Crop}, \code{PV}, \code{R1}--\code{R5}, \code{MHJ}, \code{Hindgut},
#' \code{Rectum}) from
#' \href{https://doi.org/10.1080/01677063.2020.1853722}{Lim et al. 2021}.
#' The site does not expose a machine-readable API; the search form is driven
#' by a proprietary \code{AUIquery.do} endpoint whose POST payload is signed
#' client-side.
#'
#' These helpers drive a headless Chrome session via \code{\link[chromote]{ChromoteSession}},
#' which reuses the site's own JavaScript to run the search and detail
#' look-ups. The result is the same as a user filling in the search box by
#' hand but scripted and batched.
#'
#' Call \code{kdrc_start_session()} once, then \code{kdrc_lookup_lines()}
#' with a vector of line identifiers (e.g. \code{"R20A02"}, \code{"VT012345"},
#' \code{"JK48870_1"}), and \code{kdrc_close_session()} when you are done.
#' The returned data frame has one row per (query_line, matching KDRC record)
#' pair, with the ten region Y/N flags from the KDRC detail page.
#'
#' @param session a KDRC chromote session, as returned by
#'   \code{kdrc_start_session()}. If \code{NULL}, a temporary session is
#'   created for the duration of the call.
#' @param line character; a single line identifier to search for. KDRC
#'   matches against \code{KGUTID}, \code{BDSCID}, \code{FlyLightID} and
#'   \code{AssociatedGene}.
#' @param lines character vector of line identifiers.
#' @param seq integer; the KDRC row identifier (column \code{seq} in a search
#'   result), used to fetch detail-page region flags.
#' @param chrome path to the Google Chrome / Chromium executable. If
#'   \code{NULL}, \code{chromote} uses its default detection.
#' @param timeout seconds to wait for each async query to populate a dataset
#'   before giving up on that line.
#' @param quiet logical; if \code{FALSE} (the default) print a one-line
#'   progress update per 25 lines.
#'
#' @return \code{kdrc_search_line()} returns a \code{data.frame} of zero or
#'   more KDRC hit records (columns \code{seq}, \code{KGUTID}, \code{kGutType},
#'   \code{BDSCID}, \code{FlyLightID}, \code{AssociatedGene}).
#'
#'   \code{kdrc_line_regions()} returns a named character vector of the ten
#'   region flags for a given \code{seq} (\code{"Y"}, \code{"N"} or
#'   \code{""} where not scored).
#'
#'   \code{kdrc_lookup_lines()} returns a \code{data.frame} with one row per
#'   (\code{query_line}, KDRC record) pair, columns \code{query_line},
#'   \code{kgut_hit}, \code{kGutType}, \code{KGUTID}, \code{BDSCID},
#'   \code{FlyLightID}, \code{AssociatedGene}, \code{p_seq}, and the ten
#'   region flags \code{Crop}, \code{PV}, \code{R1}..\code{R5}, \code{MHJ},
#'   \code{Hindgut}, \code{Rectum}. Lines with no hit appear as a single row
#'   with \code{kgut_hit = "N"} and everything else \code{NA}.
#'
#'   \code{kdrc_start_session()} returns an environment wrapping the
#'   chromote session; pass it back to the other functions via the
#'   \code{session} argument.
#'
#' @section Dependencies:
#' Requires the \code{chromote} R package and a local Chrome / Chromium
#' install. Install with \code{install.packages("chromote")}.
#'
#' @section Citation:
#' Lim, S.Y., et al. (2021).
#' \emph{Identification and characterization of GAL4 drivers that mark
#' distinct cell types and regions in the Drosophila adult gut.}
#' \href{https://doi.org/10.1080/01677063.2020.1853722}{J. Neurogenet. 35(1):33--44}.
#'
#' @examples
#' \donttest{
#' \dontrun{
#' # One-off lookup
#' kdrc_lookup_lines("R20A02")
#'
#' # Batch lookup with an explicit session
#' s <- kdrc_start_session()
#' res <- kdrc_lookup_lines(c("R20A02", "VT012345", "SS36097"), session = s)
#' kdrc_close_session(s)
#' }}
#' @name kdrc
NULL

KDRC_REGIONS <- c("Crop", "PV", "R1", "R2", "R3", "R4", "R5",
                  "MHJ", "Hindgut", "Rectum")

KDRC_SEARCH_URL <- "https://flyinfo.kr/FlexForm_KGutProject_Load.do?tmpl=template/E"
KDRC_DETAIL_URL <- "https://flyinfo.kr/FlexForm_KGut_Detail_Load.do?tmpl=template/E&p_seq=%s"

# hidden
# Navigate back to the search page if the session has drifted (e.g. after
# a detail-page fetch). Checks window.location rather than reloading
# unconditionally, which would slow batch lookups ~2x.
kdrc_ensure_search_page <- function(session) {
  loc <- session$chromote$Runtime$evaluate(
    "location.pathname + '|' + (typeof ds_Kgut) + '|' + (typeof fillJob_KGUT\\uc870\\ud68c)",
    returnByValue = TRUE)$result$value
  on_search <- is.character(loc) &&
    grepl("FlexForm_KGutProject_Load", loc, fixed = TRUE) &&
    grepl("function$", loc)
  if (on_search) return(invisible(NULL))
  session$chromote$Page$navigate(KDRC_SEARCH_URL)
  Sys.sleep(3)
  invisible(NULL)
}

# hidden
kdrc_require_chromote <- function() {
  if (!requireNamespace("chromote", quietly = TRUE)) {
    stop("This function requires the 'chromote' package. Install with:\n",
         "  install.packages(\"chromote\")", call. = FALSE)
  }
}

# hidden
# Poll a JS expression until it returns a non-empty JSON array, or timeout.
kdrc_wait_for <- function(session, expr, timeout = 8, interval = 0.2) {
  deadline <- Sys.time() + timeout
  repeat {
    out <- tryCatch(
      session$chromote$Runtime$evaluate(expr, returnByValue = TRUE),
      error = function(e) NULL)
    val <- if (!is.null(out)) out$result$value else NULL
    if (is.character(val) && nzchar(val)) {
      parsed <- tryCatch(jsonlite::fromJSON(val, simplifyVector = FALSE),
                         error = function(e) NULL)
      if (is.list(parsed) && length(parsed)) return(parsed)
    }
    if (Sys.time() > deadline) return(list())
    Sys.sleep(interval)
  }
}

#' @rdname kdrc
#' @export
kdrc_start_session <- function(chrome = NULL) {
  kdrc_require_chromote()
  if (!is.null(chrome)) Sys.setenv(CHROMOTE_CHROME = chrome)
  prev_args <- getOption("chromote.launch_args")
  options(chromote.launch_args = c(
    "--no-sandbox", "--disable-gpu",
    sprintf("--user-data-dir=%s", tempfile("chromote_kdrc_"))
  ))
  cs <- chromote::ChromoteSession$new()
  cs$Page$navigate(KDRC_SEARCH_URL)
  Sys.sleep(4)   # wait for DataSets/controls to register
  session <- new.env(parent = emptyenv())
  session$chromote <- cs
  session$prev_launch_args <- prev_args
  class(session) <- "kdrc_session"

  # The very first search after the page loads is occasionally flaky: the
  # DataSet initializer races with `fillJob_KGUT<U+C870><U+D68C>`. Fire a throwaway
  # query so the next real call is clean.
  try(kdrc_search_line("__warmup__", session = session, timeout = 4),
      silent = TRUE)
  session
}

#' @rdname kdrc
#' @export
kdrc_close_session <- function(session) {
  if (!inherits(session, "kdrc_session")) return(invisible(NULL))
  try(session$chromote$close(), silent = TRUE)
  options(chromote.launch_args = session$prev_launch_args)
  invisible(NULL)
}

# hidden
kdrc_ensure_session <- function(session) {
  if (is.null(session)) {
    s <- kdrc_start_session()
    attr(s, "owned") <- TRUE
    s
  } else {
    if (!inherits(session, "kdrc_session"))
      stop("`session` must be a kdrc_session created by kdrc_start_session().",
           call. = FALSE)
    session
  }
}

#' @rdname kdrc
#' @export
kdrc_search_line <- function(line, session = NULL, timeout = 8) {
  kdrc_require_chromote()
  s <- kdrc_ensure_session(session)
  on.exit(if (isTRUE(attr(s, "owned"))) kdrc_close_session(s), add = TRUE)
  line <- as.character(line)
  kdrc_ensure_search_page(s)

  # Reset and fire the site's own search function; site uses a synchronous
  # AJAX call so records appear very quickly.
  js <- sprintf('
    try { ds_Kgut.clearData && ds_Kgut.clearData(); ds_Kgut.records=[]; } catch(e){}
    try { controls["SearchText"].value=%s; } catch(e){}
    try { fillJob_KGUT\uc870\ud68c && fillJob_KGUT\uc870\ud68c(1); } catch(e){}
    "ok"', jsonlite::toJSON(line, auto_unbox = TRUE))
  s$chromote$Runtime$evaluate(js, returnByValue = TRUE)

  recs <- kdrc_wait_for(
    s,
    expr = "JSON.stringify((ds_Kgut.records||[]).map(function(r){return r.values||{};}))",
    timeout = timeout)
  if (!length(recs)) {
    return(data.frame(
      seq = integer(0), KGUTID = character(0), kGutType = character(0),
      BDSCID = character(0), FlyLightID = character(0),
      AssociatedGene = character(0),
      stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(recs, function(v) data.frame(
    seq            = as.integer(v$seq            %||% NA),
    KGUTID         = v$KGUTID         %||% NA_character_,
    kGutType       = v$kGutType       %||% NA_character_,
    BDSCID         = as.character(v$BDSCID        %||% NA),
    FlyLightID     = v$FlyLightID     %||% NA_character_,
    AssociatedGene = v$AssociatedGene %||% NA_character_,
    stringsAsFactors = FALSE)))
}

#' @rdname kdrc
#' @export
kdrc_line_regions <- function(seq, session = NULL, timeout = 8) {
  kdrc_require_chromote()
  s <- kdrc_ensure_session(session)
  on.exit(if (isTRUE(attr(s, "owned"))) kdrc_close_session(s), add = TRUE)
  stopifnot(length(seq) == 1, !is.na(seq))
  s$chromote$Page$navigate(sprintf(KDRC_DETAIL_URL, seq))
  Sys.sleep(3)

  # ds_KGut on the detail page carries the line's actual Y/N region flags.
  recs <- kdrc_wait_for(
    s,
    expr = "JSON.stringify((ds_KGut.records||[]).map(function(r){return r.values||{};}))",
    timeout = timeout)
  flags <- stats::setNames(rep(NA_character_, length(KDRC_REGIONS)),
                           KDRC_REGIONS)
  if (length(recs)) {
    v <- recs[[1]]
    for (r in KDRC_REGIONS) flags[r] <- v[[r]] %||% NA_character_
  }
  # KDRC encodes expression with "O" (observed) and "X" (not observed).
  # Normalise to "Y"/"N" so the output plays nicely with the rest of the
  # package (and with pre-existing manual lookup templates).
  kdrc_normalise_flag(flags)
}

# hidden
kdrc_normalise_flag <- function(x) {
  ifelse(is.na(x) | x == "", NA_character_,
         ifelse(x %in% c("O", "Y", "o", "y"), "Y",
                ifelse(x %in% c("X", "N", "x", "n"), "N", x)))
}

#' @rdname kdrc
#' @export
kdrc_lookup_lines <- function(lines, session = NULL, timeout = 8, quiet = FALSE) {
  kdrc_require_chromote()
  lines <- as.character(lines)
  s <- kdrc_ensure_session(session)
  on.exit(if (isTRUE(attr(s, "owned"))) kdrc_close_session(s), add = TRUE)

  # First: do every search on the search page (fast, no navigation).
  hits <- vector("list", length(lines))
  t0 <- Sys.time()
  for (i in seq_along(lines)) {
    hits[[i]] <- kdrc_search_line(lines[i], session = s, timeout = timeout)
    if (!quiet && (i %% 25 == 0 || i == length(lines))) {
      dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      message(sprintf("  [kdrc_lookup_lines] %d/%d searched  (%.0fs, %.2fs/line)",
                      i, length(lines), dt, dt / i))
    }
  }

  rows <- list()
  for (i in seq_along(lines)) {
    h <- hits[[i]]
    if (!nrow(h)) {
      rows[[length(rows) + 1L]] <- data.frame(
        query_line = lines[i], kgut_hit = "N",
        kGutType = NA_character_, KGUTID = NA_character_,
        BDSCID = NA_character_, FlyLightID = NA_character_,
        AssociatedGene = NA_character_, p_seq = NA_integer_,
        matrix(NA_character_, nrow = 1, ncol = length(KDRC_REGIONS),
               dimnames = list(NULL, KDRC_REGIONS)),
        stringsAsFactors = FALSE, check.names = FALSE)
      next
    }
    for (k in seq_len(nrow(h))) {
      # Load detail page per record for Y/N region flags.
      regions <- kdrc_line_regions(h$seq[k], session = s, timeout = timeout)
      rows[[length(rows) + 1L]] <- data.frame(
        query_line = lines[i], kgut_hit = "Y",
        kGutType = h$kGutType[k], KGUTID = h$KGUTID[k],
        BDSCID = h$BDSCID[k], FlyLightID = h$FlyLightID[k],
        AssociatedGene = h$AssociatedGene[k], p_seq = h$seq[k],
        as.data.frame(as.list(regions), stringsAsFactors = FALSE),
        stringsAsFactors = FALSE, check.names = FALSE)
    }
    # Navigating to detail pages for this line moves the tab away from the
    # search page; reload it so the next line's search works.
    s$chromote$Page$navigate(KDRC_SEARCH_URL)
    Sys.sleep(2.5)
  }
  do.call(rbind, rows)
}

# hidden
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
