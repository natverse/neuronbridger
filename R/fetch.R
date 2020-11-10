# hidden
neuronbridge_fetch <- function(path,
                               body = NULL,
                               neuronbridge_url = "https://neuronbridge.janelia.org/",
                               parse = c("json","xml","none"),
                               simplifyVector = FALSE,
                               include_headers = FALSE,
                               ...){
  parse = match.arg(parse)
  path = gsub("\\/$|^\\/","",path)
  if(grepl("^https://",path)){
    final.path = path
  }else{
    final.path = file.path(neuronbridge_url, path, fsep = "/")
  }
  req <-
    if (is.null(body)) {
      httr::GET(url = final.path,
                ...)
    }else {
      httr::POST(url = final.path,
                 body = body,
                 ...)
    }
  neuronbridge_error_check(req)
  if (parse!="none") {
    parsed = neuronbridge_parse(req, parse = parse, simplifyVector = simplifyVector, raw = FALSE, ...)
    neuronbridge_error_check(parsed)
    if (include_headers) {
      fields_to_include = c("url", "headers")
      attributes(parsed) = c(attributes(parsed), req[fields_to_include])
    }
    parsed
  }
  else req
}

# hidden
neuronbridge_parse <- function (req,
                                parse = c("json","xml","none"),
                                simplifyVector = FALSE,
                                raw = TRUE,
                                ...) {
  parse = match.arg(parse)
  if(raw){
    text <- rawToChar(req$content)
  }else{
    text <- httr::content(req, as = "text", encoding = "UTF-8")
  }
  if (identical(text, "")){
    warning("No output to parse", call. = FALSE)
    return(NULL)
  }
  if(parse=="json"){
    p = tryCatch(jsonlite::fromJSON(text, simplifyVector = simplifyVector, ...), error = function(e) NULL)
  }else{
    xmlt = XML::xmlParse(text)
    p = XML::xmlToList(xmlt)
  }
  if(is.null(p)){
    warning("error parsing ", parse)
  }
  tryCatch(nullToNA(p), error = function(e) NA )
}

# hidden
neuronbridge_error_check <- function(x){
  err_fields = c("error", "message")
  if (sum(names(x) %in% err_fields)>1) {
    stop(x$error, ": ", x$message)
  }
  NULL
}
