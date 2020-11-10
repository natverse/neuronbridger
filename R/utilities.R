# hidden
unnest_df <- function(df){
  data = as.data.frame(df, stringsAsFactors = FALSE)
  if(nrow(data)&ncol(data)){
    data = as.data.frame(df, stringsAsFactors = FALSE)
    data = apply(data,2,function(c) unlist(nullToNA(c)))
    data = as.data.frame(unlist(data), stringsAsFactors = FALSE)
    dimnames(data) = dimnames(data)
    data
  }
  data[] <- lapply(data, as.character)
  data
}

# hidden
nullToNA <- function(x) {
  x[sapply(x, is.null)] <- NA
  return(x)
}

# hidden
rep_name <- function(x, sep = "#"){
  paste0(x,sep, stats::ave(x,x,FUN= seq.int))
}
