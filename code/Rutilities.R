# Function linkMetaboliteIdentifiers
#TODO: if necessary, allow timeout option. Currently, we are probably working
# under global options timeout, getOption('timeout'), but not sure how the
# error will look like in the case of timeout (i.e., if it will be properly
# caught by try-error)

#TODO: perhaps useful to add a "bidirectional option": for each query, find each
# target, and for each target, find each of type query. This could also be done
# in subsequent runs, but having it in one function could be handy
# Do different chebi's ever return multiple KEGG? Probably not.
# But the other way around this can happen.

linkMetaboliteIdentifiers <- function(query, 
                                      target,
                                      reader = c("read.csv","fread","readLines"),
                                      database = c("pubchem_compound",
                                                   "pubchem_substance"),
                                      bidirectional = c(FALSE, TRUE),
                                      verbose = c(TRUE, FALSE)){
  
  res <- unlist(lapply(query, function(q){
    
    if(database == "pubchem_compound"){
      url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                    q,
                    "/synonyms/TXT")
    } else if(database == "pubchem_substance"){
      url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/",
                    q,
                    "/synonyms/TXT")
    }
    
    if(reader == "read.csv"){
      data <- tryCatch({
        suppressWarnings(read.csv(url,sep="\t",header = FALSE)$V1)
      },error = function(e) e)
    } else if(reader == "fread"){
      data <- tryCatch({
        suppressWarnings(fread(url,sep="\t",header = FALSE, showProgress = FALSE)$V1)
      },error = function(e) e)
    } else if(reader == "readLines"){
      data <- tryCatch({
        suppressWarnings(readLines(url))
      },error = function(e) e)
    }
    
    if (inherits(data, "error")) {
      return("Unknown query")
    }
    
    if(sum(grepl(target,data)) > 0){ # see if faster option exists
      return(paste0(data[grep(target,data)],collapse = "/")) # works also if multiple chebi's
    } else{
      return(paste0("No ", target))
    }
    
  }))
  
  if(verbose){
    warning("Number of unknown queries: ", sum(res == "Unknown query"), "\n")
    warning("Number of queries without ", target,  ": ", sum(res == paste0("No ", target)), "\n")
  }
  
  error <- rep(NA,length=length(res))
  error[res == "Unknown query"] <- "Unknown query"
  error[res == paste0("No ", target)] <- paste0("No ", target)
  res[which(res %in% c("Unknown query",
                       paste0("No ", target)))] <- NA
  
  return(data.frame(query = query,
                    target = res,
                    error = error))
}


consolidateIDs <- function(data,
                           queryID,
                           originalID,
                           link,
                           strategy = c("union","intersection","keep_original","keep_new"),
                           keepOriginal = c(TRUE, FALSE)){
  
  stopifnot("queryID not in colnames(data)" = queryID %in% colnames(data))
  #TODO make less cryptic
  stopifnot("Mismatch between link and data" = all(!is.na(match(link$query, data[[queryID]]))))
  
  newLink <- link$target[match(data[[queryID]],link$query)]
  
  consolidated <- unlist(sapply(seq_along(data[[originalID]]), 
                                function(x){
                                  if(is.na(data[[originalID]][x])){
                                    return(newLink[x])
                                  } else if(is.na(newLink[x])){
                                    return(data[[originalID]][x])
                                  } else{
                                    
                                    original_x <- unlist(str_split(data[[originalID]][x], pattern = "/"))
                                    new_x <- unlist(str_split(newLink[x], pattern = "/"))
                                    
                                    if(strategy == "union"){
                                      combined <- union(original_x,new_x)
                                    } else if(strategy == "intersect") {
                                      combined <- intersect(original_x,new_x)
                                    } else if(strategy == "original_only"){
                                      combined <- original_x
                                    } else if(strategy == "new_only"){
                                      combined <- new_x
                                    }
                                    
                                    return(paste0(combined,collapse = "/"))
                                  }
                                }))
  
  if(keepOriginal){
    data[[paste0(originalID,"_new")]] <- newLink
    data[[paste0(originalID,"_consolidated")]] <- consolidated
    colnames(data)[which(colnames(data) == originalID)] <- paste0(originalID,"_original")
    return(data)
  } else {
    data[[originalID]] <- consolidated
  }
}

read_excel_allsheets <- function(filename) {
  sheets <- openxlsx::getSheetNames(filename)
  x <- lapply(sheets, function(X) openxlsx::read.xlsx(filename, sheet = X, sep.names = " "))
  names(x) <- sheets
  x
}
