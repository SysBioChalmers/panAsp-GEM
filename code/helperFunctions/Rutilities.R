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
                                      database = c("pubchem_compound_name",
                                                   "pubchem_substance_name",
                                                   "pubchem_compound_cid"),
                                      bidirectional = c(FALSE, TRUE),
                                      verbose = c(TRUE, FALSE)){
  
  res <- unlist(lapply(query, function(q){
    
    if(database == "pubchem_compound_name"){
      url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                    q,
                    "/synonyms/TXT")
    } else if(database == "pubchem_substance_name"){
      url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/",
                    q,
                    "/synonyms/TXT")
    } else if(database == "pubchem_compound_cid"){
      url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/",
                    q,
                    "/cids/TXT")
    }
    
    url <- gsub(" ", "%20", url) # always replace spaces with %20 in URL
    
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

# Function to update column names
update_colnames <- function(x) {
  x <- sub("^[A-H][0-9]{2}\\.\\.", "", x)   # remove well prefix
  x <- gsub("\\.$", "", x)                   # remove trailing dot
  x <- gsub("\\.", " ", x)                   # replace all dots with spaces
  x <- gsub("\\s+", " ", x)                  # collapse double spaces
  x <- trimws(x)
  
  # Single-letter stereo prefixes
  prefixes <- c("L", "D", "R", "S", "a", "b", "m", "p", "o", "n", "N", "g", "d", "e")
  pattern <- paste0("\\b(", paste(prefixes, collapse="|"), ") ")
  x <- gsub(pattern, "\\1-", x)
  
  # Number prefixes
  x <- gsub("(\\d) (\\d) ", "\\1,\\2-", x)
  x <- gsub("(\\d) ([A-Z])", "\\1-\\2", x)
  
  # g Lactone -> -g-Lactone
  x <- gsub(" g Lactone", "-g-Lactone", x)
  x <- gsub(" g-Lactone", "-g-Lactone", x)
  
  # Dipeptides: two three-letter capitalised words e.g. "Gly Asp"
  x <- gsub("\\b([A-Z][a-z]{2}) ([A-Z][a-z]{2})\\b", "\\1-\\2", x)
  
  # Tween + number
  x <- gsub("(Tween) (\\d+)", "\\1-\\2", x)
  
  # myo-Inositol
  x <- gsub("myo ", "myo-", x)
  
  # Mono-, Deoxy- prefixes
  x <- gsub("\\b(Mono|Deoxy) ", "\\1-", x)
  
  # (Di)hydroxy, Keto, Bromo fuse with next word (lowercase the following word)
  x <- gsub("(Hydroxy|Dihydroxy|Keto|Bromo|Acetyl) ([A-Z])", "\\1\\L\\2", x, perl = TRUE)
  
  # But N-Acetyl- followed by stereo prefix should still hyphenate not fuse
  # so restore: N-Acetyl-[stereo] pattern
  x <- gsub("N-Acetyl([a-z])", "N-Acetyl-\\U\\1", x, perl = TRUE)  # undo fusion after N-
  
  # Acetyl, Methyl and Phtaloyl followed by stereo prefix need hyphen
  x <- gsub("(Acetyl|Methyl|Phtaloyl) ([A-Z]-)", "\\1-\\2", x)
  x <- gsub("(Acetyl|Methyl|Phtaloyl) ([a-z]-)", "\\1-\\2", x)
  
  # Methyl followed by plain word
  x <- gsub("(Methyl) ([A-Z][a-z])", "\\1-\\2", x)
  
  # Phosphate: number or word before it needs hyphen
  x <- gsub(" (\\d)-Phosphate", "-\\1-Phosphate", x)
  x <- gsub("(\\d) (Phosphate)", "\\1-\\2", x)
  x <- gsub(" Phosphate", "-Phosphate", x)
  
  # Bromo/Mono compound trailing word
  x <- gsub("(Bromo-[A-Za-z]+) ([A-Z][a-z])", "\\1-\\2", x)
  
  # Amino hyphenates to next word
  x <- gsub("(Amino) ([A-Z])", "\\1-\\2", x)
  
  # Methyl followed by plain word but NOT Ester
  x <- gsub("(Methyl) ([A-Z][a-z])(?!.*Ester)", "\\1-\\2", x, perl = TRUE)
  # Simpler: only fuse Methyl-Ester as two separate words
  x <- gsub("Methyl-Ester", "Methyl Ester", x)
  
  # Number followed by word needs hyphen (e.g. "Ribono 1,4-Lactone" -> "Ribono-1,4-Lactone")
  x <- gsub("([a-z]) (\\d)", "\\1-\\2", x)
  
  # Specific fix for 3-O-(b-D-Galactopyranosyl)-D-Arabinose
  x <- gsub("3-O-(b-D-Galactopyranosyl) D-", "3-O-(b-D-Galactopyranosyl)-D-", x)
  x <- gsub("3-O b-", "3-O-(b-", x)
  x <- gsub("Galactopyranosyl D-", "Galactopyranosyl)-D-", x)
  
  # Specific fix for sec-Butylamine
  x <- gsub("Butylamine sec", "sec-Butylamine", x)
  
  # Fix 5-Keto-d-Gluconic Acid
  x <- gsub("Ketod-", "Keto-d-", x)
  x <- gsub("KetoD-", "Keto-D-", x)
  
  # Fix Mono-Methyl Succinate
  x <- gsub("Mono-Methyl-Succinate", "Mono-Methyl Succinate", x)
  
  return(x)
}

# Function to wrangle results, i.e., move from wide to long format
wrangleResults <- function(data_wide){
  
  # split first column into Isolate and Replicate identifier
  data_wide$X <- gsub(" ", "_", data_wide$X)
  data_wide <- data_wide %>%
    separate(X, into = c("Isolate", "Replicate"), sep = "_replicate_", remove = FALSE)
  data_wide$X <- NULL
  
  data_long <- pivot_longer(data = data_wide, cols = c(3:98))
  return(data_long)
}
