save_object <- function(object, file, file_format="lz4"){

  stopifnot(file_format %in% c("zstd", "lz4", "gzip", "bzip2", "xz", "nocomp"))

  if(file_format %in% "nocomp"){
    saveRDS(object = object, file = file, compress = FALSE)
    return(invisible(NULL))
  }

  if(file_format %in% c("zstd", "lz4")){
    con <- archive::file_write(file = file, filter = file_format)
    open(con)
    saveRDS(object = object, file = con)
    close(con)
  }else{
    saveRDS(object = object, file = file, compress = file_format)
  }
}

load_object <- function(file){
  con <- archive::file_read(file = file)
  res <- readRDS(file = con)
  close(con)
  return(res)
}


