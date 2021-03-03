#' @title get_hmdb_pathway
#' @description Get compound from HMDB (SMPDB)
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param threads threads
#' @export

get_hmdb_compound <- function(
  threads = 3
){
  data("hmdbMS1Database", envir = environment())
  message(
    crayon::yellow(
      "This database is downloaded in",
      hmdbMS1Database@database.info$Version
    )
  )
  cat("\n")
  return(hmdbMS1Database)
}


#' @title get_hmdb_pathway
#' @description Get pathways from HMDB (SMPDB)
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param threads threads
#' @export

get_hmdb_pathway <- function(threads = 3) {
  data("hmdb_pathway", envir = environment())
  message(
    crayon::yellow(
      "This database is downloaded in",
      hmdb_pathway@database_info$version
    )
  )
  cat("\n")
  return(hmdb_pathway)
}
