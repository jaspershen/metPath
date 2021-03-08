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
# load("data/hmdb_pathway.rda")
# load("data/primary_pathway")
# 
# idx = match(primary_pathway, hmdb_pathway@pathway_id)
# idx = idx[!is.na(idx)]
# pathway_class =
#   hmdb_pathway@pathway_class
# 
# for(x in idx){
#   pathway_class[[x]] = paste(pathway_class[[x]], "primary_pathway", sep = ";")
# }
# 
# 
# 
# pathway_class[idx]
# 
# hmdb_pathway@pathway_class = pathway_class
# 
# save(hmdb_pathway, file = "data/hmdb_pathway.rda", )

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
