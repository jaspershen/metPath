#' An S4 class to represent pathways
#'
#' @slot database_info database_info
#' @slot pathway_id pathway_id
#' @slot pathway_name pathway_name
#' @slot describtion describtion
#' @slot pathway_class pathway_class
#' @slot gene_list gene_list
#' @slot compound_list compound_list
#' @slot protein_list protein_list
#' @slot reference_list reference_list
#' @slot related_disease related_disease
#' @slot related_module related_module
#' @exportClass pathway_database_class

setClass(
  Class = "pathway_database_class",
  representation(
    database_info = "list",
    pathway_id = "vector",
    pathway_name = "vector",
    describtion = "list",
    pathway_class = "list",
    gene_list = "list",
    compound_list = "list",
    protein_list = "list",
    reference_list = "list",
    related_disease = "list",
    related_module = "list"
  ),
  prototype = list(
    database_info = list(),
    pathway_id = list(),
    pathway_name = list(),
    describtion = list(),
    pathway_class = list(),
    gene_list = list(),
    compound_list = list(),
    protein_list = list(),
    reference_list = list(),
    related_disease = list(),
    related_module = list()
  )
)

setMethod(
  f = "show",
  signature = "pathway_database_class",
  definition = function(object) {
    version <- try(object@database_info$version, silent = TRUE)
    source <- try(object@database_info$source, silent = TRUE)
    if (class(version) != "try-error") {
      cat(crayon::green("---------Pathway source&version---------\n"))
      cat(crayon::green(source, "&", version, "\n"))
    }
    cat(crayon::green("-----------Pathway information------------\n"))
    cat(crayon::green(length(object@pathway_id), "pathways", "\n"))
    cat(
      crayon::green(
        object@gene_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        "pathways have genes",
        "\n"
      )
    )
    
    cat(
      crayon::green(
        object@protein_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        "pathways have proteins",
        "\n"
      )
    )
    
    cat(
      crayon::green(
        object@compound_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        "pathways have compounds",
        "\n"
      )
    )
    
    cat(crayon::green("Pathway class:",
                      paste(unique(
                        unlist(object@pathway_class)
                      ), collapse = ";"),
                      "\n"))
  }
)




#' @title filter_pathway
#' @description filter pathways according to pathway class
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param object pathway_database_class object.
#' @param class class pathway class you want to remain.
#' @param remain_idx Which pathways you want to remain (remain_idx).
#' @export

filter_pathway =
  function(object,
           class,
           remain_idx) {
    if (base::class(object) != "pathway_database_class"){
      stop(crayon::red('Only for pathway_database_class object.\n'))
    }
    
    if (missing(class) & missing(remain_idx)) {
      stop("Provide class or remain_idx.\n")
    }
    
    if (!missing(class) & !missing(remain_idx)) {
      message(crayon::yellow("Only use remain_idx."))
      remain_idx = remain_idx
    }
    
    if (!missing(class) & missing(remain_idx)) {
      pathway_class = object@pathway_class %>%
        unlist() %>%
        unique()
      class = class[class %in% pathway_class]
      if (length(class) == 0) {
        stop("All class are not in pathway object.\n")
      }
      
      remain_idx =
        object@pathway_class %>%
        purrr::map(function(x) {
          x %in% class
        }) %>%
        unlist() %>%
        which()
    }
    
    if (missing(class) & !missing(remain_idx)) {
      remain_idx = remain_idx
    }
    
    object@pathway_id =
      object@pathway_id[remain_idx]
    
    object@pathway_name =
      object@pathway_name[remain_idx]
    
    object@describtion =
      object@describtion[remain_idx]
    
    object@pathway_class =
      object@pathway_class[remain_idx]
    
    if (length(object@gene_list) > 0) {
      object@gene_list =
        object@gene_list[remain_idx]
    }
    
    if (length(object@compound_list) > 0) {
      object@compound_list =
        object@compound_list[remain_idx]
    }
    
    if (length(object@protein_list) > 0) {
      object@protein_list =
        object@protein_list[remain_idx]
    }
    
    if (length(object@reference_list) > 0) {
      object@reference_list =
        object@reference_list[remain_idx]
    }
    
    if (length(object@related_disease) > 0) {
      object@related_disease =
        object@related_disease[remain_idx]
    }
    
    if (length(object@related_module) > 0) {
      object@related_module =
        object@related_module[remain_idx]
    }
    return(object)
  }


#' @title get_pathway_class
#' @description Extract the class of pathways.
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param object pathway_class object.
#' @export

get_pathway_class = function(object) {
  object@pathway_class %>%
    unlist() %>%
    data.frame(class = .) %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::ungroup()
}