#' @title creatMetPathObject
#' @description Creat metPathClass object.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param ms1.data MS1 peak table name.
#' @param sample.information Sample information name.
#' @param path Work directory.
#' @return A metPathClass object.
#' @export

setGeneric(
  name = "creatMetPathObject",
  def = function(ms1.data,
                 sample.information,
                 path = ".") {
    check.result <- checkData(data = ms1.data,
                              sample.info = sample.information,
                              path = path)

    if (any(check.result$`Check result` != "Valid")) {
      stop("Error is in you data. Please check it.\n")
    }

    ms1.data.name <- ms1.data
    sample.information.name <- sample.information

    ms1.data <- pbapply::pblapply(ms1.data, function(x) {
      readr::read_csv(file = file.path(path, x),
                      col_types = readr::cols())
    })

    sample.info <-
      readr::read_csv(file.path(path, sample.information),
                      col_types = readr::cols())

    object <- new(
      Class = "metPathClass",
      ms1.data = ms1.data,
      sample.info = sample.info,
      version = "0.1.0"
    )
    invisible(object)
  }
)


#' An S4 class to represent a metPath object.
#'
#' @slot ms1.data MS1 data.
#' @slot sample.info Sample information.
#' @slot process.info MS1 data.
#' @slot version MS1 data.

##S4 class for function metIdentification
setClass(
  Class = "metPathClass",
  representation(
    ms1.data = "list",
    sample.info = "data.frame",
    process.info = "list",
    version = "character"
  )
)

setMethod(
  f = "show",
  signature = "metPathClass",
  definition = function(object) {
    # requireNamespace("magrittr")
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat(crayon::green("metPath version:", object@version, "\n"))
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat(crayon::green("MS1 data\n"))
    cat("There are",
        length(object@ms1.data),
        "peak tables in your MS1 data.\n")
    info <- lapply(object@ms1.data, dim) %>%
      do.call(rbind, .)

    colnames(info) <- c("Peak.number", "Column.number")
    rownames(info) <- paste("Batch", 1:nrow(info), sep = "")
    print(info)
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat("There are",
        nrow(object@sample.info),
        "samples in your MS1 data.\n")
    class_info <-
      as.data.frame(table(object@sample.info$class), stringsAsFactors = FALSE)
    colnames(class_info) <- c("Class", "Number")
    print(class_info)
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    group_info <-
      as.data.frame(table(object@sample.info$group), stringsAsFactors = FALSE)
    colnames(group_info) <- c("Group", "Number")
    print(group_info)
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat(crayon::green("Processing\n"))
    if (.hasSlot(object = object, name = "process.info") &
        length(object@process.info) != 0) {
      process.info <- object@process.info
      mapply(function(x, y) {
        cat(crayon::green(x, paste(rep("-", 10), collapse = ""), "\n"))
        y <- y[which(names(y) != "plot")]
        y <- data.frame(names(unlist(y)), unlist(y))
        colnames(y) <- c("Parameter", "Value")
        rownames(y) <- NULL
        print(y)
      },
      x = names(process.info),
      y = process.info)
    } else{
      cat(crayon::red("There are no processing for your data.\n"))
    }
  }
)
