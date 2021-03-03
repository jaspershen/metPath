#' #' @title check_data
#' #' @description Check data format.
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@163.com}
#' #' @param data MS1 data tables.
#' #' @param sample_info Sample information.
#' #' @return Notice of data checking.
#' #' @export
#'
#' setGeneric(
#'   name = "check_data",
#'   def = function(data,
#'                  sample_info) {
#'
#'     sample_info_record <- NULL
#'     data <- list(data)
#'
#'     sample_info <- sample_info
#'
#'     sample_info <- as.data.frame(sample_info)
#'     #####
#'
#'     ##check data, component
#'     data_record <- lapply(seq_along(data), function(x) {
#'       temp_data_record <- NULL
#'       data.col.name <- colnames(data[[x]])
#'       # if (all(data.col.name != "name")) {
#'       #   cat(crayon::red(clisymbols::symbol$cross,"Error: No name in data.\n"))
#'       #   temp_data_record <- c(temp_data_record, "Error")
#'       # } else{
#'       #   # cat("OK: The first column of data is name.\n")
#'       #   temp_data_record <- c(temp_data_record, "OK")
#'       # }
#'
#'       # if (all(data.col.name != "mz")) {
#'       #   cat(crayon::red(clisymbols::symbol$cross,"Error: No mz in data.\n"))
#'       #   temp_data_record <- c(temp_data_record, "Error")
#'       # } else{
#'       #   # cat("OK: The second column of data is mz.\n")
#'       #   temp_data_record <- c(temp_data_record, "OK")
#'       # }
#'
#'       # if (all(data.col.name != "rt")) {
#'       #   cat(crayon::red(clisymbols::symbol$cross,"Error: No rt in data.\n"))
#'       #   temp_data_record <- c(temp_data_record, "Error")
#'       # } else{
#'       #   # cat("OK: The third column of data is not rt.\n")
#'       #   temp_data_record <- c(temp_data_record, "OK")
#'       # }
#'       temp_data_record
#'     })
#'
#'     cat("--------------------------------------------------------------\n")
#'     ##check sample_info
#'     if (sum(is.na(sample_info)) > 0) {
#'       cat(crayon::red(clisymbols::symbol$cross,
#'                       "Error: There are", sum(is.na(sample_info)), "NAs in your sample_info.\n"))
#'       sample_info_record <- c(sample_info_record, "Error")
#'     } else{
#'       # cat("OK: There are no NAs in you sample_info.\n")
#'       sample_info_record <- c(sample_info_record, "OK")
#'     }
#'
#'     if (ifelse(is.na(sum(sample_info == "") > 0), FALSE, sum(sample_info == "") > 0)) {
#'       cat(crayon::red(clisymbols::symbol$cross,
#'                       "Error: There are",
#'                       sum(sample_info == ""),
#'                       "spaces in you sample_info.\n"))
#'       sample_info_record <- c(sample_info_record, "Error")
#'     } else{
#'       sample_info_record <- c(sample_info_record, "OK")
#'     }
#'
#'     if (colnames(sample_info)[1] != "sample.name") {
#'       cat(crayon::red(clisymbols::symbol$cross,
#'                       "Error: The column 1 of sample information must be sample.name\n"))
#'       sample_info_record <- c(sample_info_record, "Error")
#'     } else{
#'       # cat("OK: There are no spaces in you sample_info.\n")
#'       sample_info_record <- c(sample_info_record, "OK")
#'     }
#'
#'     if (colnames(sample_info)[2] != "injection.order") {
#'       cat(crayon::red(clisymbols::symbol$cross,
#'                       "Error: The column 2 of sample information must be injection.order\n"))
#'       sample_info_record <- c(sample_info_record, "Error")
#'     } else{
#'       # cat("OK: There are no spaces in you sample_info.\n")
#'       sample_info_record <- c(sample_info_record, "OK")
#'     }
#'
#'     if (colnames(sample_info)[3] != "class") {
#'       cat(crayon::red(clisymbols::symbol$cross,
#'                       "Error: The column 3 of sample information must be class"))
#'       sample_info_record <- c(sample_info_record, "Error")
#'     } else{
#'       # cat("OK: There are no spaces in you sample_info.\n")
#'       sample_info_record <- c(sample_info_record, "OK")
#'     }
#'
#'     if (colnames(sample_info)[4] != "batch") {
#'       cat(crayon::red(clisymbols::symbol$cross,
#'                       "Error: The column 4 of sample information must be batch"))
#'       sample_info_record <- c(sample_info_record, "Error")
#'     } else{
#'       # cat("OK: There are no spaces in you sample_info.\n")
#'       sample_info_record <- c(sample_info_record, "OK")
#'       if (any(table(sample_info[, "batch"]) <= 5)) {
#'         cat(crayon::yellow(clisymbols::symbol$warning,
#'                            "There are some batch only have less than 5 samples. Please check it.\n"))
#'         sample_info_record <-
#'           c(sample_info_record, "Error")
#'       }
#'     }
#'
#'     if (colnames(sample_info)[5] != "group") {
#'       cat(crayon::red(clisymbols::symbol$cross,
#'                       "Error: The column 5 of sample information must be group"))
#'       sample_info_record <- c(sample_info_record, "Error")
#'     } else{
#'       # cat("OK: There are no spaces in you sample_info.\n")
#'       sample_info_record <- c(sample_info_record, "OK")
#'     }
#'
#'     class <-  unique(as.character(sample_info[, 3]))
#'     if (all(c("QC", "Subject") %in% class)) {
#'       sample_info_record <- c(sample_info_record, "OK")
#'     } else{
#'       cat(crayon::red(clisymbols::symbol$cross,
#'                       "Error: The class of sample_information must be QC and Subject.\n"))
#'       sample_info_record <- c(sample_info_record, "Error")
#'     }
#'
#'     sample.idx <- match(sample_info$sample.name,
#'                         unique(unlist(lapply(data, function(x)
#'                           colnames(x)))))
#'
#'     if (any(is.na(sample.idx))) {
#'       cat(crayon::red(clisymbols::symbol$cross,
#'                       "Error: The sample names in sample_inforamtion and data are not same.\n"))
#'       sample_info_record <- c(sample_info_record, "Error")
#'     } else{
#'       sample_info_record <- c(sample_info_record, "OK")
#'     }
#'
#'     cat("--------------------------------------------------------------\n")
#'     cat("Summary:\n")
#'
#'     stat <- lapply(c(data_record, list(sample_info_record)),
#'                    function(x) {
#'                      return(c(
#'                        ifelse(all(x != "Error"), "Valid", "Invalid"),
#'                        sum(x == "OK"),
#'                        sum(x == "Warning"),
#'                        sum(x == "Error")
#'                      ))
#'                    })
#'     stat <- do.call(rbind, stat)
#'     stat <- as.data.frame(stat)
#'     colnames(stat) <-
#'       c("Check result", "OK", "Warning", "Error")
#'     rownames(stat)[nrow(stat)] <- "sample_info"
#'     rownames(stat)[1:(nrow(stat) - 1)] <-
#'       paste("batch", 1:(nrow(stat) - 1), sep = "")
#'
#'     print(stat, quote = FALSE)
#'     cat("\n")
#'     cat("data:\n")
#'     for (idx in 1:length(data_record)) {
#'       if (all(data_record[idx] != "Error")) {
#'         cat(paste("Batch", idx), "is valid.\n")
#'       } else{
#'         if (sum(data_record[idx] == "Warning") > 0) {
#'           cat(
#'             "There",
#'             ifelse(sum(data_record[idx] == "Warning") > 1, "are", "is"),
#'             sum(data_record[idx] == "Warning"),
#'             ifelse(sum(data_record[idx] == "Warning") > 1, "Warnings", "Warning"),
#'             "in your",
#'             paste("Batch", idx),
#'             ".Please check it according to the information.\n"
#'           )
#'         }
#'
#'         if (sum(data_record[idx] == "Error") > 0) {
#'           cat(
#'             "There",
#'             ifelse(sum(data_record[idx] == "Error") > 1, "are", "is"),
#'             sum(data_record[idx] == "Error"),
#'             ifelse(sum(data_record[idx] == "Error") > 1, "Errors", "Error"),
#'             "in your",
#'             paste("Batch", idx),
#'             ".Please check it according to the information.\n"
#'           )
#'         }
#'       }
#'     }
#'
#'     cat("\n")
#'     cat("sample_info:\n")
#'     if (all(sample_info_record != "Error")) {
#'       cat("sample_info is valid.\n")
#'     } else{
#'       if (sum(sample_info_record == "Warning") > 0) {
#'         cat(
#'           "There",
#'           ifelse(sum(sample_info_record == "Warning") > 1, "are", "is"),
#'           sum(sample_info_record == "Warning"),
#'           ifelse(
#'             sum(sample_info_record == "Warning") > 1,
#'             "Warnings",
#'             "Warning"
#'           ),
#'           "in your sample_info. Please check it according to the information.\n"
#'         )
#'       }
#'
#'       if (sum(sample_info_record == "Error") > 0) {
#'         cat(
#'           "There",
#'           ifelse(sum(sample_info_record == "Error") > 1, "are", "is"),
#'           sum(sample_info_record == "Error"),
#'           ifelse(sum(sample_info_record == "Error") > 1, "Errors", "Error"),
#'           "in your sample_info. Please check it according to the information.\n"
#'         )
#'       }
#'     }
#'     stat <- as.data.frame(stat)
#'     stat[, 1] <- as.character(stat[, 1])
#'     stat[, 2] <- as.numeric(stat[, 2])
#'     stat[, 3] <- as.numeric(stat[, 3])
#'     invisible(stat)
#'   })
#'
#'
