#' #' @title get_fold_change
#' #' @description Get fold changes.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param data A matarix or data frame. Each row is a metabolite (or peak) and each column is a sample.
#' #' @param sample_info A data frame. Column 1 is 'sample.name', column 2 is 'injection order', column 3 is class, column 4
#' #' is 'batch' and column 5 is 'group'.
#' #' @param control_group Control group. It should be in the sample_info 'group' column.
#' #' @param case_group Case group. It should be in the sample_info 'group' column.
#' #' @param case_group Case group. It should be in the sample_info 'group' column.
#' #' @param by 'mean' or 'median'
#' #' @export
#' #' @return A numeric vector.
#' #' @import stats
#' #' @import tidyverse
#' #' @import dplyr
#' #' @import ggplot2
#' #' @import ggfortify
#'
#' setGeneric(
#'   name = "get_fold_change",
#'   def = function(data,
#'                  sample_info,
#'                  control_group,
#'                  case_group,
#'                  by = c("mean", "median")) {
#'     ##check data
#'     check_result <-
#'       check_data(data = data, sample_info = sample_info)
#'
#'     if (any(as.numeric(as.character(check_result$Error)) != 0)) {
#'       cat("Some error in your data or sample_info, plese check them.\n")
#'     }
#'
#'     if (!control_group %in% unique(sample_info$group) |
#'         !case_group %in% unique(sample_info$group)) {
#'       stop(control_group,
#'            " or ",
#'            case_group,
#'            " are not in your sample_info group column.\n")
#'     }
#'
#'     control_idx <-
#'       sample_info %>%
#'       dplyr::filter(group == control_group) %>%
#'       dplyr::pull(sample.name) %>%
#'       match(colnames(data))
#'
#'
#'     case_idx <-
#'       sample_info %>%
#'       dplyr::filter(group == case_group) %>%
#'       dplyr::pull(sample.name) %>%
#'       match(colnames(data))
#'
#'     if (by == "mean") {
#'       fold_change <-
#'         pbapply::pbapply(data, 1, function(x) {
#'           x <- as.numeric(x)
#'           mean(x[case_idx]) / mean(x[control_idx])
#'         })
#'     } else{
#'       fold_change <-
#'         pbapply::pbapply(data, 1, function(x) {
#'           x <- as.numeric(x)
#'           median(x[case_idx]) / median(x[control_idx])
#'         })
#'     }
#'
#'     fold_change[is.na(fold_change)] <- 1
#'     fold_change[is.infinite(fold_change)] <-
#'       max(fold_change[!is.na(fold_change)])
#'
#'     fold_change[fold_change == 0] <-
#'       min(fold_change[fold_change != 0])
#'     cat(crayon::bgRed("All done.\n"))
#'     return(fold_change)
#'   }
#' )
#'
#'
#'
#'
#'
#' #' @title get_p_value
#' #' @description Get p values.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param data A matarix or data frame. Each row is a metabolite (or peak) and each column is a sample.
#' #' @param sample_info A data frame. Column 1 is 'sample.name', column 2 is 'injection order', column 3 is class, column 4
#' #' is 'batch' and column 5 is 'group'.
#' #' @param control_group Control group. It should be in the sample_info 'group' column.
#' #' @param case_group Case group. It should be in the sample_info 'group' column.
#' #' @param case_group Case group. It should be in the sample_info 'group' column.
#' #' @param by 't.test' or 'wilcox.test'
#' #' @param ... Other arguments for 't.test' or 'wilcox.test'.
#' #' @export
#' #' @return A numeric vector.
#' #' @import stats
#' #' @import tidyverse
#' #' @import dplyr
#' #' @import ggplot2
#' #' @import ggfortify
#'
#'
#'
#' setGeneric(
#'   name = "get_p_value",
#'   def = function(data,
#'                  sample_info,
#'                  control_group,
#'                  case_group,
#'                  by = c("t.test", "wilcox.test"),
#'                  ...) {
#'     ##check data
#'     check_result <-
#'       check_data(data = data, sample_info = sample_info)
#'
#'     if (any(as.numeric(as.character(check_result$Error)) != 0)) {
#'       cat("Some error in your data or sample_info, plese check them.\n")
#'     }
#'
#'     if (!control_group %in% unique(sample_info$group) |
#'         !case_group %in% unique(sample_info$group)) {
#'       stop(control_group,
#'            " or ",
#'            case_group,
#'            " are not in your sample_info group column.\n")
#'     }
#'
#'     control_idx <-
#'       sample_info %>%
#'       dplyr::filter(group == control_group) %>%
#'       dplyr::pull(sample.name) %>%
#'       match(colnames(data))
#'
#'
#'     case_idx <-
#'       sample_info %>%
#'       dplyr::filter(group == case_group) %>%
#'       dplyr::pull(sample.name) %>%
#'       match(colnames(data))
#'
#'     if (by == "t.test") {
#'       p_value <-
#'         pbapply::pbapply(data, 1, function(x) {
#'           x <- as.numeric(x)
#'           t.test(x = x[case_idx], y = x[control_idx], ...)$p.value
#'         })
#'     } else{
#'       p_value <-
#'         pbapply::pbapply(data, 1, function(x) {
#'           x <- as.numeric(x)
#'           wilcox.test(x = x[case_idx], y = x[control_idx], ...)$p.value
#'         })
#'     }
#'     p_value[is.na(p_value)] <- 1
#'     cat(crayon::bgRed("All done.\n"))
#'     return(p_value)
#'   }
#' )
#'
#'
#'
#'
