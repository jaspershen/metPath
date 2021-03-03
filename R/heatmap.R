#' # data <- object4@ms1.data[[1]]
#' # sample_info <- object4@sample.info
#' # sample_info$group
#' #
#' # tags <-
#' #   data %>%
#' #   dplyr::select(-sample_info$sample.name)
#' #
#' # data <-
#' #   data %>%
#' #   dplyr::select(sample_info$sample.name)
#' #
#' # sample_info <-
#' #   sample_info %>%
#' #   dplyr::filter(group != "QC", batch == 2)
#' #
#' # sample_info <-
#' #   sample_info[sample(1:150, 40),]
#' #
#' # data <-
#' #   data %>%
#' #   dplyr::select(one_of(sample_info$sample.name))
#' #
#' # p_value <-
#' #   get_p_value(
#' #     data = data,
#' #     sample_info = sample_info,
#' #     control_group = "0",
#' #     case_group = "1",
#' #     by = "mean"
#' #   )
#' # fold_change <-
#' #   get_fold_change(
#' #     data = data,
#' #     sample_info = sample_info,
#' #     control_group = "0",
#' #     case_group = "1",
#' #     by = "t.test"
#' #   )
#' #
#' # marker_idx <- which(p_value < 0.05 & fold_change > 1.5 |
#' #                       p_value < 0.05 & fold_change < 1 / 1.5)
#' #
#' # data <- data[marker_idx,]
#' #
#' #
#' # ##log
#' # data <- log(data + 1, 10)
#' #
#' # data <- apply(data, 1, function(x) {
#' #   (x - mean(x)) / sd(x)
#' # }) %>%
#' #   t() %>%
#' #   as.data.frame()
#'
#' #' @title draw_heatmap
#' #' @description Heatmap
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param data A matarix or data frame. Each row is a metabolite (or peak) and each column is a sample.
#' #' @param ... Other arguments for pheatmap.
#' #' @export
#' #' @return A pheatmap plot.
#' #' @import tidyverse
#' #' @import dplyr
#' #' @import ggplot2
#' #' @import pheatmap
#'
#'
#' setGeneric(
#'   name = "draw_heatmap",
#'   def = function(data,
#'                  ...) {
#'
#'     data_range <- abs(range(data))
#'     dif <- data_range[1] - data_range[2]
#'
#'     if (dif < 0) {
#'       data[data > data_range[1]] <- data_range[1]
#'     }
#'
#'     if (dif > 0) {
#'       data[data < -1 * data_range[2]] <- -1 * data_range[2]
#'     }
#'
#'     plot <-
#'     pheatmap::pheatmap(mat = data)
#'     return(plot)
#'   }
#' )
