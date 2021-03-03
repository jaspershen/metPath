#' #' @title do_pca
#' #' @description PCA analysis for MetFlowData.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param data A matarix or data frame. Each row is a metabolite (or peak) and each column is a sample.
#' #' @export
#' #' @return A pca object object.
#' #' @import stats
#' #' @import tidyverse
#' #' @import dplyr
#' #' @import ggplot2
#' #' @import ggfortify
#'
#'
#'
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
#' # ##log
#' # data <- log(data + 1, 10)
#' #
#' # data <- apply(data, 1, function(x) {
#' #   (x - mean(x)) / sd(x)
#' # }) %>%
#' #   t() %>%
#' #   as.data.frame()
#' #
#' # sample_info$batch <- as.character(sample_info$batch)
#' #
#' # object <- do_pca(
#' #   data = data,
#' #   sample_info = sample_info
#' # )
#' #
#' # draw_pca_plot(pca_object = object, sample_info = sample_info,
#' #               color_index = "batch",
#' #               loadings = TRUE,
#' #               loadings_color = "grey")
#'
#'
#' setGeneric(
#'   name = "do_pca",
#'   def = function(data,
#'                  sample_info) {
#'     if (any(check_result$Error != 0)) {
#'       stop("Some error in your data or sample_info, please check them and try again.\n")
#'     }
#'
#'     if (sum(is.na(data)) != 0) {
#'       stop("Please impute MV first.\n")
#'     }
#'
#'     data <-
#'       t(data)
#'
#'     pca_object <- prcomp(data, center = FALSE, scale. = FALSE)
#'     cat(crayon::bgRed("All done.\n"))
#'     return(pca_object)
#'   }
#' )
#'
#'
#'
#' #' @title draw_pca_plot
#' #' @description Draw pca score plot.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param pca_object A procmp object
#' #' @param sample_info A data frame. Column 1 is 'sample.name', column 2 is 'injection order', column 3 is class, column 4
#' #' is 'batch' and column 5 is 'group'.
#' #' @param color_index color_index. Which column in your sample information. Such as 'group', 'batch' or 'class'.
#' #' @param loadings Show lodings or not.
#' #' @param loadings_color Loading color.
#' #' @param loadings_label Show loading labels or not.
#' #' @param loadings_label_size loading label size.
#' #' @export
#' #' @return A ggplot object.
#' #' @import stats
#' #' @import tidyverse
#' #' @import dplyr
#' #' @import ggplot2
#' #' @import ggfortify
#'
#' setGeneric(
#'   name = "draw_pca_plot",
#'   def = function(pca_object,
#'                  sample_info,
#'                  color_index,
#'                  loadings = FALSE,
#'                  loadings_color = 'blue',
#'                  loadings_label = FALSE,
#'                  loadings_label_size = 3) {
#'
#'     if (class(pca_object) != "prcomp") {
#'       stop("Only for procomp object.\n")
#'     }
#'
#'     if(all(colnames(sample_info) != color_index)){
#'       stop("No ", color_index, " in your sample_info.\n")
#'     }
#'
#'     class <- sample_info %>%
#'       pull(color_index)
#'
#'     if (class(class) == "numeric") {
#'       plot <-
#'         ggplot2::autoplot(
#'           pca_object,
#'           data = sample_info,
#'           colour = color_index,
#'           frame = TRUE,
#'           frame.type = "norm",
#'           loadings = loadings,
#'           loadings.colour = loadings_color,
#'           loadings.label = loadings_label,
#'           loadings.label.size = loadings_label_size
#'         ) +
#'         geom_hline(yintercept = 0, color = "grey") +
#'         geom_vline(xintercept = 0, color = "grey") +
#'         scale_colour_gradientn(colours = c(
#'           alpha("#8DD3C7", 1),
#'           alpha("#8DD3C7", 0.4),
#'           alpha("#FB8072", 0.4),
#'           alpha("#FB8072", 1)
#'         ))
#'     }else{
#'       plot <-
#'         ggplot2::autoplot(
#'           pca_object,
#'           data = sample_info,
#'           colour = color_index,
#'           frame = TRUE,
#'           frame.type = "norm",
#'           loadings = loadings,
#'           loadings.colour = loadings_color,
#'           loadings.label = loadings_label,
#'           loadings.label.size = loadings_label_size
#'         ) +
#'         geom_hline(yintercept = 0, color = "grey") +
#'         geom_vline(xintercept = 0, color = "grey") +
#'         ggsci::scale_color_aaas() +
#'         ggsci::scale_fill_aaas()
#'     }
#'
#'     plot <-
#'       plot +
#'       theme_classic() +
#'       theme(
#'         axis.title = element_text(size = 15),
#'         axis.text = element_text(size = 12),
#'         legend.title = element_text(size = 15),
#'         legend.text = element_text(size = 12),
#'         strip.background = element_rect(fill = "#0099B47F"),
#'         strip.text = element_text(color = "white", size = 15)
#'       )
#'     return(plot)
#'   }
#' )
