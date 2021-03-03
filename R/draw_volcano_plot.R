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
#' # # data <- log(data + 1, 10)
#' # #
#' # # data <- apply(data, 1, function(x) {
#' # #   (x - mean(x)) / sd(x)
#' # # }) %>%
#' # #   t() %>%
#' # #   as.data.frame()
#' #
#' # ###p value
#' # p_value <- get_p_value(data = data,
#' #                        sample_info = sample_info,
#' #                        control_group = "0",
#' #                        case_group = "1",
#' #                        by = "t.test")
#' # ##fold change
#' # fold_change <- get_fold_change(data = data,
#' #                                sample_info = sample_info,
#' #                                control_group = "0",
#' #                                case_group = "1",
#' #                                by = "mean")
#' #
#' # variable_name <- rownames(data)
#' #
#' #
#' # draw_volcano_plot(p_value = p_value, fold_change = fold_change,
#' #                   variable_name = variable_name,
#' #                   fold_change_cutoff = 2, p_value_cutoff = 0.01, marker_label = TRUE)
#'
#'
#' #' @title draw_volcano_plot
#' #' @description Draw volcano plot.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param p_value A numeric vector.
#' #' @param fold_change A numeric vector.
#' #' @param variable_name A character vector.
#' #' @param fold_change_cutoff If you set as 2, the fold change cutoff will be 2 and 1/2.
#' #' @param p_value_cutoff P value cutoff.
#' #' @param xlab xlab.
#' #' @param ylab ylab.
#' #' @param marker_label Label markers or not.
#' #' @export
#' #' @return A ggplot object.
#' #' @import tidyverse
#' #' @import dplyr
#' #' @import ggplot2
#' #' @import ggrepel
#'
#'
#' setGeneric(
#'   name = "draw_volcano_plot",
#'   def = function(p_value,
#'                  fold_change,
#'                  variable_name,
#'                  fold_change_cutoff = 2,
#'                  p_value_cutoff = 0.05,
#'                  xlab = "Log(Fold change, 2)",
#'                  ylab = "-Log(P value, 10)",
#'                  marker_label = FALSE) {
#'
#'     if(length(unique(c(length(p_value), length(fold_change), length(variable_name)))) != 1){
#'       stop("p_value, fold_change and variable_name muse be same lenght.\n")
#'     }
#'
#'     if (sum(is.na(p_value)) != 0 | sum(is.na(fold_change)) != 0) {
#'       stop("NA in your p_value or fold_change.\n")
#'     }
#'
#'     temp_data <-
#'       data.frame(variable_name, p_value, fold_change, stringsAsFactors = FALSE) %>%
#'       dplyr::mutate(marker =
#'                       case_when(fold_change > fold_change_cutoff & p_value < p_value_cutoff ~ "Increase",
#'                                 fold_change < 1/fold_change_cutoff & p_value < p_value_cutoff ~ "Decrease",
#'                                 TRUE ~ "No"
#'                       ))
#'
#'   plot <-
#'     temp_data %>%
#'       ggplot(aes(log(fold_change, 2), -log(p_value, 10))) +
#'       geom_hline(yintercept = -log(p_value_cutoff, 10), color = "#D9D9D9") +
#'       geom_vline(xintercept = log(fold_change_cutoff, 2), color = "#D9D9D9") +
#'       geom_vline(xintercept = log(1/fold_change_cutoff, 2), color = "#D9D9D9") +
#'       geom_point(aes(color = marker), show.legend = FALSE) +
#'       scale_color_manual(values = c("No" = "#D9D9D9", "Increase" = "#FB8072", "Decrease" = "#80B1D3")) +
#'       theme_classic() +
#'       labs(x = xlab, y = ylab) +
#'       theme(
#'         axis.title = element_text(size = 15),
#'         axis.text = element_text(size = 12),
#'         legend.title = element_text(size = 15),
#'         legend.text = element_text(size = 12),
#'         strip.background = element_rect(fill = "#0099B47F"),
#'         strip.text = element_text(color = "white", size = 15)
#'       )
#'
#'   if(marker_label){
#'     plot <-
#'     plot +
#'       ggrepel::geom_label_repel(aes(log(fold_change, 2),
#'                                     -log(p_value, 10),
#'                                     label = variable_name),
#'                                 data = temp_data %>% dplyr::filter(marker != "No"))
#'   }
#'
#'   return(plot)
#'   }
#' )
