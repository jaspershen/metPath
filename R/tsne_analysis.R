#' #' @title do_tsne
#' #' @description tSNE analysis for metflowClass object.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@sioc.ac.cn}
#' #' @param data A matarix or data frame. Each row is a metabolite (or peak) and each column is a sample.
#' #' @param dims The number of dimensions. See ?Rtsne.
#' #' @param perplexity perplexity. See ?Rtsne.
#' #' @param verbose verbose. See ?Rtsne.
#' #' @import ggfortify
#' #' @import Rtsne
#' #' @import tidyverse
#' #' @import ggplot2
#' #' @return A tsne object.
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
#' #
#' # tsne_object <-
#' #   do_tsne(data = data)
#' #
#' # sample_info$batch <- as.factor(sample_info$batch)
#' # draw_tsne_plot(tsne_object = tsne_object, sample_info = sample_info, color_index = "batch")
#' #
#' # draw_tsne_plot(tsne_object = tsne_object, sample_info = sample_info, color_index = "injection.order")
#'
#'
#' setGeneric(name = "do_tsne",
#'            function(data,
#'                     dims = 2,
#'                     perplexity = 30,
#'                     verbose = TRUE) {
#'
#'              if (sum(is.na(data)) != 0) {
#'                stop("Please impute MV first.\n")
#'              }
#'
#'              data <-
#'                t(data)
#'
#'              tsne_object <- Rtsne::Rtsne(
#'                X = as.matrix(data),
#'                dims = dims,
#'                perplexity = perplexity,
#'                verbose = verbose
#'              )
#'              return(tsne_object)
#'            })
#'
#' #' @title draw_tsne_plot
#' #' @description Draw tsne score plot.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param tsne_object A object from do_tsne() function.
#' #' @param sample_info A data frame. Column 1 is 'sample.name', column 2 is 'injection order', column 3 is class, column 4
#' #' is 'batch' and column 5 is 'group'.
#' #' @param color_index color_index. Which column in your sample information. Such as 'group', 'batch' or 'class'.
#' #' @export
#' #' @return A ggplot object.
#' #' @import stats
#' #' @import tidyverse
#' #' @import dplyr
#' #' @import ggplot2
#' #' @import ggfortify
#'
#' setGeneric(
#'   name = "draw_tsne_plot",
#'   def = function(tsne_object,
#'                  sample_info,
#'                  color_index) {
#'
#'     if(all(colnames(sample_info) != color_index)){
#'       stop("No ", color_index, " in your sample_info.\n")
#'     }
#'
#'     class <- sample_info %>%
#'       pull(color_index)
#'
#'       Y <- tsne_object$Y
#'
#'       Y <-
#'         data.frame(Y, "Class" = class, stringsAsFactors = FALSE)
#'
#'       if(class(class) == "numeric"){
#'         plot <-
#'           ggplot(Y, aes(X1, X2, color = Class)) +
#'           geom_hline(yintercept = 0, color = "grey") +
#'           geom_vline(xintercept = 0, color = "grey") +
#'           geom_point() +
#'           labs(x = "Dimension 1",
#'                y = "Dimension 2") +
#'           theme_classic() +
#'           scale_colour_gradientn(colours = c(
#'             alpha("#8DD3C7", 1),
#'             alpha("#8DD3C7", 0.4),
#'             alpha("#FB8072", 0.4),
#'             alpha("#FB8072", 1)
#'           ))
#'       }else{
#'         plot <-
#'         ggplot(Y, aes(X1, X2, color = Class)) +
#'           geom_hline(yintercept = 0, color = "grey") +
#'           geom_vline(xintercept = 0, color = "grey") +
#'           geom_point() +
#'           labs(x = "Dimension 1",
#'                y = "Dimension 2") +
#'           theme_classic() +
#'           ggsci::scale_color_aaas()
#'
#'       }
#'
#'       plot <-
#'         plot +
#'         theme(
#'           axis.title = element_text(size = 15),
#'           axis.text = element_text(size = 12),
#'           legend.title = element_text(size = 15),
#'           legend.text = element_text(size = 12),
#'           strip.background = element_rect(fill = "#0099B47F"),
#'           strip.text = element_text(color = "white", size = 15)
#'         )
#'
#'
#'
#'
#'     return(plot)
#'   }
#' )
