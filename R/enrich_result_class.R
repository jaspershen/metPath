#' An S4 class to represent pathways
#'
#' @slot pathway_database pathway_database
#' @slot pathway_version pathway_version
#' @slot result result
#' @exportClass enrich_result_class

setClass(
  Class = "enrich_result_class",
  representation(
    pathway_database = "character",
    pathway_version = "character",
    result = "data.frame"
  ),
  prototype = list(
    pathway_database = "",
    pathway_version = "",
    result = data.frame()
  )
)

setMethod(
  f = "show",
  signature = "enrich_result_class",
  definition = function(object) {
    pathway_version <- try(object@pathway_version, silent = TRUE)
    pathway_database <- try(object@pathway_database, silent = TRUE)
    if (class(version) != "try-error") {
      cat(crayon::green("---------Pathway database&version---------\n"))
      cat(crayon::green(pathway_database, "&", pathway_version, "\n"))
    }
    cat(crayon::green("-----------Enrichment result------------\n"))
    cat(crayon::green(nrow(object@result), "pathways are enriched", "\n"))
    cat(crayon::green(nrow(object@result %>% 
                             dplyr::filter(p_value < 0.05)), "pathways p-values < 0.05", "\n"))
    if (nrow(object@result) > 0) {
      if (nrow(object@result) > 5) {
        pathway_name = object@result$pathway_name[1:5]
        cat(crayon::green(
          paste(pathway_name, collapse = ";"),
          "... (only top 5 shows)\n"
        ))
      } else{
        pathway_name = object@result$pathway_name
        cat(crayon::green(paste(pathway_name, collapse = ";")))
      }
    }
  }
)



#' @title enrich_bar_plot
#' @description Barplot for enrich_result_class
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param object enrich_result_class object.
#' @param x_axis "p_value_adjust or "p_value"
#' @param cutoff cutoff for p_value_adjust or p_value.
#' @param top top pathways to show.
#' @export

enrich_bar_plot = function(object,
                           x_axis = c("p_value_adjust", "p_value"),
                           cutoff = 0.05,
                           top = 10) {
  x_axis = match.arg(x_axis)
  if (class(object) != "enrich_result_class") {
    stop("Only for enrich_result_class")
  }
  
  temp_data =
    object@result %>%
    dplyr::rename("x_value" = x_axis) %>%
    dplyr::mutate(x = -log(x_value, 10)) %>%
    dplyr::arrange(x) %>%
    dplyr::filter(x > -log(cutoff, 10)) %>%
    dplyr::mutate(pathway_name = factor(pathway_name, levels = pathway_name))
  
  if (nrow(temp_data) > top) {
    temp_data = temp_data[1:top, ]
  }
  
  plot =
    temp_data %>%
    ggplot(aes(x = x, y = pathway_name)) +
    geom_bar(stat = "identity", aes(fill = x_value), color = "black") +
    labs(
      y = "",
      x = ifelse(
        x_axis == "p_value_adjust",
        "-log(Adjusted P-values, 10)",
        "-log(P-values, 10)"
      )
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      panel.grid = element_blank()
    ) +
    guides(fill = guide_colorbar(
      title = ifelse(x_axis == "p_value_adjust",
                     "Adjusted P-values",
                     "P-values")
    )) +
    scale_fill_gradient(low = ggsci::pal_aaas()(n = 10)[2],
                        high = alpha(ggsci::pal_aaas()(n = 10)[2], 0.1))
  plot
  
  
}



#' @title enrich_bar_plot
#' @description Barplot for enrich_result_class
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param object enrich_result_class object.
#' @param x_axis "mapped_number or "mapped_percentage"
#' @param y_axis "p_value_adjust or "p_value"
#' @param point_size "all_number" or "mapped_percentage"
#' @param x_axis_cutoff x_axis_cutoff
#' @param y_axis_cutoff y_axis_cutoff
#' @param label label
#' @param label_size label size
#' @export

# enrich_scatter_plot(object = object, x_axis = "mapped_percentage", y_axis_cutoff = 0.5)

enrich_scatter_plot = function(object,
                               x_axis = c("mapped_percentage", "mapped_number"),
                               y_axis = c("p_value_adjust", "p_value"),
                               point_size = c("mapped_percentage", "all_number"),
                               x_axis_cutoff = 0,
                               y_axis_cutoff = 0.05,
                               label = TRUE,
                               label_size = 4) {
  x_axis = match.arg(x_axis)
  y_axis = match.arg(y_axis)
  point_size = match.arg(point_size)
  
  if (class(object) != "enrich_result_class") {
    stop("Only for enrich_result_class")
  }
  
  temp_data =
    object@result %>%
    dplyr::mutate("x_axis" = object@result[, x_axis]) %>%
    dplyr::mutate("y_axis" = object@result[, y_axis]) %>%
    dplyr::mutate("point_size" = object@result[, point_size]) %>%
    # dplyr::rename("x_axis" = x_axis) %>%
    # dplyr::rename("y_axis" = y_axis) %>%
    # dplyr::rename("point_size" = point_size) %>%
    dplyr::select(pathway_name, x_axis, y_axis, point_size) %>%
    dplyr::mutate(y_axis = -log(y_axis, 10)) %>%
    dplyr::mutate(class = case_when(
      y_axis > -log(y_axis_cutoff, 10) & x_axis > x_axis_cutoff ~ "yes",
      TRUE ~ "no"
    ))
  
  plot =
    temp_data %>%
    ggplot(aes(x = x_axis, y = y_axis, size = point_size)) +
    geom_hline(
      yintercept = -log(y_axis_cutoff, 10),
      color = "grey",
      linetype = 2
    ) +
    geom_vline(xintercept = x_axis_cutoff,
               color = "grey",
               linetype = 2) +
    geom_point(aes(size = point_size, fill = class),
               show.legend = FALSE,
               shape = 21) +
    scale_fill_manual(values = c("yes" = ggsci::pal_aaas()(n = 10)[2], "no" = "grey")) +
    labs(
      y = ifelse(
        y_axis == "p_value_adjust",
        "-log(Adjusted P-values, 10)",
        "-log(P-values, 10)"
      ),
      x = ifelse(
        x_axis == "mapped_number",
        "Mapped number",
        "Mapped percentage (%)"
      )
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      panel.grid = element_blank()
    ) +
    guides(size = guide_legend(
      title = ifelse(
        point_size == "all_number",
        "Pathway size",
        "Mapped percentage (%)"
      )
    ))
  
  if (label) {
    plot =
      plot +
      ggrepel::geom_text_repel(
        data = temp_data %>%
          dplyr::filter(class == "yes"),
        aes(x = x_axis,
            y = y_axis,
            label = pathway_name),
        size = label_size
      )
  }
  
  plot
  
  
}


#' @title enrich_network
#' @description Network for enriched pathways
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param object enrich_result_class object.
#' @param point_size point_size.
#' @param label label
#' @param label_size label size
#' @param p_cutoff p_cutoff
#' @param only_significant_pathway only_significant_pathway
#' @param threads threads
#' @export

enrich_network = function(object,
                          point_size = c("p_value_adjust", "p_value"),
                          label = TRUE,
                          label_size = 4,
                          p_cutoff = 0.05,
                          only_significant_pathway = TRUE,
                          threads = 4) {
  # browser()
  if (class(object) != "enrich_result_class") {
    stop("Only for enrich_result_class")
  }

  point_size = match.arg(point_size)
  
  temp_data =
    object@result %>%
    dplyr::rename("point_size" = point_size) %>%
    dplyr::mutate(point_size = -log(point_size, 10)) %>%
    dplyr::mutate(class = case_when(
      point_size > -log(p_cutoff, 10) ~ "yes",
      TRUE ~ "no"
    ))
  
  if(only_significant_pathway){
    temp_data =
      temp_data %>%
      dplyr::filter(class == "yes")
  }
  
  node_data = 
    temp_data %>% 
    dplyr::select(pathway_name, mapped_id, point_size, class)
  
  future::plan(strategy = future::multisession, workers = threads)
  edge_data = 
    furrr::future_map(1:(length(node_data$pathway_name)-1), 
                      .f = function(idx1){
      id1 = stringr::str_split(temp_data$mapped_id[idx1], ";")[[1]]
      temp = 
      furrr::future_map(.x = (idx1+1):length(node_data$pathway_name), 
                        .f = function(idx2){
        id2 = stringr::str_split(temp_data$mapped_id[idx2], ";")[[1]]
        length(intersect(id1, idx2))/length(union(id1, id2))
      }) %>% 
        unlist()
      
      names(temp) = node_data$pathway_name[(idx1+1):length(node_data$pathway_name)]
      temp
    })
  
  names(edge_data) = node_data$pathway_name[1:(length(node_data$pathway_name)-1)]
  
  edge_data = 
  furrr::future_map2(
    .x = names(edge_data),
    .y = edge_data,
    .f = function(x, y) {
      data.frame(from = x, to = names(y), jaccard_index = y)
    }
  ) %>% 
    do.call(rbind, .) %>% 
    as.data.frame()
  
  edge_data = 
    edge_data %>% 
    dplyr::filter(jaccard_index > 0)
  
  graph_data =
    tidygraph::tbl_graph(nodes = node_data,
                         edges = edge_data,
                         directed = FALSE)
  
  plot = 
  ggraph::ggraph(graph_data,
         layout = 'kk') +
    ggraph::geom_edge_link(
                  aes(edge_width = jaccard_index),
                  alpha = 1,
                  color = "black",
                  show.legend = TRUE) +
    ggraph::geom_node_point(aes(fill = class, 
                        size = point_size), 
                    shape = 21, 
                    alpha = 1, 
                    show.legend = TRUE) +
    guides(size = guide_legend(
      title =
        ifelse(
          point_size == "p_value_adjust",
          "-log(Adjusted P-values, 10)",
          "-log(P-values, 10)"
        )
    ), fill = FALSE) +
    ggraph::scale_edge_width(range = c(0.5,2)) +
    scale_size_continuous(range = c(1, 8)) +
    scale_fill_manual(values = c("yes" = ggsci::pal_aaas()(n = 10)[2], 
                                 "no" = "grey")) +
    ggraph::theme_graph() +
    theme(plot.background = element_rect(fill = "transparent", color = NA), 
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA)) 
  
  if (label) {
    plot =
      plot +
      ggraph::geom_node_text(
        aes(
          x = x,
          y = y,
          label = ifelse(class == "yes", pathway_name, NA),
          color = class
        ),
        repel = TRUE,
        size = label_size,
        show.legend = FALSE
      ) +
      scale_color_manual(values = c("yes" = ggsci::pal_aaas()(n = 10)[2],
                                    "no" = "grey")) +
      guides(color = FALSE)
  }
  plot
}