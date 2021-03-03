# do_gsea2 <- function(feature_table,
#                      feature_set,
#                      exponent = 1,
#                      perm_num = 1000,
#                      min_size = 5,
#                      max_size = 1000,
#                      pvalue_cutoff = 0.2,
#                      p_adjust_method = "fdr",
#                      seed = FALSE,
#                      verbose = TRUE,
#                      ...) {
#   feature_set <-
#     filter_feature_set2(
#       feature_set = feature_set,
#       feature_table = feature_table,
#       min_size = min_size,
#       max_size = max_size
#     )
#
#   if (length(feature_set) == 0) {
#     return(NULL)
#   }
#
#
#   ##find where are the metabolites in the feature table
#   match_idx <-
#     purrr::map(feature_set[[1]], function(x) {
#       which(x == feature_table$Lab.ID)
#     })
#
#   ##for the metabolites in the feature set,
#   ##all the peaks should be remained, other peaks only remain one metabolite,
#   ##either one is fine
#   idx <-
#     sort(unique(unlist(match_idx)))
#
#   ##for the metabolites in feature sets
#   feature_table1 <-
#     feature_table[idx,] %>%
#     dplyr::filter(isotope == "[M]")
#
#   ##for each compound_class, if the peaks are in different fc_class,
#   ##then it should be classed into
#   ##different classes and then score again
#   ##
#   ##and if there are [M+H] or [M-H] in the compound_class, onlt remain them
#   #     feature_table1 <-
#   # feature_table1
#   #     dplyr::arrange(compound_class) %>%
#   #     dplyr::mutate(fc_class =
#   #                     dplyr::case_when(fc >= 0 ~ "increase",
#   #                                      TRUE ~ "decrease")) %>%
#   #     plyr::dlply(.variables = .(compound_class)) %>%
#   #     purrr::map(
#   #       .f = function(x) {
#   #         if (length(unique(x$fc_class)) != 1) {
#   #           new_score <-
#   #             x %>%
#   #             plyr::dlply(.variables = .(fc_class)) %>%
#   #             purrr::map(
#   #               .f = function(y) {
#   #                 score_peak_group(peak_group = y)
#   #               }
#   #             ) %>%
#   #             unlist()
#   #           x$score[x$fc_class == "decrease"] <- new_score[1]
#   #           x$score[x$fc_class == "increase"] <- new_score[1]
#   #           if (length(unique(new_score)) != 1) {
#   #             x <-
#   #               x %>%
#   #               dplyr::filter(score != min(new_score))
#   #           }
#   #         }
#   #
#   #         if(any(x$Adduct == "(M+H)+") | any(x$Adduct == "(M-H)-")){
#   #           x <-
#   #             x %>%
#   #             dplyr::filter(Adduct == "(M+H)+" | Adduct == "(M-H)-")
#   #         }
#   #       }
#   #     ) %>%
#   #     dplyr::bind_rows() %>%
#   #     dplyr::select(-fc_class) %>%
#   #     dplyr::arrange(dplyr::desc(fc)) %>%
#   #     dplyr::filter(isotope == "[M]")
#   #
#   #
#   # ##for each Lab.ID, only remain the compound_class with the max score
#   #     feature_table1 <-
#   #       feature_table1 %>%
#   #       plyr::dlply(.variables = .(Lab.ID)) %>%
#   #       purrr::map(function(x) {
#   #         x %>%
#   #           dplyr::filter(score == max(score))
#   #       }) %>%
#   #       dplyr::bind_rows() %>%
#   #       dplyr::arrange(desc(fc))
#
#
#   ##for the metabolites are not in the feature set
#   feature_table2 <-
#     feature_table[-idx,] %>%
#     dplyr::filter(isotope == "[M]") %>%
#     dplyr::filter(!name %in% feature_table1$name) %>%
#     dplyr::distinct(name, .keep_all = TRUE)
#
#   feature_table <-
#     rbind(feature_table1,
#           feature_table2) %>%
#     dplyr::arrange(dplyr::desc(fc))
#
#   if (verbose)
#     message("calculating observed enrichment scores...")
#   # browser()
#
#   observed_info <- lapply(feature_set, function(x)
#     get_gsea_score2(
#       feature_set = x,
#       feature_table = feature_table,
#       exponent = exponent
#     ))
#
#
#   observed_score <- sapply(observed_info, function(x)
#     x$ES)
#
#   if (verbose) {
#     message("calculating permutation scores...")
#   }
#
#   if (seed) {
#     seeds <- sample.int(perm_num)
#   }
#
#   perm_scores <- BiocParallel::bplapply(1:perm_num, function(i) {
#     if (seed) {
#       set.seed(seeds[i])
#     }
#
#     perm_gsea_es2(
#       feature_table = feature_table,
#       feature_set = feature_set,
#       exponent = exponent
#     )
#   })
#
#
#   perm_scores <- do.call(what = "cbind",
#                          args = perm_scores)
#
#   rownames(perm_scores) <- names(feature_set)
#
#   pos.m <- apply(perm_scores, 1, function(x)
#     mean(x[x >= 0]))
#
#   neg.m <- apply(perm_scores, 1, function(x)
#     abs(mean(x[x < 0])))
#
#   normalize_es <- function(ES, pos.m, neg.m) {
#     s <- sign(ES)
#     m <- numeric(length(ES))
#     m[s == 1] <- pos.m[s == 1]
#     m[s == -1] <- neg.m[s == -1]
#     ES / m
#   }
#
#   NES <- normalize_es(observed_score, pos.m, neg.m)
#
#   perm_scores <-
#     apply(perm_scores,
#           2,
#           normalize_es,
#           pos.m = pos.m,
#           neg.m = neg.m)
#
#   if (class(perm_scores) == "numeric") {
#     perm_scores <- matrix(perm_scores, nrow = 1)
#     rownames(perm_scores) <- names(feature_set)
#   }
#
#   if (verbose)
#     message("calculating p values...")
#
#   p_values <- sapply(seq_along(observed_score), function(i) {
#     if (is.na(NES[i])) {
#       NA
#     } else if (NES[i] >= 0) {
#       (sum(perm_scores[i, ] >= NES[i]) + 1) / (sum(perm_scores[i,] >= 0) + 1)
#     } else {
#       # NES[i] < 0
#       (sum(perm_scores[i, ] <= NES[i]) + 1) / (sum(perm_scores[i,] < 0) +
#                                                  1)
#     }
#
#   })
#
#   p_adjust <- p.adjust(p_values, method = p_adjust_method)
#
#   qvalues <- calculate_qvalue(pvals = p_values)
#
#   gs.name <- names(feature_set)
#
#   description <- stringr::str_split(gs.name, pattern = ";") %>%
#     lapply(function(x) {
#       x[2]
#     }) %>%
#     unlist()
#
#   id <- stringr::str_split(gs.name, pattern = ";") %>%
#     lapply(function(x) {
#       x[1]
#     }) %>%
#     unlist()
#
#   params <- list(
#     pvalue_cutoff = pvalue_cutoff,
#     perm_num = perm_num,
#     p_adjust_method = p_adjust_method,
#     exponent = exponent,
#     min_size = min_size,
#     max_size = max_size
#   )
#
#
#   res <- data.frame(
#     ID = id,
#     Description = description,
#     setSize = sapply(feature_set, length),
#     enrichmentScore = observed_score,
#     NES = NES,
#     pvalue = p_values,
#     p_adjust = p_adjust,
#     qvalues = qvalues,
#     stringsAsFactors = FALSE
#   )
#
#
#   res <-
#     res %>%
#     dplyr::filter(!is.na(pvalue)) %>%
#     dplyr::filter(pvalue <= pvalue_cutoff) %>%
#     dplyr::filter(p_adjust <= pvalue_cutoff) %>%
#     dplyr::arrange(pvalue)
#
#   if (nrow(res) == 0) {
#     message("no term enriched under specific pvalue_cutoff...")
#     return(
#       new(
#         "metPathGSEA",
#         result     = res,
#         feature_set   = feature_set,
#         feature_list   = feature_list,
#         params     = params
#       )
#     )
#   }
#
#   row.names(res) <- res$ID
#   observed_info <-
#     observed_info[sapply(res$ID, function(x)
#       grep(x, names(observed_info)))]
#
#   if (verbose)
#     message("leading edge analysis...")
#
#   ledge <- get_leading_edge(observed_info)
#
#   res$rank <- ledge$rank
#   res$leading_edge <- ledge$leading_edge
#   res$core_enrichment <-
#     sapply(ledge$core_enrichment, paste0, collapse = '/')
#
#
#   if (verbose)
#     message("done...")
#
#   feature_list <- feature_table$fc
#   names(feature_list) <- feature_table$Lab.ID
#   new(
#     "metPathGSEA",
#     result     = res,
#     feature_set   = feature_set,
#     feature_list   = feature_list,
#     perm_scores = perm_scores,
#     params     = params
#   )
# }
#
#
#
# ##' @name metPathGSEA-class
# ##' @aliases metPathGSEA-class
# ##'   show,metPathGSEA-method summary,metPathGSEA-method
# ##'
# ##' @docType class
# ##' @slot result GSEA anaysis result.
# ##' @slot geneSets geneSets
# ##' @slot feature_list order rank feature_list
# ##' @slot permScores permutation scores
# ##' @slot params parameters
# ##' @exportClass metPathGSEA
# ##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
# ##' @keywords classes
# setClass(
#   "metPathGSEA",
#   representation   = representation(
#     result          = "data.frame",
#     feature_set        = "list",
#     feature_list        = "numeric",
#     perm_scores      = "matrix",
#     params          = "list"
#   )
# )
#
#
#
#
#
#
#
#
#
#
#
#
#
# ###----------------------------------------------------------------------------
# # feature_set <- kegg_hsa_compound_pathway
# # feature_list <- feature_list
#
# filter_feature_set2 <- function(feature_set,
#                                 feature_table,
#                                 min_size = 5,
#                                 max_size = 1000) {
#   len <- lapply(feature_set, function(x) {
#     length(intersect(x, feature_table$Lab.ID))
#   }) %>%
#     unlist() %>%
#     unname()
#
#   remain_idx <-
#     which(len > min_size & len < max_size)
#
#   feature_set <- feature_set[remain_idx]
#   return(feature_set)
# }
#
#
#
#
#
#
# get_gsea_score2 <- function(feature_table,
#                             feature_set,
#                             exponent = 1,
#                             fortify = FALSE) {
#   ###################################################################
#   ##    feature_list                                               ##
#   ##                                                               ##
#   ## 1. Rank order the N genes in D to form L = { g_1, ... , g_N}  ##
#   ##    according to the correlation, r(g_j)=r_j,                  ##
#   ##    of their expression profiles with C.                       ##
#   ##                                                               ##
#   ###################################################################
#
#   ###################################################################
#   ##    exponent                                                   ##
#   ##                                                               ##
#   ## An exponent p to control the weight of the step.              ##
#   ##   When p = 0, Enrichment Score ( ES(S) ) reduces to           ##
#   ##   the standard Kolmogorov-Smirnov statistic.                  ##
#   ##   When p = 1, we are weighting the genes in S                 ##
#   ##   by their correlation with C normalized                      ##
#   ##   by the sum of the correlations over all of the genes in S.  ##
#   ##                                                               ##
#   ###################################################################
#
#   ## genes defined in feature_set should appear in feature_table.
#   feature_set <- intersect(feature_set, feature_table$Lab.ID)
#
#   ##N the length of feature_list
#   num_feature_list <- nrow(feature_table)
#   # num_feature_list <-
#   # length(unique(feature_table$Lab.ID))
#
#   ##num_feature_set the length of feature set
#   num_feature_set <- length(feature_set)
#
#   Phit <- Pmiss <- numeric(num_feature_list)
#   hits <- feature_table$Lab.ID %in% feature_set ## logical
#
#   Phit[hits] <-
#     (abs(feature_table$fc[hits]) * (feature_table$score[hits] / 100)) ^ exponent
#   #   score <- feature_table$score[hits]
#   #   compound_class <- feature_table$compound_class[hits]
#   #
#   #   temp_data <- data.frame(
#   #     index = 1:sum(hits),
#   #     phit = Phit[hits],
#   #     feature_table[hits, ],
#   #     stringsAsFactors = FALSE
#   #   )
#   #
#   # idx <- 17
#   #   temp_data %>%
#   #     ggplot(aes(index, phit)) +
#   #     geom_point() +
#   #     geom_point(
#   #       aes(x = index, y = phit, color = compound_class),
#   #       show.legend = FALSE,
#   #       data = temp_data %>% filter(Lab.ID == unique(Lab.ID)[idx])
#   #     ) +
#   #     ggrepel::geom_label_repel(
#   #       aes(
#   #         x = index,
#   #         y = phit,
#   #         label = paste(compound_class, score, Adduct, isotope, sep = ":"),
#   #         color = compound_class
#   #       ),
#   #       show.legend = FALSE,
#   #       data = temp_data %>% filter(Lab.ID == unique(Lab.ID)[idx])
#   #     ) +
#   #     ggsci::scale_color_aaas() +
#   #     theme_bw()
#
#   # write.csv(feature_table[hits,], "test.csv", row.names = FALSE)
#
#   NR <- sum(Phit)
#   Phit <- cumsum(Phit / NR)
#
#   Pmiss[!hits] <- 1 / (num_feature_list - num_feature_set)
#   Pmiss <- cumsum(Pmiss)
#
#   runningES <- Phit - Pmiss
#
#   # data.frame(index = 1:length(Phit),
#   #            Phit, Pmiss, runningES,
#   #            stringsAsFactors = FALSE) %>%
#   #   tidyr::pivot_longer(cols = -index, names_to = "class", values_to = "value") %>%
#   #   ggplot(aes(x = index, y = value, color = class)) +
#   #   geom_point()
#
#   ## ES is the maximum deviation from zero of Phit-Pmiss
#   max.ES <- max(runningES)
#   min.ES <- min(runningES)
#
#   if (abs(max.ES) > abs(min.ES)) {
#     ES <- max.ES
#   } else {
#     ES <- min.ES
#   }
#
#   df <- data.frame(
#     x = seq_along(runningES),
#     runningScore = runningES,
#     position = as.integer(hits)
#   )
#
#   # df %>%
#   #   ggplot(aes(x, runningScore)) +
#   #   geom_point(aes(color = as.character(position)))
#
#   if (fortify == TRUE) {
#     return(df)
#   }
#
#   df$gene = feature_table$Lab.ID
#   res <- list(ES = ES, runningES = df)
#   return(res)
# }
#
#
#
# perm_feature_table <- function(feature_table) {
#   perm.idx <- sample.int(nrow(feature_table))
#   perm_feature_table <- feature_table
#   perm_feature_table$Lab.ID <- perm_feature_table$Lab.ID[perm.idx]
#   return(perm_feature_table)
# }
#
#
# perm_gsea_es2 <- function(feature_table,
#                           feature_set,
#                           exponent = 1) {
#   feature_table <- perm_feature_table(feature_table)
#   res <- sapply(1:length(feature_set), function(i)
#     get_gsea_score2(
#       feature_table = feature_table,
#       feature_set = feature_set[[i]],
#       exponent = exponent
#     )$ES)
#   return(res)
# }
#
#
#
#
#
#
#
#
#
#
# calculate_qvalue <- function(pvals) {
#   if (length(pvals) == 0)
#     return(numeric(0))
#
#   qobj <- tryCatch(
#     qvalue::qvalue(pvals, lambda = 0.05,
#                    pi0.method = "bootstrap"),
#     error = function(e)
#       NULL
#   )
#
#   if (class(qobj) == "qvalue") {
#     qvalues <- qobj$qvalues
#   } else {
#     qvalues <- NA
#   }
#   return(qvalues)
# }
#
#
# ###----------------------------------------------------------------------------
# get_leading_edge <- function(observed_info) {
#   core_enrichment <- lapply(observed_info, function(x) {
#     runningES <- x$runningES
#     runningES <- runningES[runningES$position == 1, ]
#     ES <- x$ES
#     if (ES >= 0) {
#       i <- which.max(runningES$runningScore)
#       leading_gene <- runningES$gene[1:i]
#     } else {
#       i <- which.min(runningES$runningScore)
#       leading_gene <- runningES$gene[-c(1:(i - 1))]
#     }
#     return(leading_gene)
#   })
#
#   rank <- sapply(observed_info, function(x) {
#     runningES <- x$runningES
#     ES <- x$ES
#     if (ES >= 0) {
#       rr <- which.max(runningES$runningScore)
#     } else {
#       i <- which.min(runningES$runningScore)
#       rr <- nrow(runningES) - i + 1
#     }
#     return(rr)
#   })
#
#   tags <- sapply(observed_info, function(x) {
#     runningES <- x$runningES
#     runningES <- runningES[runningES$position == 1, ]
#     ES <- x$ES
#     if (ES >= 0) {
#       i <- which.max(runningES$runningScore)
#       res <- i / nrow(runningES)
#     } else {
#       i <- which.min(runningES$runningScore)
#       res <- (nrow(runningES) - i + 1) / nrow(runningES)
#     }
#     return(res)
#   })
#
#   ll <- sapply(observed_info, function(x) {
#     runningES <- x$runningES
#     ES <- x$ES
#     if (ES >= 0) {
#       i <- which.max(runningES$runningScore)
#       res <- i / nrow(runningES)
#     } else {
#       i <- which.min(runningES$runningScore)
#       res <- (nrow(runningES) - i + 1) / nrow(runningES)
#     }
#     return(res)
#   })
#
#   N <- nrow(observed_info[[1]]$runningES)
#   setSize <-
#     sapply(observed_info, function(x)
#       sum(x$runningES$position))
#   signal <- tags * (1 - ll) * (N / (N - setSize))
#
#   tags <- paste0(round(tags * 100), "%")
#   ll <- paste0(round(ll * 100), "%")
#   signal <- paste0(round(signal * 100), "%")
#   leading_edge <-
#     paste0('tags=', tags, ", list=", ll, ", signal=", signal)
#
#   res <- list(
#     rank = rank,
#     tags = tags,
#     list = ll,
#     signal = signal,
#     leading_edge = leading_edge,
#     core_enrichment = core_enrichment
#   )
#   return(res)
# }
#
#
#
# ##' show method for \code{metPathGSEA} instance
# ##'
# ##' @name show
# ##' @docType methods
# ##' @rdname show-methods
# ##'
# ##' @title show method
# ##' @return message
# ##' @importFrom methods show
# ##' @exportMethod show
# ##' @usage show(object)
# ##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
# setMethod("show", signature(object = "metPathGSEA"),
#           function (object) {
#             params <- object@params
#             cat("#\n# Gene Set Enrichment Analysis\n#\n")
#             cat("#...@feature_list", "\t")
#             str(object@feature_list)
#             cat("#...nPerm", "\t", params$perm_num, "\n")
#             cat(
#               "#...pvalues adjusted by",
#               paste0("'", params$p_adjust_method, "'"),
#               paste0("with cutoff <", params$pvalueCutoff),
#               "\n"
#             )
#             cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
#             str(object@result)
#           })
#
#
#
# ##' visualize analyzing result of GSEA
# ##'
# ##' plotting function for gseaResult
# ##' @title gsea_plot
# ##' @rdname gseaplot
# ##' @param x object of gsea result
# ##' @param geneSetID geneSet ID
# ##' @param by one of "runningScore" or "position"
# ##' @param title plot title
# ##' @param ... additional parameters
# ##' @return ggplot2 object
# ##' @export
# ##' @examples
# ##' library(DOSE)
# ##' data(geneList)
# ##' x <- gseDO(geneList)
# ##' gseaplot(x, geneSetID=1)
# setGeneric(name = "gsea_plot",
#            function(x,
#                     feature_set_idx = 1,
#                     by = "all",
#                     title = "",
#                     color = 'black',
#                     color.line = "#8DD3C7",
#                     color.vline = "#FB8072",
#                     ...) {
#              standardGeneric("gsea_plot")
#            })
#
#
# ##' @rdname gsea_plot
# ##' @exportMethod gsea_plot
# setMethod(f = "gsea_plot", signature(x = "metPathGSEA"),
#           function (x,
#                     feature_set_idx = 1,
#                     by = "all",
#                     title = "",
#                     color = 'black',
#                     color.line = "#8DD3C7",
#                     color.vline = "#FB8072",
#                     ...) {
#             gsea_plot.metPathGSEA(
#               x,
#               feature_set_idx = feature_set_idx,
#               by = by,
#               title = title,
#               color = color,
#               color.line = color.line,
#               color.vline = color.vline,
#               ...
#             )
#           })
#
# ##' @rdname gsea_plot
# ##' @param color color of line segments
# ##' @param color.line color of running enrichment score line
# ##' @param color.vline color of vertical line which indicating the maximum/minimal running enrichment score
# ##' @return ggplot2 object
# ##' @importFrom ggplot2 ggplot
# ##' @importFrom ggplot2 geom_linerange
# ##' @importFrom ggplot2 geom_line
# ##' @importFrom ggplot2 geom_vline
# ##' @importFrom ggplot2 geom_hline
# ##' @importFrom ggplot2 xlab
# ##' @importFrom ggplot2 ylab
# ##' @importFrom ggplot2 xlim
# ##' @importFrom ggplot2 aes
# ##' @importFrom ggplot2 ggplotGrob
# ##' @importFrom ggplot2 geom_segment
# ##' @importFrom ggplot2 ggplot_gtable
# ##' @importFrom ggplot2 ggplot_build
# ##' @importFrom ggplot2 ggtitle
# ##' @importFrom ggplot2 element_text
# ##' @importFrom ggplot2 rel
# ##' @importFrom cowplot plot_grid
# ##' @author Guangchuang Yu
#
# gsea_plot.metPathGSEA <-
#   function (x,
#             feature_set_idx = 1,
#             by = "all",
#             title = "",
#             color = 'black',
#             color.line = "green",
#             color.vline = "#FA5860",
#             ...) {
#     if (is.null(x)) {
#       return(NULL)
#     }
#     by <- match.arg(by, c("runningScore", "preranked", "all"))
#     gs_data <- get_gs_info(x, feature_set_idx)
#
#     p <- ggplot(gs_data, aes_(x = ~ x)) +
#       theme_dose() +
#       xlab("Position in the ranked list of features")
#
#     if (by == "runningScore" || by == "all") {
#       p.res <-
#         p + geom_linerange(aes_(ymin =  ~ ymin, ymax =  ~ ymax), color = color)
#       p.res <-
#         p.res + geom_line(aes_(y = ~ runningScore), color = color.line, size =
#                             1)
#       enrichmentScore <-
#         x@result[feature_set_idx, "enrichmentScore"]
#
#       es.df <-
#         data.frame(es = which.min(abs(
#           p$data$runningScore - enrichmentScore
#         )))
#
#       p.res <-
#         p.res + geom_vline(
#           data = es.df,
#           aes_(xintercept = ~ es),
#           colour = color.vline,
#           linetype = "dashed"
#         )
#       p.res <- p.res + ylab("Running enrichment score")
#       p.res <- p.res + geom_hline(yintercept = 0)
#     }
#
#     if (by == "preranked" || by == "all") {
#       df2 <- data.frame(x = which(p$data$position == 1))
#       df2$y <- p$data$feature_list[df2$x]
#       p.pos <-
#         p + geom_segment(data = df2,
#                          aes_(
#                            x =  ~ x,
#                            xend =  ~ x,
#                            y =  ~ y,
#                            yend = 0
#                          ),
#                          color = color)
#       p.pos <-
#         p.pos + ylab("Ranked list metric") + xlim(0, length(p$data$feature_list))
#     }
#     if (by == "runningScore")
#       return(p.res + ggtitle(title))
#     if (by == "preranked")
#       return(p.pos + ggtitle(title))
#
#     p.pos <-
#       p.pos + xlab(NULL) + theme(axis.text.x = element_blank(),
#                                  axis.ticks.x = element_blank())
#     p.pos <- p.pos + ggtitle(title) +
#       theme(plot.title = element_text(hjust = 0.5, size = rel(2)))
#     cowplot::plot_grid(p.pos, p.res, ncol = 1, align = "v")
#   }
#
#
#
# get_gs_info <- function(object,
#                         feature_set_idx = 1) {
#   feature_list <- object@feature_list
#   # if (is.numeric(feature_set_idx))
#   #   feature_set_idx <- object@result[feature_set_idx, "ID"]
#   feature_set <- object@feature_set[[feature_set_idx]]
#   exponent <- object@params[["exponent"]]
#   df <-
#     get_gsea_score(feature_list, feature_set, exponent, fortify = TRUE)
#   df$ymin = 0
#   df$ymax = 0
#   pos <- df$position == 1
#   h <- diff(range(df$runningScore)) / 20
#   df$ymin[pos] <- -h
#   df$ymax[pos] <- h
#   df$feature_list <- feature_list
#   return(df)
# }
#
#
# ##' ggplot theme of DOSE
# ##'
# ##' @title theme_dose
# ##' @param font.size font size
# ##' @return ggplot theme
# ##' @importFrom ggplot2 theme_bw
# ##' @importFrom ggplot2 theme
# ##' @importFrom ggplot2 element_text
# ##' @importFrom ggplot2 margin
# ##' @examples
# ##' library(ggplot2)
# ##' qplot(1:10) + theme_dose()
# ##' @export
# theme_dose <- function(font.size = 13) {
#   theme_bw() +
#     theme(
#       axis.text.x = element_text(
#         colour = "black",
#         size = 12,
#         vjust = 1
#       ),
#       axis.text.y = element_text(
#         colour = "black",
#         size = 12,
#         hjust = 1
#       ),
#       axis.title = element_text(
#         margin = margin(10, 5, 0, 0),
#         color = "black",
#         size = 13
#       ),
#       axis.title.y = element_text(angle = 90)
#     )
# }
