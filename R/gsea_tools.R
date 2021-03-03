# # feature_list <- feature_list
# # feature_set <- kegg_hsa_compound_pathway
# # exponent = 1
# # perm_num = 1000
# # min_size = 5
# # max_size = 1000
# # pvalue_cutoff = 0.2
# # p_adjust_method = "fdr"
# # seed = FALSE
# # verbose = TRUE
# # library(tidyverse)
# # result <-do_gsea(feature_list = feature_list,
# #                  feature_set = feature_set,
# #                  exponent = 1, perm_num = 1000,
# #                  min_size = 5, max_size = 1000,
# #                  pvalue_cutoff = 0.2,
# #                  p_adjust_method = "fdr",
# #                  seed = TRUE,
# #                  verbose = TRUE)
#
# do_gsea <- function(feature_list,
#                     feature_set,
#                     exponent = 1,
#                     perm_num = 1000,
#                     min_size = 5,
#                     max_size = 1000,
#                     pvalue_cutoff = 0.2,
#                     p_adjust_method = "fdr",
#                     seed = FALSE,
#                     verbose = TRUE,
#                     ...) {
#   if (verbose)
#     message("preparing feature_set collections...")
#
#   feature_set <-
#     filter_feature_set(feature_set, feature_list, min_size, max_size)
#
#   if(length(feature_set) == 0){
#     return(NULL)
#   }
#
#   if (verbose)
#     message("calculating observed enrichment scores...")
# # browser()
#   observed_info <- lapply(feature_set, function(x)
#     get_gsea_score(
#       feature_set = x,
#       feature_list = feature_list,
#       exponent = exponent
#     )
#     )
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
#     if (seed){
#       set.seed(seeds[i])
#     }
#
#     perm_gsea_es(feature_list = feature_list,
#                  feature_set = feature_set,
#                  exponent = exponent)
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
#   if(class(perm_scores) == "numeric"){
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
#       (sum(perm_scores[i,] >= NES[i]) + 1) / (sum(perm_scores[i, ] >= 0) + 1)
#     } else {
#       # NES[i] < 0
#       (sum(perm_scores[i,] <= NES[i]) + 1) / (sum(perm_scores[i, ] < 0) +
#                                                1)
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
#     lapply(function(x){
#       x[2]
#     }) %>%
#     unlist()
#
#   id <- stringr::str_split(gs.name, pattern = ";") %>%
#     lapply(function(x){
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
#   observed_info <- observed_info[sapply(res$ID, function(x) grep(x, names(observed_info)))]
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
# setClass("metPathGSEA",
#          representation   = representation(
#            result          = "data.frame",
#            feature_set        = "list",
#            feature_list        = "numeric",
#            perm_scores      = "matrix",
#            params          = "list")
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
# filter_feature_set <- function(feature_set,
#                                feature_list,
#                                min_size = 5,
#                                max_size = 1000) {
#   len <- lapply(feature_set, function(x) {
#     length(intersect(x, names(feature_list)))
#   }) %>%
#     unlist() %>%
#     unname()
#
#   # len <- lapply(feature_set, length) %>%
#   #   unlist()
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
# get_gsea_score <- function(feature_list,
#                            feature_set,
#                            exponent=1,
#                            fortify=FALSE) {
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
#   ## genes defined in feature_set should appear in feature_list.
#   ## this is a must, see https://github.com/GuangchuangYu/DOSE/issues/23
#   feature_set <- intersect(feature_set, names(feature_list))
#
#   ##N the length of feature_list
#   N <- length(feature_list)
#   ##Nh the length of feature set
#   Nh <- length(feature_set)
#
#   Phit <- Pmiss <- numeric(N)
#   hits <- names(feature_list) %in% feature_set ## logical
#
#   Phit[hits] <- abs(feature_list[hits])^exponent
#   NR <- sum(Phit)
#   Phit <- cumsum(Phit/NR)
#
#   Pmiss[!hits] <- 1/(N-Nh)
#   Pmiss <- cumsum(Pmiss)
#
#   runningES <- Phit - Pmiss
#
# # data.frame(index = 1:length(Phit),
# #            Phit, Pmiss, runningES,
# #            stringsAsFactors = FALSE) %>%
# #   tidyr::pivot_longer(cols = -index, names_to = "class", values_to = "value") %>%
# #   ggplot(aes(x = index, y = value, color = class)) +
# #   geom_point()
#
#   ## ES is the maximum deviation from zero of Phit-Pmiss
#   max.ES <- max(runningES)
#   min.ES <- min(runningES)
#
#   if( abs(max.ES) > abs(min.ES) ) {
#     ES <- max.ES
#   } else {
#     ES <- min.ES
#   }
#
#   df <- data.frame(x = seq_along(runningES),
#                    runningScore = runningES,
#                    position = as.integer(hits)
#   )
#
#   # df %>%
#   #   ggplot(aes(x, runningScore)) +
#   #   geom_point(aes(color = as.character(position)))
#
#   if(fortify==TRUE) {
#     return(df)
#   }
#
#   df$gene = names(feature_list)
#   res <- list(ES=ES, runningES = df)
#   return(res)
# }
#
#
#
# perm_feature_list <- function(feature_list) {
#   perm.idx <- sample.int(length(feature_list))
#   perm_feature_list <- feature_list
#   names(perm_feature_list) <- names(feature_list)[perm.idx]
#   return(perm_feature_list)
# }
#
#
# perm_gsea_es <- function(feature_list,
#                          feature_set,
#                          exponent=1) {
#   feature_list <- perm_feature_list(feature_list)
#   res <- sapply(1:length(feature_set), function(i)
#     get_gsea_score(feature_list = feature_list,
#                    feature_set = feature_set[[i]],
#                    exponent = exponent)$ES
#   )
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
#   qobj <- tryCatch(qvalue::qvalue(pvals, lambda=0.05,
#                           pi0.method="bootstrap"),
#                    error=function(e) NULL)
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
#     runningES <- runningES[runningES$position == 1,]
#     ES <- x$ES
#     if (ES >= 0) {
#       i <- which.max(runningES$runningScore)
#       leading_gene <- runningES$gene[1:i]
#     } else {
#       i <- which.min(runningES$runningScore)
#       leading_gene <- runningES$gene[-c(1:(i-1))]
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
#     runningES <- runningES[runningES$position == 1,]
#     ES <- x$ES
#     if (ES >= 0) {
#       i <- which.max(runningES$runningScore)
#       res <- i/nrow(runningES)
#     } else {
#       i <- which.min(runningES$runningScore)
#       res <- (nrow(runningES) - i + 1)/nrow(runningES)
#     }
#     return(res)
#   })
#
#   ll <- sapply(observed_info, function(x) {
#     runningES <- x$runningES
#     ES <- x$ES
#     if (ES >= 0) {
#       i <- which.max(runningES$runningScore)
#       res <- i/nrow(runningES)
#     } else {
#       i <- which.min(runningES$runningScore)
#       res <- (nrow(runningES) - i + 1)/nrow(runningES)
#     }
#     return(res)
#   })
#
#   N <- nrow(observed_info[[1]]$runningES)
#   setSize <- sapply(observed_info, function(x) sum(x$runningES$position))
#   signal <- tags * (1-ll) * (N / (N - setSize))
#
#   tags <- paste0(round(tags * 100), "%")
#   ll <- paste0(round(ll * 100), "%")
#   signal <- paste0(round(signal * 100), "%")
#   leading_edge <- paste0('tags=', tags, ", list=", ll, ", signal=", signal)
#
#   res <- list(rank = rank,
#               tags = tags,
#               list = ll,
#               signal = signal,
#               leading_edge = leading_edge,
#               core_enrichment = core_enrichment)
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
# setMethod("show", signature(object="metPathGSEA"),
#           function (object){
#             params <- object@params
#             cat("#\n# Gene Set Enrichment Analysis\n#\n")
#             cat("#...@feature_list", "\t")
#             str(object@feature_list)
#             cat("#...nPerm", "\t", params$perm_num, "\n")
#             cat("#...pvalues adjusted by", paste0("'", params$p_adjust_method, "'"),
#                 paste0("with cutoff <", params$pvalueCutoff), "\n")
#             cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
#             str(object@result)
#           }
# )
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
#     if(is.null(x)){
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
#       enrichmentScore <- x@result[feature_set_idx, "enrichmentScore"]
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
#     p.pos <- p.pos + xlab(NULL) + theme(axis.text.x = element_blank(),
#                                         axis.ticks.x = element_blank())
#     p.pos <- p.pos + ggtitle(title) +
#       theme(plot.title = element_text(hjust = 0.5, size = rel(2)))
#     cowplot::plot_grid(p.pos, p.res, ncol = 1, align = "v")
#   }
#
#
#
# get_gs_info <- function(object,
#                    feature_set_idx = 1) {
#   feature_list <- object@feature_list
#   # if (is.numeric(feature_set_idx))
#   #   feature_set_idx <- object@result[feature_set_idx, "ID"]
#   feature_set <- object@feature_set[[feature_set_idx]]
#   exponent <- object@params[["exponent"]]
#   df <- get_gsea_score(feature_list, feature_set, exponent, fortify = TRUE)
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
# theme_dose <- function(font.size=13) {
#   theme_bw() +
#     theme(axis.text.x = element_text(colour = "black",
#                                      size = 12, vjust =1 ),
#           axis.text.y = element_text(colour = "black",
#                                      size = 12, hjust =1 ),
#           axis.title = element_text(margin=margin(10, 5, 0, 0),
#                                     color = "black",
#                                     size = 13),
#           axis.title.y = element_text(angle=90)
#     )
# }
