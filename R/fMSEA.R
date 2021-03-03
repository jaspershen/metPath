# #to avoind source
# no_exist_function()
#
# sxtTools::setwd_project()
# setwd("demo_data/")
# library(tidyverse)
# # load("fc_p_value")
# # load("variable_info")
# # load("kegg_hsa_compound_pathway")
# #
# #
# # ms1_table <-
# #   data.frame(variable_info[,c(1:3),], fc_p_value[[6]], stringsAsFactors = FALSE)
# #
# # ms1_table <-
# #   ms1_table %>%
# #   dplyr::mutate(polarity = case_when(
# #     stringr::str_detect(name, "POS") ~ "positive",
# #     TRUE ~ "negative"
# #   )) %>%
# #   dplyr::select(name, mz, rt, polarity, everything())
# #
# #
# # write.csv(ms1_table, "ms1_table.csv", row.names = FALSE)
#
# ###annotation
# library(metID)
# ms1_data <- readr::read_csv("ms1_table.csv")
#
# ms1_data_pos <-
#   ms1_data %>%
#   dplyr::filter(polarity == "positive") %>%
#   dplyr::select(name, mz, rt)
#
# ms1_data_neg <-
#   ms1_data %>%
#   dplyr::filter(polarity == "negative") %>%
#   dplyr::select(name, mz, rt)
#
#
# # if(nrow(ms1_data_pos) == 0){
# #   annotation_table_pos <- NULL
# # }else{
# #   write.csv(ms1_data_pos,
# #             file = "ms1_data_pos.csv",
# #             row.names = FALSE)
# #   annotation_table_pos <-
# #     metID::identify_metabolites(ms1.data = "ms1_data_pos.csv",
# #                                 ms1.match.ppm = 15,
# #                                 column = "rp",
# #                                 polarity = "positive",
# #                                 database = "kegg_ms1_database",
# #                                 candidate.num = 10000)
# #
# #   annotation_table_pos <-
# #     metID::get_identification_table(annotation_table_pos,
# #                                     candidate.num = 10000,
# #                                     type = "new")
# #
# #
# #   annotation_table_pos <-
# #     annotation_table_pos %>%
# #     # dplyr::select(name, KEGG.ID, Adduct, mz.error) %>%
# #     dplyr::left_join(ms1_data %>% dplyr::select(-c(mz, rt)), by = c("name"))
# # }
# #
# #
# # save(annotation_table_pos, file = "annotation_table_pos")
#
# load("annotation_table_pos")
#
# # if(nrow(ms1_data_neg) == 0){
# #   annotation_table_neg <- NULL
# # }else{
# #   write.csv(ms1_data_neg,
# #             file = "ms1_data_neg.csv", row.names = FALSE)
# #   annotation_table_neg <-
# #     metID::identify_metabolites(ms1.data = "ms1_data_neg.csv",
# #                                 ms1.match.ppm = 15,
# #                                 column = "rp",
# #                                 polarity = "negative",
# #                                 database = "kegg_ms1_database",
# #                                 candidate.num = 10000)
# #   annotation_table_neg <-
# #     metID::get_identification_table(annotation_table_neg,
# #                                     candidate.num = 10000,
# #                                     type = "new")
# #
# #   annotation_table_neg <-
# #     annotation_table_neg %>%
# #     # dplyr::select(name, KEGG.ID, Adduct, mz.error) %>%
# #     dplyr::left_join(ms1_data %>% dplyr::select(-c(mz, rt)), by = c("name"))
# # }
# #
# # save(annotation_table_neg, file = "annotation_table_neg")
#
# load("annotation_table_neg")
#
#
# annotation_table <-
#   rbind(annotation_table_pos, annotation_table_neg) %>%
#   dplyr::select(name, polarity, mz, rt, Compound.name, Lab.ID, Adduct,
#                 mz.error, mz.match.score, RT.error, RT.match.score, p_value, fc) %>%
#   dplyr::mutate(fc = -log(fc, 2)) %>%
#   dplyr::arrange(desc(fc))
#
# load("kegg_ms1_database")
#
# annotation_table <-
#   annotation_table %>%
#   dplyr::left_join(kegg_ms1_databasee@spectra.info[,c("Lab.ID", "Formula")], by = "Lab.ID")
#
# annotation_table <-
#   annotation_table %>%
#   dplyr::filter(!is.na(Lab.ID))
#
#
# annotation_table <-
#   annotation_table %>%
#   dplyr::group_by(name) %>%
#   dplyr::mutate(mz = mz[1], rt = rt[1]) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(mz = as.numeric(mz), rt = as.numeric(rt))
#
#
#
# ####isotope annotation
# ###positive
# annotation_table_pos <-
#   annotation_table %>%
#   dplyr::filter(polarity == "positive")
#
#
# # system.time(
# #   isotope_pos <-
# #     purrr::map(
# #       as.data.frame(t(annotation_table_pos)),
# #       .f = function(x) {
# #         adduct <-
# #           stringr::str_extract(x[7], "\\(.+\\)") %>%
# #           stringr::str_replace("\\(", "") %>%
# #           stringr::str_replace("\\)", "")
# #
# #         temp_iso <- try(
# #           annotate_isotope(
# #             formula = x[14],
# #             adduct = adduct,
# #             mz = as.numeric(x[3]),
# #             rt = as.numeric(x[4]),
# #             peak.mz = ms1_data_pos$mz,
# #             peak.rt = ms1_data_pos$rt,
# #             rt.tol = 5,
# #             mz.tol = 15,
# #             max.isotope = 3
# #           ), silent = TRUE
# #         )
# #         if(class(temp_iso) == "try-error"){
# #           return(NULL)
# #         }
# #         temp_iso
# #       }
# #     )
# # )
# #
# #
# # ###
# # isotope_pos <-
# #   purrr::map2(
# #     isotope_pos,
# #     data.frame(t(annotation_table_pos)),
# #     .f = function(x, y) {
# #      if(is.null(x)) return(NULL)
# #       x <- cbind(ms1_data_pos[x$peakIndex, ], x) %>%
# #         dplyr::select(-c(peakIndex))
# #       colnames(x) <- c("name", "mz", "rt", "mz.error", "isotope", "RT.error")
# #
# #       y <- matrix(y, nrow = 1) %>% as.data.frame()
# #       colnames(y) <- colnames(annotation_table_pos)
# #
# #       x$polarity <- y$polarity
# #       x$Compound.name <- y$Compound.name
# #       x$Lab.ID <- y$Lab.ID
# #       x$Adduct <- y$Adduct
# #       x$p_value <- y$p_value
# #       x$fc <- y$fc
# #       x$Formula <- y$Formula
# #       x$rt <- y$rt
# #
# #       x$mz.match.score <- (15 - x$mz.error)/15
# #       x$RT.match.score <- (15 - x$RT.error)/15
# #       x
# #     }
# #   )
# #
# # isotope_pos <- isotope_pos %>%
# #   do.call(rbind, .)
# #
# #
# #
# # isotope_pos <-
# #   isotope_pos %>%
# #   dplyr::select(colnames(annotation_table_pos))
# #
# #
# # save(isotope_pos, file = "isotope_pos")
#
# load("isotope_pos")
#
# annotation_table_pos$isotope <-
#   "[M]"
#
# annotation_table_pos <-
#   rbind(annotation_table_pos, isotope_pos) %>%
#   dplyr::arrange(name)
#
#
#
# ##negative
# annotation_table_neg <-
#   annotation_table %>%
#   dplyr::filter(polarity == "negative")
#
#
# # system.time(
# #   isotope_neg <-
# #     purrr::map(
# #       as.data.frame(t(annotation_table_neg)),
# #       .f = function(x) {
# #         adduct <-
# #           stringr::str_extract(x[7], "\\(.+\\)") %>%
# #           stringr::str_replace("\\(", "") %>%
# #           stringr::str_replace("\\)", "")
# #
# #         temp_iso <- try(
# #           annotate_isotope(
# #             formula = x[14],
# #             adduct = adduct,
# #             mz = as.numeric(x[3]),
# #             rt = as.numeric(x[4]),
# #             peak.mz = ms1_data_neg$mz,
# #             peak.rt = ms1_data_neg$rt,
# #             rt.tol = 5,
# #             mz.tol = 15,
# #             max.isotope = 3
# #           ), silent = TRUE
# #         )
# #         if(class(temp_iso) == "try-error"){
# #           return(NULL)
# #         }
# #         temp_iso
# #       }
# #     )
# # )
# #
# #
# # ###
# # isotope_neg <-
# #   purrr::map2(
# #     isotope_neg,
# #     data.frame(t(annotation_table_neg)),
# #     .f = function(x, y) {
# #       if(is.null(x)) return(NULL)
# #       x <- cbind(ms1_data_neg[x$peakIndex, ], x) %>%
# #         dplyr::select(-c(peakIndex))
# #       colnames(x) <- c("name", "mz", "rt", "mz.error", "isotope", "RT.error")
# #
# #       y <- matrix(y, nrow = 1) %>% as.data.frame()
# #       colnames(y) <- colnames(annotation_table_neg)
# #
# #       x$polarity <- y$polarity
# #       x$Compound.name <- y$Compound.name
# #       x$Lab.ID <- y$Lab.ID
# #       x$Adduct <- y$Adduct
# #       x$p_value <- y$p_value
# #       x$fc <- y$fc
# #       x$Formula <- y$Formula
# #       x$rt <- y$rt
# #
# #       x$mz.match.score <- (15 - x$mz.error)/15
# #       x$RT.match.score <- (15 - x$RT.error)/15
# #       x
# #     }
# #   )
# #
# # isotope_neg <- isotope_neg %>%
# #   do.call(rbind, .)
# #
# # isotope_neg <-
# #   isotope_neg %>%
# #   dplyr::select(colnames(annotation_table_neg))
# #
# # save(isotope_neg, file = "isotope_neg")
#
# load("isotope_neg")
#
# annotation_table_neg$isotope <-
#   "[M]"
#
# annotation_table_neg <-
#   rbind(annotation_table_neg, isotope_neg) %>%
#   dplyr::arrange(name)
#
#
# ####group peaks groups according to RT and ID
# annotation_table <-
#   rbind(annotation_table_pos, annotation_table_neg)
#
#
# library(plyr)
#
# annotation_table <-
#   annotation_table %>%
#   dplyr::arrange(rt) %>%
#   plyr::dlply(.variables = .(Lab.ID)) %>%
#   purrr::map(.f = function(x){
#     x$rt <- as.numeric(x$rt)
#     x <- x %>%
#       dplyr::arrange(rt)
#     rt_class <- group_peaks_rt(rt = x$rt, rt.tol = 10) %>%
#       dplyr::arrange(rt)
#
#     rt_class <- paste(x$Lab.ID[1], rt_class$class, sep = "_")
#     data.frame(x, compound_class = rt_class, stringsAsFactors = FALSE)
#   })
#
#
# ###score for each compound class
# annotation_table <-
#   annotation_table %>%
#   purrr::map(.f = function(x){
#     x <-
#       x %>%
#       plyr::dlply(.(compound_class)) %>%
#       purrr::map(.f = function(y){
#         score <- score_peak_group(y)
#         data.frame(y, score, stringsAsFactors = FALSE)
#       }) %>%
#       do.call(rbind, .)
#   }
#   )
#
#
# ####remove redundant
# # annotation_table <-
# #   annotation_table %>%
# #   do.call(rbind, .)
# # temp_data <-
# # annotation_table$score %>%
# #   table() %>%
# #   as.data.frame()
# # colnames(temp_data) <- c("Score", "Freq")
# #
# #
# # plot <-
# # temp_data %>%
# #   ggplot(aes(Score, Freq)) +
# #   geom_bar(stat = "identity", aes(fill = Score), show.legend = FALSE) +
# #   scale_fill_manual(values = c("20" = "#3F4041FF",
# #                      "30" = "#3F4041FF",
# #                      "40" = "#84D7E1FF",
# #                      "50" = "#84D7E1FF",
# #                      "60" = "#84D7E1FF",
# #                      "70" = "#84D7E1FF",
# #                      "80" = "#84D7E1FF",
# #                      "90" = "#84D7E1FF",
# #                      "100" = "#008EA0FF",
# #                      "110" = "#008EA0FF",
# #                      "120" = "#008EA0FF",
# #                      "130" = "#008EA0FF",
# #                      "140" = "#008EA0FF",
# #                      "150" = "#008EA0FF",
# #                      "160" = "#008EA0FF",
# #                      "170" = "#008EA0FF",
# #                      "180" = "#008EA0FF",
# #                      "190" = "#008EA0FF",
# #                      "200" = "#C71000FF")) +
# #   theme_classic() +
# #   scale_y_continuous(expand = c(0, 10)) +
# #   theme(axis.title = element_text(size = 13),
# #         axis.text = element_text(size = 12),
# #         plot.background  = element_rect(fill = "transparent",
# #                                       color = NA),
# #         panel.background = element_rect(fill = "transparent",
# #                                       color = NA),
# #         legend.background = element_rect(fill = "transparent",
# #                                       color = NA),
# #         axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1))
# #
# #
# #
# # plot
# #
# #
# # ggsave(plot, filename = "score_distributation.png",
# #        width = 7, height = 7, bg = "transparent")
# #
# # ggsave(plot, filename = "score_distributation.pdf",
# #        width = 7, height = 7, bg = "transparent")
# #
# #
# # idx <- which(annotation_table$score == 200)
# #
# # library(plyr)
# # temp_data <-
# # annotation_table[idx,] %>%
# #   plyr::dlply(.variables = .(compound_class)) %>%
# #   purrr::map(.f = function(x){
# #     dplyr::select(x, name, polarity, mz, rt, Compound.name, Adduct, isotope, score, compound_class)
# #       # dplyr::distinct(polarity, isotope)
# #   })
# #
# # len <-
# # purrr::map(temp_data, function(x){
# #   nrow(x)
# # }) %>%
# #   unlist()
# #
# # table(len)
# #
# # test <-
# # variable_info %>%
# #   dplyr::filter(!is.na(Level)) %>%
# #   dplyr::filter(Level == 1 | Level == 2) %>%
# #   # dplyr::filter(name %in% temp_data[[1]]$name) %>%
# #   dplyr::select(name, Compound.name, KEGG.ID, Adduct, Database)
# #
# # temp_data[[2]] %>%
# #   dplyr::mutate(Compound.name = stringr::str_split(Compound.name, ";") %>% lapply(function(x)x[1]) %>% unlist) %>%
# #   dplyr::select(name, Compound.name, Adduct, isotope, compound_class) %>%
# #   left_join(test, by = "name") %>%
# #   as.data.frame()
#
#
#
# # annotation_table %>%
# #   dplyr::filter(compound_class == "C00191_32")
#
#
#
# ##score distributation
# # score_rule <-
# #   data.frame(
# #     rule =
# #       c(
# #         "Adduct M+H",
# #         "Isotope (Adduct M+H)",
# #         "Other positive adduct",
# #         "Isotope (Other positive adduct)",
# #         "Adduct M-H",
# #         "Isotope (Adduct M-H)",
# #         "Other negative adduct",
# #         "Isotope (Other negative adduct)"
# #       ),
# #     score = c(50, 20, 20, 10, 50, 20, 20, 10),
# #     stringsAsFactors = FALSE
# #   )
# #
# # 3*3*3*3
# #
# # test <- list(c(1, 1), c(1, 0), c(0, 0))
# #
# # # comb <- vector(mode = "list", length = 3^4)
# # comb <- NULL
# # for(i in 1:3){
# #
# #   for(j in 1:3){
# #
# #     for(k in 1:3){
# #
# #       for(z in 1:3){
# #         comb <- c(comb, list(c(test[[i]], test[[j]], test[[k]], test[[z]])))
# #       }
# #     }
# #   }
# # }
# #
# #
# # comb <-
# #   comb %>%
# #   do.call(cbind, .)
# #
# # remove_idx <-
# # apply(comb, 2, function(x){
# #   all(x == 0)
# # }) %>%
# #   which()
# #
# # comb <- comb[,-remove_idx]
# #
# # score <- score_rule$score * comb
# # score <- score %>% colSums()
# #
# # idx <- order(score, decreasing = TRUE)
# #
# # score <- score[idx]
# #
# # comb <- comb[,idx]
# #
# # score_rule <-
# #   data.frame(score_rule, comb, stringsAsFactors = FALSE)
# #
# #
# # plot1 <-
# #   data.frame(index = factor(paste('X', 1:length(score), sep = ""),
# #                             levels = paste('X', 1:length(score), sep = "")),
# #              score, stringsAsFactors = TRUE) %>%
# #   ggplot(aes(index, score)) +
# #   geom_point(stat = "identity",
# #              aes(color = as.character(score)), shape = 16, show.legend = FALSE) +
# #   geom_segment(aes(x = index, y = 0, xend = index, yend = score,
# #                    color = as.character(score)),
# #                show.legend = FALSE) +
# #   theme_classic() +
# #   scale_y_continuous(expand = expansion(mult = c(0, .05))) +
# #   labs(x = "", y = "Score") +
# #   theme(axis.title.x = element_blank(),
# #         axis.text.x = element_blank(),
# #         axis.ticks.x = element_blank(),
# #         axis.title.y = element_text(size = 10),
# #         axis.text.y = element_text(size = 10),
# #         plot.background = element_rect(fill = "transparent", color = NA),
# #         panel.background = element_rect(fill = "transparent", color = NA),
# #         legend.background = element_rect(fill = "transparent", color = NA),
# #         plot.margin = unit(c(0, 0, 0, 0), "pt")
# #         )
# #
# # plot1
# #
# # plot2 <-
# #   score_rule %>%
# #   dplyr::select(-score) %>%
# #   tidyr::pivot_longer(cols = -rule, names_to = "index",
# #                       values_to = "value") %>%
# #   dplyr::mutate(rule = factor(rule,
# #                               levels =
# #                                 c(
# #                                   "Adduct M+H",
# #                                   "Isotope (Adduct M+H)",
# #                                   "Other positive adduct",
# #                                   "Isotope (Other positive adduct)",
# #                                   "Adduct M-H",
# #                                   "Isotope (Adduct M-H)",
# #                                   "Other negative adduct",
# #                                   "Isotope (Other negative adduct)"
# #                                 ) %>% rev()
# #                               )) %>%
# #   ggplot(aes(index, rule)) +
# #   geom_tile(aes(fill = as.character(value)),
# #             color = "#FF6F00FF",
# #             show.legend = FALSE) +
# #   theme_bw() +
# #   labs(x = "", y = "") +
# #   scale_fill_manual(values = c("0" = "white", "1" = "#008EA0FF")) +
# #   theme(panel.grid = element_blank(), axis.title.x = element_blank(),
# #         axis.text.x = element_blank(), axis.ticks.x = element_blank(),
# #         axis.text.y = element_text(size = 10),
# #         plot.margin = unit(c(0, 0, 0, 0), "pt"))
# #
# # plot2
# #
# #
# # library(patchwork)
# #
# # plot <-
# # plot1 + plot2 + patchwork::plot_layout(ncol = 1, heights = c(1,3))
# #
# #
# # plot
# #
# # ggsave(plot, filename = "score_rule.pdf", width = 12, height = 7,
# #        bg = "transparent")
# #
# # ggsave(plot, filename = "score_rule.png", width = 12, height = 7,
# #        bg = "transparent")
#
#
#
# ###remove redunt annotaiton according to compound
# library(plyr)
#
# # annotation_table <-
# #   annotation_table %>%
# #   dplyr::arrange(Lab.ID, compound_class) %>%
# #   plyr::dlply(.(Lab.ID))
#
# save(annotation_table, file = "annotation_table")
#
#
#
# ##-----------------------------------------------------------------------------
# load("annotation_table")
# class(annotation_table)
#
# annotation_table <-
#   annotation_table %>%
#   dplyr::bind_rows()
#
# annotation_redundance <- NULL
#
# annotation_redundance <- c(annotation_redundance,
#                            list(calculate_redundance(x = annotation_table)))
#
# redundance_diff <- c(-1, -1)
#
# round <- 1
#
# while(round < 5 & all(redundance_diff < 0)){
#   cat(round, " ")
#   round <- round + 1
#   ##for one compound, if it has one compound class with score > =100,
#   ##we remove the compound with 20
#   annotation_table <-
#     annotation_table %>%
#     dplyr::group_by(Lab.ID) %>%
#     dplyr::filter(if(any(score >= 100)) score > 20 else score > 0) %>%
#     dplyr::ungroup()
#
#   ##for one peak, if it has one annotation score > =100,
#   ##we remove the annotation <= 30
#   annotation_table <-
#     annotation_table %>%
#     dplyr::group_by(name) %>%
#     dplyr::filter(if(any(score >= 100)) score > 20 else score > 0) %>%
#     dplyr::ungroup()
#
#   annotation_redundance<-
#     c(annotation_redundance, list(calculate_redundance(x = annotation_table)))
#
#   # cat(annotation_redundance[[length(annotation_redundance)]], "\n")
#   redundance_diff <- annotation_redundance[[length(annotation_redundance)]] -
#     annotation_redundance[[length(annotation_redundance) - 1]]
#
# }
#
#
# ####Metabolite set enrichment
# annotation_table <-
#   annotation_table %>%
#   dplyr::mutate(fc = as.numeric(fc)) %>%
#   dplyr::arrange(dplyr::desc(fc))
#
# temp_annotation <-
#   annotation_table %>%
#   dplyr::arrange(dplyr::desc(fc))
#
# load("kegg_hsa_compound_pathway")
#
# intersect(unique(unlist(kegg_hsa_compound_pathway)),
#           unique(annotation_table$Lab.ID))
#
# set1 <- kegg_hsa_compound_pathway[15]
#
# feature_table = temp_annotation
# feature_set = set1
# exponent = 1
# perm_num = 1000
# min_size = 5
# max_size = 1000
# pvalue_cutoff = 0.2
# p_adjust_method = "fdr"
# seed = FALSE
# verbose = TRUE
#
#
# result <- do_gsea2(feature_table = feature_table,
#                    feature_set = kegg_hsa_compound_pathway[15],
#                    exponent = 1,
#                    perm_num = 1000,
#                    min_size = 5,
#                    max_size = 1000,
#                    pvalue_cutoff = 0.25,
#                    p_adjust_method = "fdr",
#                    seed = TRUE, verbose = TRUE)
#
# gsea_plot(x = result, title = result@result$Description[1])
#
# gsea_plot(x = head_result[[idx]],
#           title = paste(head_result[[idx]]@result$p_adjust))
#
#
# gsea_plot(x = tail_result[[idx]],
#           title = paste(tail_result[[idx]]@result$p_adjust))
