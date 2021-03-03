# table_for_gsea
#
# plot(table_for_gsea$fc, type = "h")
#
# nrow(table_for_gsea) / length(unique(table_for_gsea$name))
#
#
# feature_set <- kegg_hsa_compound_pathway
#
# intersect(unique(unlist(kegg_hsa_compound_pathway)),
#           unique(table_for_gsea$KEGG.ID))
#
# feature_list <- table_for_gsea$fc
# names(feature_list) <- table_for_gsea$KEGG.ID
#
# feature_set <-
#   filter_feature_set(
#     feature_set = kegg_hsa_compound_pathway,
#     feature_list = feature_list,
#     min_size = 5,
#     max_size = 1000
#   )
#
# table_for_gsea
#
# set1 <- feature_set[[1]]
#
# match_idx <-
#   lapply(set1, function(x) {
#     which(x == table_for_gsea$KEGG.ID)
#   })
#
# remain_idx <-
#   lapply(match_idx, length) %>%
#   unlist() %>%
#   `>`(0) %>%
#   which()
#
# set1 <- set1[remain_idx]
# match_idx <- match_idx[remain_idx]
#
# lapply(match_idx, function(idx) {
#   table_for_gsea[idx, ]
# })
#
# ###First we think the pathways are enriched in the head of ranked list
# table_for_gsea1 <-
#   table_for_gsea %>%
#   dplyr::arrange(desc(fc))
#
# table_for_gsea1 <-
#   table_for_gsea1[!duplicated(table_for_gsea1$KEGG.ID), ] %>%
#   dplyr::arrange(desc(fc))
#
#
# ###First we think the pathways are enriched in the tail of ranked list
# table_for_gsea2 <-
#   table_for_gsea %>%
#   dplyr::arrange(fc)
#
# table_for_gsea2 <-
#   table_for_gsea2[!duplicated(table_for_gsea2$KEGG.ID), ] %>%
#   dplyr::arrange(desc(fc))
#
# # match_idx <-
# #   lapply(set1, function(x){
# #     which(x == table_for_gsea$KEGG.ID)
# #   })
# #
# # match_idx1 <-
# #   lapply(set1, function(x){
# #     which(x == table_for_gsea1$KEGG.ID)
# #   })
#
#
#
# # feature_list1 <-
# # table_for_gsea1$fc
# # names(feature_list1) <- table_for_gsea1$KEGG.ID
# #
# # feature_list2 <-
# #   table_for_gsea2$fc
# # names(feature_list2) <- table_for_gsea2$KEGG.ID
#
# head_result <-
#   vector(mode = "list",
#          length = length(kegg_hsa_compound_pathway))
# tail_result <-
#   vector(mode = "list",
#          length = length(kegg_hsa_compound_pathway))
#
# names(head_result) <- names(tail_result) <-
#   names(kegg_hsa_compound_pathway)
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
# for (i in seq_along(kegg_hsa_compound_pathway)) {
#   cat(i, " ")
#   gene_set <- kegg_hsa_compound_pathway[i]
#
#   ##head
#   ##remove the redunt annotation for each metabolites in gene_set
#   peak_name <-
#     lapply(gene_set[[1]], function(x) {
#       which(table_for_gsea1$KEGG.ID == x)
#     }) %>%
#     unlist() %>%
#     `[`(table_for_gsea1, ., ) %>%
#     pull(name)
#
#   temp_data1 <-
#     table_for_gsea1 %>%
#     dplyr::filter(name %in% peak_name)
#
#   temp_data2 <-
#     table_for_gsea1 %>%
#     dplyr::filter(!name %in% peak_name)
#
#   temp_data1 <-
#     temp_data1 %>%
#     dplyr::filter(KEGG.ID %in% gene_set[[1]])
#
#   temp_data2 <-
#     temp_data2 %>%
#     dplyr::group_by(name) %>%
#     dplyr::filter(mz != "") %>%
#     dplyr::ungroup()
#
#   temp_data <-
#     rbind(temp_data1, temp_data2) %>%
#     dplyr::arrange(desc(fc))
#
#   feature_list1 <- temp_data$fc
#   names(feature_list1) <- temp_data$KEGG.ID
#
#   head_result[[i]] <- do_gsea(
#     feature_list = feature_list1,
#     feature_set = gene_set,
#     exponent = 1,
#     perm_num = 1000,
#     min_size = 5,
#     max_size = 1000,
#     pvalue_cutoff = 0.25,
#     p_adjust_method = "fdr",
#     seed = TRUE,
#     verbose = TRUE
#   )
#
#
#
#   ##tail
#   ##remove the redunt annotation for each metabolites in gene_set
#   peak_name <-
#     lapply(gene_set[[1]], function(x) {
#       which(table_for_gsea2$KEGG.ID == x)
#     }) %>%
#     unlist() %>%
#     `[`(table_for_gsea2, ., ) %>%
#     pull(name)
#
#   temp_data1 <-
#     table_for_gsea2 %>%
#     dplyr::filter(name %in% peak_name)
#
#   temp_data2 <-
#     table_for_gsea2 %>%
#     dplyr::filter(!name %in% peak_name)
#
#   temp_data1 <-
#     temp_data1 %>%
#     dplyr::filter(KEGG.ID %in% gene_set[[1]])
#
#   temp_data2 <-
#     temp_data2 %>%
#     dplyr::group_by(name) %>%
#     dplyr::filter(mz != "") %>%
#     dplyr::ungroup()
#
#   temp_data <-
#     rbind(temp_data1, temp_data2) %>%
#     dplyr::arrange(desc(fc))
#
#   feature_list2 <- temp_data$fc
#   names(feature_list2) <- temp_data$KEGG.ID
#
#   tail_result[[i]] <- do_gsea(
#     feature_list = feature_list2,
#     feature_set = gene_set,
#     exponent = 1,
#     perm_num = 1000,
#     min_size = 5,
#     max_size = 1000,
#     pvalue_cutoff = 0.25,
#     p_adjust_method = "fdr",
#     seed = TRUE,
#     verbose = TRUE
#   )
#
# }
#
#
# #####output plots
# dir.create("head_result")
# lapply(head_result, function(x) {
#   if (is.null(x)) {
#     return(NULL)
#   }
#
#   if (nrow(x@result) == 0) {
#     return(NULL)
#   } else{
#     plot <-
#       gsea_plot(
#         x = x,
#         feature_set_idx = 1,
#         title = x@result$p_adjust
#       )
#     name <- stringr::str_split(x@feature_set %>% names(), ";")[[1]][1]
#     ggsave(
#       plot,
#       filename = file.path("head_result", paste(name, ".png", sep = "")),
#       width = 7,
#       height = 7
#     )
#   }
#
#
# })
#
#
# #####output plots
# dir.create("tail_result")
# lapply(tail_result, function(x) {
#   if (is.null(x)) {
#     return(NULL)
#   }
#
#   if (nrow(x@result) == 0) {
#     return(NULL)
#   } else{
#     plot <-
#       gsea_plot(
#         x = x,
#         feature_set_idx = 1,
#         title = x@result$p_adjust
#       )
#     name <-
#       stringr::str_split(x@feature_set %>% names(), ";")[[1]][1]
#     ggsave(
#       plot,
#       filename = file.path("tail_result", paste(name, ".png", sep = "")),
#       width = 7,
#       height = 7
#     )
#   }
#
#
# })
#
#
#
# lapply(head_result, function(x) {
#   if (is.null(x)) {
#     return(0)
#   }
#   nrow(x@result)
# }) %>% unlist() %>%
#   unname() %>%
#   plot()
#
# idx <- 22
#
# gsea_plot(x = head_result[[idx]],
#           title = paste(head_result[[idx]]@result$p_adjust))
#
#
# gsea_plot(x = tail_result[[idx]],
#           title = paste(tail_result[[idx]]@result$p_adjust))
