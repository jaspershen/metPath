# # ###database construction
# # ##KEGG
# # library(tidyverse)
# # setwd("/Users/shenxt/projects/getDatabase/KEGG")
# # load("compound/kegg_compound_database")
# #
# # Lab.ID <-
# #   lapply(kegg_compound_database, function(x){
# #    x$ENTRY
# #   }) %>%
# #   unlist() %>%
# #   unname()
# #
# # Compound.name <-
# #   lapply(kegg_compound_database, function(x){
# #     paste(x$NAME, collapse = "")
# #   }) %>%
# #   unlist() %>%
# #   unname()
# #
# # mz <-
# #   lapply(kegg_compound_database, function(x){
# #     paste(x$EXACT_MASS, collapse = "")
# #   }) %>%
# #   unlist() %>%
# #   unname() %>%
# #   as.numeric()
# #
# # Formula <-
# #   lapply(kegg_compound_database, function(x){
# #     paste(x$FORMULA, collapse = "")
# #   }) %>%
# #   unlist() %>%
# #   unname()
# #
# # CAS.ID <-
# #   lapply(kegg_compound_database, function(x){
# #     id <-
# #       try(expr =  grep(pattern = "CAS", x$DBLINKS, value = TRUE) %>%
# #             stringr::str_replace("CAS: ", ""), silent = TRUE)
# #     if(class(id) == "try-error"){
# #     return(NA)
# #     }
# #     if(length(id) == 0){
# #       return(NA)
# #     }
# #     stringr::str_split(id, pattern = " ")[[1]] %>% paste(collapse = ";")
# #   }) %>%
# #   unlist() %>%
# #   unname()
# #
# # HMDB.ID <-
# #   lapply(kegg_compound_database, function(x){
# #     id <-
# #       try(expr =  grep(pattern = "HMDB", x$DBLINKS, value = TRUE) %>%
# #             stringr::str_replace("HMDB: ", ""), silent = TRUE)
# #     if(class(id) == "try-error"){
# #       return(NA)
# #     }
# #     if(length(id) == 0){
# #       return(NA)
# #     }
# #     stringr::str_split(id, pattern = " ")[[1]] %>% paste(collapse = ";")
# #   }) %>%
# #   unlist() %>%
# #   unname()
# #
# # KEGG.ID <- Lab.ID
# #
# # keggMS1database <-
# #   data.frame(Lab.ID,
# #              Compound.name,
# #              mz,
# #              RT = NA,
# #              CAS.ID,
# #              HMDB.ID,
# #              KEGG.ID,
# #              Formula,
# #              mz.pos = NA,
# #              mz.neg = NA,
# #              Submitter = "KEGG",
# #              stringsAsFactors = FALSE)
# #
# # keggMS1database <-
# #   keggMS1database %>%
# #   dplyr::filter(!is.na(mz))
# #
# # ##remove Mn H20 and other items
# # remove_idx <-
# #   keggMS1database$KEGG.ID[which(nchar(keggMS1database$Formula) <= 3)]
# #
# # keggMS1database <-
# #   keggMS1database %>%
# #   dplyr::filter(!KEGG.ID %in% remove_idx)
# #
# # ##remove something like that (C11H9NO2)2. H2SO4
# # keggMS1database %>%
# #   dplyr::filter(Formula == "(C11H9NO2)2. H2SO4")
# #
# # remove_idx <-
# #   keggMS1database %>%
# #   dplyr::filter(stringr::str_detect(Formula, "\\.")) %>%
# #   pull(KEGG.ID)
# #
# # keggMS1database <-
# #   keggMS1database %>%
# #   dplyr::filter(!KEGG.ID %in% remove_idx)
# #
# #
# # database.info <- list(
# #   "Version" = "0.0.1",
# #   "Source" = "Xiaotao Shen",
# #   "Link" = "shenxt.me",
# #   "Creater" = "Xiaotao Shen",
# #   "Email" = "shenxt1990@163.com",
# #   "RT" = FALSE
# # )
# #
# # spectra.info <- keggMS1database
# # Spectra <- list()
# #
# #
# # setClass(
# #   Class = "databaseClass",
# #   representation(
# #     database.info = "list",
# #     spectra.info = "data.frame",
# #     spectra.data = "list"
# #   ),
# #   prototype = list(
# #     database.info = list(),
# #     spectra.info = data.frame(matrix(nrow = 0, ncol = 0), stringsAsFactors = FALSE),
# #     spectra.data = list()
# #   )
# # )
# #
# #
# # kegg_ms1_database <-
# #   new(
# #     Class = "databaseClass",
# #     database.info = database.info,
# #     spectra.info = spectra.info,
# #     spectra.data = Spectra
# #   )
# #
# #
# # save(kegg_ms1_database, file = "kegg_ms1_database")
# #
# #
# # ##Pathway
# # rm(list = ls())
# # load("pathway/kegg_hsa_pathway_database")
# # ENTRY <-
# #   lapply(kegg_hsa_pathway_database, function(x){
# #    x$ENTRY
# #   }) %>%
# #   unlist() %>%
# #   unname()
# #
# # NAME <-
# #   lapply(kegg_hsa_pathway_database, function(x){
# #     x$NAME
# #   }) %>%
# #   unlist() %>%
# #   unname() %>%
# #   stringr::str_replace(" - Homo sapiens \\(human\\)", "")
# #
# # kegg_hsa_compound_pathway <-
# #   lapply(kegg_hsa_pathway_database, function(x){
# #    names(x$COMPOUND) %>% unique()
# #   })
# #
# # names(kegg_hsa_compound_pathway) <-
# #   paste(ENTRY, NAME, sep = ";")
# #
# # remain_idx <- lapply(kegg_hsa_compound_pathway, length) %>% unlist() %>% `>`(0) %>% which()
# #
# # kegg_hsa_compound_pathway <- kegg_hsa_compound_pathway[remain_idx]
# #
# # save(kegg_hsa_compound_pathway, file = "kegg_hsa_compound_pathway")
# #
# #
# # ###KEGG module
# # load("module/kegg_module_database")
# # kegg_module_database[[1]]$ENTRY
# # kegg_module_database[[1]]$NAME
# #
# #
# # ENTRY <-
# #   lapply(kegg_module_database, function(x){
# #     x$ENTRY
# #   }) %>%
# #   unlist() %>%
# #   unname()
# #
# # NAME <-
# #   lapply(kegg_module_database, function(x){
# #     x$NAME
# #   }) %>%
# #   unlist() %>%
# #   unname() %>%
# #   stringr::str_replace(" - Homo sapiens \\(human\\)", "")
# #
# # kegg_hsa_compound_module <-
# #   lapply(kegg_module_database, function(x){
# #     names(x$COMPOUND) %>% unique()
# #   })
# #
# # names(kegg_hsa_compound_module) <-
# #   paste(ENTRY, NAME, sep = ";")
# #
# # remain_idx <-
# #   lapply(kegg_hsa_compound_module, length) %>% unlist() %>% `>`(0) %>% which()
# #
# # kegg_hsa_compound_module <- kegg_hsa_compound_module[remain_idx]
# #
# # save(kegg_hsa_compound_module, file = "kegg_hsa_compound_module")
#
#
#
#
# ##According to PIUMet and GSEA method, we can just from untargeted metabolomics to
# setwd("/Users/shenxt/projects/demark_multiomics/denmark_project/data_analysis_2020_02_13/metabolome_DEG_analysis")
# load("fc_p_value")
#
# setwd("/Users/shenxt/projects/demark_multiomics/denmark_project/data_analysis_2020_02_13/metabolome_data_preparation/peaks/")
# load("variable_info")
#
#
# sxtTools::setwd_project()
# setwd("demo_data/")
# # save(fc_p_value, file = "fc_p_value")
# # save(variable_info, file = "variable_info")
#
# load("fc_p_value")
# load("variable_info")
# load("kegg_hsa_compound_pathway")
# # dim(fc_p_value[[6]])
# # dim(variable_info)
# #
# # kegg_hsa_compound_pathway %>%
# #   lapply(length) %>%
# #   unlist() %>%
# #   `==`(0) %>%
# #   which()
# #
# #
# # table <-
# #   data.frame(variable_info[,c(1:3),], fc_p_value[[6]], stringsAsFactors = FALSE)
# #
# # fc_p_value[[6]]
# #
# # ms1_data <- table %>%
# #   dplyr::select(1:3) %>%
# #   dplyr::mutate(polarity = case_when(
# #    stringr::str_detect(name, "POS") ~ "positive",
# #    TRUE ~ "negative"
# #   ))
# #
# #
# # write.csv(ms1_data, file = "ms1.data.csv", row.names = FALSE)
# #
# #
# # library(metID)
# #
# # ms1_data <- readr::read_csv("ms1.data.csv")
# #
# # ms1_data_pos <-
# #   ms1_data %>%
# #   dplyr::filter(polarity == "positive")
# #
# # ms1_data_neg <-
# #   ms1_data %>%
# #   dplyr::filter(polarity == "negative")
# #
# #
# # if(nrow(ms1_data_pos) == 0){
# #   annotation_table_pos <- NULL
# # }else{
# #   write.csv(ms1_data_pos[,-4], file = "ms1_data_pos.csv", row.names = FALSE)
# #   annotation_table_pos <-
# #     metID::identify_metabolites(ms1.data = "ms1_data_pos.csv",
# #                                 ms1.match.ppm = 15,
# #                                 column = "rp",
# #                                 polarity = "positive",
# #                                 database = "kegg_ms1_database",
# #                                 candidate.num = 1000)
# #   annotation_table_pos <-
# #     metID::get_identification_table(annotation_table_pos, candidate.num = 10000, type = "new")
# #
# #
# #   annotation_table_pos <-
# #     annotation_table_pos %>%
# #     # dplyr::select(name, KEGG.ID, Adduct, mz.error) %>%
# #     dplyr::left_join(table, by = c("name"))
# # }
# #
# #
# # if(nrow(ms1_data_neg) == 0){
# #   annotation_table_neg <- NULL
# # }else{
# #   write.csv(ms1_data_neg[,-4], file = "ms1_data_neg.csv", row.names = FALSE)
# #   annotation_table_neg <-
# #     metID::identify_metabolites(ms1.data = "ms1_data_neg.csv",
# #                                 ms1.match.ppm = 15,
# #                                 column = "rp",
# #                                 polarity = "negative",
# #                                 database = "kegg_ms1_database",
# #                                 candidate.num = 1000)
# #   annotation_table_neg <-
# #     metID::get_identification_table(annotation_table_neg,
# #                                     candidate.num = 10000,
# #                                     type = "new")
# #
# #
# #   annotation_table_neg <-
# #     annotation_table_neg %>%
# #     # dplyr::select(name, KEGG.ID, Adduct, mz.error) %>%
# #     dplyr::left_join(table, by = c("name"))
# # }
# #
# # annotation_table_pos$polarity <- "positive"
# #
# # annotation_table_neg$polarity <- "negative"
# #
# # annotation_table <- rbind(annotation_table_pos, annotation_table_neg)
# #
# # save(annotation_table, file = "annotation_table")
# load("annotation_table")
#
# table_for_gsea <-
# annotation_table %>%
#   dplyr::select(name, polarity, mz.x, rt.x, KEGG.ID, Adduct, mz.error, mz.match.score, p_value, fc) %>%
#   dplyr::rename(mz = mz.x, rt = rt.x) %>%
#   dplyr::mutate(fc = log(fc, 2)) %>%
#   dplyr::arrange(desc(fc))
#
# load("kegg_ms1_database")
# table_for_gsea <-
#   table_for_gsea %>%
#   dplyr::left_join(kegg_ms1_databasee@spectra.info[,c("KEGG.ID", "Formula")], by = "KEGG.ID")
#
# table_for_gsea <-
#   table_for_gsea %>%
#   dplyr::filter(!is.na(KEGG.ID))
#
#
# table_for_gsea <-
#   table_for_gsea %>%
#   dplyr::group_by(name) %>%
#   dplyr::mutate(mz = mz[1], rt = rt[1]) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(mz = as.numeric(mz), rt = as.numeric(rt))
#
# ####isotope annotation
# ms1_table_pos <- readr::read_csv("ms1_data_pos.csv")
#
# ms1_table_neg <- readr::read_csv("ms1_data_neg.csv")
#
# ###positive
#
# system.time(
#   isotope_pos <-
#     purrr::map(
#       as.data.frame(t(table_for_gsea_pos)),
#       .f = function(x) {
#         adduct <-
#           stringr::str_extract(x[6], "\\(.+\\)") %>%
#           stringr::str_replace("\\(", "") %>%
#           stringr::str_replace("\\)", "")
#
#         temp_iso <- try(
#           annotate_isotope(
#             formula = x[11],
#             adduct = adduct,
#             mz = as.numeric(x[3]),
#             rt = as.numeric(x[4]),
#             peak.mz = ms1_table_pos$mz,
#             peak.rt = ms1_data_pos$rt,
#             rt.tol = 5,
#             mz.tol = 15,
#             max.isotope = 3
#           ), silent = TRUE
#         )
#         if(class(temp_iso) == "try-error"){
#           return(NULL)
#         }
#         temp_iso
#       }
#     )
# )
#
#
# length(isotope_pos)
#
# dim(table_for_gsea_pos)
#
#
# # isotope_pos <-
# # purrr::map2
#
#
#
#
# #####for each metabolite, just check them level
# length(unique(table_for_gsea$name))
# length(unique(table_for_gsea$KEGG.ID))
#
#
#
# library(plyr)
#
# table_for_gsea1 <-
#   table_for_gsea %>%
#   plyr::dlply(.variables = .(KEGG.ID)) %>%
#   lapply(function(x){
#   x %>%
#       arrange(rt)
#   })
#
#
# ####combine peaks as group according to rt
# table_for_gsea2 <-
#   table_for_gsea1 %>%
#   purrr::map(.f = function(x){
#   rt_class <- group_peaks_rt(rt = x$rt, rt.tol = 10)
#   rt_class <- paste(x$KEGG.ID[1], rt_class$class, sep = "_")
#   data.frame(x, compound_class = rt_class, stringsAsFactors = FALSE)
#   }) %>%
#   do.call(rbind, .) %>%
#   dplyr::arrange(compound_class)
#
