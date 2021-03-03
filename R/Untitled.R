# geneList <- geneList
# geneSet <- kegg_hsa_compound_pathway[[2]]
#
# lapply(kegg_hsa_compound_pathway, function(x){
#   length(intersect(x, names(geneList)))
# }) %>%
#   unlist() %>%
#   unname()
#
# gseaScores <- function(geneList, geneSet, exponent=1, fortify=FALSE) {
#   ###################################################################
#   ##    geneList                                                   ##
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
#   ## genes defined in geneSet should appear in geneList.
#   ## this is a must, see https://github.com/GuangchuangYu/DOSE/issues/23
#   geneSet <- intersect(geneSet, names(geneList))
#
#   N <- length(geneList)
#   Nh <- length(geneSet)
#
#   Phit <- Pmiss <- numeric(N)
#   hits <- names(geneList) %in% geneSet ## logical
#
#   Phit[hits] <- abs(geneList[hits])^exponent
#   NR <- sum(Phit)
#   Phit <- cumsum(Phit/NR)
#
#   Pmiss[!hits] <-  1/(N-Nh)
#   Pmiss <- cumsum(Pmiss)
#
#   runningES <- Phit - Pmiss
#
#   ## ES is the maximum deviation from zero of Phit-Pmiss
#   max.ES <- max(runningES)
#   min.ES <- min(runningES)
#   if( abs(max.ES) > abs(min.ES) ) {
#     ES <- max.ES
#   } else {
#     ES <- min.ES
#   }
#
#   df <- data.frame(x=seq_along(runningES),
#                    runningScore=runningES,
#                    position=as.integer(hits)
#   )
#
#   if(fortify==TRUE) {
#     return(df)
#   }
#
#   df$gene = names(geneList)
#   res <- list(ES=ES, runningES = df)
#   return(res)
# }
