#' @title enrich_pathway
#' @description Pathway enrichment for metabolomics.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param metabolite_id The vector of gene IDs.
#' @param pathway_database KEGG or other metabolomics pathway database.
#' @param method Hypergeometric or fisher test.
#' @return  The MSE analysis result.
#' @export

setGeneric(
  name = "enrich_pathway",
  def = function(metabolite_id,
                 pathway_database,
                 method = c("hypergeometric", "fisher")) {

    id <- unique(id)
    method <- match.arg(method)

    all_id <-
      database %>%
      unlist() %>%
      unique() %>%
      unname() %>%
      as.character()

    ##remove the ID which is not in the all_id
    id <- id[which(is.element(id, all_id))]

    if (length(id) == 0) {
      return(NULL)
    }

    sig_id <- as.character(id)

    num_all <- length(all_id)
    num_sig <- length(sig_id)

    # --------------------------------------------------------------
    #Generating label matrix for detected metabolites
    # --------------------------------------------------------------

    all_matrix <- setlabel(all_id, database)

    # delete metabolite set
    all_matrix2 <-
      all_matrix[, colSums(all_matrix) != 0, drop = FALSE]

    # error handling
    if (ncol(all_matrix2) < 2) {
      # stop function
      return(NULL)
    }

    # --------------------------------------------------------------
    #Generating label matrix for significant metabolites
    # --------------------------------------------------------------

    sig_matrix <- setlabel(sig_id, database)
    sig_matrix <-
      sig_matrix[, colSums(all_matrix) != 0, drop = FALSE]

    # -------------------------------
    #Calculating  ORA
    # -------------------------------

    p_value <- rep(NA, ncol(all_matrix2))

    #for each pathway

    for (i in 1:ncol(all_matrix2)) {
      # ------------------------------------
      #Generating 2~2 table
      # -------------------------------------
      a1 <- sum(sig_matrix[, i])# significant and including pathway
      a2 <-
        sum(all_matrix2[, i]) - sum(sig_matrix[, i])# not significant and including pathway
      a3 <-
        length(sig_id) - a1# significant and not including pathway
      a4 <-
        (length(all_id) - length(sig_id)) - a2# not significant and not including pathway

      tab <- t(matrix(c(a1, a2, a3, a4), 2))

      if (method == "hypergeometric") {
        p_value[i] <-
          phyper(
            q = a1 - 1,
            m = sum(all_matrix2[, i]),
            n = num_all - sum(all_matrix2[, i]),
            k = num_sig,
            lower.tail = FALSE
          )
      } else{
        # ----------------------------------
        # Fisher's exact test
        # ----------------------------------
        check <- tryCatch({
          resfish <- fisher.test(tab, alternative = "greater")
        }, error = function(e) {
          NA
        })

        if (class(check) != "htest") {
          p_value[i] <- 1
        } else{
          resfish <- fisher.test(tab, alternative = "greater")
          p_value[i] <- resfish$p.value
        }
      }

    }

    # -----------------------
    #q-value
    # -----------------------
    p_value_bh <- p.adjust(p_value, method = "BH")
    p_value_fdr <- p.adjust(p_value, method = "fdr")

    # ----------------------
    #Result
    # ----------------------
    PQ <- data.frame(p_value, p_value_bh, p_value_fdr,
                     stringsAsFactors = FALSE)
    rownames(PQ) <- colnames(sig_matrix)
    colnames(PQ) <- c("p.value", "p.value.bh", "p.value.fdr")

    PQ <-
      PQ %>%
      rownames_to_column(., "Pathway.name.ID")

    ##calculate the impact of pathway
    info <- lapply(database, function(module) {
      overlap.number <- length(intersect(module, id))
      pathway.number <- length(module)
      c(pathway.number, overlap.number)
    })

    info <- do.call(rbind, info)
    colnames(info) <- c("Pathway.length", "Overlap")

    info <-
      info %>%
      as_tibble() %>%
      mutate(Overlap.frac = Overlap / Pathway.length)

    rownames(info) <-
      names(database)

    info <-
      info %>%
      rownames_to_column(var = 'Pathway.name.ID')

    info <-
      info %>%
      left_join(PQ, info, by = 'Pathway.name.ID')

    info <-
      info %>%
      mutate(Pathway.name = unlist(lapply(stringr::str_split(info$Pathway.name.ID, ";"), function(x)
        x[1])),
        Pathway.ID = unlist(lapply(stringr::str_split(info$Pathway.name.ID, ";"), function(x)
          x[2]))) %>%
      select(., Pathway.name, Pathway.ID, everything(), -Pathway.name.ID) %>%
      filter(., Overlap != 0) %>%
      arrange(., p.value)

    info <- info

  }
)



#' @title setlabel
#' @description setlabel
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param M_ID M_ID.
#' @param M M
#' @return  The MSE analysis result.

setlabel <- function (M_ID, M) {
  temp <- lapply(M_ID, function(x) {
    unlist(lapply(M, function(y) {
      sum(is.element(x, y))
    }))
  })
  temp <- do.call(rbind, temp)
  rownames(temp) <- M_ID
  temp <- temp
}
