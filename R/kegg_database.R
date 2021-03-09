####KEGG database
#' @title get_kegg_compound
#' @description Get compound from KEGG
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param local Use have downloaded database or from online.
#' @param threads threads
#' @importFrom KEGGREST keggList keggGet
#' @export

get_kegg_compound <-
  function(local = TRUE,
           threads = 3) {
    if (local) {
      data("keggMS1database", envir = environment())
      message(
        crayon::yellow(
          "This database is downloaded in",
          keggMS1database@database.info$Version
        )
      )
      cat("\n")
      return(keggMS1database)
    } else{
      message(crayon::yellow("It may take a while...\n"))
      compound_ID <-
        KEGGREST::keggList(database = "compound") %>%
        names() %>%
        unique() %>%
        stringr::str_replace_all(., "cpd:", "")

      temp_fun = function(x) {
        KEGGREST::keggGet(dbentries = x)[[1]]
      }
      kegg_compound_database <-
        pbapply::pblapply(compound_ID, function(x) {
          KEGGREST::keggGet(dbentries = x)[[1]]
        })

      future::plan(future::multisession, workers = threads)
      kegg =
        kegg_compound_database %>%
        furrr::future_map(
          .f = function(x) {
            KEGG.ID = x$ENTRY
            Compound.name = paste(x$NAME, collapse = "{}")
            Formula = x$FORMULA
            if (is.null(x$FORMULA)) {
              Formula = NA
            }
            mz = as.numeric(x$EXACT_MASS)
            if (is.null(x$EXACT_MASS)) {
              mz = NA
            }
            CAS.ID = stringr::str_replace(grep("CAS", x$DBLINKS, value = TRUE), "CAS: ", "") %>%
              stringr::str_trim(side = "both")

            PubChem.ID = stringr::str_replace(grep("PubChem", x$DBLINKS, value = TRUE), "PubChem: ", "") %>%
              stringr::str_trim(side = "both")

            if (length(CAS.ID) == 0) {
              CAS.ID = NA
            }

            if (length(PubChem.ID) == 0) {
              PubChem.ID = NA
            }

            data.frame(
              Lab.ID = KEGG.ID,
              Compound.name,
              Formula,
              mz,
              CAS.ID,
              HMDB.ID = NA,
              KEGG.ID,
              PubChem.ID
            )
          },
          .progress = TRUE
        ) %>%
        do.call(rbind, .) %>%
        as.data.frame()

      kegg =
        kegg %>%
        dplyr::filter(!is.na(mz)) %>%
        dplyr::mutate(synonym = Compound.name) %>%
        dplyr::mutate(
          RT = NA,
          mz.pos = NA,
          mz.neg = NA,
          Submitter = "KEGG"
        ) %>%
        dplyr::select(
          Lab.ID,
          Compound.name,
          mz,
          RT,
          CAS.ID,
          HMDB.ID,
          KEGG.ID,
          Formula,
          mz.pos,
          mz.neg,
          Submitter,
          dplyr::everything()
        )

      kegg$Compound.name =
        kegg$Compound.name %>%
        stringr::str_split(pattern = "\\{\\}") %>%
        purrr::map(function(x) {
          x[1]
        }) %>%
        unlist() %>%
        stringr::str_replace(pattern = ";", "")

      openxlsx::write.xlsx(kegg, file = "kegg.xlsx", asTable = TRUE)

      keggMS1datbase =
        metID::construct_database(
          path = ".",
          version = as.character(Sys.Date()),
          metabolite.info.name = "kegg.xlsx",
          source = "KEGG",
          link = "https://www.genome.jp/kegg/compound/",
          creater = "metPath",
          email = "shenxt@stanford.edu",
          rt = FALSE,
          threads = threads
        )
      return(keggMS1datbase)
      unlink(x = "kegg.xlsx")
    }
  }



#' @title get_kegg_metabolite_pathway
#' @description Get pathway from KEGG
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param local Use have downloaded database or from online.
#' @param organism organism
#' @param threads threads
#' @importFrom KEGGREST keggList keggGet
#' @export

get_kegg_pathway <- function(local = TRUE,
                             organism = "hsa",
                             threads = 3) {
  organism = match.arg(organism)
  if (local) {
    if (organism == "hsa") {
      data("kegg_hsa_pathway", envir = environment())
      message(
        crayon::yellow(
          "This database is downloaded in",
          kegg_hsa_pathway@database_info$version
        )
      )
      cat("\n")
      return(kegg_hsa_pathway)
    }
  } else{
    message(crayon::yellow("It may take a while...\n"))
    organism = match.arg(organism)
    pathway_ID <-
      KEGGREST::keggList(database = "pathway", organism = organism) %>%
      names() %>%
      unique() %>%
      stringr::str_replace_all(., "path:", "")

    kegg_hsa_pathway_database <-
      pbapply::pblapply(pathway_ID, function(x) {
        KEGGREST::keggGet(dbentries = x)[[1]]
      })

    pathway_id =
      kegg_hsa_pathway_database %>%
      purrr::map(function(x) {
        unname(x$ENTRY)
      }) %>%
      unlist()

    pathway_name =
      kegg_hsa_pathway_database %>%
      purrr::map(function(x) {
        unname(x$PATHWAY_MAP)
        # stringr::str_split(pattern = " - ") %>%
        # `[[`(1) %>%
        # `[`(1)
      }) %>%
      unlist()

    pathway_name =
      kegg_hsa_pathway_database %>%
      purrr::map(function(x) {
        unname(x$PATHWAY_MAP)
      }) %>%
      unlist()

    describtion =
      kegg_hsa_pathway_database %>%
      purrr::map(function(x) {
        unname(x$DESCRIPTION)
      })

    pathway_class =
      kegg_hsa_pathway_database %>%
      purrr::map(function(x) {
        unname(x$CLASS)
      })

    gene_list =
      kegg_hsa_pathway_database %>%
      purrr::map(function(x) {
        gene = x$GENE
        if (is.null(gene)) {
          return(data.frame())
        }
        data.frame(
          KEGG.ID = gene[seq(1, length(gene) - 1, by = 2)],
          Gene.name = gene[seq(2, length(gene), by = 2)],
          stringsAsFactors = FALSE
        )
      })

    compound_list =
      kegg_hsa_pathway_database %>%
      purrr::map(function(x) {
        data.frame(
          KEGG.ID = names(x$COMPOUND),
          Compound.name = x$COMPOUND,
          stringsAsFactors = FALSE
        )
      })

    reference_list =
      kegg_hsa_pathway_database %>%
      purrr::map(function(x) {
        purrr::map(
          x$REFERENCE,
          .f = function(y) {
            y = lapply(y, function(z) {
              if (length(z) > 1) {
                paste(z, collapse = "{}")
              } else{
                z
              }
            })
            y = unlist(y)
            if (any(names(y) == "JOURNAL")) {
              names(y)[names(y) == "JOURNAL"] = "JOURNAL1"
              c(y, JOURNAL2 = "")
            }
          }
        ) %>%
          do.call(rbind, .) %>%
          as.data.frame()
      })

    related_disease =
      kegg_hsa_pathway_database %>%
      purrr::map(function(x) {
        data.frame(
          Disease.ID = names(x$DISEASE),
          Disease.name = x$DISEASE,
          stringsAsFactors = FALSE
        )
      })


    related_module =
      kegg_hsa_pathway_database %>%
      purrr::map(function(x) {
        data.frame(
          Module.ID = names(x$MODULE),
          Module.name = x$MODULE,
          stringsAsFactors = FALSE
        )
      })

    pathway =
      new(
        Class = "pathway_database_class",
        database_info = list(source = "KEGG",
                             version = as.character(Sys.Date())),
        pathway_id = pathway_id,
        pathway_name = pathway_name,
        describtion = describtion,
        pathway_class = pathway_class,
        gene_list = gene_list,
        compound_list = compound_list,
        protein_list = list(),
        reference_list = reference_list,
        related_disease = related_disease,
        related_module = related_module
      )
    return(pathway)
  }
}
