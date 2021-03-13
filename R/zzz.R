.onAttach <- function(libname, pkgname){

  needed <- core[!is_attached(core)]
  if (length(needed) == 0)
    return()

  crayon::num_colors(TRUE)
  metPath_attach()


  packageStartupMessage(crayon::green(
    "metPath,
More information can be found at https://jaspershen.github.io/metPath/
Authors: Xiaotao Shen (shenxt@stanford.edu)
Maintainer: Xiaotao Shen.
Version 0.0.9 (20200406)"
  )
  )
}




is_attached <- function(x) {
  paste0("package:", x) %in% search()
}



globalVariables(
  names = c(
    "CAS.ID",
    "Compound.name",
    "Formula",
    "HMDB.ID",
    "KEGG.ID",
    "Lab.ID",
    "RT",
    "Submitter",
    "compound_class",
    "data",
    "hmdb_pathway",
    "keggMS1database",
    "kegg_hsa_pathway",
    "mz",
    "mz.neg",
    "mz.pos",
    "n",
    "name",
    ".",
    "hmdbMS1Database",
    "x_value",
    "x",
    "pathway_name",
    "mapped_id",
    "jaccard_index",
    "geom_node_text",
    "y",
    "p_value",
    "mapped_number",
    "all_number"
  )
)
