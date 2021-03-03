core <- c("metID")

core_unloaded <- function() {
  search <- paste0("package:", core)
  core[!search %in% search()]
}

# Attach the package from the same package library it was
same_library <- function(pkg) {
  loc <-
    if (pkg %in% loadedNamespaces())
      dirname(getNamespaceInfo(pkg, "path"))
  do.call("library",
          list(
            pkg,
            lib.loc = loc,
            character.only = TRUE,
            warn.conflicts = FALSE
          ))
}

metPath_attach <- function() {
  to_load <- core_unloaded()
  if (length(to_load) == 0)
    return(invisible())

  # msg(
  #   cli::rule(
  #     left = crayon::bold("Attaching packages"),
  #     right = paste0("metPath ", package_version("metPath"))
  #   ),
  #   startup = TRUE
  # )

  versions <- vapply(to_load, package_version, character(1))
  packages <- paste0(
    crayon::green(cli::symbol$tick),
    " ",
    crayon::blue(format(to_load)),
    " ",
    crayon::col_align(versions, max(crayon::col_nchar(versions)))
  )

  if (length(packages) %% 2 == 1) {
    packages <- append(packages, "")
  }
  col1 <- seq_len(length(packages) / 2)
  info <- paste0(packages[col1], "     ", packages[-col1])

  msg(paste(info, collapse = "\n"), startup = TRUE)

  suppressPackageStartupMessages(lapply(to_load, same_library))

  invisible()
}

package_version <- function(x) {
  version <- as.character(unclass(utils::packageVersion(x))[[1]])

  if (length(version) > 3) {
    version[4:length(version)] <-
      crayon::red(as.character(version[4:length(version)]))
  }
  paste0(version, collapse = ".")
}


msg <- function(..., startup = FALSE) {
  if (startup) {
    if (!isTRUE(getOption("metPath.quiet"))) {
      packageStartupMessage(text_col(...))
    }
  } else {
    message(text_col(...))
  }
}

text_col <- function(x) {
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }

  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }

  theme <- rstudioapi::getThemeInfo()

  if (isTRUE(theme$dark))
    crayon::white(x)
  else
    crayon::black(x)

}

#' List all packages in the metPath
#'
#' @param include_self Include metPath in the list?
#' @export
#' @examples
#' metPath_packages()
metPath_packages <- function(include_self = TRUE) {
  raw <- utils::packageDescription("metPath")$Imports
  imports <- strsplit(raw, ",")[[1]]
  parsed <- gsub("^\\s+|\\s+$", "", imports)
  names <-
    vapply(strsplit(parsed, "\\s+"), "[[", 1, FUN.VALUE = character(1))

  if (include_self) {
    names <- c(names, "metPath")
  }

  names
}

invert <- function(x) {
  if (length(x) == 0)
    return()
  stacked <- utils::stack(x)
  tapply(as.character(stacked$ind), stacked$values, list)
}


style_grey <- function(level, ...) {
  crayon::style(paste0(...),
                crayon::make_style(grDevices::grey(level), grey = TRUE))
}
