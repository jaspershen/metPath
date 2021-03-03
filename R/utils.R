#' @title group peaks according to RT
#' @description group peaks according to RT
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param rt rt
#' @param rt.tol rt.tol
#' @export
## group peaks according to RT
group_peaks_rt <- function(rt, rt.tol = 10) {
  rt_class <-
    lapply(rt, function(x) {
      which(rt >= x & rt < x + rt.tol)
    })

  rt_class <-
    lapply(seq_along(rt_class)[-1], function(i) {
      setdiff(rt_class[[i]], unlist(rt_class[1:(i - 1)]) %>% unique())
    }) %>%
    `c`(rt_class[1], .)

  rt_class <-
    rt_class[which(lapply(rt_class, length) != 0)]

  names(rt_class) <- seq_along(rt_class)

  rt_class <-
    purrr::map2(
      .x = rt_class,
      .y = names(rt_class),
      .f = function(x, y) {
        data.frame(rt = rt[x],
                   class = y,
                   stringsAsFactors = FALSE)
      }
    ) %>%
    do.call(rbind, .)
  rownames(rt_class) <- NULL
  return(rt_class)
}














###positive
### [M+H] 50
### [M+H] isotope 20
### other adduct 20
### other adduct 10

###negative
### [M-H] 50
### [M-H] isotope 20
### other adduct 20
### other adduct 10


#' @title score_peak_group
#' @description score_peak_group
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param peak_group peak_group
#' @export
score_peak_group <- function(peak_group){
score <- 0

if(any(peak_group$Adduct == "(M+H)+")){
  score <- score + 50
}

if(any(peak_group$Adduct == "(M+H)+" & peak_group$isotope != "[M]")){
  score <- score + 20
}

if(any(peak_group$polarity == "positive" & peak_group$Adduct != "(M+H)+")){
  score <- score + 20
}

if(any(peak_group$polarity == "positive" & peak_group$Adduct != "(M+H)+" &
       peak_group$isotope != "[M]")){
  score <- score + 10
}

if(any(peak_group$Adduct == "(M-H)-")){
  score <- score + 50
}

if(any(peak_group$Adduct == "(M-H)-" & peak_group$isotope != "[M]")){
  score <- score + 20
}

if(any(peak_group$polarity == "negative" & peak_group$Adduct != "(M-H)-")){
  score <- score + 20
}


if(any(peak_group$polarity == "negative" & peak_group$Adduct != "(M-H)-" &
       peak_group$isotope != "[M]")){
  score <- score + 10
}

return(score)
}


#' @title calculate_redundance
#' @description calculate_redundance
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param x x
#' @export
calculate_redundance <- function(x){
x <-
  dplyr::bind_rows(x)
x1 <-
  x %>%
  dplyr::group_by(Lab.ID) %>%
  dplyr::distinct(compound_class) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::pull(n) %>%
  mean()
x2 <-
  x %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::pull(n) %>%
  mean()
c(x1, x2)
}
