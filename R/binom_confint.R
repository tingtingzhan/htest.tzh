
#' @title Clopper-Pearson Exact Binomial Confidence Interval
#' 
#' @description 
#' Clopper-Pearson exact binomial confidence interval.
#' 
#' @param x positive \link[base]{integer} scalar or \link[base]{vector}, counts
#' 
#' @param n positive \link[base]{integer} scalar or \link[base]{vector}, sample sizes \eqn{n}
#' 
#' @param conf.level,alternative,... additional parameters of function \link[stats]{binom.test}
#' 
#' @returns 
#' Function [binom_confint()] returns an S3 class `'binom_confint'`, 
#' which is a \link[base]{matrix} with additional \link[base]{attributes}
#' `'conf.level'`, `'alternative'`, `'x'` and `'n'`.
#' 
#' @examples 
#' library(flextable)
#' binom_confint(0:10, 10L) |> as_flextable()
#' binom_confint(0:10, 10L, alternative = 'less') |> as_flextable()
#' binom_confint(0:10, 10L, alternative = 'greater') |> as_flextable()
#' @keywords internal
#' @name binom_confint
#' @importFrom stats binom.test
#' @export
binom_confint <- function(x, ...) UseMethod(generic = 'binom_confint')
  
#' @rdname binom_confint
#' @examples
#' state.region |> binom_confint() |> as_flextable()
#' @export binom_confint.default
#' @export
binom_confint.factor <- function(x, ...) {
  binom_confint.default(x = c(table(x)), n = length(x), ...)
}
  
#' @rdname binom_confint
#' @export binom_confint.default
#' @export
binom_confint.default <- function(x, n, conf.level = .95, alternative = c('two.sided', 'less', 'greater'), ...) {
  
  if (!is.integer(x) || !length(x) || anyNA(x)) stop('x must be integer')
  if (!is.integer(n) || !length(n) || anyNA(n) || any(n <= 0L)) stop('n must be positive integer')
  if (!is.double(conf.level) || length(conf.level) != 1L || anyNA(conf.level) || conf.level < 0 || conf.level > 1) stop('illegal level')
  
  alternative <- match.arg(alternative)
  
  ht <- mapply(
    FUN = binom.test, 
    x = x, n = n, # this is not compute intensive; no need to fully vectorize
    MoreArgs = list(
      conf.level = conf.level,
      alternative = alternative
    ), SIMPLIFY = FALSE) 
  
  cint <- ht |> 
    lapply(FUN = `[[`, 'conf.int') |>
    do.call(what = rbind)
  
  # to match behavior of ?stats::t.test
  attr(cint, which = 'conf.level') <- conf.level
  
  # additional attributes
  
  attr(cint, which = 'alternative') <- alternative
   
  attr(cint, which = 'x') <- ht |>
    vapply(FUN = `[[`, 'statistic', FUN.VALUE = NA_real_) |>
    as.integer()
  
  attr(cint, which = 'n') <- ht |> 
    vapply(FUN = `[[`, 'parameter', FUN.VALUE = NA_real_) |>
    as.integer()
  
  class(cint) <- c('binom_confint', class(cint)) |>
    unique.default()
  return(cint)
  
}



#' @title Convert [binom_confint] to \link[flextable]{flextable}
#' 
#' @param x a [binom_confint]
#' 
#' @param ... ..
#' 
#' @keywords internal
#' @importFrom flextable as_flextable flextable autofit align vline
#' @export as_flextable.binom_confint
#' @export
as_flextable.binom_confint <- function(x, ...) {
  
  obj <- x
  
  x <- attr(obj, which = 'x', exact = TRUE)
  n <- attr(obj, which = 'n', exact = TRUE)
  conf.level <- attr(obj, which = 'conf.level', exact = TRUE)
  alternative <- attr(obj, which = 'alternative', exact = TRUE)
  
  d <- data.frame(
    Count = sprintf(fmt = '%d / %d', x, n),
    sprintf(fmt = '%.1f%% (%.1f%%, %.1f%%)', 1e2*(x/n), 1e2*obj[,1L], 1e2*obj[,2L])
  )
  names(d)[2L] <- sprintf(fmt = 'Percentage\n(%.f%% %s-Sided Exact CI)', 1e2*conf.level, switch(alternative, two.sided = 'Two', 'One'))
  
  nm <- rownames(obj)
  if (length(nm)) d <- data.frame(Name = nm, d, check.names = FALSE)
  
  d |>
    flextable() |>
    autofit(part = 'all') |>
    align(j = if (length(nm)) 2:3 else 1:2, align = 'right', part = 'all') |>
    vline(j = if (length(nm)) 1L)
  
}






#' @title View Binomial Confidence Interval
#' 
#' @description ..
#' 
#' @param x a \link[base]{logical} \link[base]{matrix}
#'
#' @examples 
#' swiss |> is.na() |> viewBinomCI() # no missing
#' airquality |> is.na() |> viewBinomCI()
#' @keywords internal
#' @export
viewBinomCI <- function(x) {
  
  obj <- x; x <- NULL
  if (!is.matrix(obj) || !is.logical(obj)) stop('input must be `logical` `matrix`')
  
  x <- obj |> colSums() |> as.integer()
  id <- (x > 0L)
  if (!any(id)) return(invisible())
  
  n <- obj |> nrow()
  nm <- obj |> colnames()
  
  x_ <- x[id]
  o <- order(x_, decreasing = TRUE)
  .x <- x_[o]
  names(.x) <- nm[id][o]
  
  binom_confint(
    x = .x, 
    n = n
  ) |>
    as_flextable.binom_confint()
  
}






#' @title md_.binom_confint
#' 
#' @description ..
#' 
#' @param x a [binom_confint]
#' 
#' @param xnm ..
#'  
#' @param ... ..
#' 
#' @examples
#' list(
#'  'State Region' = state.region |> binom_confint()
#' ) |> rmd.tzh::render_(file = 'binom_confint')
#' @keywords internal
#' @importFrom methods new
#' @importClassesFrom rmd.tzh md_lines  
#' @importFrom rmd.tzh md_
#' @export md_.binom_confint
#' @export
md_.binom_confint <- function(x, xnm, ...) {
  
  #z1 <- x$call$formula[[2L]] |> 
  #  deparse1() |> 
  #  sprintf(fmt = '@KaplanMeier58 estimates and curves of time-to-event endpoint **`%s`** are obtained using <u>**`R`**</u> package <u>**`survival`**</u>.') |>
  #  new(Class = 'md_lines', package = 'survival', bibentry = KaplanMeier58())
  
  z2 <- c(
    '```{r}',
    '#| echo: false', 
    xnm |> sprintf(fmt = 'as_flextable.binom_confint(%s)'),
    '```'
  ) |> 
    new(Class = 'md_lines')
  
  #c(z1, z2) # ?rmd.tzh::c.md_lines
  return(z2)
  
}




