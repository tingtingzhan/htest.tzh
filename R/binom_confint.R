
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
#' The `S3` methods of the generic function [binom_confint()] returns an `S3` class `'binom_confint'`, 
#' which is a \link[base]{matrix} with additional \link[base]{attributes}
#' `'conf.level'`, `'alternative'`, `'x'` and `'n'`.
#' 
#' @examples 
#' binom_confint(0:10, 10L)
#' binom_confint(0:10, 10L, alternative = 'less')
#' binom_confint(0:10, 10L, alternative = 'greater')
#' @keywords internal
#' @name binom_confint
#' @importFrom stats binom.test
#' @export
binom_confint <- function(x, ...) UseMethod(generic = 'binom_confint')
  
#' @rdname binom_confint
#' @examples
#' state.region |> binom_confint()
#' @export binom_confint.default
#' @export
binom_confint.factor <- function(x, ...) {
  binom_confint.default(x = c(table(x)), n = length(x), ...)
}


#' @rdname binom_confint
#' @examples 
#' swiss |> is.na() |> binom_confint() # no missing
#' airquality |> is.na() |> binom_confint() 
#' @keywords internal
#' @export binom_confint.matrix
#' @export
binom_confint.matrix <- function(x, n = c('col', 'row'), ...) {
  
  obj <- x; x <- NULL
  if (!is.matrix(obj) || !is.logical(obj)) stop('input must be `logical` `matrix`')
  
  switch(EXPR = match.arg(n), col = {
    obj |> 
      colSums() |> # percentages *by column*
      as.integer() |> 
      setNames(nm = colnames(obj)) |>
      binom_confint.default(n = nrow(obj), ...)
  }, row = {
    obj |> 
      rowSums() |> # percentages *by row*
      as.integer() |> 
      setNames(nm = rownames(obj)) |>
      binom_confint.default(n = ncol(obj), ...)
  })

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
    x = x, n = n, # not compute intensive; no need to fully vectorize
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



#' @title Subset [binom_confint]
#' 
#' @description
#' Subset a [binom_confint].
#' 
#' @param x [binom_confint]
#' 
#' @param subset \link[base]{language}
#' 
#' @param ... additional parameters, currently of no use
#' 
#' @examples
#' (cint = penguins |> is.na() |> binom_confint())
#' cint |> 
#'   subset(subset = x > 0L) |>
#'   sort()
#' 
#' @keywords internal
#' @export subset.binom_confint
#' @export
subset.binom_confint <- function(x, subset, ...) {
  
  if (missing(subset)) return(x) # exception handling
  
  obj <- x; x <- NULL
  x <- obj |>
    attr(which = 'x', exact = TRUE)
  
  e <- substitute(subset)
  id <- eval(e) # only allow `x` in `subset`
  return(obj[id]) # `[.binom_confint`
  
}


#' @export
`[.binom_confint` <- function(x, i) {
  z <- unclass(x)[i, , drop = FALSE]
  attr(z, which = 'conf.level') <- attr(x, which = 'conf.level', exact = TRUE)
  attr(z, which = 'alternative') <- attr(x, which = 'alternative', exact = TRUE)
  attr(z, which = 'x') <- attr(x, which = 'x', exact = TRUE)[i]
  attr(z, which = 'n') <- attr(x, which = 'n', exact = TRUE)[i]
  class(z) <- c('binom_confint', class(z)) |>
    unique.default()
  return(z)
}









#' @title Sort [binom_confint]
#' 
#' @description
#' Sort a [binom_confint] by count `x`, if the sample size `n` is the same.
#' 
#' @param x [binom_confint]
#' 
#' @param decreasing,... parameters of function \link[base]{order}
#' 
#' @examples
#' penguins |> is.na() |> binom_confint()
#' penguins |> is.na() |> binom_confint() |> sort()
#' 
#' @keywords internal
#' @export sort.binom_confint
#' @export
sort.binom_confint <- function(x, decreasing = TRUE, ...) {
  
  obj <- x; x <- NULL
  
  n <- obj |>
    attr(which = 'n', exact = TRUE)
  if (!all(duplicated.default(n)[-1L])) {
    warning('cannot sort if not all `n`s are equal')
    return(x) # exception handling
  }
  
  o <- obj |>
    attr(which = 'x', exact = TRUE) |>
    order(decreasing = decreasing, ...)
  return(obj[o]) # `[.binom_confint`
  
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
  
  obj0 <- unclass(obj)
  
  d <- data.frame(
    Count = sprintf(fmt = '%d / %d', x, n),
    Percentage = sprintf(fmt = '%.1f%%', 1e2*(x/n)),
    sprintf(fmt = '(%.1f%%, %.1f%%)', 1e2*obj0[,1L], 1e2*obj0[,2L])
  )
  names(d)[3L] <- sprintf(fmt = '%.f%% %s-Sided Exact\nConfidence Interval', 1e2*conf.level, switch(alternative, two.sided = 'Two', 'One'))
  
  nm <- rownames(obj)
  if (length(nm)) d <- data.frame(Name = nm, d, check.names = FALSE)
  
  d |>
    flextable() |>
    autofit(part = 'all') |>
    align(j = if (length(nm)) 2:4 else 1:3, align = 'right', part = 'all') |>
    vline(j = if (length(nm)) 1L else integer())
  
}



#' @export
print.binom_confint <- function(x, ...) {
  x |>
    as_flextable.binom_confint(...) |>
    print() # ?flextable:::print.flextable
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
#'  'State Region' = state.region |> binom_confint() |> sort()
#' ) |> fastmd::render_(file = 'binom_confint')
#' @keywords internal
#' @importFrom methods new
#' @importClassesFrom fastmd md_lines  
#' @importFrom fastmd md_
#' @export md_.binom_confint
#' @export
md_.binom_confint <- function(x, xnm, ...) {
  
  z1 <- c(
    '```{r}',
    '#| echo: false', 
    xnm |> sprintf(fmt = 'as_flextable.binom_confint(%s)'),
    '```'
  ) |> 
    new(Class = 'md_lines')
  
  return(z1)
  
}




