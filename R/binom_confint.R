
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
#' @importFrom stats binom.test
#' @export
binom_confint <- function(x, n, conf.level = .95, alternative = c('two.sided', 'less', 'greater'), ...) {
  
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
  
  class(cint) <- 'binom_confint'
  return(cint)
  
}



#' @title Convert [binom_confint] to \link[flextable]{flextable}
#' 
#' @param x a [binom_confint]
#' 
#' @param ... ..
#' 
#' @keywords internal
#' @importFrom flextable as_flextable flextable autofit
#' @export as_flextable.binom_confint
#' @export
as_flextable.binom_confint <- function(x, ...) {
  
  obj <- x
  
  x <- attr(obj, which = 'x', exact = TRUE)
  n <- attr(obj, which = 'n', exact = TRUE)
  nm <- attr(obj, which = 'nm', exact = TRUE)
  conf.level <- attr(obj, which = 'conf.level', exact = TRUE)
  alternative <- attr(obj, which = 'alternative', exact = TRUE)
  
  d <- data.frame(
    Count = sprintf(fmt = '%d / %d', x, n),
    sprintf(fmt = '%.1f%% (%.1f%%, %.1f%%)', 1e2*(x/n), 1e2*obj[,1L], 1e2*obj[,2L])
  )
  names(d)[2L] <- sprintf(fmt = 'Percentage\n(%.f%% %s-Sided Exact CI)', 1e2*conf.level, switch(alternative, two.sided = 'Two', 'One'))
  
  if (length(nm)) d <- data.frame(Name = nm, d, check.names = FALSE)
  
  d |>
    flextable() |>
    autofit(part = 'all')
  
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
  
  ret <- binom_confint(
    x = x_[o], 
    n = n
  )
  attr(ret, which = 'nm') <- nm[id][o]
  
  as_flextable.binom_confint(ret)
  
}




