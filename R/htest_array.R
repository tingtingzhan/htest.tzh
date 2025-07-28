


#' @title Array of `htest`
#' 
#' @description ..
#' 
#' @param X,Y two \link[base]{matrix}es with same number of rows
#' 
#' @param ... additional parameters of \link[stats]{cor.test}
#' 
#' @details .. 
#' 
#' @returns 
#' Function [outer.cor.test()] returns a [htest_array] object.
#' 
#' @note
#' This is a factorial structure of `xc`-by-`yc`
#' *not* a pairwise *combination* of `x`.
#' 
#' @examples 
#' set.seed(100)
#' x = matrix(rnorm(50), ncol = 5, dimnames = list(NULL, LETTERS[1:5]))
#' y = matrix(rnorm(30), ncol = 3, dimnames = list(NULL, letters[1:3]))
#' m = outer.cor.test(x, y) # do no have print method yet
#' 
#' library(rmd.tzh); list(
#'  '`htest_array`' = m
#' ) |> render_(file = 'htest_array')
#' 
#' @keywords internal
#' @importFrom stats cor.test 
#' @name htest_array
#' @export
outer.cor.test <- function(X, Y = X, ...) {
  
  DNAME <- paste(deparse1(substitute(X)), 'and', deparse1(substitute(Y)))
  
  x <- as.matrix(X) # use S3 generic
  y <- as.matrix(Y) # use S3 generic
  
  if (!(nx <- length(xc <- dimnames(x)[[2L]]))) stop('`x` must colnames')
  if (!(ny <- length(yc <- dimnames(y)[[2L]]))) stop('`y` must colnames')
  
  # ?stats:::cor.test.default errs on different lengths of `x` and `y`
  if (nrow(x) != nrow(y)) stop('')
  
  statistic <- estimate <- p.value <- array(0, dim = c(nx, ny), dimnames = list(xc, yc))
  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      tmp <- cor.test(x = x[,i], y = y[,j], ...)
      statistic[i,j] <- tmp$statistic
      estimate[i,j] <- tmp$estimate
      p.value[i,j] <- tmp$p.value
    }
  }
  
  # just grab the last `tmp`
  # see inside ?stats:::cor.test.default
  names(dimnames(statistic)) <- names(dimnames(p.value)) <- names(dimnames(estimate)) <- c(tmp$method, '')
  tmp$statistic <- statistic
  tmp$p.value <- p.value
  tmp$estimate <- estimate
  tmp$data.name <- DNAME
  class(tmp) <- c('htest_array', class(tmp))
  return(tmp)
  
}


#' @title Markdown Script of [htest_array]
#' 
#' @description ..
#' 
#' @param x an [htest_array] object
#' 
#' @param xnm \link[base]{character} scalar, call of `x`
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns 
#' Function [md_.htest_array()] returns a \link[base]{character} \link[base]{vector}.
#' 
#' @keywords internal
#' @importFrom methods new
#' @importFrom rmd.tzh md_
#' @importClassesFrom rmd.tzh md_lines
#' @export md_.htest_array
#' @export
md_.htest_array <- function(x, xnm, ...) {
  c(
    '```{r}', 
    '#| results: asis',
    sprintf(fmt = '%s$estimate |> as_flextable.matrix() |> set_caption(caption = \'Correlation Coefficients\')', xnm),
    sprintf(fmt = '%s$p.value |> label_pvalue_sym()() |> as_flextable.matrix() |> set_caption(caption = \'p-values\')', xnm), 
    sprintf(fmt = '%s |> p_adjust_.htest_array() |> label_pvalue_sym()() |> as_flextable.matrix() |> set_caption(caption = \'Multiple Testing Adjusted p-values\')', xnm), 
    '```'
  ) |> 
    new(Class = 'md_lines')
}



#' @title Adjusted \eqn{p}-values for [htest_array]
#' 
#' @param x a [htest_array] object
#' 
#' @keywords internal
#' @importFrom flextable.tzh p_adjust_ p_adjust_.numeric
#' @export p_adjust_.htest_array
#' @export
p_adjust_.htest_array <- function(x) {
  pv0 <- x$p.value
  dnm <- dimnames(pv0)
  pv <- c(pv0)
  names(pv) <- c(outer(dnm[[1L]], dnm[[2L]], FUN = \(...) paste(..., sep = ' & ')))
  ret <- p_adjust_.numeric(pv) # 'matrix'
  names(dimnames(ret)) <- c(x$method, '') # for ?as_flextable.matrix
  return(ret)
}



  




