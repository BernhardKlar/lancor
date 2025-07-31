#' Lancaster correlation test
#' 
#' @description
#' Lancaster correlation test of bivariate independence. Lancaster correlation is a bivariate measures of dependence.
#' 
#' @details
#' For more details on the testing procedure see \eqn{Remark \, 2} in Holzmann, Klar (2024).
#' 
#' 
#' @param x a numeric vector, or a matrix or data frame with two columns.
#' @param y NULL (default) or a vector with same length as x
#' @param type a character string indicating which lancaster correlation is to be computed. One of "rank" (default), or "linear": can be abbreviated.
#' @param nperm number of permutations.
#' @param method a character string indicating how the p-value is computed if type ="linear". One of "permutation" (default), "asymptotic" or "symmetric": can be abbreviated.
#' 
#' @return A list containing the following components:
#' \item{lcor}{the value of the test statistic}
#' \item{pval}{the p-value of the test}
#' 
#' @author Hajo Holzmann, Bernhard Klar
#' 
#' @references
#' Holzmann, Klar (2024). "Lancester correlation - a new dependence measure linked to maximum correlation". \doi{https://doi.org/10.1111/sjos.12733} 
#' 
#' @seealso \code{\link{lcor}, \link{lcor.comp}, \link{lcor.ci}} and for performing an ACE permutation test of independence see \code{acepack} (\url{https://cran.r-project.org/package=acepack}). 
#' 
#' @examples 
#' n <- 200
#' x <- matrix(rnorm(n*2), n)
#' nu <- 2
#' y <- x / sqrt(rchisq(n, nu)/nu)
#' cor.test(y[,1], y[,2], method = "spearman")
#' lcor.test(y, type = "rank")
#'
#' @export 
lcor.test = function(x, y = NULL, type = c("rank", "linear"), nperm = 999,
                     method = c("permutation", "asymptotic", "symmetric")) {
  if (is.data.frame(x)) x = as.matrix(x)
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (is.matrix(x)) {
    y = x[,2]
    x = x[,1]
  }
  type = match.arg(type)
  if (!(type %in% c("rank", "linear")))
    stop("type for lcor can only be rank or linear")
  #
  if (type == "rank") {
    lc = lcor(x, y, type="rank")
    n = length(x)
    method = match.arg(method)
    if (!(method %in% c("asymptotic", "permutation")))
      stop("method for rank test can only be permutation or asymptotic")
    if (method == "asymptotic") { #asymptotic theory
      ts = sqrt(n) * lc
      pval = 1 - (2*pnorm(ts) - 1)^2
      return(list(lcor = lc, pval = pval))
    }
    if (method == "permutation") {
      ts = lc
      if (n <= 6) {
        nperm = arrangements::npermutations(n) #use all permutations
        perm = arrangements::permutations(x)
      } else {
        perm = arrangements::permutations(x, nsample = nperm)
      }
      tp = rep(0, nperm)
      for (i in 1:nperm) {
        tp[i] = lcor( perm[i,], y, type="rank")
      }
      pval = (sum(tp > ts) + 1) / (nperm + 1)
      return(list(lcor = lc, pval = pval))
    }
  }
  #
  if (type == "linear") {
    lc = lcor(x, y, type="linear")
    n = length(x)
    if (!(method %in% c("permutation", "asymptotic", "symmetric")))
      stop("method for test can only be permutation, asymptotic or symmetric")
    if (method == "symmetric") { #asymptotic theory with vanishing third moments
      ts = sqrt(n) * lc
      pval = 1 - (2*pnorm(ts) - 1)^2
      return(list(lcor = lc, pval = pval))
    }
    if (method == "asymptotic") { #asymptotic theory with arbitrary third moments
      ts = sqrt(n) * lc
      x = sqrt(n/(n-1)) * scale(x) # factor ensures that mean(x^2)=1, hence e40>1
      y = sqrt(n/(n-1)) * scale(y)
      e30 = mean(x^3)
      e40 = mean(x^4)
      e03 = mean(y^3)
      e04 = mean(y^4)
      tau = e30*e03 / sqrt( (e40-1) * (e04-1) )
      alpha1 = sqrt( (1-tau) / (1+tau) )
      alpha2 = sqrt( (1+tau) / (1-tau) )
      pval = 1 - 2 * ( sn::psn(ts,0,1,alpha1) - sn::psn(0,0,1,alpha1)
                       + sn::psn(-ts,0,1,alpha2) - sn::psn(0,0,1,alpha2) )
      return(list(lcor = lc, pval = pval))
    }
    if (method == "permutation") {
      ts = lc
      if (n <= 6) {
        nperm = arrangements::npermutations(n) #use all permutations
        perm = arrangements::permutations(x)
      } else {
        perm = arrangements::permutations(x, nsample = nperm)
      }
      tp = rep(0, nperm)
      for (i in 1:nperm) {
        tp[i] = lcor( perm[i,], y, type="linear")
      }
      pval = (sum(tp > ts) + 1) / (nperm + 1)
      return(list(lcor = lc, pval = pval))
    }
  }
}
