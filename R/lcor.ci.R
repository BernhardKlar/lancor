#' Confidence intervals for the Lancaster correlation coefficient
#' 
#' @description
#' Computes confidence intervals for the Lancaster correlation coefficient. Lancaster correlation is a bivariate measures of dependence.
#' 
#' @details
#' Computes asymptotic and bootstrap-based confidence intervals for the (linear) Lancaster correlation coefficient \eqn{\rho_L} (\eqn{\rho_{L,1}}). For more details see \code{\link{lcor}}.
#'
#' Asymptotic confidence intervals are derived under two cases (analogously for \eqn{\rho_{L}}; see Holzmann and Klar (2024)):
#'
#' \strong{Case 1:} If \eqn{|\rho_{L1}|\neq|\rho_{L2}|}, the \eqn{1-\alpha} asymptotic interval is
#' \deqn{ \left[ \max\{\hat\rho_{L,1} - z_{1-\alpha/2}\,s/\sqrt{n}, 0\},\ \min\{\hat\rho_{L,1} + z_{1-\alpha/2}\,s/\sqrt{n}, 1\} \right], }
#' where \eqn{z_{1-\alpha/2}} is the standard normal quantile and \eqn{s} is an estimator of the corresponding standard deviation.
#'
#' \strong{Case 2:} If \eqn{|\rho_{L1}|=|\rho_{L2}|=a>0}, let \eqn{\tau} denote the correlation between the two components and let \eqn{q_{1-\alpha/2}} be the \eqn{1-\alpha/2} quantile of the asymptotic distribution of \eqn{\sqrt{n}(\hat\rho_{L,1} - a)}. A conservative asymptotic interval is
#' \deqn{ \left[ \max\{\hat\rho_{L,1} - q_{1-\alpha/2}/\sqrt{n}, 0\},\ \min\{\hat\rho_{L,1} + z_{1-\alpha/2}\,s/\sqrt{n}, 1\} \right]. }
#'
#' Additionally, bootstrap-based intervals can be obtained by resampling and estimating the covariance matrix of the rank or linear correlation components.
#' 
#' @param x a numeric vector, or a matrix or data frame with two columns.
#' @param y NULL (default) or a vector with same length as x.
#' @param conf.level confidence level of the interval.
#' @param type a character string indicating which lancaster correlation is to be computed. One of "rank" (default), or "linear": can be abbreviated.
#' @param con logical; if TRUE (default), conservative asymptotic confidence intervals are computed.
#' @param R number of bootstrap replications. 
#' @param method a character string indicating how the asymptotic covariance matrix is computed if type ="linear". One of "plugin" (default), "boot" or "symmetric": can be abbreviated.
#' 
#' @return a vector containing the lower and upper limits of the confidence interval.
#' 
#'
#' @author Hajo Holzmann, Bernhard Klar
#' 
#' 
#' @references
#' Holzmann, Klar (2024). "Lancester correlation - a new dependence measure linked to maximum correlation". \doi{https://doi.org/10.1111/sjos.12733}
#' 
#' @seealso \code{\link{lcor}, \link{lcor.comp}, \link{lcor.test}}
#' 
#' @examples
#' n <- 1000
#' x <- matrix(rnorm(n*2), n)
#' nu <- 2
#' y <- x / sqrt(rchisq(n, nu)/nu) # multivariate t
#' lcor(y, type = "rank")
#' lcor.ci(y, type = "rank")
#' 
#' @export 
lcor.ci = function(x, y = NULL, conf.level = 0.95, type = c("rank", "linear"),
                   con = TRUE, R = 1000, method = c("plugin", "boot", "pretest")) {
  if (is.data.frame(x)) x = as.matrix(x)
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (!is.matrix(x)) x = cbind(x, y)
  xx = x
  type = match.arg(type)
  if (!(type %in% c("rank", "linear")))
    stop("type for lcor can only be rank or linear")
  requireNamespace("boot", quietly = TRUE)
  #
  if (type == "rank") {
    lc.res = lcor.comp(xx, type="rank")
    rho1 = lc.res[1]
    rho2 = lc.res[2]
    lc = lc.res[3]
    n = dim(xx)[1]
    alpha = 1 - conf.level
    #compute bootstrap covariance matrix of (rho1,rho2)
    Sigma =  n * cov( boot::boot(xx, statistic = function(xx, id) lcor.comp( xx[id, ], type="rank"), R = R)$t[,1:2] )
    #case |rho1| != |rho2|
    if (abs(rho1) >= abs(rho2)) {
      se = sqrt(Sigma[1,1])
    } else {
      se = sqrt(Sigma[2,2])
    }
    z = qnorm(1 - alpha/2)
    l1 = max( lc - z * se / sqrt(n), 0)
    u1 = min( lc + z * se / sqrt(n), 1)
    #case |rho1|=|rho2|>0
    F = function(x, Sigma, p) {
      sigma1 = sqrt(Sigma[1,1])
      sigma2 = sqrt(Sigma[2,2])
      tau = Sigma[1,2] / ( sigma1 * sigma2 )
      alpha1 = (sigma1/sigma2-tau) / sqrt(1-tau^2)
      alpha2 = (sigma2/sigma1-tau) / sqrt(1-tau^2)
      0.5 * ( sn::psn(x,0,sigma1,alpha1) + sn::psn(x,0,sigma2,alpha2) ) - p
    }
    si = 5 * max(Sigma[1,1],Sigma[2,2]) #guess for search interval
    c.u = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = alpha/2)$root
    c.o = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = 1-alpha/2)$root
    l2 = max( lc - c.o / sqrt(n), 0)
    u2 = min( lc - c.u / sqrt(n), 1)
    if (con==TRUE) return(setNames(c(l2, u1), c("lower", "upper"))) else
      return(setNames(c(l1, u1), c("lower", "upper"))) #by construction: l2<l1, u2<u1
  }
  #
  if (type == "linear") {
    method = match.arg(method)
    if (!(method %in% c("plugin", "boot", "pretest")))
      stop("method for test can only be plugin, boot or pretest")
    if (method == "plugin" | method == "boot") {
      lc.res = lcor.comp(xx, type="linear")
      rho1 = lc.res[1]
      rho2 = lc.res[2]
      lc = lc.res[3]
      n = dim(xx)[1]
      alpha = 1 - conf.level
      if (method == "plugin") Sigma = Sigma.est(xx) #general asymptotic theory
      if (method == "boot") { #bootstrap estimate of covariance matrix
        Sigma = n * cov( boot::boot(xx, statistic = function(xx, id) lcor.comp(xx[id, ], type="linear"), R = R)$t[,1:2] )
      }
      #case |rho1| != |rho2|
      if (abs(rho1) >= abs(rho2)) {
        se = sqrt(Sigma[1,1])
      } else {
        se = sqrt(Sigma[2,2])
      }
      z = qnorm(1 - alpha/2)
      l1 = max( lc - z * se / sqrt(n), 0)
      u1 = min( lc + z * se / sqrt(n), 1)
      #case |rho1|=|rho2|>0
      F = function(x, Sigma, p) {
        sigma1 = sqrt(Sigma[1,1])
        sigma2 = sqrt(Sigma[2,2])
        tau = Sigma[1,2] / ( sigma1 * sigma2 )
        alpha1 = (sigma1/sigma2-tau) / sqrt(1-tau^2)
        alpha2 = (sigma2/sigma1-tau) / sqrt(1-tau^2)
        0.5 * ( sn::psn(x,0,sigma1,alpha1) + sn::psn(x,0,sigma2,alpha2) ) - p
      }
      si = 5 * max(Sigma[1,1],Sigma[2,2]) #guess for search interval
      c.u = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = alpha/2)$root
      c.o = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = 1-alpha/2)$root
      l2 = max( lc - c.o / sqrt(n), 0)
      u2 = min( lc - c.u / sqrt(n), 1)
      if (con==TRUE) return(setNames(c(l2, u1), c("lower", "upper"))) else
        return(setNames(c(l1, u1), c("lower", "upper"))) #by construction: l2<l1, u2<u1
    }
    #
    if (method == "pretest") { #asymptotic theory with pretest
      lc.res = lcor.comp(xx, type="linear")
      rho1 = lc.res[1]
      rho2 = lc.res[2]
      lc = lc.res[3]
      n = dim(xx)[1]
      alpha = 1 - conf.level
      Sigma = Sigma.est(xx) #general asymptotic theory
      #perform pretest for equality of |rho1| and |rho2|
      sigma12 = Sigma[1,2]
      #sign change in asymptotic covariance if -rho1=rho2!=0 (case 2)
      if (sign( rho1*rho2 ) < 0) sigma12 = -sigma12
      S = sqrt(n) * (abs(rho1) - abs(rho2)) / sqrt(Sigma[1,1] - 2*sigma12 + Sigma[2,2])
      if (abs(S) < qnorm(1-alpha/2)) { #assume abs(rho1) = abs(rho2) (!=0)
        F = function(x, Sigma, p) {
          sigma1 = sqrt(Sigma[1,1])
          sigma2 = sqrt(Sigma[2,2])
          tau = Sigma[1,2] / ( sigma1 * sigma2 )
          alpha1 = (sigma1/sigma2-tau) / sqrt(1-tau^2)
          alpha2 = (sigma2/sigma1-tau) / sqrt(1-tau^2)
          0.5 * ( sn::psn(x,0,sigma1,alpha1) + sn::psn(x,0,sigma2,alpha2) ) - p
        }
        si = 5 * max(Sigma[1,1],Sigma[2,2]) #guess for search interval
        c.u = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = alpha/2)$root
        c.o = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = 1-alpha/2)$root
        l = max( lc - c.o / sqrt(n), 0)
        u = min( lc - c.u / sqrt(n), 1)
        return( setNames(c(l, u), c("lower", "upper")) )
      }
      #if tests decides for abs(rho1) != abs(rho2)
      if (abs(rho1) >= abs(rho2)) {
        rho = abs(rho1)
        se = sqrt(Sigma[1,1])
      } else {
        rho = abs(rho2)
        se = sqrt(Sigma[2,2])
      }
      z = qnorm(1 - alpha/2)
      l = max( rho - z * se / sqrt(n), 0)
      u = min( rho + z * se / sqrt(n), 1)
      return( setNames(c(l, u), c("lower", "upper")) )
    }
  }
}
