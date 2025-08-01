#' Covariance matrix of components of Lancaster correlation coefficient
#' 
#' @description 
#' Estimate of covariance matrix of the two components of Lancaster correlation. Lancaster correlation is a bivariate measures of dependence.
#'
#' @details 
#' For more details see the Appendix in Holzmann, Klar (2024).
#' 
#' @param xx a matrix or data frame with two columns.
#' 
#' @return the estimated covariance matrix.
#'
#' @author Hajo Holzmann, Bernhard Klar
#' 
#' @references 
#' Holzmann, Klar (2024). "Lancester correlation - a new dependence measure linked to maximum correlation". \doi{https://doi.org/10.1111/sjos.12733}
#' 
#' @seealso \code{\link{lcor.ci}}
#' 
#' @examples 
#' Sigma <- matrix(c(1,0.1,0.1,1), ncol=2)
#' R <- chol(Sigma)
#' n <- 1000
#' x <- matrix(rnorm(n*2), n)
#' nu <- 8
#' y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t
#' Sigma.est(y)
#'
#' @export 
Sigma.est = function(xx) {
  x = xx[,1]
  n = length(x)
  x = sqrt(n/(n-1)) * scale( as.matrix(x) )  # mean(x^2)=1
  y = xx[,2]
  y = sqrt(n/(n-1)) * scale( as.matrix(y) )
  e = matrix(nrow = 9, ncol = 9)
  for (i in 1:9){ for (j in 1:9) e[i,j] = mean( x^{i-1}*y^{j-1} ) }
  rho = e[1+1,1+1]

  EUV = matrix( c(
    1,         rho,       e[3+1,0+1],e[1+1,2+1],e[2+1,1+1],e[4+1,0+1],e[1+1,3+1],e[3+1,1+1],e[2+1,2+1],e[5+1,0+1],e[1+1,4+1],e[3+1,2+1],
    rho,       1,         e[2+1,1+1],e[0+1,3+1],e[1+1,2+1],e[3+1,1+1],e[0+1,4+1],e[2+1,2+1],e[1+1,3+1],e[4+1,1+1],e[0+1,5+1],e[2+1,3+1],
    e[3+1,0+1],e[2+1,1+1],e[4+1,0+1],e[2+1,2+1],e[3+1,1+1],e[5+1,0+1],e[2+1,3+1],e[4+1,1+1],e[3+1,2+1],e[6+1,0+1],e[2+1,4+1],e[4+1,2+1],
    e[1+1,2+1],e[0+1,3+1],e[2+1,2+1],e[0+1,4+1],e[1+1,3+1],e[3+1,2+1],e[0+1,5+1],e[2+1,3+1],e[1+1,4+1],e[4+1,2+1],e[0+1,6+1],e[2+1,4+1],
    e[2+1,1+1],e[1+1,2+1],e[3+1,1+1],e[1+1,3+1],e[2+1,2+1],e[4+1,1+1],e[1+1,4+1],e[3+1,2+1],e[2+1,3+1],e[5+1,1+1],e[1+1,5+1],e[3+1,3+1],
    e[4+1,0+1],e[3+1,1+1],e[5+1,0+1],e[3+1,2+1],e[4+1,1+1],e[6+1,0+1],e[3+1,3+1],e[5+1,1+1],e[4+1,2+1],e[7+1,0+1],e[3+1,4+1],e[5+1,2+1],
    e[1+1,3+1],e[0+1,4+1],e[2+1,3+1],e[0+1,5+1],e[1+1,4+1],e[3+1,3+1],e[0+1,6+1],e[2+1,4+1],e[1+1,5+1],e[4+1,3+1],e[0+1,7+1],e[2+1,5+1],
    e[3+1,1+1],e[2+1,2+1],e[4+1,1+1],e[2+1,3+1],e[3+1,2+1],e[5+1,1+1],e[2+1,4+1],e[4+1,2+1],e[3+1,3+1],e[6+1,1+1],e[2+1,5+1],e[4+1,3+1],
    e[2+1,2+1],e[1+1,3+1],e[3+1,2+1],e[1+1,4+1],e[2+1,3+1],e[4+1,2+1],e[1+1,5+1],e[3+1,3+1],e[2+1,4+1],e[5+1,2+1],e[1+1,6+1],e[3+1,4+1],
    e[5+1,0+1],e[4+1,1+1],e[6+1,0+1],e[4+1,2+1],e[5+1,1+1],e[7+1,0+1],e[4+1,3+1],e[6+1,1+1],e[5+1,2+1],e[8+1,0+1],e[4+1,4+1],e[6+1,2+1],
    e[1+1,4+1],e[0+1,5+1],e[2+1,4+1],e[0+1,6+1],e[1+1,5+1],e[3+1,4+1],e[0+1,7+1],e[2+1,5+1],e[1+1,6+1],e[4+1,4+1],e[0+1,8+1],e[2+1,6+1],
    e[3+1,2+1],e[2+1,3+1],e[4+1,2+1],e[2+1,4+1],e[3+1,3+1],e[5+1,2+1],e[2+1,5+1],e[4+1,3+1],e[3+1,4+1],e[6+1,2+1],e[2+1,6+1],e[4+1,4+1]
  ), ncol=12, byrow=TRUE)

  EU = c(0,0,1,1,rho,e[3+1,0+1],e[0+1,3+1],e[2+1,1+1],e[1+1,2+1],e[4+1,0+1],e[0+1,4+1],e[2+1,2+1])
  Sigma = EUV - EU %*% t(EU)

  A = matrix( c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                -4*e[3+1,0+1], 0, -2*e[4+1,0+1], 0, 0, 0, 0, 0, 0, 1, 0, 0,
                0, -4*e[0+1,3+1], 0, -2*e[0+1,4+1], 0, 0, 0, 0, 0, 0, 1, 0,
                -2*e[1+1,2+1], -2*e[2+1,1+1], -e[2+1,2+1], -e[2+1,2+1], 0, 0, 0, 0, 0, 0, 0, 1), ncol=12, byrow=TRUE)

  B = matrix( c(-rho/2, -rho/2, 1, 0, 0, 0,
                0, 0, 0, -(e[2+1,2+1]-1)/(2*(e[4+1,0+1]-1)*((e[4+1,0+1]-1)*(e[0+1,4+1]-1))^(1/2)),
                -(e[2+1,2+1]-1)/(2*(e[0+1,4+1]-1)*((e[4+1,0+1]-1)*(e[0+1,4+1]-1))^(1/2)),
                1/((e[4+1,0+1]-1)*(e[0+1,4+1]-1))^(1/2) ), ncol=6, byrow=TRUE)

  Sigma = B %*% A %*% Sigma %*% t(A) %*% t(B)
  delta = 1e-6
  if (Sigma[1,1] < delta | Sigma[2,2] < delta) {
    Sigma[1,2] = 0.0
    if (Sigma[1,1] < delta) Sigma[1,1] = delta
    if (Sigma[2,2] < delta) Sigma[2,2] = delta
  }
  return(Sigma)
}
