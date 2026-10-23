## (OLD) Internal function to compute the log likelihood for paired data Dij
## @param xpar vector of log(s2w), log(s2b), mu0
## @param Dij  n by m matrix of paired measures
## @param rho  null value of sqrt(s2w+s2b+mu0^2)
llkfun <- function(xpar, Dij, REML=TRUE){
   n = dim(Dij)[1]; m = dim(Dij)[2]
   s2w = exp(xpar[1])
   s2b = exp(xpar[2])
   mu0 = xpar[3]
   ## 
   D1 = rowMeans(Dij)
   E = Dij-D1
   s2s = m*s2b+s2w
   llk0 = log(s2s) + (m-1)*xpar[1]
   llk1 = sum(E^2)/s2w + m*sum((D1-mu0)^2)/s2s
   allk = 0
   if(REML) allk = -log(s2s)
   n*llk0+llk1+allk
}
## (OLD) Internal constraint function
## @param xpar vector of log(s2w), log(s2b), mu0
## @param rho  null value of sqrt(s2w+s2b+mu0^2)
Cfun <- function(xpar, rho=3){
   s2w = exp(xpar[1])
   s2b = exp(xpar[2])
   mu0 = xpar[3]
   s2w+s2b+mu0^2-rho^2
}

## Internal P-value calc func.
Tpval0 <- function(xq, n, m, s2w, s2b, mu0){
  lam = c(s2w+m*s2b, s2w);   h = c(n,n*(m-1))
  eta0 = n*m*mu0^2/lam[1]
  dta = c(eta0, 0)
  1-pchisum(xq,lam,h,dta)
} 


#' Conduct estimation and RMS test for a paired repeated measures (PRM) design study
#'
#' The problem can be reduced to a mixed effects model estimation under a quadratic constraint. Due to the special structure of the
#' mixed effects model for the PRM designs, we can compute the likelihood analytically and hence efficiently calculate
#' MLE/REML. See the reference of Bai et. al (2026). The score and Wald Z-tests are also implemented.
#' This function works for balanced designs, and is mainly used for the purpose of sample size/power calculation. 
#'
#' @param X n by m matrix of paired measures
#' @param rho  null threshold of acceptable RMS value
#' @param REML using REML instead of MLE. Default to TRUE.
#' @return
#' \describe{
#'   \item{p.value}{ test p-values for: QMS test, score Z-test, Wald Z-test }
#'   \item{pars}{ estimated null parameter values }
#' }
#' @export
#' @references
#' Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2026) An Analytical Approach to Sample Size and Power Calculation for Equivalence Testing in Medical Device Comparison Studies. tech report
#' @examples
#' n=20; m=20; s2w=4.44; s2b=1.2; mu0=0.94
#' e = matrix(rnorm(n*m), n,m)*sqrt(s2w)
#' u = rnorm(n)*sqrt(s2b)
#' X = mu0 + u + e
#' RMStest(X,3)
#' RMStest(X, sqrt(s2w+s2b+mu0^2))
RMStest <- function(X, rho, REML=TRUE){
  n = dim(X)[1]; m = dim(X)[2]
  Y = as.vector(X); subj = rep(1:n, m)
  aaa = PBrms(Y, subj, rho, REML)
  bbb = Zrms(Y,subj, rho, REML)
  pval = c(aaa$p.val, bbb$p.val)
  names(pval) = c('QMS', 'Z-score', 'Z-Wald')
  return(list(p.value=pval, pars=aaa$pars) )
}

#' Analytical power calculation for RMS test problem
#'
#' Due to its nature of composite null hypothesis, it is not trivial to derive the analytical power calculation. We develop an average log likelihood
#' based approach to quickly compute the power.
#' 
#' @param alpha desired significance level. Default to 0.05
#' @param n  number of subjects
#' @param m  number of repeated measures for each subject
#' @param s2w  within subject variation
#' @param s2b  between subject variation
#' @param mu   mean measure difference
#' @param rho   null threshold of acceptable RMS value
#' @param REML using REML instead of MLE. Default to TRUE.
#' @return
#' \describe{
#'   \item{xq}{ computed quantile for the QMS test statistic }
#'   \item{size}{ computed actual type I error }
#'   \item{pwr}{ computed power }
#'   \item{par0}{ estimated null parameter values }
#' }
#' @export
#' @references
#' Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2026) An Analytical Approach to Sample Size and Power Calculation for Equivalence Testing in Medical Device Comparison Studies. tech report
RMSap <- function(alpha=0.05, n, m,s2w, s2b, mu, rho, REML=TRUE){
  ## est
  lfn = function(xpar){
    r0 = xpar[1]; a1 = exp(xpar[2]); a2 = exp(xpar[3])
    b1 = 4*a1/(a1+1)^2; b2 = 1/(1+a2)
    eu0 = r0*(a1-1)/(a1+1); e2w = r0^2*b1*b2; e2b = r0^2*b1*(1-b2)
    e2s = m*e2b+e2w
    llk0 = log(e2s) + (m-1)*log(e2w)
    llk1 = (m-1)*s2w/e2w + (m*(eu0-mu)^2 + m*s2b+s2w)/e2s
    allk = 0
    if(REML) allk = -log(e2s)/n
    llk0+llk1+allk
  }
  ## est
  xpar = optim(c(rho+0.001,0,0), fn=lfn, method="L-BFGS-B", lower=c(rho,-20,-20), upper=c(20,20,20) )$par
  r0 = xpar[1]; a1 = exp(xpar[2]); a2 = exp(xpar[3])
  b1 = 4*a1/(a1+1)^2; b2 = 1/(1+a2)
  eu0 = r0*(a1-1)/(a1+1); e2w = r0^2*b1*b2; e2b = r0^2*b1*(1-b2)
  ##  
  A0 = round(1/alpha)
  f0 = function(x) ( Tpval0(x, n,m,e2w,e2b,eu0)/alpha - 1 )*A0
  xq = uniroot(f0, c(1,1e5))$root
  size = Tpval0(xq, n,m,e2w,e2b,eu0)
  pwr1 = Tpval0(xq, n,m,s2w,s2b,mu)
  return( list(xq=xq,size=size,pwr=pwr1, par0=c(s2w=e2w,s2b=e2b,mu=eu0)) )
}


#' Monte Carlo power calculation for RMS test problem
#'
#' Due to its nature of composite null hypothesis, it is not trivial to derive the analytical power calculation. 
#' Monte Carlo simulation offers a straightforward approach by generating p-value random numbers.
#' 
#' @param alpha desired significance level. Default to 0.05
#' @param n  number of subjects
#' @param m  number of repeated measures for each subject
#' @param s2w  within subject variation
#' @param s2b  between subject variation
#' @param mu   mean measure difference
#' @param rho   null threshold of acceptable RMS value
#' @param REML using REML instead of MLE. Default to TRUE.
#' @param B   total number of Monte Carlo simulations. Default to 1000.
#' @return
#' \describe{
#'   \item{pwr}{ computed power }
#'   \item{p.value}{ p-values from Monte Carlo simulations }
#' }
#' @export
#' @references
#' Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2026) An Analytical Approach to Sample Size and Power Calculation for Equivalence Testing in Medical Device Comparison Studies. tech report
RMSmcp <- function(alpha=0.05,n, m, s2w, s2b, mu, rho, REML=TRUE, B=1e3){
  pval = rep(1, B)
  for(it in 1:B){
    e = matrix(rnorm(n*m), n,m)*sqrt(s2w)
    u = rnorm(n)*sqrt(s2b)
    X = mu + u + e
    pval[it] = RMStest(X, rho, REML)$p.value[1]
  }
  pwr = mean(pval<=alpha)
  return( list( pwr=pwr, p.value=pval) )
}


#' Analytical sample size calculation for RMS test problem
#'
#' Due to discrete nature of sample size, the power is not a continuous function. This function utilizes the analytical power calculation function RMSap()
#' to compute an approximate sample size n under a given m (number of measures per unit), which has power closest to the desired power.
#' 
#' @param pwr desired power. Default to 0.8
#' @param alpha desired significance level. Default to 0.05
#' @param m  number of repeated measures for each subject
#' @param s2w  within subject variation
#' @param s2b  between subject variation
#' @param mu   mean measure difference
#' @param rho   null threshold of acceptable RMS value
#' @param REML using REML instead of MLE. Default to TRUE.
#' @return
#' \describe{
#'   \item{n}{ computed number of units }
#'   \item{pwr}{ actual power under the computed sample size }
#' }
#' @export
#' @references
#' Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2026) An Analytical Approach to Sample Size and Power Calculation for Equivalence Testing in Medical Device Comparison Studies. tech report
RMSas <- function(pwr=0.8,alpha=0.05, m,s2w, s2b, mu, rho, REML=TRUE){
  ## est
  fn = function(n) RMSap(alpha, round(n), m,s2w, s2b, mu, rho, REML)$pwr - pwr
  a1 = 5; while(fn(a1)>0) a1 = a1/2
  a2 = 30; while(fn(a2)<0) a2 = a2*2
  n = round( uniroot(fn, c(a1,a2))$root )
  pwr1 = RMSap(alpha, n, m,s2w, s2b, mu, rho, REML)$pwr
  return( list(n=n, pwr=pwr1) )
}


