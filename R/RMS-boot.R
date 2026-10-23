#' Conduct estimation and RMS test for an one-way random-effects ANOVA model
#'
#' We are conducting hypothesis test for a composite parameter, the RMS, defined as \eqn{\sqrt{\mu^2+\sigma_b^2+\sigma_w^2}:=\rho }, where
#' \eqn{\mu} is the overall mean, and \eqn{(\sigma_b^2,\sigma_w^2)} are the between/within-subject variances in the
#'  one-way random-effects ANOVA model, \eqn{y_{ij}=\mu+u_i+\epsilon_{ij}}, where \eqn{u_i\sim N(0,\sigma_b^2) } and
#'  \eqn{\epsilon_{ij}\sim N(0,\sigma_w^2) }. We want to test \eqn{H_0: \rho\ge \rho_0}.
#' We implement a parametric Bootstrap based test with ``exact'' p-value calculation, voiding the need for Bootstrap Monte Carlo simulation.
#' See the reference of Bai et. al (2026).  
#'
#' @param Y vector of outcomes
#' @param subj subject id (factors). Observations with the sam id are coming from the same individual.
#' @param rho0  null threshold of acceptable RMS value.
#' @param REML using REML instead of MLE. Default to TRUE.
#' @return
#' \describe{
#'   \item{p.value}{ QMS test p-value }
#'   \item{pars}{ estimated null parameter values }
#' }
#' @export
#' @references
#' Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2026) An Analytical Approach to Sample Size and Power Calculation for Equivalence Testing in Medical Device Comparison Studies. tech report.
#' @examples
#' s2w=1.4^2; s2b=1.7^2; mu0=-0.4
#' ng = c(10,2,10,10,5,7,9,10)
#' A = rep(1:8, times=ng)
#' Y = mu0 + (rnorm(8)*sqrt(s2b))[A] + rnorm(sum(ng))*sqrt(s2w)
#' PBrms(Y,A)
PBrms <- function(Y, subj, rho0=3, REML=TRUE){
  ng = as.vector( table(subj) ); N = sum(ng); K = length(ng)
  mus = tapply(Y, subj, mean)
  sse = sum( tapply(Y, subj, function(yi) sum((yi-mean(yi))^2) ) )
  ## est
  lfn = function(xpar){
    r0 = xpar[1]; a1 = exp(xpar[2]); a2 = exp(xpar[3])
    b1 = 4*a1/(a1+1)^2; b2 = 1/(1+a2)
    mu = r0*(a1-1)/(a1+1); s2w = r0^2*b1*b2; s2b = r0^2*b1*(1-b2)
    ll1 = sum( log(s2w+ng*s2b) ) + (N-K)*log(s2w) + sse/s2w + sum( ng/(s2w+ng*s2b)*(mus-mu)^2 )
    ans = ll1 + REML*log( sum(ng/(s2w+ng*s2b)) )
    return(ans)
  }
  xpar = optim(c(rho0+0.001,0,0), fn=lfn, method="L-BFGS-B", lower=c(rho0,-20,-20), upper=c(20,20,20) )$par
  r0 = xpar[1]; a1 = exp(xpar[2]); a2 = exp(xpar[3])
  b1 = 4*a1/(a1+1)^2; b2 = 1/(1+a2)
  mu = r0*(a1-1)/(a1+1); s2w = r0^2*b1*b2; s2b = r0^2*b1*(1-b2)
  ## RMS test
  lam = h = dta = NULL
  for(i in 1:K){
    lam = c(lam, s2w+ng[i]*s2b, s2w);   h = c(h, 1,ng[i]-1)
    eta0 = ng[i]*mu^2/(s2w+ng[i]*s2b)
    dta = c(dta,eta0, 0)
  }
  Qy = sum(Y^2)
  pvalt = 1 - pchisum(Qy,lam,h,dta)
  return( list(p.value=pvalt, pars=c(s2w=s2w,s2b=s2b,mu=mu)) )
}





#' Large-sample Z-tests for the total absolute variation in a random-effects ANOVA model
#'
#' Given a random-effects ANOVA model, \eqn{Y_{ij}=\mu+u_i+\epsilon_{ij}}
#' with \eqn{u_i\sim N(0,\sigma_b^2), \epsilon_{ij}\sim N(0,\sigma_w^2)}, we
#' want to make inference on the total absolute variation parameter, \eqn{\rho^2=\mu^2+\sigma_b^2+\sigma_w^2}.
#' Large-sample Z-tests can be conducted based on approximating \eqn{\sum_{i,j}Y_{ij}^2} with a normal distribution.
#'
#' @param Y vector of outcomes
#' @param subj subject id (factors). Observations with the sam id are coming from the same individual.
#' @param rho0 the null threshold value to test \eqn{H_0: \rho\ge \rho_0}
#' @param REML whether we use REML or MLE, default to TRUE
#' @return
#' \describe{
#'   \item{p.value}{ Z-score and Z-Wald test p-values }
#'   \item{par0}{ estimated parameters, test statistic and its variance, under the null hypothesis }
#'   \item{pars}{ estimated parameters, test statistic and its variance }
#' }
#'
#' @export
#' @references
#' Ndikintum, N. K. and Rao, M. (2016) A Special Inference Problem in Repeated Measures Design Test of Statistical Hypothesis on Accuracy Root Mean Square Application to Pulse Oximetry Studies. Statistics in Biopharmaceutical Research 8, 60–76.
#'
#' Pennello, G. A. (2002) Sample size for paired repeated measures designs of method comparison studies: Application to pulse oximetry. ASA Proceedings of the Joint Statistical Meetings pages 2671–2675.
#'
#' Pennello, G. A. (2003) Comparing monitoring devices when a gold standard in unavailable: Application to pulse oximeters. ASA Proceedings of the Joint Statistical Meetings pages 3256–3263.
#' 
#' Bai,Y., Wang,Z., Lystig,T., and Wu,B. (2026) An Analytical Approach to Sample Size and Power Calculation for Equivalence Testing in Medical Device Comparison Studies. tech report
#' @examples
#' s2w=1.4^2; s2b=1.7^2; mu0=-0.4
#' ng = c(10,2,10,10,5,7,9,10)
#' A = rep(1:8, times=ng)
#' Y = mu0 + (rnorm(8)*sqrt(s2b))[A] + rnorm(sum(ng))*sqrt(s2w)
#' Zrms(Y,A, rho0=2.5)
Zrms <- function(Y,subj, rho0=3, REML=TRUE){
  ng = as.vector( table(subj) ); N = sum(ng); K = length(ng)
  mus = tapply(Y, subj, mean)
  sse = sum( tapply(Y, subj, function(yi) sum((yi-mean(yi))^2) ) )
  rho2 = rho0^2
  Qy = sum(Y^2) ## sse + sum(ng*mus^2)
  ## score Z
  lfn = function(xpar){
    r0 = xpar[1]; a1 = exp(xpar[2]); a2 = exp(xpar[3])
    b1 = 4*a1/(a1+1)^2; b2 = 1/(1+a2)
    mu = r0*(a1-1)/(a1+1); s2w = r0^2*b1*b2; s2b = r0^2*b1*(1-b2)
    ll1 = sum( log(s2w+ng*s2b) ) + (N-K)*log(s2w) + sse/s2w + sum( ng/(s2w+ng*s2b)*(mus-mu)^2 )
    ans = ll1 + REML*log( sum(ng/(s2w+ng*s2b)) )
    return(ans)
  }
  ## xpar = optim(c(rho0+0.001,0,0), fn=lfn, method="L-BFGS-B", lower=c(rho0,-Inf,-Inf) )$par
  xpar = optim(c(rho0+0.001,0,0), fn=lfn, method="L-BFGS-B", lower=c(rho0,-20,-20), upper=c(20,20,20) )$par
  r0 = xpar[1]; a1 = exp(xpar[2]); a2 = exp(xpar[3])
  b1 = 4*a1/(a1+1)^2; b2 = 1/(1+a2)
  mu = r0*(a1-1)/(a1+1); s2w = r0^2*b1*b2; s2b = r0^2*b1*(1-b2)
  tau2 = ( 2*(s2w+ng*s2b)^2+2*(ng-1)*s2w^2 + 4*ng*(s2w+ng*s2b)*mu^2 )/ng^2
  V0 = sum(tau2*ng^2/N^2)  
  Zs = (Qy/N-rho2)/sqrt(V0)
  pvals = pnorm(Zs)
  ## Wald Z
  lfn = function(xpar){
    s2b = exp(xpar[1]); s2w = exp(xpar[2])
    aa1 = ng/(s2w+ng*s2b); aas = sum(aa1)
    mu = sum(aa1*mus)/aas
    ll1 = sum( log(s2w+ng*s2b) ) + (N-K)*xpar[2] + sse/s2w + sum( ng/(s2w+ng*s2b)*(mus-mu)^2 )
    ll1 + REML*log(aas)
  }
  xpar = optim(c(0,0), lfn, method='BFGS')$par
  s2b1 = exp(xpar[1]); s2w1 = exp(xpar[2])
  aa1 = ng/(s2w1+ng*s2b1); mu1 = sum(aa1*mus)/sum(aa1)
  tau2 = ( 2*(s2w1+ng*s2b1)^2+2*(ng-1)*s2w1^2 + 4*ng*(s2w1+ng*s2b1)*mu1^2 )/ng^2
  V1 = sum(tau2*ng^2/N^2)
  Zw = (Qy/N-rho2)/sqrt(V1)
  pvalw = pnorm(Zw)
  ##
  pval = c(pvals,pvalw)
  names(pval) = c('Z-score', 'Z-Wald')
  return( list(p.value=pval, pars0=c(s2w=s2w,s2b=s2b,mu=mu,R=Qy/N,V=V0),
               pars=c(s2w=s2w1,s2b=s2b1,mu=mu1,R=Qy/N,V=V1)) )
}
