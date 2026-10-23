## Internal func.
# ' Compute the tail probability of weighted sum of non-central chi-square random variables
# '
# ' Use Davies' method and apply a relative error bound (avoid redundant computations): efficient and accurate.
# ' @param Q.all  test statistics
# ' @param lambda  mixing coefficients
# ' @param h  DFs of chi-square dist
# ' @param delta  non-centrality parameters of chi-square dist
# ' @param acc  relative error bound
# ' @param lim  maximum number of integration terms
# ' @return tail probability
# ' @export
# ' @references
# ' Wu,B., Guan,W., Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. \emph{Annals of human genetics}, 80(2), 123-135.
# '
# ' Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2026) An Analytical Approach to Sample Size and Power Calculation for Equivalence Testing in Medical Device Comparison Studies. tech report.
pchisum <- function(Q.all, lambda, h = rep(1, length(lambda)), delta = rep(0, length(lambda)), acc=1e-2,lim=1e7){
  pval = pnchisum(Q.all,lambda,h,delta)
  i1 = which(is.finite(Q.all))
  for(i in i1){
    ## tmp = Wdavies(Q.all[i],lambda,h,delta,acc=acc*pval[i],lim=lim)
    tmp = Wdavies(Q.all[i],lambda,h,delta,acc=acc*min(pval[i],1-pval[i]),lim=lim)
    if((tmp$ifault<=0)&(tmp$pval[i]>=0)&(tmp$pval[i]<=1))  pval[i] = tmp$pval
  }
  return(pval)
}





## internal func.
# ' Compute the tail prob for weighted sum of 1-DF chi-square random variables
# '
# ' Imported from CompQuadForm R package (version 1.4.3; 10/4/2017). Slight simplifications for sum of central chi-squares.
# ' @param qx quantile to compute tail prob
# ' @param lambda weights for 1-DF chi-square rvs
# ' @param delta non-centrality parameters for 1-DF chi-square rvs. Default to zeros.
# ' @param lim total number of integration terms
# ' @param acc accuracy for computed probability
# ' 
# ' @return
# ' \describe{
# '  \item{pval}{ tail prob }
# '  \item{trace}{vector, indicating performance of procedure, with the following components: 1: absolute value sum, 2: total number of integration terms, 3: number of integrations, 4: integration interval in main integration, 5: truncation point in initial integration, 6: standard deviation of convergence factor term, 7: number of cycles to locate integration parameters}
# '  \item{ifault}{ fault indicator: 0: no error, 1: requested accuracy could not be obtained, 2: round-off error possibly significant, 3: invalid parameters, 4: unable to locate integration parameters }
# ' }
# ' @export
Wdavies <- function(q, lambda, h = rep(1, length(lambda)), delta = rep(0, length(lambda)), sigma = 0, lim = 1e6, acc = 1e-9){
  r = length(lambda)
  if (any(delta < 0)) stop("All non centrality parameters in 'delta' should be positive!") 
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  out <- .C("KATqfc", lambdas = as.double(lambda), noncentral = as.double(delta), df = as.integer(h), r = as.integer(r),
            sigma = as.double(sigma), q = as.double(q), lim = as.integer(lim), acc = as.double(acc),
            trace = as.double(rep(0, 7)), ifault = as.integer(0), res = as.double(0), PACKAGE = "RAMpb")
  return(list(trace = out$trace, ifault = out$ifault, pval = 1-out$res))
}


#' @useDynLib RAMpb
NULL
#> NULL

.onUnload <- function (libpath) {
  library.dynam.unload("RAMpb", libpath)
}



## sum of non-central chi_1^2
pnchisum <- function(q, lambda, h = rep(1, length(lambda)), delta = rep(0, length(lambda))){
  r = length(lambda)
  ## if (any(delta < 0)) stop("All non centrality parameters in 'delta' should be positive!")
  ## if (length(h) != r) stop("lambda and h should have the same length!")
  ## if (length(delta) != r) stop("lambda and delta should have the same length!")
  c1 <- sum(lambda * h) + sum(lambda * delta)
  c2 <- sum(lambda^2 * h) + 2 * sum(lambda^2 * delta)
  c3 <- sum(lambda^3 * h) + 3 * sum(lambda^3 * delta)
  c4 <- sum(lambda^4 * h) + 4 * sum(lambda^4 * delta)
  s1 = c3/(c2^(3/2));  s2 = c4/c2^2
  muQ = c1;  sigmaQ = sqrt(2 * c2)
  tstar = (q - muQ)/sigmaQ
  if (s1^2 > s2) {
    a = 1/(s1 - sqrt(s1^2 - s2))
    delta = s1 * a^3 - a^2
    l = a^2 - 2 * delta
  } else {
    a = 1/s1
    delta = 0
    l = c2^3/c3^2
  }
  muX = l + delta; sigmaX = sqrt(2)*a
  pval = pchisq(tstar * sigmaX + muX, df = l, ncp = delta, lower.tail = FALSE)
  return(pval)
}
