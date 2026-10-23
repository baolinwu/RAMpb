# RAMpb
 - An R package for Statistical Test and Sample Size/Power Calculation for Equivalency Assessment Using Absolute Deviation 
  from Repeated Measures (as commonly applied in medical comparison studies)
 - Reference
    - Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2026) An Analytical Approach to Sample Size and Power Calculation for Equivalence Testing in Medical Device Comparison Studies. *tech rep*.
 - Sample R codes
```R
 ## install the package
 library(RAMpb)
 ## simulate RAM data
 n = 20;  m = 10
 A = rep(1:m, each=n)
 s2w = s2b = 4;  mu = 1
 Y = mu + rnorm(n*m)*sqrt(s2w) + rnorm(n)[A]*sqrt(s2b)
 ## QMS test
 rho = sqrt(s2w+s2b+mu^2)
 PBrms(Y,A, rho)$p.val
 PBrms(Y,A, rho+0.75)$p.val
 Zrms(Y,A, rho0=rho)$p.val
 Zrms(Y,A, rho0=rho+0.75)$p.val
 ## power calculation
 ## analytical calc
 RMSap(alpha=0.05,n, m, s2w, s2b, mu, rho+0.75)$pwr
 ## Monte Carlo sim
 RMSmcp(alpha=0.05,n, m, s2w, s2b, mu, rho+0.75)$pwr
 ## sample size calc
 aa = RMSas(pwr=0.8,alpha=0.05, m,s2w, s2b, mu, rho+0.75)
 aa
 ## RMSap(alpha=0.05,n=aa$n, m, s2w, s2b, mu, rho+0.75)$pwr
 RMSmcp(alpha=0.05,n=aa$n, m, s2w, s2b, mu, rho+0.75)$pwr

 ## analyze the oxymetry comparison data
 data(poData)
 PBrms(poData$Y,poData$A,rho0=3)
 Zrms(poData$Y,poData$A,rho0=3)

```
