# --------------------------------------------------------------------------
# power (or alpha) of 2-stage studies according to Potvin et. al. 
# method "B" modified to include a futility criterion for PE or CI and
# modified to use PE of stage 1 in sample size estimation
#
# author D.L.
# --------------------------------------------------------------------------
# require(PowerTOST)
# source("C:/Users/dlabes/workspace/Power2stage/R/sampsiz.R")
# source("C:/Users/dlabes/workspace/Power2stage/R/power.R")


power.2stage.Bf <- function(alpha=c(0.0294,0.0294), n1, CV, GMR, targetpower=0.8, 
                            pmethod=c("nct","exact"), usePE=FALSE, powerstep=TRUE, 
                            fPElower, fPEupper, theta0, theta1, theta2, 
                            npct=c(0.05, 0.5, 0.95), 
                            nsims=1e5, setseed=TRUE, print=TRUE, details=TRUE)
{
  message("This function is depreciated! Use power.2stage.fC instead.")
  
  if (missing(CV)) stop("CV must be given!")
  if (CV<=0)       stop("CV must be >0!")
  
  if (missing(n1)) stop("Number of subjects in stage 1 (n1) must be given!")
  if (n1<=2)       stop("Number of subjects in stage 1 (n1) must be >2!")
  
  if (missing(GMR)) GMR <- 0.95
  
  if (missing(theta1)  & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2))  theta2 <- 1/theta1
  if (missing(theta1)  & !missing(theta2)) theta1 <- 1/theta2
  
  if (GMR<=theta1 | GMR>=theta2) stop("GMR must be within acceptance range!")
  
  if (missing(theta0)) theta0 <- GMR
  
  if (missing(fPElower) & missing(fPEupper))  fPElower <- 0.825
  if (missing(fPElower) & !missing(fPEupper)) fPElower <- 1/fPEupper
  if (!missing(fPElower) & missing(fPEupper)) fPEupper <- 1/fPElower
  
  res <- power.2stage.fC(method="B", alpha=alpha, n1=n1, CV=CV, 
                         GMR=GMR, targetpower=targetpower, pmethod=pmethod, 
                         usePE=usePE, powerstep=powerstep, fCrit="PE", 
                         fClower=fPElower, fCupper=fPEupper, theta0=theta0, 
                         theta1=theta1, theta2=theta2, npct=npct, nsims=nsims, 
                         setseed=setseed, print=print, details=details)
  
  if (print) return(invisible(res)) else return(res)
  
} #end function
