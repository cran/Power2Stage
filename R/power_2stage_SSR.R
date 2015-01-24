# helper function to get CV from 'required sample size'
reqN2CV <- function(alpha=0.05, targetpower=0.8, theta0=1., theta1=0.8, n)
{
  z1 <- qnorm(1-alpha)
  if (log(theta0)!=0) z2 <- qnorm(targetpower) else z2 <- qnorm(1-(1-targetpower)/2) 

  s2e <- n*(log(theta0)-log(theta1))^2/2/(z1+z2)^2
  return(mse2CV(s2e))
}
  
# --------------------------------------------------------------------------
# power (or alpha) of 2-stage studies with (blinded) sample size re-estimation
# see Golkowski et al.
#
# author D.L.
# --------------------------------------------------------------------------
# require(PowerTOST)
# source("C:/Users/dlabes/workspace/Power2Stage/R/sampsiz.R")
# source("C:/Users/dlabes/workspace/Power2Stage/R/power.R")

power.2stage.ssr <- function(alpha=0.05, n1, GMR, CV, targetpower=0.8, 
                             pmethod=c("nct","exact", "shifted","ls"),
                             blind=FALSE, min.n=0, max.n=Inf, 
                             theta0, theta1, theta2, npct=c(0.05, 0.5, 0.95), 
                             nsims=1e5, setseed=TRUE, print=TRUE, details=TRUE)
{ # seems to give errors if alpha is a vector
  alpha <- alpha[1L]
  
  if (missing(CV)) stop("CV must be given!")
  if (CV<=0)       stop("CV must be >0!")
  
  if (missing(n1)) stop("Number of subjects in stage 1 must be given!")
  if (n1<=0)       stop("Number of subjects in stage 1 must be >0!")
  
  if (missing(GMR)) GMR <- 0.95
  
  if (missing(theta1) & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) & !missing(theta2)) theta1 <- 1/theta2
  
  if (GMR<=theta1 | GMR>=theta2) stop("GMR must be within acceptance range!")
  
  if (missing(theta0)) theta0 <- GMR
  
  if (n1>max.n) stop("max.n<n1 doestn't make sense!")
  if (min.n!=0 & min.n<n1) stop("min.n<n1 doestn't make sense!")
  
  # check if power calculation method is nct or exact
  pmethod <- match.arg(pmethod)
  
  if(print & details){
    cat(nsims,"sims. Stage 1")
  }
  # start timer
  ptm  <- proc.time()
  
  if (setseed) set.seed(1234567)

  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  lGMR    <- log(GMR)
  mlog    <- log(theta0)
  mse     <- CV2mse(CV)
  bk      <- 2   # 2x2x2 crossover design const
  # reserve memory
  BE      <- rep.int(NA, times=nsims)
  
# ----- interim ----------------------------------------------------------
# simulate blinded est. of variance
  Cfact <- bk/n1
  df    <- n1-2
  tval  <- qt(1-alpha[1], df)
  sdm   <- sqrt(mse*Cfact)
  # simulate point est. via normal distribution
  pes   <- rnorm(n=nsims, mean=mlog, sd=sdm)
  # simulate mse via chi-squared distribution
  mses  <- mse*rchisq(n=nsims, df=df)/df
   
  #s2os  <- mses + pes^2/2   # is this correct? example shows: it's not correct.
                             # s2os is smaller than mse! I would expect 
                             # the other way round!
  # SS of tmt formula is only valid if sequence balanced
  #        res. SS + SS of tmt 
  s2os <- (df*mses + n1*pes^2/2)/(n1-1)
  # unblind
  if (!blind) s2os <- mses

  if(print & details){
    # time for stage 1 sims
    cat(" - Time consumed (secs):\n")
    print(round((proc.time()-ptm),2))
    cat("Keep calm. ")
  }

  # ------ recalculate sample size -----------------------------------------
  ptms <- proc.time()
  # first calc. power to avoid unnecessary time consuming sample size estimation
  # this is not really part of the method (?) but a run-time goody
  # will give some boost for small CV's and/or high n1
  if(pmethod!="ls"){
    pwr <- .calc.power(alpha=alpha, ltheta1=ltheta1, ltheta2=ltheta2, 
                       diffm=lGMR, sem=sqrt(bk*s2os/n1), df=df, method=pmethod)
    if(print & details){
      cat("Sample sizes (", sum(pwr<targetpower),
          " studies) will be re-estimated.\n", sep="")
      cat("May need some time.\n")
    }
  } else {
    # don't have power for large sample approx.
    # moreover ssr is so fast that we don't need the power step
    pwr  <- rep(0, times=nsims)
    if(print & details){
      cat("Sample sizes will be re-estimated.\n")
      cat("May need some time.\n")
    }
  }
  # total sample size
  ntot <- rep(n1, times=nsims)
  mse_tmp <- s2os[pwr<targetpower]
  # large sample approx. via normal distribution
  if(pmethod=="ls"){
    nt <- .sampleN00(alpha=alpha, targetpower=targetpower, se=sqrt(mse_tmp), 
                     diffm=lGMR, ltheta1=ltheta1, ltheta2=ltheta2, bk=2, 
                     steps=2, diffmthreshold=0.0)
  } else {
#     nt <- mapply(FUN=.sampleN, mse=mse_tmp, 
#                  MoreArgs=list(alpha=alpha, targetpower=targetpower, 
#                                ltheta0=lGMR, ltheta1=ltheta1, ltheta2=ltheta2,
#                                method=pmethod, bk=2))
     nt <- .sampleN2(alpha=alpha, targetpower=targetpower, ltheta0=lGMR,
                     mse=mse_tmp, ltheta1=ltheta1, ltheta2=ltheta2, 
                     method=pmethod, bk=2)
  }
  # maybe we have enough power in all cases and thus no re-estimated sample size
  if(length(nt)>0) ntot[pwr<targetpower] <- nt
  # take care of memory
  rm(nt, pwr, mse_tmp)
  # sample size returns Inf if pe outside acceptance range, then stay with n1 
  # but this should not occure here since pe1 is not used
  ntot <- ifelse(is.finite(ntot), ntot, n1)
  # use max.n if nt > max.n
  ntot <- ifelse(ntot>max.n, max.n, ntot)
  # use min.n if nt < min.n
  ntot <- ifelse(ntot<min.n, min.n, ntot)
  
  n2   <- ntot-n1
  # do not fall below n1
  n2 <- ifelse(n2<0, 0, n2)

  if(print & details){
    cat("Time consumed (secs):\n")
    print(round((proc.time()-ptms),1))
  }

  # ---------- stage 2 evaluation --------------------------------------
  m1    <- pes
  SS1   <- (n1-2)*mses
  # to avoid warnings for n2=0 in rnorm() and rchisq()
  ow    <- options("warn")
  options(warn=-1)
  m2    <- ifelse(n2>0, rnorm(n=nsims, mean=mlog, sd=sqrt(mse*bk/n2)), 0)
  # ??? (n2-2) cancels out! 
  SS2   <- ifelse(n2>2, (n2-2)*mse*rchisq(n=nsims, df=n2-2)/(n2-2), 0)
  # reset options
  options(ow)

  ntot <- n1+n2

  # do we need this 'stage' contribution?
  # can't discover such in Golkowski et al., but ... 
  withStage <- TRUE
  if(withStage){
    SSmean <- ifelse(n2>0, (m1-m2)^2/(2/n1+2/n2), 0)
    df2    <- ifelse(n2>0, ntot-3, n1-2)
  } else {
    SSmean <- 0
    df2    <- ntot-2
  }
  pe2    <- ifelse(n2>0, (n1*m1+n2*m2)/ntot, pes)
  mse2   <- ifelse(n2>0, (SS1+SSmean+SS2)/df2, mses)
  # take care of memory
  rm(m1, m2, SS1, SS2, SSmean)
  # calculate CI and decide BE 
  hw    <- qt(1-alpha, df2)*sqrt(mse2*bk/ntot)
  lower <- pe2 - hw
  upper <- pe2 + hw
  BE    <- (lower>=ltheta1) & (upper<=ltheta2)
  # take care of memory
  rm(lower, upper, hw)

  # the return list
  res <- list(alpha=alpha, CV=CV, n1=n1, GMR=GMR, targetpower=targetpower, 
              pmethod=pmethod, theta0=theta0, theta1=theta1, theta2=theta2, 
              max.n=max.n, min.n=min.n, nsims=nsims,
              # results 
              pBE=sum(BE)/nsims, 
              #pBE_s1 ?
              pct_s2=100*length(ntot[ntot>n1])/nsims, 
              nmean=mean(ntot), nrange=range(ntot), nperc=quantile(ntot, p=npct))
  # output
  if (print) {
    if (details){
      cat("Total time consumed (secs):\n")
      print(round((proc.time()-ptm),1))
      cat("\n")
    }
    if(blind) blinded <- "blinded " else blinded <- ""
    cat("2-stage design with ", blinded, "sample size re-estimation\n", sep="")
    cat("Nominal alpha =", alpha, "\n")
    cat("Sample size based on power calculated via", pmethod, "method\n")
    cat("with GMR =", GMR,"and targetpower =", targetpower,"\n")
    if(is.finite(max.n)){
      cat("Maximum sample size max.n = ",max.n,"\n", sep="")
    }
    cat("BE margins = ", theta1," ... ", theta2,"\n", sep="")
    cat("CV= ", CV,"; n(stage 1) = ", n1,"\n", sep="")
    
    cat("\n",nsims," sims at theta0 = ", theta0, sep="")
    if(theta0<=theta1 | theta0>=theta2) cat(" (p(BE)='alpha').\n") else { 
       cat(" (p(BE)='power').\n")}
    cat("p(BE)   = ", res$pBE,"\n", sep="")
    cat("Studies in stage 2 = ", round(res$pct_s2,2), "%\n", sep="")
    cat("\nDistribution of n(total)\n")
    cat("- mean (range) = ", round(res$nmean,1)," (", res$nrange[1]," ... ",
        res$nrange[2],")\n", sep="")
    cat("- percentiles\n")
    print(res$nperc)
    cat("\n")
  } 
  
  if (print) return(invisible(res)) else return(res)
  
} #end function
