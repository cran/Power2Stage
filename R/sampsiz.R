# ----- helper function ----------------------------------------------------
# Sample size for a desired power, large sample approx.
# 
# author D. Labes
# bk = design constant, see known.designs()
.sampleN0 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, 
                          se, steps=2, bk=4, diffmthreshold=0.04)
{
  
  z1 <- qnorm(1-alpha)
  # value diffmthreshold=0.04 corresponds roughly to log(0.96)
  # with lower values there are many steps around between 0.95 and 1
  # in sampleN.TOST
  if (abs(diffm)>diffmthreshold) z2 <- qnorm(targetpower) else {
    z2 <- qnorm(1-(1-targetpower)/2) # #1-beta/2 for diffm ~0 (log(theta0=1))
    diffm <- 0
  }
  n01<-(bk/2)*((z1+z2)*(se*sqrt(2)/(diffm-ltheta1)))^2;
  n02<-(bk/2)*((z1+z2)*(se*sqrt(2)/(diffm-ltheta2)))^2;
  n0 <- ceiling(max(n01,n02))
  
  if (n0>2){
  # make another step with t-distri
    z1 <- qt(1-alpha, df=n0-2)
    if (abs(diffm)>diffmthreshold) z2 <- qt(targetpower, df=n0-2) else {
      z2 <- qt(1-(1-targetpower)/2, df=n0-2)  
      diffm <- 0
    }
    n01<-(bk/2)*((z1+z2)*(se*sqrt(2)/(diffm-ltheta1)))^2;
    n02<-(bk/2)*((z1+z2)*(se*sqrt(2)/(diffm-ltheta2)))^2;
    n0 <- ceiling(max(n01,n02))
  }
  # make an even multiple of step (=2 in case of 2x2 cross-over)
  n0 <- steps*trunc(n0/steps)
  
  # minimum sample size will be checked outside
  return(n0)
}

# -------------------------------------------------------------------------
# sample size function without all the overhead
# for 2x2 crossover or 2-group parallel
# power via nct approximation or exact
#
# author D. Labes
# -------------------------------------------------------------------------

.sampleN <- function(alpha=0.05, targetpower=0.8, ltheta0, 
                     ltheta1=log(0.8), ltheta2=log(1.25), mse, bk=2,
                     method=c("nct","exact"))
{
  # return 'Inf' if ltheta0 not between or very near to ltheta1, ltheta2
  if ((ltheta0-ltheta1)<1.25e-5 | (ltheta2-ltheta0)<1.25e-5) {
    # debug
    # cat ("Inf returned for",exp(ltheta0),"\n")
    return(Inf)
  }

  # design characteristics for 2-group parallel and 2x2 crossover design
  # df for the design as an unevaluated expression
  dfe   <- parse(text="n-2", srcfile=NULL)
  steps <- 2     # stepsize for sample size search
  nmin  <- 4     # minimum n
  
  # log transformation assumed
  se     <- sqrt(mse)
  diffm  <- ltheta0
  
  # start value from large sample approx. (hidden func.)
  n  <- .sampleN0(alpha, targetpower, ltheta1, ltheta2, diffm, se, steps, bk)
  if (n<nmin) n <- nmin
  df <- eval(dfe)
  pow <- .calc.power(alpha, ltheta1, ltheta2, diffm, se, n, df, bk, method)
  
  iter <- 0; imax <- 50
  # iter>50 is emergency brake
  # this is eventually not necessary, depends on quality of sampleN0
  # in experimentation I have seen max of six steps
  # reformulation with only one loop does not shorten the code considerable
  # --- loop until power <= target power, step-down
  down <- FALSE; up <- FALSE
  while (pow>targetpower) {
    if (n<=nmin) { 
      break
    }
    down <- TRUE
    n    <- n-steps     # step down if start power is to high
    iter <- iter+1
    df   <- eval(dfe)
    pow  <- .calc.power(alpha, ltheta1, ltheta2, diffm, se, n, df, bk, method)
    
    if (iter>imax) break  
    # loop results in n with power too low
    # must step one up again. is done in the next loop
  }
  # --- loop until power >= target power
  while (pow<targetpower) {
    up   <- TRUE; down <- FALSE
    n    <- n+steps
    iter <- iter+1
    df   <- eval(dfe)
    pow  <- .calc.power(alpha, ltheta1, ltheta2, diffm, se, n, df, bk, method)
    if (iter>imax) break 
  }
  nlast <- n
  if ((up & pow<targetpower) | (down & pow>targetpower) ) {
    n <- NA
  }
  
  #return results as data.frame
#   res <- data.frame(design="2x2", alpha=alpha, CV=mse2CV(mse), 
#                     theta0=exp(ltheta0), theta1=exp(ltheta1), 
#                     theta2=exp(ltheta2), n=n, power=pow, 
#                     targetpower=targetpower)
#   names(res) <-c("Design","alpha","CV","theta0","theta1","theta2",
#                  "Sample size", "Achieved power", "Target power")
#   
#   return(res)

  #  return only n
  return(n)
  
} # end of function
