#------------------------------------------------------------------------------
# Author: dlabes
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# helper function to allow partial match of the method
.powerMethod <- function(method){
  meth <- tolower(method[1])
  if (method=="") method <- "exact"
  methods <- c("exact","owenq","noncentral","nct","shifted","central")
  #                         ^ = match at start
  meth <- methods[grep(paste("^",meth,sep=""), methods)]
  if (length(meth)==0) meth <- tolower(method)
  if (length(meth)>1)  meth <- "nct" # happens if only "n" is given
  
  return(meth)
}
#------------------------------------------------------------------------------
# 'raw' power function without any error checks,
# does not vectorize propperly!
# to be used in sampleN.TOST avoiding overhead of redundant calculations
# in case of multiplicative model:
# diffm=log(null ratio), theta1=log(lower BE limit), theta2=log(upper BE limit)
# in case of additive model:
# diffm=1-null ratio, theta1=lower BE limit-1, theta2=upper BE limit -1
.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, se, n, df, bk=2)
{
  tval   <- qt(1 - alpha, df, lower.tail = TRUE)
  # if alpha>0.5 (very unusual) then R is negative 
  # in the application of OwensQ the upper integration limit 
  # is lower then the lower integration limit!
  # SAS OwenQ gives missings if b or a are negative!
  
  # 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2
  # and se=0!
  delta1 <- (diffm-ltheta1)/(se*sqrt(bk/n))
  delta2 <- (diffm-ltheta2)/(se*sqrt(bk/n))
  # is this correct?
  delta1[is.nan(delta1)] <- 0
  delta2[is.nan(delta2)] <- 0
  # R is infinite in case of alpha=0.5
  R <- (delta1-delta2)*sqrt(df)/(2.*abs(tval))
  # in case of se=0 it results: delta1=Inf, delta2=inf if diffm>ltheta2
  # Inf - Inf is NaN
  R[is.nan(R)] <- 0
  
  # to avoid numerical errors in OwensQ implementation
  if (min(df)>10000) { 
    # Joulious formula (57) or (67), large sample normal approximation
    p1 <- pnorm( (abs(delta1)-tval), lower.tail = TRUE)
    p2 <- pnorm( (abs(delta2)-tval), lower.tail = TRUE)
    return(p1 + p2 - 1.)
  }
  if (min(df)>=5000 & min(df)<=10000) {
    # approximation via non-central t-distribution
    return(.approx.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk))
  }
  
  # attempt to vectorize (it vectorizes properly if diffm is a vector
  # OR se OR n,df are vectors) 
  nel <- length(delta1)
  dl <- length(tval)
  p1 <- c(1:nel)	
  p2 <- p1
  for (i in seq_along(delta1)) {
    if (dl>1) {
      ddf <- df[i]; ttt <- tval[i]
    } else {
      ddf <- df[1]; ttt <- tval[1]
    }
    p1[i] <- OwensQ(ddf,  ttt, delta1[i], 0, R[i])
    p2[i] <- OwensQ(ddf, -ttt, delta2[i], 0, R[i])
  }
  pow <- p2-p1
  # due to numeric inaccuracies
  pow[pow<0] <- 0
  return( pow )
}
#------------------------------------------------------------------------------
# 'raw' approximate power function without any error checks, 
# approximation based on non-central t
# this vectorizes ok
.approx.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, 
		                           se, n, df, bk=2)
{
  tval <- qt(1 - alpha, df, lower.tail = TRUE, log.p = FALSE)
  
  # 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2
  # and se=0!
  delta1 <- (diffm-ltheta1)/(se*sqrt(bk/n))
  delta2 <- (diffm-ltheta2)/(se*sqrt(bk/n))
  # is this correct?
  delta1[is.nan(delta1)] <- 0
  delta2[is.nan(delta2)] <- 0
  
  pow <- pt(-tval, df, ncp=delta2)-pt(tval, df, ncp=delta1)
  pow[pow<0] <- 0 # this is to avoid neg. power due to approx. (vector form)
  
  return(pow)
}
#------------------------------------------------------------------------------
# 'raw' power function without any error checks, 
# approximation based on central 'shifted' t distribution
# according to Chow, Liu "Design and Analysis of Bioavailability ..."
# Chapter 9.6 and implemented in PASS 2008
# where does this all came from?
.approx2.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, 
                                se, n, df, bk=2)
{
	tval   <- qt(1 - alpha, df, lower.tail = TRUE)
	# 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2
	# and se=0!
	delta1 <- (diffm-ltheta1)/(se*sqrt(bk/n))
	delta2 <- (diffm-ltheta2)/(se*sqrt(bk/n))
	# is this correct?
	delta1[is.nan(delta1)] <- 0
	delta2[is.nan(delta2)] <- 0
	
	pow <- pt(-tval-delta2,df) - pt(tval-delta1,df)
	pow[pow<0] <- 0 # this is to avoid neg. power due to approx. (vector form)
	
	return(pow)
}
#------------------------------------------------------------------------------
# function for merging the various power calculations
.calc.power <- function(alpha=0.05, ltheta1, ltheta2, diffm, se, n, df, bk, 
                        method="exact")
{ 
  pow <- switch(
      method,
      exact=.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk),
      owenq=.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk),
      nct=  .approx.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk),
      noncentral=.approx.power.TOST(alpha, ltheta1, ltheta2, diffm, 
                                    se, n, df, bk),
      shifted=.approx2.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk),
      central=.approx2.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk),
      stop("Method '", method, "' unknown!\n", call.=TRUE)
  ) 
  return(pow)
}
