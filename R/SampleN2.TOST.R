# ----------------------------------------------------------------
# function to estimate the sample size of the second stage (in the
# pooled analysis one degree of freedom less than in the common
# model due to the 'stage'-term).
#
# Author: H. Schuetz, modified by D. Labes
# ----------------------------------------------------------------
sampleN2.TOST <- function(alpha=0.0294, CV, n1, theta0=0.95,
                          theta1=0.8, theta2=1.25, targetpower=0.8,
                          design="2x2", method="exact", imax=100)
  {
  if (missing(CV) || CV <=0) stop("CV must be given and positive.")
  if (missing(n1) || n1 < 4) stop("n1 must be given and >=4.")
  if (targetpower <=0 || targetpower >= 1) stop("targetpower must be >0 and <1.")
  design <- tolower(design)
  if (!design %in% c("2x2", "2x2x2", "parallel"))
     stop("design must be any of \"2x2\", \"2x2x2\", or \"parallel\"")
  if(theta0 < theta1 || theta0 > theta2) stop("theta0 not within BE limits.")
  # design constant in terms of ntot
  ifelse (design == "parallel", bk <- 4, bk <- 2)
  dfe   <- expression(n - 3) # pooled evaluation with additional 'stage' effect
  steps <- 2 # stepsize for sample size search
  nmin  <- 4 # minimum n
  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  diffm   <- log(theta0)
  se      <- CV2se(CV)
  # start value by Zhang's method
  n       <- .sampleN0_3(alpha, targetpower, ltheta1, ltheta2, diffm, se, steps, bk)
  if (n < nmin) n <- nmin
  df  <- eval(dfe)
  pow <- .calc.power(alpha, ltheta1, ltheta2, diffm, sem=se*sqrt(bk/n), df, method)
  iter <- 0
  down <- FALSE; up <- FALSE
  while (pow > targetpower) {
    if (n <= nmin) break
    down <- TRUE
    n    <- n - steps
    iter <- iter + 1
    df   <- eval(dfe)
    pow  <- .calc.power(alpha, ltheta1, ltheta2, diffm, sem=se*sqrt(bk/n), df, method)
    if (iter > imax) break
  }
  while (pow < targetpower) {
    up   <- TRUE; down <- FALSE
    n    <- n + steps
    iter <- iter + 1
    df   <- eval(dfe)
    pow  <- .calc.power(alpha, ltheta1, ltheta2, diffm, sem=se*sqrt(bk/n), df, method)
    if (iter > imax) break
  }
  # total sample size to stage 2 sample size
  ifelse (n - n1 > 0, n2 <- n - n1, n2 <- 0)
  if (n2 > 0) {
    pow <- .calc.power(alpha, ltheta1, ltheta2, diffm, sem=se*sqrt(bk/(n1+n2)),
                       df=n1+n2-3, method)
  } else {
    # usual df if study stops in stage 1
    pow <- power.TOST(alpha=alpha, theta1=theta1, theta2=theta2,
                      theta0=theta0, CV=CV, n=n1, method=method)
  }
  res <- data.frame(Design=design, alpha=alpha, CV=CV, theta0=theta0,
                    theta1=theta1, theta2=theta2, n1=n1, n2=n2, pow=pow,
                    target=targetpower)
  names(res) <- c("Design", "alpha", "CV", "theta0", "theta1", "theta2",
                  "n1", "Sample size", "Achieved power", "Target power")
  return(res)
}
