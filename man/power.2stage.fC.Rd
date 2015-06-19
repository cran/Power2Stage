\name{power.2stage.fC}
\alias{power.2stage.fC}
\title{
Power calculation of adaptive 2-stage BE studies (2x2 crossover) with a 
futility criterion for the point estimator or CI
}
\description{
This function calculates the empirical power of 2-stage BE studies according
to Potvin et al. Method B/C via simulations. The Potvin methods are modified to 
include a futility criterion for the point estimator or for the CI and to 
allow the sample size estimation step to be done with point estimator & mse 
of stage 1.
}
\usage{
power.2stage.fC(method = c("B", "C"), alpha0 = 0.05, alpha = c(0.0294, 0.0294), 
                n1, CV, GMR, targetpower = 0.8, pmethod = c("nct", "exact", "shifted"), 
                usePE = FALSE, powerstep = TRUE, min.n2=0, fCrit=c("PE","CI"), 
                fClower, fCupper, theta0, theta1, theta2, npct = c(0.05, 0.5, 0.95), 
                nsims, setseed = TRUE, print = TRUE, details = FALSE)
}
\arguments{
  \item{method}{
Decision schemes according to Potvin et.al. Default is "B". \cr
Potvin method D can be obtained by choosing "C" but setting 
\code{alpha=c(0.028,0.028)}.
}
  \item{alpha0}{
Alpha value for the first step(s) in Potvin C, the power inspection
and BE decision if power > targetpower.\cr
Defaults to 0.05.
}
  \item{alpha}{
Vector of the nominal alpha's for the two stages.\cr
Defaults to the Pocock's alpha setting \code{alpha=c(0.0294,0.0294)}.
}
  \item{n1}{
Sample size of stage 1.
}
  \item{CV}{
Coefficient of variation of the intra-subject variability as ratio.
}
  \item{GMR}{
Ratio T/R to be used in Potvin "B" decision scheme (power calculations 
in stage 1 and sample size estimation for stage 2).
}
  \item{targetpower}{
Power threshold for the power to achieve in the sample size estimation step.
}
  \item{pmethod}{
Power calculation method, also to be used in the sample size estimation for stage 2.\cr 
Implemented are "nct" = approximate calculations via non-central t-distribution, 
"exact" = exact calculations via Owen's Q and "shifted" = approximate calculation
via shifted central t-distribution like in the paper of Potvin et al.\cr
Defaults to "nct" as reasonable compromise between speed and accuracy in the 
sample size estimation step.
}
  \item{usePE}{
If \code{TRUE} the sample size estimation step is done with mse 
\bold{and} PE of stage 1.\cr
Defaults to \code{FALSE} in wich case the sample size is estimated with 
(constant) GMR given as argument and mse of stage 1 (analogous to Potvin et. al.).
}
  \item{powerstep}{
If \code{TRUE} (the default) the interim power monitoring step in the stage 1 
evaluation of method "B" will be done as described in Potvin et.al.\cr
Setting this argument to \code{FALSE} will omit this step.\cr
Has no effect if \code{method="C"} is choosen.
}
  \item{min.n2}{
Minimum sample size of stage 2. Defaults to zero.\cr
If the sample size estimation step gives N < n1 the sample size for stage 2 
will be set to \code{min.n2}, i.e. the total sample size to \code{n1+min.n2}.
}
  \item{fCrit}{
Futility criterion to use for PE or CI.
}
  \item{fClower}{
Lower futility limit for the PE or CI of stage 1.\cr
If the PE or CI is outside \code{fClower} ... \code{fCupper} the study is stopped
with the result FAIL (not BE).\cr
May be missing. Defaults then to 0.825 if \code{fCrit="PE"} or 0.9 if \code{fCrit="CI"}.
}
  \item{fCupper}{
Upper futility limit for the PE or CI of stage 1.\cr
If missing it will be set to 1/fClower.
}
  \item{theta0}{
True ratio of T/R for simulating. Defaults to the \code{GMR} argument.
}
  \item{theta1}{
Lower bioequivalence limit. Defaults to 0.8.
}
  \item{theta2}{
Upper bioequivalence limit. Defaults to 1.25.
}
  \item{npct}{
p-values to be used for the percentiles of the distribution of n(total)=n1+n2.\cr
Defaults to \code{c(0.05, 0.5, 0.95)} to obtain the 5\% and 95\% percentiles
and the median.
}
  \item{nsims}{
Number of studies to simulate.\cr 
If missing \code{nsims} is set to 1E+05 = 100 000 or to 1E+06 = 1 Mio if you are 
calculating 'alpha', i.e. with \code{theta0} at border or outside acceptance range 
\code{theta1} ... \code{theta2}.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power for different runs a 
\code{set.seed(1234567)} is issued if \code{setseed=TRUE}, the default.\cr
Set this argument to \code{FALSE} to view the variation in power between 
different runs.
}
  \item{print}{
This argument is now defunct since printing is done via S3 method print.\cr  
It is retained for compatibility but will be removed in the next version. 
}
  \item{details}{
If set to \code{TRUE} the function prints the results of time measurements
of the simulation steps. Default is \code{FALSE}.
}
}
\details{
The calculations follow in principle the simulations as described in Potvin 
et al.\cr
The underlying subject data are assumed to be evaluated after log-transformation.
But instead of simulating subject data the statistics pe1, mse1 and pe2, SS2 are
simulated via their associated distributions (normal and chi-squared 
distri's).
}
\value{
Returns an object of class "pwrtsd" with all the input arguments and results 
as components.\cr
The class "pwrtsd" has an S3 print method.\cr
The results are in the components:
\item{pBE}{Contains the ratio of studies found BE.}
\item{pBE_s1}{Ratio of studies found BE in stage 1.}
\item{pct_s2}{Percentage of studies continuing to stage 2.}
\item{nmean}{Mean of n(total).}
\item{nrange}{Range (min, max) of n(total).}
\item{nperc}{Percentiles of the distribution of n(total).}
\item{ntable}{Object of class "table" summarizing the discrete distribution of 
n(total) via its distinct values and counts of occurences of these values.\cr
This component is only given back if \code{usePE==FALSE} or  
\code{usePE==TRUE & fClower>0 & is.finite(fCupper)}, i.e. a futility range is used.}
}
\references{
Potvin D et al.\cr
"Sequential design approaches for bioequivalence studies with crossover designs"\cr
Pharmaceut. Statist. 7(4), 245-62 (2008)

Schuetz H.\cr
"Two-stage designs in bioequivalence trials"\cr
Eur J Clin Pharmacol 71(3):271-81 (2015)\cr
published inline 22 Jan 2015

Kieser M, Rauch G\cr
"Two-stage designs for cross-over bioequivalence trials"\cr
Stat Med 2015 Mar 24\cr
epub ahead of print

Zheng Ch, Zhao L, Wang J\cr
"Modifications of sequential designs in bioequivalence trials"\cr
Pharm Stat 14 (3), 180-188, May/June 2015\cr
published online: 9 FEB 2015
}
\author{
D. Labes
}
\seealso{
\code{\link{power.2stage}}
}
\examples{
# using all the defaults
power.2stage.fC(CV=0.25, n1=24)
# run-time ~1 sec
\dontrun{
# as above but storing the results
res <- power.2stage.fC(CV=0.25, n1=24)
# making a plot of the discrete distribution of ntotal
plot(res$ntable/sum(res$ntable), ylab="Density", main="Distribution of n(total)")
}
}
