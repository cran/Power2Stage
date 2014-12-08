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
                n1, CV, GMR, targetpower = 0.8, pmethod = c("nct", "exact"), 
                usePE = FALSE, powerstep = TRUE, min.n2=0, fCrit=c("PE","CI"), 
                fClower, fCupper, theta0, theta1, theta2, npct = c(0.05, 0.5, 0.95), 
                nsims = 1e+05, setseed = TRUE, print = TRUE, details = TRUE)
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
Defaults to the Pocock alpha setting \code{alpha=c(0.0294,0.0294)}.
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
Implemented are "nct" = approximate calculations via non-central t-distribution 
and "exact" = exact calculations via Owen's Q.\cr
Defaults to "nct" for speed reasons in the sample size estimation step.
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
Number of studies to simulate. Defaults to 1E+05 = 100 000.\cr
Set this to 1E+06 = 1 mio if you are calculating 'alpha'.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power for different runs a 
\code{set.seed(1234567)} is issued if \code{setseed=TRUE}, the default.\cr
Set this argument to \code{FALSE} to view the variation in power between 
different runs.
}
  \item{print}{
If \code{TRUE} (default) the function prints its results.\cr 
If \code{FALSE} only the list with the results will be returned (see value). 
}
  \item{details}{
If \code{TRUE} (default) the function prints the results of time measurements
of the simulation steps. 
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
Returns a list with all the input arguments and results as components.\cr
The results are in the components:
\item{pBE}{Contains the ratio of studies found BE.}
\item{pBE_s1}{Ratio of studies found BE in stage 1.}
\item{pct_s2}{Percentage of studies continuing to stage 2.}
\item{nmean}{Mean of n(total).}
\item{nrange}{Range (min, max) of n(total).}
\item{nperc}{Percentiles of the distribution of n(total).}
}
\references{
Potvin D et al.\cr
"Sequential design approaches for bioequivalence studies with 
crossover designs"\cr
Pharmaceut. Statist. 7(4), 245-62 (2008)
}
\author{
D. Labes
}
\note{
The computation time is between ~ 0.01 - 0.4 min for 100 000 sims on my machine
(Intel core i3 2.93 GHz, 4GB RAM). Thus be patient and take a pot of coffee if 
you simulate for 'alpha' with 1 Mio (1E6) sims.
}
\seealso{
\code{\link{power.2stage}}
}
\examples{
# using all the defaults
# but too low number of sims to complain with CRAN policy 
# "examples check time a few seconds" (a few means <5 sec)
power.2stage.fC(CV=0.25, n1=24, nsims=1E4)
}
