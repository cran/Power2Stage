\name{power.2stage.Bf}
\alias{power.2stage.Bf}
\title{
Power calculation of adaptive 2-stage BE studies with a futility criterion
for the point estimator
}
\description{
This function calculates the empirical power of 2-stage BE studies according
to Potvin et.al. Method B via simulations. The Potvin methods are modified to 
include a futility criterion for the point estimator and to allow the sample size 
estimation step to be done with point estimator & mse of stage 1.
}
\usage{
power.2stage.Bf(alpha = c(0.0294, 0.0294), n1, CV, GMR, targetpower = 0.8, 
                pmethod = c("nct", "exact"), usePE = FALSE, powerstep = TRUE, 
                fPElower, fPEupper, theta0, theta1, theta2, 
                npct = c(0.05, 0.5, 0.95), nsims = 1e+05, setseed = TRUE, 
                print = TRUE, details = TRUE)
}
\arguments{
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
If \code{TRUE} (the default) the interim power analysis step in the stage 1 
evaluation will be done as described in Potvin et.al.\cr
Setting this argument to \code{FALSE} will omit this step.
}
  \item{fPElower}{
Lower futility limit for the PE of stage 1.\cr
If the PE is outside \code{fPElower} ... \code{fPEupper} the study is stopped
with the result FAIL (not BE).\cr
May be missing. Defaults then to 0.825.
}
  \item{fPEupper}{
Upper futility limit for the PE of stage 1.\cr
If missing it will be set to 1/fPElower (default 1/0.825=1.212121...).
}
  \item{theta0}{
True ratio of T/R for simulating. Defaults to the \code{GMR} argument.
}
  \item{theta1}{
Lower bioequivalence limit. Defaults to 0.8.
}
  \item{theta2}{
Lower bioequivalence limit. Defaults to 1.25.
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
The calculations follow in principle the simulations as described in Potvin et al.\cr
But instead of simulating subject data the statistics pe1, mse1 and pe2, SS2 are
simulated via their associated distributions (normal and chi-squared distri's).
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
Potvin D et.al.\cr
"Sequential design approaches for bioequivalence studies with 
crossover designs"\cr
Pharmaceut. Statist. 7(4), 245-62 (2008)

Montague TH et al.\cr
"Additional results for 'Sequential design approaches for bio-equivalence studies 
with crossover designs'"\cr
Pharmaceut. Statist. 11(1), 8-13 (2011)
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
# "examples check time a few seconds"
power.2stage.Bf(CV=0.25, n1=24, nsims=1E4)
# slight alpha inflation using the O'Brian-Fleming nominal alpha settings
# and small 'hidden pilot'
# commented out to complain with CRAN policy "examples check time a few seconds"
#power.2stage.Bf(alpha=c(0.005, 0.048), CV=0.2, n1=12, theta0=1.25, nsims=1E6)
}
