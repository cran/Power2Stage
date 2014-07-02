\name{power.2stage.p}
\alias{power.2stage.p}
\title{
Power calculation of adaptive 2-stage BE studies with 2-group parallel design
}
\description{
This function calculates the 'empirical' power of 2-stage BE studies according
to A. Fuglsang via simulations. The Potvin decision schemes are modified to 
include a futility criterion Nmax and to allow the sample size estimation 
step to be done with point estimator & variabilities from stage 1.
}
\usage{
power.2stage.p(method = c("B", "C"), alpha0 = 0.05, alpha = c(0.0294, 0.0294), 
              n1, GMR, CV, targetpower = 0.8, pmethod = c("nct", "exact"), 
              usePE = FALSE, Nmax = Inf, test = c("welch", "t-test", "anova"), 
              theta0, theta1, theta2, npct = c(0.05, 0.5, 0.95), 
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
Vector (length=2) of the nominal alpha's for the two stages.\cr
Defaults to the Pocock alpha setting \code{alpha=c(0.0294,0.0294)}.
}
  \item{n1}{
(Total) Sample size for stage 1.
}
  \item{GMR}{
Ratio T/R to be used in the Potvin decision schemes (power calculations 
in stage 1 and sample size estimation for stage 2).
}
  \item{CV}{
Coefficient of variation of the total variability as ratio.\cr
May be a vector (length=2). In that case CV[1] is for Test, CV[2] for Reference 
treatment.
}
  \item{targetpower}{
Power threshold in the interim power inspection steps of Potvin "C" and "B" and
power to achieve in the sample size estimation step.
}
  \item{pmethod}{
Power calculation method, also to be used in the sample size estimation for 
stage 2.\cr 
Implemented are "nct" = approximate calculations via non-central t-distribution 
and "exact" = exact calculations via Owen's Q.\cr
Defaults to "nct" for speed reasons in the sample size estimation step.
}
  \item{usePE}{
If \code{TRUE} the sample size estimation step is done with PE and mse of 
stage 1.\cr
Defaults to \code{FALSE} in wich case the sample size is estimated with GMR and 
mse of stage 1 (analogous to Potvin et. al.).\cr
The power inspection steps in the Potvin methods are always done with the
GMR argument and mse (CV) of stage 1.
}
  \item{Nmax}{
Futility criterion. If set to a finite value all studies simulated in which a 
sample size >Nmax is obtained will be regarded as BE=FAIL.\cr
Set this argument to \code{Inf}, the default, to work without that futility
criterion.
}
  \item{test}{
Test on which the CI calculations are based on.\cr
Defaults to \code{"welch"} = Welch's t-test accounting for heteroscedasticity in the
variabilities of Test and Reference, but neglecting stage effects.\cr
\code{"anova"} calculates the 1-2*alpha CI's based on an ANOVA with treatment and stage
in the model.\cr
\code{"t-test"} calculates the 1-2*alpha CI's based on t-test assuming equal 
variabilities of Test and Reference and neglecting stage effects.
}
  \item{theta0}{
True ratio of T/R for simulating. Defaults to the \code{GMR} argument if missing.
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
generator. To avoid differences in power for different runs a \code{set.seed(1234567)} 
is issued if \code{setseed=TRUE}, the default.\cr
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
The calculations follow in principle the simulations as described in Fuglsang A.\cr
The underlying subject data are assumed to be evaluated after log-transformation.
But instead of simulating subject data the statistics (mean and variance of Test
and Reference of stage 1 and stage 2) are simulated via their associated 
distributions (normal and chi-squared distri's).
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
Fuglsang A.\cr
"Sequential Bioequivalence Approaches for Parallel Design"\cr
AAPS J. 2014 Feb 14. [Epub ahead of print]

Potvin D et al.\cr
"Sequential design approaches for bioequivalence studies with crossover designs"\cr
Pharmaceut. Statist. 7(4), 245-62 (2008)
}
\author{
D. Labes
}
\note{
The computation time is between ~ 0.01 - 1 min for 100 000 sims on my machine
(Intel core i3 2.93 GHz, 4GB RAM). Thus be patient and take a pot of coffee if 
you simulate for 'alpha' with 1 Mio sims.
}
\seealso{
\code{\link{power.2stage}} for analogous calculations for the 2x2 crossover.
}
\examples{
# using all the defaults
# but with too low number of sims to complain with CRAN policy: 
# "Examples should run for no more than a few seconds each"
# minimum number of sims should be 1E5 for 'power', 1E6 sims for 'alpha'!
power.2stage.p(n1=48, CV=0.25, nsims=1E4)
# computation time ~ 8 sec if 1E5 sims were used.
}
