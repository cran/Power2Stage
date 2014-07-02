\name{power.2stage}
\alias{power.2stage}
\title{
Power calculation of adaptive 2-stage BE studies with 2x2 crossover design
}
\description{
This function calculates the empirical power of 2-stage BE studies according
to Potvin et.al. via simulations. The Potvin methods are modified to include
a futility criterion Nmax and to allow the sample size estimation step to be 
done with point estimator & mse of stage 1.
}
\usage{
power.2stage(method = c("B", "C"), alpha0 = 0.05, alpha = c(0.0294, 0.0294), 
             n1, GMR, CV, targetpower = 0.8, 
             pmethod = c("nct", "exact", "shifted"),
             usePE = FALSE, Nmax = Inf, theta0, theta1, theta2,  
             npct = c(0.05, 0.5, 0.95), nsims = 1e+05, setseed = TRUE, 
             print = TRUE, details = TRUE)
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
Sample size for stage 1.
}
  \item{GMR}{
Ratio T/R to be used in the Potvin decision schemes (power calculations 
in stage 1 and sample size estimation for stage 2).
}
  \item{CV}{
Coefficient of variation of the intra-subject variability as ratio.
}
  \item{targetpower}{
Power threshold in the first step of Potvin "C" and 
power to achieve in the sample size estimation step.
}
\item{pmethod}{
Power calculation method, also to be used in the sample size estimation for 
stage 2.\cr 
Implemented are \code{"nct"} = approximate calculations via non-central 
t-distribution and \code{"exact"} = exact calculations via Owen's Q as well as
the approximation via shifted central t-distribution, as used by Potvin et al.\cr
Defaults to \code{"nct"} for speed reasons in the sample size estimation step.
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
The calculations follow in principle the simulations as described in 
Potvin et al.\cr
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
"Sequential design approaches for bioequivalence studies with crossover designs"\cr
Pharmaceut. Statist. 7(4), 245-62 (2008)

Montague TH et al.\cr
"Additional results for 'Sequential design approaches for bio-equivalence studies 
with crossover designs'"\cr
Pharm Stat. 2012 Jan-Feb;11(1):8-13. Epub 2011 Feb 10.

Fuglsang A\cr
"Controlling type I errors for two-stage bioequivalence study designs"\cr
Clinical Research and Regulatory Affairs 28(4), 100-105 (2011)

Fuglsang A\cr
"Sequential Bioequivalence Trial Designs with Increased Power and 
Controlled Type I Error Rates"\cr
AAPS J. 2013 Jul;15(3):659-61. Epub 2013 Mar 30.

Fuglsang A\cr
"Futility Rules in Bioequivalence Trials with Sequential Designs"\cr
AAPS J. 2014 Jan;16(1):79-82. Epub 2013 Nov 12.
}
\author{
D. Labes
}
\note{
The computation time is between ~ 0.01 - 0.4 min for 100 000 sims on my machine
(Intel core i3 2.93 GHz, 4GB RAM). Thus be patient and take a pot of coffee if 
you simulate for 'alpha' with 1 Mio sims.
}
\seealso{
\code{\link{power.2stage.p}} for analogous calculations for 2-group parallel 
design.
}
\examples{
# using all the defaults
# but too low number of sims to complain with CRAN policy 
# "Examples should run for no more than a few seconds each"
# minimum number of sims should be 1E5 for 'power', 1E6 sims for 'alpha'
power.2stage(n1=24, CV=0.25, nsims=1E4)
# computation time ~ 16 sec if 1E5 sims were used.
}
