\name{power.2stage.GS}
\alias{power.2stage.GS}
\title{
Power calculation for non-adaptive group sequential (2-stage) BE studies
}
\description{
This function calculates the 'empirical' power of group sequential 2-stage BE 
with 2x2 crossover designs via simulations.\cr
The number of subjects in both stages have to be prespecified (non-daptive).
}
\usage{
power.2stage.GS(alpha = c(0.0294, 0.0294), n, CV, theta0, theta1, theta2, 
                fCrit = c("PE", "CI"), fClower, fCupper, nsims = 1e+05, 
                setseed = TRUE, print = TRUE, details = TRUE)
}
\arguments{
  \item{alpha}{
Vector of the two nominal alpha values to be used in the 1-2*alpha confidence
intervals calculations.
}
  \item{n}{
Vector of the two sample sizes in stage 1 and stage 2.\cr
n(total is n1+n2).
}
  \item{CV}{
Coefficient of variation of the intra-subject variability as ratio.
}
  \item{theta0}{
True ratio of T/R for simulating. Defaults to 0.95 if missing.
}
  \item{theta1}{
Lower bioequivalence limit. Defaults to 0.8.
}
  \item{theta2}{
Upper bioequivalence limit. Defaults to 1.25.
}
  \item{fCrit}{
Futility criterion.\cr
If set to \code{"PE"} the study stops after stage 1 if not BE and if the point 
estimator of stage 1 evaluation is outside the range defined in the next two 
argumenmts \code{"fClower"} and \code{"fCupper"}.\cr
If set to \code{"CI"} the study stops after stage 1 if not BE and if the 
confidence interval of stage 1 evaluation is outside the range defined in the 
next two argumenmts.\cr
Defaults to \code{"PE"}.
}
  \item{fClower}{
Lower limit of the futility criterion.\cr
Defaults to \code{theta1} if missing, i.e. no futility criterion.
}
  \item{fCupper}{
Upper limit of the futility criterion.\cr
Defaults to \code{1/fClower} if missing.
}
  \item{nsims}{
Number of studies to simulate. Defaults to 1E+05 = 100 000.\cr
Set this at least to 1E+06 = 1 mio if you are calculating 'alpha'.
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
The calculations follow in principle the simulations as described in 
Potvin et al. for adaptive designs, but with no interim power inspection and 
pre-specified (fixed) number of subjects in stage 2.\cr
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
}
\references{
Gould, A.L. \cr
"Group sequential extensions of a standard bioequivalence testing procedure."\cr
J. Pharmacokinetics and Biopharmaceutics, 23, 57-86. (1995)

Patterson S and Jones B\cr
"BIOEQUIVALENCE and STATISTICS in CLINICAL PHARMACOLOGY"\cr
Chapman & Hall/CRC\cr
Boca Raton 2006\cr
Chapter 5.6 Optional Designs

Potvin D et.al.\cr
"Sequential design approaches for bioequivalence studies with crossover designs"\cr
Pharmaceut. Statist. 7(4), 245-62 (2008)
}
\author{
D. Labes
}
\note{
The code is reasonable fast. 1E6 sims take ~ 1 sec on my machine.
Even 1E7 sims are possible without meanwhile too much beer :-)).\cr\cr
To calculate without a futility criterion choose \code{fCrit="PE"} and set 
\code{fClower=0, fCupper=Inf}.
}
\section{Warning }{
This function has to be considered as 'experimental' since it was not thoroughly tested.
}
\seealso{
\code{\link{power.2stage}} and \code{\link{power.2stage.p}} for adptive sequential designs.
}
\examples{
# using the Pocock alpha settings and the point estimator futility rule
# (defaults) and 12 subjects in both stages
power.2stage.GS(CV=0.2, n=c(12,12))
}
