\name{approx_entropy}
\alias{approx_entropy}
\title{
  Approximate Entropy
}
\description{
  Calculates the approximate entropy of a time series.
}
\usage{
approx_entropy(ts, edim = 2, r = 0.2*sd(ts), elag = 1)
}
\arguments{
  \item{ts}{a time series.}
  \item{edim}{the embedding dimension, as for chaotic time series;
              a preferred value is 2.}
  \item{r}{filter factor; work on heart rate variability has suggested
           setting r to be 0.2 times the standard deviation of the data.}
  \item{elag}{embedding lag; defaults to 1, more appropriately it should be
              set to the smallest lag at which the autocorrelation function
              of the time series is close to zero.
              (At the moment it cannot be changed by the user.)}
}
\details{
  Approximate entropy was introduced to quantify the the amount of
  regularity and the unpredictability of fluctuations in a time series.
  A low value of the entropy indicates that the time series is deterministic;
  a high value indicates randomness.
}
\value{
  The approximate entropy, a scalar value.
}
\note{
  There exists a translation of the Kaplan code to R by Ben Bolker, see
  \url{http://www.macalester.edu/~kaplan/hrv/doc/funs/apen.html}.\cr
  This code here derives from a Matlab version at Mathwork's File Exchange,
  ``Fast Approximate Entropy'', by Kijoon Lee under BSD license.
}
\author{
  HwB  email: <hwborchers@googlemail.com>
}
\references{
  Pincus, S.M. (1991). Approximate entropy as a measure of system complexity.
  Proc. Natl. Acad. Sci. USA, Vol. 88, pp. 2297--2301.

  Kaplan, D., M. I. Furman, S. M. Pincus, S. M. Ryan, L. A. Lipsitz, and
  A. L. Goldberger (1991). Aging and the complexity of cardiovascular
  dynamics, Biophysics Journal, Vol. 59, pp. 945--949.
}
\seealso{
  \code{RHRV::CalculateApEn}
}
\examples{
ts <- rep(61:65, 10)
approx_entropy(ts, edim = 2)                      # -0.000936195

set.seed(8237)
approx_entropy(rnorm(500), edim = 2)              # 1.48944  high, random
approx_entropy(sin(seq(1,100,by=0.2)), edim = 2)  # 0.22831  low,  deterministic

\dontrun{(Careful: This will take several minutes.)
# generate simulated data
N <- 1000; t <- 0.001*(1:N)
sint   <- sin(2*pi*10*t);    sd1 <- sd(sint)    # sine curve
chirpt <- sint + 0.1*whitet; sd2 <- sd(chirpt)  # chirp signal
whitet <- rnorm(N);          sd3 <- sd(whitet)  # white noise

# calculate approximate entropy
rnum <- 30; result <- zeros(3, rnum)
for (i in 1:rnum) {
    r <- 0.02 * i
    result[1, i] <- approx_entropy(sint,   2, r*sd1)
    result[2, i] <- approx_entropy(chirpt, 2, r*sd2)
    result[3, i] <- approx_entropy(whitet, 2, r*sd3)
}

# plot curves
r <- 0.02 * (1:rnum)
plot(c(0, 0.6), c(0, 2), type="n",
     xlab = "", ylab = "", main = "Approximate Entropy")
points(r, result[1, ], col="red");    lines(r, result[1, ], col="red")
points(r, result[2, ], col="green");  lines(r, result[2, ], col="green")
points(r, result[3, ], col="blue");   lines(r, result[3, ], col="blue")
grid()}
}
\keyword{ timeseries }