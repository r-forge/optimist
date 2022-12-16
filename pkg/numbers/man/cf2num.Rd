\name{cf2num}
\alias{cf2num}
\alias{num2cf}
\title{
  Generalized Continous Fractions
}
\description{
  Evaluate a generalized continuous fraction as an alternating sum.
}
\usage{
cf2num(b0, b, a = 1, scaled = FALSE, tol = 1e-12)

num2cf(x, nterms=20)
}
\arguments{
  \item{b0}{absolute term, integer part of the continuous fraction.}
  \item{b}{numeric vector of length greater than 2.}
  \item{a}{numeric vector of length 1 or the same length as a.}
  \item{scaled}{logical; shall the convergents be scaled.}
  \item{tol}{relative tolerance.}

  \item{x}{real number.}
  \item{nterms}{number of terms.}
}
\details{
  \code{cf2num} calculates the numerical value of (simple or generalized) 
  continued fractions of the form
  \deqn{ b_0 + \frac{a1}{b1+} \frac{a2}{b2+} \frac{a3}{b3+...} }
  by converting its convergents into an alternating sum.
  The argument \eqn{a} is by default set to \eqn{a = (1, 1, ...)},
  that is the continued fraction is treated in its simple form.

  The convergents may grow and become too big, especially for large and
  growing coefficients \code{a}. Then \code{NA} will probably be
  returned. In this case use \code{scaled = TRUE} and the convergents
  will be scaled in every iteration, thus avoiding FP overflow
  .

  \code{num2cf} applies the direct calculation of the continued fraction.
  Attention: The input is not checked -- this allows for applying 'mpfr'
  numbers for longer, correct continued fractions, see the examples.

  The relation between the number of terms \code{nterms} and the
  precision of the 'mpfr' numbers should be about
  \code{nterm = 0.30..35 precBits}.
}
\value{
  Returns a numerical value, an approximation of the continued fraction.
}
\references{
  Press, Teukolsky, Vetterling, Flannery. NUMERICAL RECIPES: The Art of
  Scientific Computing. 3rd Edition, Cambridge University Press 2007.
}
\note{
  This function is \emph{not} vectorized.
}
\seealso{
  \code{\link{contfrac}}
}
\examples{
#-- Continued fraction of sqrt(2) is [1; 2, 2, 2, ...] -----------------
b0 <- 1; b <- rep(2, 20)                # sqrt(2)
cf2num(b0, b)                           # 1.414213562, error = eps()

#-- Approximate an analytic function -----------------------------------
# tan(x) = 0 + x / (1 + (-x^2) / 3 + (-x^2) / (5 + ...))
x  <- 0.5                           # tan(0.5) = 0.546302489844
n  <- 10                            # CF of length 20
b0 <- 0                             # b0 = 0
b  <- seq(1, by=2, length=n)        # b = c(1, 3, 5, ...)
a  <- c(x, rep(-x^2, times=n-1))    # a = c(x, -x^2, -x^2, ...)

cf2num(b0, b, a)                # 0.546302489844, error: eps()

\dontrun{
#-- Continued fraction of 1/(exp(0.5)-1) -------------------------------
e0 = 1/(exp(1/mpfr(2, precBits=128)) - 1)  #=> 1.54149408253679828413...
x0 <- num2cf(e0, nterms=20)
## [1]  1  1  1  5  1  1  9  1  1 13  1  1 17  1  1 21  1  1 25  1  1

# Determine e0 from its continued fraction
b0 <- x0[1]; b <- x0[2:19]
cf2num(b0, b)                       # 1.541494082536798, error = eps()

#-- pi as arctan(1.0) --------------------------------------------------
n = 24
b0 = 0
b = seq(1, by=2, length=n)          # b = c(1, 3, 5, 7, ...)
a = c(1, seq(1,by=1, length=n-1)^2) # a = c(1, 1, 4, 9, ...)
4 * cf2num(b0, b, a) - pi           # error: 0

#-- Leibniz-Wallis continued fraction for pi ---------------------------
# 4/pi = 1 + 1^2/(2 + 3^2 / (2 + 5^2 / 2 + ...))
pi4 = 4 / pi                        # 1.27323954474

n  = 20
b0 = 1
b  = rep(2, times=n)
a  = seq(1, by = 2, length=n)^2
cf2num(b0, b, a)                    # NA, i.e. convergents overflow
cf2num(b0, b, a, scaled = TRUE)
# 1.273272 with estimated precision 0.004032851 
# [1] 1.273272                      # actual error: 3.3e-05
  }
}
