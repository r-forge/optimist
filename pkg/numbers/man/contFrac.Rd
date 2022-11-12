\name{contFrac}
\alias{contFrac}
\alias{convergents}
\title{
  Continued Fractions
}
\description{
  Evaluate a continued fraction or generate one.
}
\usage{
contFrac(x, tol = 1e-12)

convergents(x)
}
\arguments{
  \item{x}{a numeric scalar or vector.}
  \item{tol}{tolerance; default \code{1e-12}.}
}
\details{
  If \code{x} is a scalar its continued fraction will be generated up to
  the accuracy prescribed in \code{tol}. If it is of length greater 1, the
  function assumes this is a continued fraction and computes its value.

  For implementation \code{contfrac} uses the representation of continued
  fractions through 2-by-2 matrices, i.e. the recursion formula.

  \code{convergents} will explicitly compute and return the convergents of
  the continued fraction through this recursion.
}
\value{
  Either a numeric value, or a list with components \code{cf}, numeric vector
  representing the continued fraction \eqn{[b_0; b_1, \ldots, b_{n-1}]};
  \code{rat}, the rational number as a vector with (numerator, denumerator);
  and \code{prec}, the difference between \code{x} and the value of the
  contimued fraction.

  \code{convergents} returns a list with components \code{p} and \code{q}.
}
\references{
  Hardy, G. H., and E. M. Wright (1979). An Introduction to the Theory of
  Numbers. Fifth Edition, Oxford University Press, New York.
}
\note{
  This function are \emph{not} vectorized.
}
\seealso{
\code{\link{cf2num}}, \code{\link{ratFarey}}
}
\examples{
contFrac(pi)
contFrac(c(3, 7, 15, 1))        # rational Approx: 355/113

contFrac(0.555)                 #  0  1  1  4 22
contFrac(c(1, rep(2, 25)))      #  1.414213562373095, sqrt(2)

convergents(c(3, 7, 15, 1))     # $p : 3  22 333 355
                                # $q : 1   7 106 113
}
