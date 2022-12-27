\name{modNthroot}
\alias{modNthroot}
\title{
  N-th root modulo p
}
\description{
  Find all n-th roots r^n = a mod p of an integer a for p prime.
}
\usage{
modNthroot(a, n, p)
}
\arguments{
  \item{a}{an integer.}
  \item{n}{exponent, an integer.}
  \item{p}{a prime number.}
}
\details{
  Computes the n-th root of an integer modulo a prime number p, i.e., 
  solves the equation \eqn{r^n = a \mod p} with \eqn{p} a prime number.
}
\value{
  Returns a unique solution an integer, the solution \code{r} of 
  \code{r^n = a mod p} when \code{coprime(n, p-1)} -- or is empty
  when there is no solution. Returns an array of integer solutions else.

  In the first case the code is very efficient. In the second case, the
  search is exhaustive, but still quite fast for not too big numbers.
}
\references{
  E. Bach and J. Shallit. Algorithmic Number Theory.
  Vol 1: Efficient Algorithms.The MIT Press, Cambridge, MA, 1996.
}
\note{
  There is a more efficient algorithm if \code{n} and \code{p-1} have
  common prime divisors. This may be implemented in a future version.
}
\seealso{
  \code{\link{modpower}}
}
\examples{
a = 10; n = 5; p = 13           # the best case
modNthroot(a, n, p)             # 4
a = 10; n = 17; p = 13          # n greater than p-1
modNthroot(a, n, p)             # 4
a = 9; n = 4; p = 13            # n and p-1 not coprime
modNthroot(a, n, p)             # 4 6 7 9
a = 17; n = 35; p = 101         # 7 is a prime divisor of n and p-1
modNthroot(a, n, p)             # 9 21 47 49 76
a = 5001; n = 5; p = 100003     # some bigger numbers
modNthroot(a, n, p)             # 47768
}
