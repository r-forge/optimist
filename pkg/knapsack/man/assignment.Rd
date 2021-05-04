\name{assignment}
\alias{assignment}
\title{
 Linear Sum Assignment Problem
}
\description{
  Linear (sum) assignment problem, or LSAP.
}
\usage{
assignment(cmat)
}
\arguments{
  \item{cmat}{quadratic integer matrix, the cost matrix.}
}
\details{
  Solves the linear (sum) assignment problem for quadratic matrices with
  integer entries.
}
\value{
  List with components \code{perm}, the permutation that defines the minimum 
  solution, \code{min}, the minimum value, and \code{err}, which is \code{-1} 
  if an integer overflow occured.
}
\references{
  Burkard, R., M. Dell'Amico, and S. Martello (2009). Assignment Problems.
  Society for Industrial and Applied Mathematics (SIAM).

  Martello, S., and P. Toth (1990). Knapsack Problems: Algorithms and
  Computer Implementations. John Wiley & Sons, Ltd.
}
\author{
  Copyright(c) 1993 A. H. Morris, Jr., Naval Surface Warfare Center, using
  Fortran routines written by S. Martello and P. Toth, University of Bologna.
  Released for free and general use, now under GPL license, wrapped for R by
  Hans W Borchers <hwborchers@googlemail.com>.
}
\note{
  Faster than the Hungarian algorithm, but only applicable to quadratic cost
  matrices with integer components.
}
\seealso{
  \code{clue::solve_LSAP}
}
\examples{
##  Example similar to clue::solve_LSAP
set.seed(8237)
x <- matrix(sample(1:100), nrow = 10)
y <- assignment(x)
# show permutation and check minimum sum
y$perm; y$min                   # 7  6  8  5  1  9  2  4 10  3
z <- cbind(1:10, y$perm)        # 173
x[z]                            # 16  9 15  6 34  3  4 44 35  7
y$min == sum(x[z])              # TRUE

\dontrun{
##  Example: minimize sum of distances of complex points
n <- 100
x <- rt(n, df=3) + 1i * rt(n, df=3)
y <- runif(n) + 1i * runif(n)
cmat <- round(outer(x, y, FUN = function(x,y) Mod(x - y)), 2)
dmat <- round(100*cmat, 2)
system.time(T1 <- assignment(dmat))   # elapsed: 0.002
T1$min / 100                          # 145.75

library("clue")
system.time(T2 <- solve_LSAP(cmat))     # elapsed: 0.016
sum(cmat[cbind(1:n, T2)])               # 145.75
}
}