\name{mknapsack}
\alias{mknapsack}
\title{
  Multiple 0-1 Knapsack Problem
}
\description{
  Solves the 0-1 (binary) multiple knapsack problem.
}
\usage{
mknapsack(w, p, cap)
}
\arguments{
  \item{w}{vector of (positive) weights.}
  \item{p}{vector of (positive) profits.}
  \item{cap}{vector of capacities of different knapsacks.}
}
\details{
  Solves the 0-1 multiple knapsack problem for a set of profits and weights.\cr
  A multiple 0-1 knapsack problem can be formulated as:

  \code{
  maximize  vstar = p(1)*(x(1,1) + ... + x(m,1)) + ...
                        ... + p(n)*(x(1,n) + ... + x(m,n))
  subject to
      w(1)*x(i,1) + ... + w(n)*x(i,n) <= cap(i)   for  i=1,...,m
      x(1,j) + ... + x(m,j) <= 1   for  j=1,...,n
      x(i,j) = 0 or 1   for  i=1,...,m ,  j=1,...,n ,
  }

  The multiple knapsack problem is reformulated as a linear program and
  solved with the help of package \code{lpSolve}.

  This function can be used for the single knapsack problem as well,
  but the 'dynamic programming' version in the \code{knapsack} function
  is faster (but: allows only integer values).

  Knapsack procedures can even solve subset sum problems,
  see the examples 4 and 5 below.

  The solution found is most often not unique and may not be the most compact
  one. In the future, we will attempt to 'compactify' through backtracking.
  The number of backtracks will be returned in list element \code{bs}.
}
\value{
  A list with components, \code{ksack} the knapsack numbers the items are
  assigned to, \code{value} the total value/profit of the solution found, and 
  \code{bs} the number of backtracks used.
}
\note{
  Contrary to earlier versions, the sequence of profits and weights has been
  interchanged: first the weights, then profits.
}
\references{
  Kellerer, H., U. Pferschy, and D. Pisinger (2004). Knapsack Problems.
  Springer-Verlag, Berlin Heidelberg.

  Martello, S., and P. Toth (1990). Knapsack Problems: Algorithms and
  Computer Implementations. John Wiley & Sons, Ltd.
}
\seealso{
  Other packages implementing knapsack routines.
}
\examples{
## Example 1: single knapsack
w <- c( 2,  20, 20, 30, 40, 30, 60, 10)
p <- c(15, 100, 90, 60, 40, 15, 10,  1)
cap <- 102
(is <- mknapsack(w, p, cap))
which(is$ksack == 1)
# [1] 1 2 3 4 6 , capacity 102 and total profit 280

## Example 2: multiple knapsack
w <- c( 40,  60, 30, 40, 20, 5)
p <- c(110, 150, 70, 80, 30, 5)
cap <- c(85, 65)
is <- mknapsack(w, p, cap)
# kps 1: 1,4;  kps 2: 2,6;  value: 345

## Example 3: multiple knapsack
p <- c(78, 35, 89, 36, 94, 75, 74, 79, 80, 16)
w <- c(18,  9, 23, 20, 59, 61, 70, 75, 76, 30)
cap <- c(103, 156)
is <- mknapsack(w, p, cap)
# kps 1: 3,4,5;  kps 2: 1,6,9; value: 452 

\dontrun{
##  Example 4: subset sum
p <- seq(2, 44, by = 2)^2
w <- p
is <- mknapsack(w, p, 2012)
sum((2 * which(is$ksack == 1))^2)  # 4 + 36 + 64 + 144 + 1764 = 2012

##  Example 5: maximize number of items
w <- seq(2, 44, by = 2)^2
p <- numeric(22) + 1
is <- mknapsack(w, p, 2012)  # 4+36+64+100+144+196+256+324+400+484 = 2012

## Example from Rosetta Code:
w = c(  9, 13, 153,  50, 15, 68, 27, 39, 23, 52, 11,
       32, 24,  48,  73, 42, 43, 22,  7, 18,  4, 30)
p = c(150, 35, 200, 160, 60, 45, 60, 40, 30, 10, 70,
       30, 15,  10,  40, 70, 75, 80, 20, 12, 50, 10)
cap = 400
system.time(is <- mknapsack(w, p, cap))  # 0.001 sec
}
}

