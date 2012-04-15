##
##  s i m p l e D E . R  Simple Differential Evolution Algorithm
##


simpleDE <-
function(fun, lower, upper, N = 60, nmax = 300, r = 0.4, log = FALSE)
{
    n <- length(lower)
    if (length(upper) != n)
        stop("'lower' and 'upper' must of the same length.")

    G <- H <- matrix(runif(N*n), nrow = N, ncol = n)
    F <- apply(G, 1, fun)

    for (g in 1:nmax) {
        for (i in 1:N) {
            ii <- sample(1:N, 3)
            ci <- G[ii[1], ] + r * (G[ii[2], ] - G[ii[3], ])
            fi <- fun(ci)
            if (fi < F[i]) {
                H[i, ] <- ci
                F[i] <- fi
            }
        }
        G <- H

        if (log && (g %% 10 == 0))
            cat(g, "  ", "\t", min(F), "\n", sep = "")
    }

    i0 <- which.min(F)
    return( list(fmin = F[i0], xmin = G[i0, ]) )
}