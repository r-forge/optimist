##  --------------------------------------------------------------------
##  modNthroot
##  --------------------------------------------------------------------


# library(numbers)

modNthroot = function(a, n, p) {
    ## Solve r^n = a mod p
    stopifnot(is.numeric(a), is.numeric(n), is.numeric(p),
              length(a) == 1, length(n) == 1, length(p) == 1)
    if (! isPrime(p)) stop("Argument 'p' must be a prime number")
    if (floor(n) != ceiling(n) || n <= 0)
        stop("Argument 'r' mut be an integer greater 0.")
    if (a == 0 || n == 1) return(a)
    if (n >= p-1) n = n %% (p-1)

    if (coprime(n, p-1)) {
        e = modinv(n, p-1)
        r = modpower(a, e, p)  #  a^e %% p

    } else {
        # f = c(1:(p-1))^n %% p
        f = numeric(p-1)
        for (i in 1:(p-1)) f[i] = modpower(i, n, p)
        r = which(f == a)
    }
    return(r)
}
