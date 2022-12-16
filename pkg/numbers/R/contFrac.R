##
##  c o n t F r a c . R  Continued Fractions
##


contfrac <- function(x, tol = 1e-12) {
    if (!is.numeric(x) || is.matrix(x))
        stop("Argument 'x' must be a numeric scalar or vector.")

    if (length(x) > 1) {
        # Compute value of a continuous fraction
        n <- length(x)
        p <- q <- numeric(n)
        B <- diag(1, 2)
        for (i in seq(along=x)) {
            B <- B %*% matrix(c(x[i], 1, 1, 0), 2, 2)
            p[i] <- B[1, 1]; q[i] <- B[2, 1]
        }
        return(list(f = B[1,1]/B[2,1], p = p, q = q, prec = 1/B[2, 1]^2))

    } else {
        # Generate the continuous fraction of a value
        sgnx <- sign(x)
        x <- abs(x)

        b <- floor(x)
        k <- b
        r <- x - b
        B <- matrix(c(b, 1, 1, 0), 2, 2)
        while ( abs(x - B[1,1]/B[2,1])  > tol) {
            b <- floor(1/r)
            k <- c(k, b)
            r <- 1/r - b
            B <- B %*% matrix(c(b, 1, 1, 0), 2, 2)
        }
        return(list(cf = sgnx * k, rat = c(sgnx*B[1,1], B[2,1]),
                    prec = abs(x - B[1,1]/B[2,1])))
    }
}


cf2num <- function(b0, b, a = 1, scaled = FALSE, tol = 1e-12) {
    stopifnot(is.numeric(a), is.numeric(b), is.numeric(b0))
    if (length(b0) != 1)
        stop("Argument 'b0' must be a single real value.")
    n <- length(b)
    if (n == 0) return(b0)
    if (length(a) != n) {
        if (length(a) == 1) a <- rep(a, n)
    } else if (length(b) != length(a)) {
        stop("length(a)==length(b) or length(a)==1 required.")
    }
    
    # Very short continued fractions
    if (n == 1) return (b0 + a[1]/b[1])
    else if (n == 2) return (b0 + b[2]*a[1]/(b[1]*b[2] + a[2]))
    
    if (! scaled) {
        # Calculate CF as an alternating sum
        q <- numeric(n)  # q_{-1} = 0; q_0 = 1
        q[1] <- b[1]; q[2] <- b[2]*b[1] + a[2]*1
        for (j in 3:n) {
            q[j] <-b[j]*q[j-1] + a[j]*q[j-2]
        }
        qq <- c(1, q[1:(n-1)]) * q
        pp <- (-1)^(0:(n-1)) * cumprod(a)
        aa <- pp / qq
        ss = sum(aa)
        return(b0 + ss)
        
    } else {
        # Calculate the scaled convergents
        n = length(b)
        pm1 = 1;  qm1 = 0
        p0  = b0; q0  = 1
        P = matrix(c(p0, q0, pm1, qm1), 2, 2)
        
        f1 = b0; f2 = Inf
        j = 0
        while (j < n && abs(f1-f2) > tol) {
            j = j + 1
            P = P %*% matrix(c(b[j], a[j], 1, 0), 2, 2)
            if (P[2,1] != 0) P = P / P[2, 1]
            f2 = f1; f1 = P[1,1]
        }
        f = (f1 + f2)/2.0
        cat(f, "with estimated precision", abs(f1-f2)/2.0, "\n")
        return(f)
    }
}


num2cf <- function(x, nterms=20) {
    # allow numerical or 'mpfr' type numbers
    # `no_entries = 0.3..0.35*precBits` for 'mpfr' numbers
    p = as.integer(floor(x))
    x = 1 / (x - floor(x))
    for (i in 1:nterms) {
        xf = floor(x)
        p = c(p, as.integer(xf))
        x = 1 / (x - xf)
    }
    return(p)
}


# cf2num <- function(a, b = 1, a0 = 0, finite = FALSE) {
#     stopifnot(is.numeric(a), is.numeric(b), is.numeric(a0))
#     n <- length(a)
#     if (length(b) != n) {
#         if (length(b) == 1) b <- rep(b, n)
#     } else if (length(a) != length(b)) {
#         stop("length(a)==length(b) or length(b)==1 required.")
#     }
#     
#     # Calculate CF as an alternating sum
#     q <- numeric(n)  # q_{-1} = 0; q_0 = 1
#     q[1] <- a[1]; q[2] <- a[2]*a[1] + b[2]*1
#     for (j in 3:n) {
#         q[j] <-a[j]*q[j-1] + b[j]*q[j-2]
#     }
#     qq <- c(1, q[1:(n-1)]) * q
#     pp <- (-1)^(0:(n-1)) * cumprod(b)
#     aa <- pp / qq
#     
#     if (finite) {
#         ss <- sum(aa)
#     } else {
#     # Apply Algorithm 1 from Cohen et al. (2000)
#         bb <- 2^(2 * n - 1)
#         cc <- bb
#         ss <- 0
#         for (k in (n-1):0) {
#             tt <- aa[k+1]
#             ss <- ss + cc * tt
#             bb <- bb * (2 * k + 1) * (k + 1)/(2 * (n - k) * (n + k))
#             cc <- cc + bb
#         }
#         ss <- ss / cc
#     }
#     # Don't forget the absolute term
#     return(a0 + ss)
# }

