
 library(aster)

 ifam <- fam.truncated.poisson(truncation = 0)

 # change to include all cases in case splitting
 # see ../src/astfam.c lines 338 to 373
 theta <- seq(-6, 6, 0.5)
 mu <- exp(theta)

 zeroth <- double(length(theta))
 first <- double(length(theta))
 second <- double(length(theta))

 for (i in seq(along = theta)) {
    zeroth[i] <- famfun(ifam, 0, theta[i])
    first[i] <- famfun(ifam, 1, theta[i])
    second[i] <- famfun(ifam, 2, theta[i])
 }

 all.equal(zeroth, log(exp(mu) - 1))
 tau <- mu / (1 - exp(- mu))
 all.equal(first, tau)
 all.equal(second, tau * (1 - tau * exp(- mu)))

