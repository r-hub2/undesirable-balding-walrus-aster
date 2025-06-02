
 library(aster)

 data(echinacea)

 vars <- c("ld02", "ld03", "ld04", "fl02", "fl03", "fl04",
     "hdct02", "hdct03", "hdct04")
 redata <- reshape(echinacea, varying = list(vars), direction = "long",
     timevar = "varb", times = as.factor(vars), v.names = "resp")
 redata <- data.frame(redata, root = 1)
 pred <- c(0, 1, 2, 1, 2, 3, 4, 5, 6)
 fam <- c(1, 1, 1, 1, 1, 1, 3, 3, 3)
 hdct <- grepl("hdct", as.character(redata$varb))
 redata <- data.frame(redata, hdct = as.integer(hdct))
 level <- gsub("[0-9]", "", as.character(redata$varb))
 redata <- data.frame(redata, level = as.factor(level))
 aout <- aster(resp ~ varb + level : (nsloc + ewloc) + hdct : pop,
     pred, fam, varb, id, root, data = redata)
 sout <- summary(aout)
 vout <- vcov(aout)
 all.equal(sqrt(diag(vout)), sout$coefficients[ , "Std. Error"])

 rm(list = ls())

 data(radish2)

 pred <- c(0,1,2)
 fam <- c(1,3,2)

 rout <- reaster(resp ~ varb + fit : (Site * Region),
     list(block = ~ 0 + fit : Block, pop = ~ 0 + fit : Pop),
     pred, fam, varb, id, root, data = radish2)

 # re.too == FALSE and complete == FALSE and standard.deviation == FALSE

 sout <- summary(rout, standard.deviation = FALSE)
 vout <- vcov(rout, standard.deviation = FALSE, complete = FALSE)

 foo <- sout$alpha[ , "Std. Error"]
 bar <- sqrt(diag(vout)[attr(vout, "is.alpha")])
 all.equal(foo, bar, tolerance = 1e-5)

 foo <- sout$nu[ , "Std. Error"]
 foo <- foo[! is.na(foo)]
 bar <- sqrt(diag(vout)[attr(vout, "is.nu")])
 all.equal(foo, bar, tolerance = 1e-4)

 izzy <- attributes(vout)
 izzy <- izzy[grepl("^is\\.", names(izzy))]
 all(sapply(izzy, length) == nrow(vout))
 all(Reduce("+", izzy) == 1)

 # re.too == FALSE and complete == FALSE and standard.deviation == TRUE

 sout <- summary(rout, standard.deviation = TRUE)
 vout <- vcov(rout, standard.deviation = TRUE, complete = FALSE)

 foo <- sout$alpha[ , "Std. Error"]
 bar <- sqrt(diag(vout)[attr(vout, "is.alpha")])
 all.equal(foo, bar, tolerance = 1e-5)

 foo <- sout$sigma[ , "Std. Error"]
 foo <- foo[! is.na(foo)]
 bar <- sqrt(diag(vout)[attr(vout, "is.sigma")])
 all.equal(foo, bar, tolerance = 1e-4)

 izzy <- attributes(vout)
 izzy <- izzy[grepl("^is\\.", names(izzy))]
 all(sapply(izzy, length) == nrow(vout))
 all(Reduce("+", izzy) == 1)

 # re.too == FALSE and complete == TRUE and standard.deviation == FALSE

 sout <- summary(rout, standard.deviation = FALSE)
 vout <- vcov(rout, standard.deviation = FALSE, complete = TRUE)

 foo <- sout$alpha[ , "Std. Error"]
 bar <- sqrt(diag(vout)[attr(vout, "is.alpha")])
 all.equal(foo, bar, tolerance = 1e-5)

 foo <- sout$nu[ , "Std. Error"]
 bar <- sqrt(diag(vout)[attr(vout, "is.nu")])
 foo[is.na(foo)] <- 0
 all.equal(foo, bar, tolerance = 1e-4)

 izzy <- attributes(vout)
 izzy <- izzy[grepl("^is\\.", names(izzy))]
 all(sapply(izzy, length) == nrow(vout))
 all(Reduce("+", izzy) == 1)

 # re.too == FALSE and complete == TRUE and standard.deviation == TRUE

 sout <- summary(rout, standard.deviation = TRUE)
 vout <- vcov(rout, standard.deviation = TRUE, complete = TRUE)

 foo <- sout$alpha[ , "Std. Error"]
 bar <- sqrt(diag(vout)[attr(vout, "is.alpha")])
 all.equal(foo, bar, tolerance = 1e-5)

 foo <- sout$sigma[ , "Std. Error"]
 bar <- sqrt(diag(vout)[attr(vout, "is.sigma")])
 foo[is.na(foo)] <- 0
 all.equal(foo, bar, tolerance = 1e-4)

 izzy <- attributes(vout)
 izzy <- izzy[grepl("^is\\.", names(izzy))]
 all(sapply(izzy, length) == nrow(vout))
 all(Reduce("+", izzy) == 1)

 # re.too == TRUE
 # see equations (4.16) and (4.17) and (4.21) of the Aster Theory book
 # make as function to do two cases with one bit of code

 doit <- function(standard.deviation) {

 objfun <- objfun.factory(rout$fixed, rout$random, rout$response,
     rout$obj$pred, rout$obj$fam, as.vector(rout$obj$root),
     rout$zwz, deriv = 2, standard.deviation = standard.deviation)
 if (standard.deviation) {
     theta.hat <- c(rout$alpha, rout$c, rout$sigma)
 } else {
     theta.hat <- c(rout$alpha, rout$b, rout$nu)
 }
 oout <- objfun(theta.hat)

 is.alpha <- seq_along(theta.hat) <= length(rout$alpha)
 is.nu <- seq_along(theta.hat) > length(rout$alpha) + length(rout$b)
 is.bee <- (! (is.alpha | is.nu))
 is.zero.nu <- rout$nu == 0.0
 nrand <- sapply(rout$random, ncol)
 is.zero.bee <- rep(is.zero.nu, times = nrand)
 is.zero.alpha <- rep(FALSE, sum(is.alpha))
 is.zero <- c(is.zero.alpha, is.zero.bee, is.zero.nu)

 is.nonzero.alpha <- (! is.zero) & is.alpha
 is.nonzero.bee <- (! is.zero) & is.bee
 is.nonzero.nu <- (! is.zero) & is.nu

 pW_alpha_alpha <- oout$hessian[is.nonzero.alpha, is.nonzero.alpha]
 pW_alpha_bee <- oout$hessian[is.nonzero.alpha, is.nonzero.bee]
 pW_alpha_nu <- oout$hessian[is.nonzero.alpha, is.nonzero.nu]
 pW_bee_alpha <- oout$hessian[is.nonzero.bee, is.nonzero.alpha]
 pW_bee_bee <- oout$hessian[is.nonzero.bee, is.nonzero.bee]
 pW_bee_nu <- oout$hessian[is.nonzero.bee, is.nonzero.nu]
 pW_nu_alpha <- oout$hessian[is.nonzero.nu, is.nonzero.alpha]
 pW_nu_bee <- oout$hessian[is.nonzero.nu, is.nonzero.bee]
 pW_nu_nu <- oout$hessian[is.nonzero.nu, is.nonzero.nu]

 pW_bee_bee_inv <- solve(pW_bee_bee)
 # equation (4.16) of the Aster Theory book
 b_star_alpha <- (-1) * pW_bee_bee_inv %*% pW_bee_alpha
 # equation (4.17) of the Aster Theory book
 b_star_nu <- (-1) * pW_bee_bee_inv %*% pW_bee_nu
 
 # equation (4.19a) of the Aster Theory book
 qW_alpha_alpha <- pW_alpha_alpha -
     pW_alpha_bee %*% pW_bee_bee_inv %*% pW_bee_alpha
 # equation (4.19b) of the Aster Theory book
 qW_alpha_nu <- pW_alpha_nu -
     pW_alpha_bee %*% pW_bee_bee_inv %*% pW_bee_nu
 # equation (4.19c) of the Aster Theory book
 qW_nu_nu <- pW_nu_nu -
     pW_nu_bee %*% pW_bee_bee_inv %*% pW_bee_nu
 qW_nu_alpha <- t(qW_alpha_nu)

 qW_psi_psi <- rbind(
     cbind(qW_alpha_alpha, qW_alpha_nu),
     cbind(qW_nu_alpha, qW_nu_nu)
 )
 qW_psi_psi_inv <- solve(qW_psi_psi)

 # equation (4.21) of the Aster Theory book
 vee <- matrix(0, nrow = length(theta.hat),
     ncol = sum((! is.zero) & (! is.bee)))
 stopifnot(ncol(vee) == ncol(qW_psi_psi_inv))
 vee.col.is.alpha <- seq(1, ncol(vee)) <= sum(is.alpha)
 vee.col.is.nu <- (! vee.col.is.alpha)
 vee[is.alpha, vee.col.is.alpha] <- diag(sum(is.alpha))
 vee[is.nonzero.bee, vee.col.is.alpha] <- b_star_alpha
 vee[is.nonzero.bee, vee.col.is.nu] <- b_star_nu
 vee[is.nonzero.nu, vee.col.is.nu] <- diag(sum(vee.col.is.nu))
 # text following equation (4.21) of the Aster Theory book
 foo <- vee %*% qW_psi_psi_inv %*% t(vee)

 # re.too == TRUE and complete == TRUE and standard.deviation == FALSE

 vout <- vcov(rout, re.too = TRUE, complete = TRUE,
     standard.deviation = standard.deviation)
 print(all.equal(vout, foo, check.attributes = FALSE))

 izzy <- attributes(vout)
 izzy <- izzy[grepl("^is\\.", names(izzy))]
 print(all(sapply(izzy, length) == nrow(vout)))
 print(all(Reduce("+", izzy) == 1))

 vout <- vcov(rout, re.too = TRUE, complete = FALSE,
     standard.deviation = standard.deviation)
 foo <- foo[! is.zero, ! is.zero]
 print(all.equal(vout, foo, check.attributes = FALSE))

 izzy <- attributes(vout)
 izzy <- izzy[grepl("^is\\.", names(izzy))]
 print(all(sapply(izzy, length) == nrow(vout)))
 print(all(Reduce("+", izzy) == 1))

 }

 doit(standard.deviation = FALSE)
 doit(standard.deviation = TRUE)

