
 set.seed(42)

 library(aster)

 a <- matrix(rnorm(1000), 5)
 a <- apply(a, 1, cumsum)
 a <- var(a)

 b <- .Call(aster:::C_pos_def_mat_inv, a)

 all.equal(b %*% a, diag(nrow(a)))

