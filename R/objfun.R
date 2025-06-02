
# objective function factory for reaster

objfun.factory <- function(fixed, random, response, pred, fam, root,
    zwz, famlist = fam.default(), offset, standard.deviation = FALSE,
    deriv = 0)
{
    stopifnot(is.matrix(fixed))
    stopifnot(is.matrix(random) || is.list(random))
    if (is.matrix(random)) random <- list(random)
    stopifnot(sapply(random, is.matrix))
    stopifnot(sapply(random, nrow) == nrow(fixed))
    stopifnot(sapply(random, ncol) > 0)
    nfix <- ncol(fixed)
    nrand <- sapply(random, ncol)
    modmat <- Reduce(cbind, random, init = fixed)
    stopifnot(is.numeric(modmat))
    stopifnot(is.finite(modmat))
    storage.mode(modmat) <- "double"
    stopifnot(is.numeric(response))
    stopifnot(is.finite(response))
    stopifnot(length(response) == nrow(fixed))
    storage.mode(response) <- "double"
    stopifnot(is.numeric(pred))
    stopifnot(pred == round(pred))
    stopifnot(pred >= 0)
    stopifnot(pred <= seq_along(pred))
    storage.mode(pred) <- "integer"
    nnode <- length(pred)
    nind <- length(response) %/% nnode
    stopifnot(length(response) == nind * nnode)
    stopifnot(length(fam) == length(pred))
    stopifnot(fam %in% seq_along(famlist))
    storage.mode(fam) <- "integer"
    stopifnot(is.numeric(root))
    stopifnot(is.finite(root))
    stopifnot(length(root) == length(response))
    stopifnot(root == round(root))
    stopifnot(root > 0)
    if (missing(offset)) {
        # copy from R function aster
        setfam(famlist)
        offset <- .C(C_aster_default_origin,
            nind = as.integer(nind),
            nnode = as.integer(nnode),
            fam = as.integer(fam),
            theta = matrix(as.double(0), nind, nnode))$theta
        offset <- .C(C_aster_theta2phi,
            nind = as.integer(nind),
            nnode = as.integer(nnode),
            pred = as.integer(pred),
            fam = as.integer(fam),
            theta = offset,
            phi = matrix(as.double(0), nind, nnode))$phi
        clearfam()
    } else {
        stopifnot(is.numeric(offset))
        stopifnot(is.finite(offset))
        stopifnot(length(offset) == length(response))
    }
    storage.mode(offset) <- "double"
    stopifnot(is.logical(standard.deviation))
    stopifnot(length(standard.deviation) == 1)

    stopifnot(is.matrix(zwz))
    if (any(dim(zwz) != sum(nrand)))
        stop("dimension of zwz must be number of random effects")
    stopifnot(isSymmetric(zwz))
    stopifnot(is.matrix(zwz))
    stopifnot(is.numeric(zwz))
    stopifnot(is.finite(zwz))
    storage.mode(zwz) <- "double"
    # should also check positive definite and symmetric, but do not

    storage.mode(deriv) <- "integer"
    stopifnot(deriv %in% 0:2)

    function(theta) {
        setfam(famlist)
        on.exit(clearfam())
        .Call(C_objfun, theta, modmat, nfix, nrand,
            response, pred, fam, root, zwz,
            offset, standard.deviation, deriv)
    }
}

