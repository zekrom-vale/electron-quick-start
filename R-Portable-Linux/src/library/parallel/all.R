#  File src/library/parallel/R/RngStream.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2017 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

nextRNGStream <- function(seed)
{
    if(!is.integer(seed) || seed[1L] %% 100L != 7L)
	stop(gettextf("invalid value of %s", "'seed'"), domain = NA)
    .Call(C_nextStream, seed)
}

nextRNGSubStream <- function(seed)
{
    if(!is.integer(seed) || seed[1L] %% 100L != 7L)
	stop(gettextf("invalid value of %s", "'seed'"), domain = NA)
    .Call(C_nextSubStream, seed)
}

## Different from snow's RNG code
clusterSetRNGStream <- function(cl = NULL, iseed = NULL)
{
    cl <- defaultCluster(cl)
    oldseed <-
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        else NULL
    RNGkind("L'Ecuyer-CMRG")
    if(!is.null(iseed)) set.seed(iseed)
    nc <- length(cl)
    seeds <- vector("list", nc)
    seeds[[1L]] <- .Random.seed
    for(i in seq_len(nc-1L)) seeds[[i+1L]] <- nextRNGStream(seeds[[i]])
    ## Reset the random seed in the master.
    if(!is.null(oldseed))
        assign(".Random.seed", oldseed, envir = .GlobalEnv)
    else rm(.Random.seed, envir = .GlobalEnv)
    for (i in seq_along(cl)) {
        expr <- substitute(assign(".Random.seed", seed, envir = .GlobalEnv),
                           list(seed = seeds[[i]]))
        sendCall(cl[[i]], eval, list(expr))
    }
    checkForRemoteErrors(lapply(cl, recvResult))
    invisible()
}

RNGenv <- new.env()

mc.reset.stream <- function() {
    if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            sample.int(1L)
        assign("LEcuyer.seed",
               get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
               envir = RNGenv)
    }
}

## For use in the master before forking
mc.advance.stream <- function(reset = FALSE)
{
    if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
        if (reset ||
            !exists("LEcuyer.seed", envir = RNGenv, inherits = FALSE)) {
            if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                sample.int(1L)
            assign("LEcuyer.seed",
                   get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
                   envir = RNGenv)
        } else {
            assign("LEcuyer.seed",
                   nextRNGStream(get("LEcuyer.seed", envir = RNGenv)),
                   envir = RNGenv)
        }
    }
}

## For use in the child
mc.set.stream <- function()
{
    if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
            assign(".Random.seed", get("LEcuyer.seed", envir = RNGenv),
                   envir = .GlobalEnv)
    } else {
        ## It is random to simply unset the seed
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
}
#  File src/library/parallel/R/clusterApply.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2020 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

## Derived from snow 0.3-6 by Luke Tierney

staticClusterApply <- function(cl = NULL, fun, n, argfun) {
    cl <- defaultCluster(cl)
    p <- length(cl)
    if (n > 0L && p) {
        val <- vector("list", n)
        start <- 1L
        while (start <= n) {
            end <- min(n, start + p - 1L)
	    jobs <- end - start + 1L
            for (i in 1:jobs)
                sendCall(cl[[i]], fun, argfun(start + i - 1L))
            val[start:end] <- lapply(cl[1:jobs], recvResult)
            start <- start + jobs
        }
        checkForRemoteErrors(val)
    }
}

dynamicClusterApply <- function(cl = NULL, fun, n, argfun) {
    cl <- defaultCluster(cl)
    p <- length(cl)
    if (n > 0L && p) {
        submit <- function(node, job)
            sendCall(cl[[node]], fun, argfun(job), tag = job)
        for (i in 1:min(n, p)) submit(i, i)
        val <- vector("list", n)
        for (i in 1:n) {
            d <- recvOneResult(cl)
            j <- i + min(n, p)
            if (j <= n) submit(d$node, j)
            val[d$tag] <- list(d$value)
        }
        checkForRemoteErrors(val)
    }
}

## exported and documented from here down unless otherwise stated.

clusterCall  <- function(cl = NULL, fun, ...)
{
    cl <- defaultCluster(cl)
    for (i in seq_along(cl)) sendCall(cl[[i]], fun, list(...))
    checkForRemoteErrors(lapply(cl, recvResult))
}


clusterEvalQ <- function(cl = NULL, expr)
    clusterCall(cl, eval, substitute(expr), env=.GlobalEnv)

clusterExport <- local({
    gets <- function(n, v) { assign(n, v, envir = .GlobalEnv); NULL }
    function(cl = NULL, varlist, envir = .GlobalEnv) {
        ## do this with only one clusterCall--loop on workers?
        for (name in varlist) {
            clusterCall(cl, gets, name, get(name, envir = envir))
        }
    }
})

clusterApply <- function(cl = NULL, x, fun, ...)
{
    ## **** this closure is sending all of x to all nodes
    argfun <- function(i) c(list(x[[i]]), list(...))
    staticClusterApply(cl, fun, length(x), argfun)
}

clusterApplyLB <- function(cl = NULL, x, fun, ...)
{
    ## **** this closure is sending all of x to all nodes
    argfun <- function(i) c(list(x[[i]]), list(...))
    dynamicClusterApply(cl, fun, length(x), argfun)
}

clusterMap <- function (cl = NULL, fun, ..., MoreArgs = NULL, RECYCLE = TRUE,
                        SIMPLIFY = FALSE, USE.NAMES = TRUE,
                        .scheduling = c("static", "dynamic"))
{
    cl <- defaultCluster(cl)
    args <- list(...)
    if (length(args) == 0) stop("need at least one argument")
    .scheduling <- match.arg(.scheduling)
    n <- lengths(args)
    if (RECYCLE) {
        vlen <- max(n)
        if(vlen && min(n) == 0L)
            stop("zero-length inputs cannot be mixed with those of non-zero length")
        if (!all(n == vlen))
            for (i in seq_along(args)) # why not lapply?
                args[[i]] <- rep(args[[i]], length.out = vlen)
    }
    else vlen <- min(n)
    ## **** this closure is sending all of ... to all nodes
    argfun <- function(i) c(lapply(args, function(x) x[[i]]), MoreArgs)
    answer <-
        if(.scheduling == "dynamic") dynamicClusterApply(cl, fun, vlen, argfun)
    else staticClusterApply(cl, fun, vlen, argfun)
    ## rest matches mapply(): with a different default for SIMPLIFY
    if (USE.NAMES && length(args)) {
        if (is.null(names1 <- names(args[[1L]])) && is.character(args[[1L]]))
            names(answer) <- args[[1L]]
        else if (!is.null(names1))
            names(answer) <- names1
    }
    if (!isFALSE(SIMPLIFY))
        simplify2array(answer, higher = (SIMPLIFY == "array"))
    else answer
}

## splitIndices <- function(nx, ncl)
## {
##     i <- seq_len(nx)
##     if (ncl == 1L) i
##     else structure(split(i, cut(i, ncl)), names = NULL)
## }

# The fuzz used by cut() is too small when nx and ncl are both large
# and causes some groups to be empty. The definition below avoids that
# while minimizing changes from the results produced by the definition
# above.
splitIndices <- function(nx, ncl) {
    i <- seq_len(nx)
    if (ncl == 0L) list()
    else if (ncl == 1L || nx == 1L) list(i)
    else {
        fuzz <- min((nx - 1L) / 1000, 0.4 * nx / ncl)
        breaks <- seq(1 - fuzz, nx + fuzz, length.out = ncl + 1L)
        structure(split(i, cut(i, breaks)), names = NULL)
    }
}

clusterSplit <- function(cl = NULL, seq) {
    cl <- defaultCluster(cl)
    lapply(splitIndices(length(seq), length(cl)), function(i) seq[i])
}

#internal
splitList <- function(x, ncl)
    lapply(splitIndices(length(x), ncl), function(i) x[i])

#internal
splitRows <- function(x, ncl)
    lapply(splitIndices(nrow(x), ncl), function(i) x[i, , drop=FALSE])

#internal
splitCols <- function(x, ncl)
    lapply(splitIndices(ncol(x), ncl), function(i) x[, i, drop=FALSE])

#internal
staticNChunks <- function(nelems, nnodes, chunk.size) {
    if (is.null(chunk.size) || chunk.size <= 0)
        nnodes
    else
        max(1, ceiling(nelems / chunk.size))
}

#internal
dynamicNChunks <- function(nelems, nnodes, chunk.size) {
    if (is.null(chunk.size))
        2 * nnodes
    else if (chunk.size <= 0)
        nelems
    else
        max(1, ceiling(nelems / chunk.size))
}

parLapply <- function(cl = NULL, X, fun, ..., chunk.size = NULL)
{
    cl <- defaultCluster(cl)
    nchunks <- staticNChunks(length(X), length(cl), chunk.size)
    do.call(c,
            clusterApply(cl = cl, x = splitList(X, nchunks),
                         fun = lapply, FUN = fun, ...),
            quote = TRUE)
}

parLapplyLB <- function(cl = NULL, X, fun, ..., chunk.size = NULL)
{
    cl <- defaultCluster(cl)
    nchunks <- dynamicNChunks(length(X), length(cl), chunk.size)
    do.call(c,
            clusterApplyLB(cl = cl, x = splitList(X, nchunks),
                           fun = lapply, FUN = fun, ...),
            quote = TRUE)
}

parRapply <- function(cl = NULL, x, FUN, ..., chunk.size = NULL)
{
    cl <- defaultCluster(cl)
    nchunks <- staticNChunks(nrow(x), length(cl), chunk.size)
    do.call(c,
            clusterApply(cl = cl, x = splitRows(x, nchunks),
                         fun = apply, MARGIN = 1L, FUN = FUN, ...),
            quote = TRUE)
}

parCapply <- function(cl = NULL, x, FUN, ..., chunk.size = NULL) {
    cl <- defaultCluster(cl)
    nchunks <- staticNChunks(ncol(x), length(cl), chunk.size)
    do.call(c,
            clusterApply(cl = cl, x = splitCols(x, nchunks),
                         fun = apply, MARGIN = 2L, FUN = FUN, ...),
            quote = TRUE)
}


parSapply <-
    function (cl = NULL, X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE,
              chunk.size = NULL)
{
    FUN <- match.fun(FUN) # should this be done on worker?
    answer <- parLapply(cl = cl, X = as.list(X), fun = FUN, ...,
                        chunk.size = chunk.size)
    if(USE.NAMES && is.character(X) && is.null(names(answer)))
	names(answer) <- X
    if(!isFALSE(simplify))
	simplify2array(answer, higher = (simplify == "array"))
    else answer
}

parSapplyLB <-
    function (cl = NULL, X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE,
              chunk.size = NULL)
{
    FUN <- match.fun(FUN) # should this be done on worker?
    answer <- parLapplyLB(cl = cl, X = as.list(X), fun = FUN, ...,
                          chunk.size = chunk.size)
    if(USE.NAMES && is.character(X) && is.null(names(answer)))
	names(answer) <- X
    if(!isFALSE(simplify))
	simplify2array(answer, higher = (simplify == "array"))
    else answer
}


parApply <- function(cl = NULL, X, MARGIN, FUN, ..., chunk.size = NULL)
{
    cl <- defaultCluster(cl) # initial sanity check
    FUN <- match.fun(FUN) # should this be done on worker?

    ## Ensure that X is an array object
    dl <- length(dim(X))
    if(!dl) stop("dim(X) must have a positive length")
    if(is.object(X))
	X <- if(dl == 2L) as.matrix(X) else as.array(X)
    ## now record dim as coercion can change it
    ## (e.g. when a data frame contains a matrix).
    d <- dim(X)
    dn <- dimnames(X)
    ds <- seq_len(dl)

    ## Extract the margins and associated dimnames

    if (is.character(MARGIN)) {
        if(is.null(dnn <- names(dn))) # names(NULL) is NULL
           stop("'X' must have named dimnames")
        MARGIN <- match(MARGIN, dnn)
        if (anyNA(MARGIN))
            stop("not all elements of 'MARGIN' are names of dimensions")
    }
    s.call <- ds[-MARGIN]
    s.ans  <- ds[MARGIN]
    d.call <- d[-MARGIN]
    d.ans  <- d[MARGIN]
    dn.call <- dn[-MARGIN]
    dn.ans <- dn[MARGIN]
    ## dimnames(X) <- NULL

    ## do the calls

    d2 <- prod(d.ans)
    if(d2 == 0L) {
        ## arrays with some 0 extents: return ``empty result'' trying
        ## to use proper mode and dimension:
        ## The following is still a bit `hackish': use non-empty X
        newX <- array(vector(typeof(X), 1L), dim = c(prod(d.call), 1L))
        ans <- FUN(if(length(d.call) < 2L) newX[,1] else
                   array(newX[, 1L], d.call, dn.call), ...)
        return(if(is.null(ans)) ans else if(length(d.ans) < 2L) ans[1L][-1L]
               else array(ans, d.ans, dn.ans))
    }
    ## else
    newX <- aperm(X, c(s.call, s.ans))
    dim(newX) <- c(prod(d.call), d2)
    ans <- vector("list", d2)
    arglist <- if(length(d.call) < 2L) {# vector
        if (length(dn.call)) dimnames(newX) <- c(dn.call, list(NULL))
        lapply(seq_len(d2), function(i) newX[,i])
    } else
        lapply(seq_len(d2), function(i) array(newX[,i], d.call, dn.call))
    ans <- parLapply(cl = cl, X = arglist, fun = FUN, ...,
                     chunk.size = chunk.size)

    ## answer dims and dimnames

    ans.list <- is.recursive(ans[[1L]])
    l.ans <- length(ans[[1L]])

    ans.names <- names(ans[[1L]])
    if(!ans.list)
	ans.list <- any(lengths(ans) != l.ans)
    if(!ans.list && length(ans.names)) {
        all.same <- vapply(ans, function(x) identical(names(x), ans.names), NA)
        if (!all(all.same)) ans.names <- NULL
    }
    len.a <- if(ans.list) d2 else length(ans <- unlist(ans, recursive = FALSE))
    if(length(MARGIN) == 1L && len.a == d2) {
	names(ans) <- if(length(dn.ans[[1L]])) dn.ans[[1L]] # else NULL
	return(ans)
    }
    if(len.a == d2)
	return(array(ans, d.ans, dn.ans))
    if(len.a && len.a %% d2 == 0L) {
        if(is.null(dn.ans)) dn.ans <- vector(mode="list", length(d.ans))
        dn.ans <- c(list(ans.names), dn.ans)
	return(array(ans, c(len.a %/% d2, d.ans),
		     if(!all(vapply(dn.ans, is.null, NA))) dn.ans))
    }
    return(ans)
}

#  File src/library/parallel/R/detectCores.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2020 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

## In part based on code in package multicore 0.1-6 by Simon Urbanek

detectCores <-
    if(.Platform$OS.type == "windows") {
        function(all.tests = FALSE, logical = TRUE) {
            ## result is # cores, logical processors.
            res <- .Call(C_ncpus, FALSE)
	    res[if(logical) 2L else 1L]
        }
    } else {
        function(all.tests = FALSE, logical = TRUE) {
            ## Commoner OSes first
            ## for Linux systems, physical id is 1 for second hyperthread
            ## Irix support removed in R 4.1.0
            systems <-
                ## quoting needed for a Bourne shell
                list(linux = 'grep "^processor" /proc/cpuinfo 2>/dev/null | wc -l',
                     ## hw.physicalcpu is not documented for 10.9, but works
                     darwin = if(logical) "/usr/sbin/sysctl -n hw.logicalcpu 2>/dev/null" else "/usr/sbin/sysctl -n hw.physicalcpu 2>/dev/null",
                     solaris = if(logical) "/usr/sbin/psrinfo -v | grep 'Status of.*processor' | wc -l" else "/bin/kstat -p -m cpu_info | grep :core_id | cut -f2 | uniq | wc -l",
                     freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null",
                     openbsd = "/sbin/sysctl -n hw.ncpuonline 2>/dev/null")
            nm <- names(systems)
            m <- pmatch(nm, R.version$os); m <- nm[!is.na(m)]
            if (length(m)) {
                cmd <- systems[[m]]
                if(!is.null(a <- tryCatch(suppressWarnings(system(cmd, TRUE)),
                                         error = function(e) NULL))) {
                    a <- gsub("^ +","", a[1])
                    if (grepl("^[1-9]", a)) return(as.integer(a))
                }
            }
            if (all.tests) {
                for (i in seq(systems))
                    for (cmd in systems[i]) { # Irix had two commands
			if(is.null(a <- tryCatch(suppressWarnings(system(cmd, TRUE)),
						 error = function(e) NULL)))
			    next
                        a <- gsub("^ +","", a[1])
                        if (grepl("^[1-9]", a)) return(as.integer(a))
                    }
            }
            NA_integer_
        }
    }
## added in R 3.0.3
.check_ncores <- function(nc)
{
    chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
    if (nzchar(chk) && (chk != "false") && nc > 2L) {
        msg <- sprintf("%d simultaneous processes spawned", nc)
        if(chk == "warn") warning(msg, call. = FALSE, immediate. = TRUE)
        else stop(msg, call. = TRUE)
    }
}
#  File src/library/parallel/R/snow.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2021 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

## Derived from snow 0.3-6 by Luke Tierney

.reg <-  new.env()
assign("default", NULL, envir = .reg)

defaultCluster <- function(cl = NULL)
{
    if(is.null(cl)) cl <- get("default", envir = .reg)
    if(is.null(cl)) stop("no cluster 'cl' supplied and none is registered")
    checkCluster(cl)
    cl
}

setDefaultCluster <- function(cl = NULL)
{
    if(!is.null(cl)) checkCluster(cl)
    assign("default", cl, envir = .reg)
}

getDefaultCluster <-
function()
    get("default", envir = .reg)

#
# Checking and subsetting
#

checkCluster <- function(cl)
    if (!inherits(cl, "cluster")) stop("not a valid cluster");

`[.cluster` <- function(cl, ...) {
    v <- NextMethod()
    class(v) <- class(cl)
    v
}


#
# Higher-Level Node Functions
#

closeNode <- function(node) UseMethod("closeNode")
closeNode.default <- function(node) {}

## These have SOCK methods
sendData <- function(node, data) UseMethod("sendData")
recvData <- function(node) UseMethod("recvData")
recvOneData <- function(cl) UseMethod("recvOneData")

postNode <- function(con, type, value = NULL, tag = NULL)
    sendData(con, list(type = type, data = value, tag = tag))

stopNode <- function(n) {
    postNode(n, "DONE")
    closeNode(n)
}



#
#  Cluster Creation and Destruction
#

defaultClusterOptions <- NULL

#**** check valid cluster option

initDefaultClusterOptions <- function(libname)
{
    rscript <- file.path(R.home("bin"), "Rscript")
    port <- Sys.getenv("R_PARALLEL_PORT")
    port <- if (identical(port, "random")) NA else as.integer(port)
    if (is.na(port)) {
	seed <- .GlobalEnv$.Random.seed
        ran1 <- sample.int(.Machine$integer.max - 1L, 1L) / .Machine$integer.max
        port <- 11000 + 1000 * ((ran1 + unclass(Sys.time()) / 300) %% 1)
	if(is.null(seed)) ## there was none, initially
	    rm(    ".Random.seed",       envir = .GlobalEnv, inherits = FALSE)
	else # reset
	    assign(".Random.seed", seed, envir = .GlobalEnv, inherits = FALSE)
    }
    Sys.i <- Sys.info()
    options <- list(port = as.integer(port),
                    setup_timeout = 60 * 2,      # 2 minutes
                    timeout = 60 * 60 * 24 * 30, # 30 days
                    setup_strategy = "parallel",
                    master = Sys.i[["nodename"]],
                    homogeneous = TRUE,
                    type = "PSOCK",
                    outfile = "/dev/null",
                    rscript = rscript,
                    rscript_args = character(),
                    user = NULL,
                    rshcmd = "ssh",
                    manual = FALSE,
                    methods = TRUE,
                    renice = NA_integer_,
                    ## rest are unused in parallel
                    rhome = R.home(),
                    rlibs = Sys.getenv("R_LIBS"),
                    scriptdir = file.path(libname, "parallel"),
                    rprog = file.path(R.home("bin"), "R"),
                    snowlib = .libPaths()[1],
                    useRscript = TRUE, # for use by snow clusters
                    useXDR = TRUE)
    defaultClusterOptions <<- addClusterOptions(emptyenv(), options)
}

addClusterOptions <- function(options, new) {
    if (!is.null(new)) {
        options <- new.env(parent = options)
        names <- names(new)
        for (i in seq_along(new))
            assign(names[i], new[[i]], envir = options)
    }
    options
}

getClusterOption <- function(name, options = defaultClusterOptions)
    get(name, envir = options)

setDefaultClusterOptions <- function(...) {
    list <- list(...)
    names <- names(list)
    for (i in seq_along(list))
        assign(names[i], list[[i]], envir = defaultClusterOptions)
}


makeCluster <-
    function (spec, type = getClusterOption("type"), ...)
{
    switch(type,
           PSOCK = makePSOCKcluster(names = spec, ...),
           FORK = makeForkCluster(nnodes = spec, ...),
           SOCK = snow::makeSOCKcluster(names = spec, ...),
           MPI = snow::makeMPIcluster(count = spec, ...),
           ## NWS = snow::makeNWScluster(names = spec, ...),
           stop("unknown cluster type"))
}


stopCluster <- function(cl = NULL)
{
    cl <- defaultCluster(cl)
    if(identical(cl, get("default", envir = .reg)))
        assign("default", NULL, envir = .reg)
    UseMethod("stopCluster")
}

stopCluster.default <- function(cl) for (n in cl) stopNode(n)


#
# Cluster Functions
#

sendCall <- function (con, fun, args, return = TRUE, tag = NULL)
{
    timing <-  .snowTimingData$running()
    if (timing)
        start <- proc.time()[3L]
    postNode(con, "EXEC",
             list(fun = fun, args = args, return = return, tag = tag))
    if (timing)
        .snowTimingData$enterSend(con$rank, start, proc.time()[3L])
    NULL
}

recvResult <- function(con)
{
    if (.snowTimingData$running()) {
        start <- proc.time()[3L]
        r <- recvData(con)
        end <- proc.time()[3L]
        .snowTimingData$enterRecv(con$rank, start, end, r$time[3L])
    }
    else r <- recvData(con)
    r$value
}

checkForRemoteErrors <- function(val)
{
    count <- 0
    firstmsg <- NULL
    for (v in val) {
        if (inherits(v, "try-error")) {
            count <- count + 1
            if (count == 1) firstmsg <- v
        }
    }
    ## These will not translate
    if (count == 1)
        stop("one node produced an error: ", firstmsg, domain = NA)
    else if (count > 1)
        stop(count, " nodes produced errors; first error: ", firstmsg, domain = NA)
    val
}

recvOneResult <- function (cl) {
    if (.snowTimingData$running()) {
        start <- proc.time()[3]
        v <- recvOneData(cl)
        end <- proc.time()[3]
        .snowTimingData$enterRecv(v$node, start, end, v$value$time[3])
    }
    else v <- recvOneData(cl)
    list(value = v$value$value, node = v$node, tag = v$value$tag)
}

findRecvOneTag <- function(cl, anytag) {
    rtag <- NULL
    for (node in cl) {
        if (is.null(rtag))
            rtag <- node$RECVTAG
        else if (rtag != node$RECVTAG) {
            rtag <- anytag
            break;
        }
    }
    rtag
}

### ========== snow support ===========

## place holder for now.
.snowTimingData <-
    list(running = function() FALSE,
         enterSend = function(...) {},
         enterRecv = function(...) {})

### NWS support no longer exists in snow

##closeNode.NWSnode <- function(node) snow::closeNode.NWSnode(node)

recvData.MPInode <- function(node) snow::recvData.MPInode(node)
##recvData.NWSnode <- function(node) snow::recvData.NWSnode(node)

recvOneData.MPIcluster <- function(cl) snow::recvOneData.MPIcluster(cl)
##recvOneData.NWScluster <- function(cl) snow::recvOneData.NWScluster(cl)

sendData.MPInode <- function(node, data) snow::sendData.MPInode(node, data)
##sendData.NWSnode <- function(node, data) snow::sendData.NWSnode(node, data)

## these use NextMethod() so need copies.
stopCluster.MPIcluster <- function(cl) {
    NextMethod()
    snow::setMPIcluster(NULL)
}

stopCluster.spawnedMPIcluster <- function(cl) {
    comm <- 1
    NextMethod()
    Rmpi::mpi.comm.disconnect(comm)
}

## stopCluster.NWScluster <- function(cl) {
##     NextMethod()
##     nws::nwsDeleteWs(cl[[1]]$wsServer, nws::nwsWsName(cl[[1]]$ws))
##     close(cl[[1]]$wsServer)
## }


#  File src/library/parallel/R/snowSOCK.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2020 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

## Derived from snow 0.3-6 by Luke Tierney
## Uses solely Rscript, and a function in the package rather than scripts.

## NB: there is also workCommand in worker.R: this one is run on the master
workerCommand <- function(machine, options, setup_strategy = "sequential")
{
    outfile <- getClusterOption("outfile", options)
    master <- if (machine == "localhost") "localhost"
              else getClusterOption("master", options)
    port <- getClusterOption("port", options)
    setup_timeout <- getClusterOption("setup_timeout", options)
    manual <- getClusterOption("manual", options)
    timeout <- getClusterOption("timeout", options)
    methods <- getClusterOption("methods", options)
    useXDR <- getClusterOption("useXDR", options)
    homogeneous <- getClusterOption("homogeneous", options)

    ## build the local command for starting the worker
    env <- paste0("MASTER=", master,
                 " PORT=", port,
                 " OUT=", shQuote(outfile),
                 " SETUPTIMEOUT=", setup_timeout,
                 " TIMEOUT=", timeout,
                 " XDR=", useXDR,
                 " SETUPSTRATEGY=", setup_strategy)
    ## Should cmd be run on a worker with R <= 4.0.2,
    ## .workRSOCK will not exist, so fallback to .slaveRSOCK
    arg <- "tryCatch(parallel:::.workRSOCK,error=function(e)parallel:::.slaveRSOCK)()"
    ## option rscript got set by initDefaultClusterOptions to the full path
    ## on the master, but can be overridden in the makePSOCKcluster call.
    rscript <-
        if (homogeneous) shQuote(getClusterOption("rscript", options)) else "Rscript"
    rscript_args <- getClusterOption("rscript_args", options)
    if(methods)
        rscript_args <-c("--default-packages=datasets,utils,grDevices,graphics,stats,methods",
                         rscript_args)

    ## in principle we should quote these,
    ## but the current possible values do not need quoting
    cmd <- paste(rscript,
                 if(length(rscript_args)) paste(rscript_args, collapse = " "),
                 "-e", shQuote(arg), env)

    ## We do redirection of connections at R level once the process is
    ## running.  We could instead do it at C level here, at least on
    ## a Unix-alike.
    renice <- getClusterOption("renice", options)
    if(!is.na(renice) && renice) ## ignore 0
        cmd <- sprintf("nice +%d %s", as.integer(renice), cmd)

    ## add the remote shell command if needed
    if (!manual && machine != "localhost") {
        ## This assumes an ssh-like command
        rshcmd <- getClusterOption("rshcmd", options)
        user <- getClusterOption("user", options)
        ## this assume that rshcmd will use a shell, and that is
        ## the same shell as on the master.
        cmd <- paste(rshcmd,
                     if(length(user) == 1L) paste("-l", user),
                     machine, shQuote(cmd))
    }
    cmd
}

newPSOCKnode <- function(machine = "localhost", ...,
                         options = defaultClusterOptions, rank)
{
    options <- addClusterOptions(options, list(...))
    if (is.list(machine)) {
        options <- addClusterOptions(options, machine)
        machine <- machine$host
    }
    port <- getClusterOption("port", options)
    manual <- getClusterOption("manual", options)
    timeout <- getClusterOption("timeout", options)
    useXDR <- getClusterOption("useXDR", options)
    cmd <- workerCommand(machine, options)

    if (manual) {
        cat("Manually start worker on", machine, "with\n    ", cmd, "\n")
        utils::flush.console()
    } else {
        if (.Platform$OS.type == "windows") {
            ## snow said:
            ## On Windows using input = something seems needed to
            ## disconnect standard input of an ssh process when run
            ## from Rterm (at least using putty's plink).  In
            ## principle this could also be used for supplying a
            ## password, but that is probably a bad idea. So, for now
            ## at least, on Windows password-less authentication is
            ## necessary.
            ##
            ## (Not clear if that is the current behaviour: works for me)
            system(cmd, wait = FALSE, input = "")
        }
        else {
            ## If workers are running a different R version, avoid a WARNING
            cmd <- paste("R_HOME=", cmd)
            system(cmd, wait = FALSE)
        }
    }

    con <- socketConnection("localhost", port = port, server = TRUE,
                            blocking = TRUE, open = "a+b", timeout = timeout)
    structure(list(con = con, host = machine, rank = rank),
              class = if(useXDR) "SOCKnode" else "SOCK0node")
}

closeNode.SOCKnode <- closeNode.SOCK0node <- function(node) close(node$con)

sendData.SOCKnode <- function(node, data) serialize(data, node$con)
sendData.SOCK0node <- function(node, data) serialize(data, node$con, xdr = FALSE)

recvData.SOCKnode <- recvData.SOCK0node <- function(node) unserialize(node$con)

recvOneData.SOCKcluster <- function(cl)
{
    socklist <- lapply(cl, function(x) x$con)
    repeat {
        ready <- socketSelect(socklist)
        if (length(ready) > 0) break;
    }
    n <- which.max(ready) # may need rotation or some such for fairness
    list(node = n, value = unserialize(socklist[[n]]))
}

makePSOCKcluster <- function(names, ...)
{
    options <- addClusterOptions(defaultClusterOptions, list(...))
    manual <- getClusterOption("manual", options)
    homogeneous <- getClusterOption("homogeneous", options)
    setup_strategy <- match.arg(getClusterOption("setup_strategy",
                                                  options),
                                c("sequential", "parallel"))
    setup_timeout <- getClusterOption("setup_timeout", options)
    local <- is.numeric(names) || (is.character(names) &&
                 identical(names, rep('localhost', length(names))))
    if (is.numeric(names)) {
        names <- as.integer(names[1L])
        if(is.na(names) || names < 1L) stop("numeric 'names' must be >= 1")
        names <- rep('localhost', names)
    }
    .check_ncores(length(names))

    cl <- vector("list", length(names))
    if (!manual && homogeneous && local && setup_strategy == "parallel") {
        port <- getClusterOption("port", options)
        timeout <- getClusterOption("timeout", options)
        useXDR <- getClusterOption("useXDR", options)
        cmd <- workerCommand("localhost", options,
                             setup_strategy = "parallel" )

        ## Start listening and start workers.
        socket <- serverSocket(port = port)
        on.exit(close(socket), add = TRUE)
        if (.Platform$OS.type == "windows") {
            for(i in seq_along(cl))
                ## see newPSOCKnode for the input = ""
                system(cmd, wait = FALSE, input = "")
        } else {
            ## Asynchronous lists are defined by POSIX
            cmd <- paste(rep(cmd, length(cl)), collapse = " & ")
            system(cmd, wait = FALSE)
        }

        ## Accept connections and send the first command as initial
        ## handshake.  The handshake makes TCP synchronization detect and
        ## err on half-opened connections, which arise during parallel setup
        ## of client-server connections (due to internal timeouts, limited
        ## length of the listen backlog queue, race in timing out on
        ## creating a connection and probably more).
        ##
        ## The handshake looks like a regular server command followed by
        ## client response, which is compatible with older versions of R.

        cls <- if(useXDR) "SOCKnode" else "SOCK0node"
        ready <- 0
        pending <- list()
        on.exit(lapply(pending, function(x) close(x$con)), add = TRUE)
        t0 <- Sys.time()
        while (ready < length(cl)) {
            cons <- lapply(pending, function(x) x$con)

            if (difftime(Sys.time(), t0, units="secs") > setup_timeout + 5) {
                ## The workers will give up after setup_timeout, so there is
                ## no point waiting for them much longer.
                failed <- length(cl) - ready
                msg <- sprintf(ngettext(failed,
                           "Cluster setup failed. %d worker of %d failed to connect.",
                           "Cluster setup failed. %d of %d workers failed to connect."),
                               failed, length(cl))
                stop(msg)
            }
            a <- socketSelect(append(list(socket), cons), FALSE,
                              timeout = setup_timeout)
            canAccept <- a[1]
            canReceive <- seq_along(pending)[a[-1]]

            if (canAccept) {
                con <- socketAccept(socket = socket, blocking = TRUE,
                                    open = "a+b", timeout = timeout)
                scon <- structure(list(con = con, host = "localhost",
                                  rank = ready), class = cls)
                tryCatch({ sendCall(scon, eval, list(quote(Sys.getpid()))) },
                         error = identity)
                pending <- append(pending, list(scon))
            }
            for (scon in pending[canReceive]) {
                pid <- tryCatch({ recvResult(scon) }, error = identity)
                if (is.integer(pid)) {
                    ready <- ready + 1
                    cl[[ready]] <- scon
                } else
                    close(scon$con)
            }
            if (length(canReceive) > 0)
                pending <- pending[-canReceive]
        }
    } else {
        ## We could also use socketAccept() with sequential setup, but only
        ## if all workers connect to the same port.
        for (i in seq_along(cl))
            cl[[i]] <- newPSOCKnode(names[[i]], options = options, rank = i)
    }
    class(cl) <- c("SOCKcluster", "cluster")
    cl
}

print.SOCKcluster <- function(x, ...)
{
    nc <- length(x)
    hosts <- unique(sapply(x, `[[`, "host"))
    msg <- sprintf(ngettext(length(hosts),
                            "socket cluster with %d nodes on host %s",
                            "socket cluster with %d nodes on hosts %s"),
                   nc, paste(sQuote(hosts), collapse = ", "))
    cat(msg, "\n", sep = "")
    invisible(x)
}

print.SOCKnode <- print.SOCK0node <- function(x, ...)
{
    sendCall(x, eval, list(quote(Sys.getpid())))
    pid <- recvResult(x)

    msg <- gettextf("node of a socket cluster on host %s with pid %d",
                    sQuote(x[["host"]]), pid)
    cat(msg, "\n", sep = "")
    invisible(x)
}

.workRSOCK <- function()
{
    makeSOCKmaster <- function(master, port, setup_timeout, timeout, useXDR,
                               setup_strategy)
    {
        port <- as.integer(port)
        timeout <- as.integer(timeout)
        stopifnot(setup_timeout >= 0)
        cls <- if(useXDR) "SOCKnode" else "SOCK0node"

        ## Retry scheme parameters (do these need to be customizable?)
        retryDelay <- 0.1     # 0.1 second initial delay before retrying
        retryScale <- 1.5     # 50% increase of delay at each retry

        ## Retry multiple times in case the master is not yet ready
        t0 <- Sys.time()

        scon_timeout <- 1
        repeat {
            ## Set up a short timeout for the connection with parallel
            ## setup, which is needed to deal with half-opened connections
            ## (opened on client, closed on server).  The final connection
            ## timeout defaults to a large number, updated after the setup.
            if (setup_strategy == "parallel")
                scon_timeout <- scon_timeout + 0.2
            else
                ## Using "timeout" makes socketConnection() essentially
                ## blocking, which has been the practice for many years.
                ## Perhaps we could now use values similar to those for
                ## parallel setup.
                scon_timeout <- timeout

            con <- tryCatch({
                socketConnection(master, port = port, blocking = TRUE,
                                 open = "a+b",
                                 timeout = as.integer(scon_timeout))
            }, error = identity)

            hres <- NULL
            if (inherits(con, "connection")) {
                scon <- structure(list(con = con), class = cls)
                if (setup_strategy == "sequential")
                    return(scon)

                ## Serve the first command as a handshake during connection
                ## setup.  This is to get rid of half-opened connections.
                hres <- tryCatch({ workCommand(scon) }, error = identity)
                if (identical(hres, TRUE)) {
                    if (setup_strategy == "parallel")
                        socketTimeout(socket = con, timeout = timeout)
                    return(scon)
                } else if (identical(hres, FALSE)) {
                    ## "Done" command from server.  Could happen with server
                    ## accidentally not performing parallel setup.
                    return(NULL)
                } else
                    ## Error, possibly half-opened connection.
                    close(con)
            }

            if (difftime(Sys.time(), t0, units="secs") > setup_timeout) {
                if (inherits(hres, "error"))
                    stop(hres)
                if (inherits(con, "error"))
                    stop(con)
                stop("Connection setup failed or timed out.")
            }
            Sys.sleep(retryDelay)
            retryDelay <- retryScale * retryDelay
        }
    }

    ## set defaults in case run manually without args.
    master <- "localhost" # hostname of master process
    port <- NA_integer_   # no point in getting option on worker
    outfile <- Sys.getenv("R_SNOW_OUTFILE") # defaults to ""
    setup_timeout <- 120  # retry setup for 2 minutes before failing
    timeout <- 2592000L   # wait 30 days for new cmds before failing
    useXDR <- TRUE        # binary serialization
    setup_strategy <- "sequential"

    for (a in commandArgs(TRUE)) {
        ## Or use strsplit?
        pos <- regexpr("=", a)
        name <- substr(a, 1L, pos - 1L)
        value <- substr(a, pos + 1L, nchar(a))
        switch(name,
               MASTER = {master <- value},
               PORT = {port <- value},
               OUT = {outfile <- value},
               SETUPTIMEOUT = {setup_timeout <- as.numeric(value)},
               TIMEOUT = {timeout <- value},
               XDR = {useXDR <- as.logical(value)},
               SETUPSTRATEGY = {
                   setup_strategy <- match.arg(value,
                                               c("sequential", "parallel"))
               })
    }
    if (is.na(port)) stop("PORT must be specified")

    ## We should not need to attach parallel, as we are running in the namespace.

    sinkWorkerOutput(outfile)
    msg <- sprintf("starting worker pid=%d on %s at %s\n",
                   Sys.getpid(), paste(master, port, sep = ":"),
                   format(Sys.time(), "%H:%M:%OS3"))
    cat(msg)
    workLoop(makeSOCKmaster(master, port, setup_timeout, timeout, useXDR,
                            setup_strategy))
}
#  File src/library/parallel/R/unix/forkCluster.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2020 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

makeForkCluster <- function(nnodes = getOption("mc.cores", 2L), ...)
{
    options <- addClusterOptions(defaultClusterOptions, list(...))
    port <- getClusterOption("port", options)

    nnodes <- as.integer(nnodes)
    if(is.na(nnodes) || nnodes < 1L) stop("'nnodes' must be >= 1")
    .check_ncores(nnodes)
    cl <- vector("list", nnodes)
    socket <- serverSocket(port = port)
    on.exit(close(socket))
    for (i in seq_along(cl)) {
        node <- tryCatch( newForkNode(..., rank = i, socket = socket,
                                      server_socket = socket),
                          error = identity )
        if (inherits(node, "forknode"))
            cl[[i]] <- node
        else {
            if (i > 1)
                for(j in 1:i) stopNode(cl[[j]])
            stop("Cluster setup failed.")
        }
    }
    class(cl) <- c("SOCKcluster", "cluster")
    cl
}


newForkNode <- function(..., options = defaultClusterOptions, rank, socket,
                        setup_timeout, server_socket)
{
    options <- addClusterOptions(options, list(...))
    outfile <- getClusterOption("outfile", options)
    port <- getClusterOption("port", options)
    timeout <- getClusterOption("timeout", options)
    renice <- getClusterOption("renice", options)
    setup_timeout <- 10

    f <- mcfork(TRUE)
    if (inherits(f, "masterProcess")) { # the worker
        on.exit(mcexit(1L, structure("fatal error in wrapper code",
                                  class = "try-error")))
        # closeStdout()
        close(server_socket)
        master <- "localhost"
        makeSOCKmaster <- function(master, port, timeout)
        {
            port <- as.integer(port)
            con <- socketConnection(master, port = port, blocking = TRUE,
                                    open = "a+b", timeout = setup_timeout)
            socketTimeout(socket = con, timeout = timeout)
            structure(list(con = con), class = "SOCK0node")
        }
        sinkWorkerOutput(outfile)
        msg <- sprintf("starting worker pid=%d on %s at %s\n",
                       Sys.getpid(), paste(master, port, sep = ":"),
                       format(Sys.time(), "%H:%M:%OS3"))
        cat(msg)
        if(!is.na(renice) && renice) ## ignore 0
            tools::psnice(Sys.getpid(), renice)
        workLoop(makeSOCKmaster(master, port, timeout))
        mcexit(0L)
    }


    con <- socketAccept(socket = socket, blocking = TRUE, open = "a+b",
                        timeout = setup_timeout + 5)
    socketTimeout(socket = con, timeout = timeout)
    structure(list(con = con, host = "localhost", rank = rank),
              class = c("forknode", "SOCK0node"))
}
#  File src/library/parallel/R/unix/mcfork.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2018 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

### Derived from multicore version 0.1-6 by Simon Urbanek

### --- multicore --- low-level functions ---

## all not exported in parallel.

## registered as finalizer in .onLoad() to kill all child processes
clean_pids <- function(e)
    cleanup(kill = tools::SIGKILL, detach = TRUE, shutdown = TRUE)

## used in mclapply, mcparallel, newWorkNode
mcfork <- function(estranged = FALSE) {
    r <- .Call(C_mc_fork, estranged)

    # Disable JIT in the child process because it could lead to repeated
    # compilation of the same functions in each forked R process. Ideally
    # the compiled code would propagate to other processes, but it is not
    # currently possible.
    processClass <- if (!r[1L]) { compiler::enableJIT(0) ; "masterProcess" }  else
    		    if (is.na(r[2L])) "estrangedProcess" else "childProcess"
    structure(list(pid = r[1L], fd = r[2:3]), class = c(processClass, "process"))
}

## not used
readChildren <- function(timeout = 0)
    .Call(C_mc_read_children, as.double(timeout))

## used by mccollect, mclapply
readChild <- function(child)
{
    if (inherits(child, "process")) child <- processID(child)
    if (!is.numeric(child)) stop("invalid 'child' argument")
    .Call(C_mc_read_child, as.integer(child))
}

## used by mccollect, mclapply
selectChildren <- function(children = NULL, timeout = 0)
{
    if (!length(children)) children <- integer()
    if (inherits(children, "process")) children <- processID(children)
    if (is.list(children))
        children <- unlist(lapply(children, function(x)
	    if (inherits(x, "process")) processID(x)
            else stop("'children' must be a list of processes or a single process")
        ))
    if (!is.numeric(children))
        stop("'children' must be a list of processes or a single process")
    .Call(C_mc_select_children, as.double(timeout), as.integer(children))
}

## used by mccollect
rmChild <- function(child)
{
    if (inherits(child, "process")) child <- processID(child)
    if (!is.numeric(child)) stop("invalid 'child' argument")
    .Call(C_mc_rm_child, as.integer(child))
}

## not used
mckill <- function(process, signal = 2L)
{
    process <- processID(process)
    ## or simply tools::pskill(process, signal)
    unlist(lapply(process, function(p)
                  .Call(C_mc_kill, as.integer(p), as.integer(signal))))
}

## used by mcparallel, mclapply
sendMaster <- function(what, raw.asis=TRUE)
{
    # This is talking to the same machine, so no point in using xdr.
    if (!raw.asis || !is.raw(what)) what <- serialize(what, NULL, xdr = FALSE)
    .Call(C_mc_send_master, what)
}

## used widely, not exported
processID <- function(process) {
    if (inherits(process, "process")) process$pid
    else if (is.list(process)) unlist(lapply(process, processID))
    else stop(gettextf("'process' must be of class %s", dQuote("process")),
              domain = NA)
}

# not used
sendChildStdin <- function(child, what)
{
    if (inherits(child, "process") || is.list(child)) child <- processID(child)
    if (!is.numeric(child) || !length(child))
        stop("'child' must be a valid child process")
    child <- as.integer(child)
    if (is.character(what)) what <- charToRaw(paste(what, collapse='\n'))
    if (!is.raw(what)) stop("'what' must be a character or raw vector")
    invisible(unlist(lapply(child, function(p)
                            .Call(C_mc_send_child_stdin, p, what))))
}

## used by mcparallel, mclapply, newForkNode
mcexit <- function(exit.code = 0L, send = NULL)
{
    if (!is.null(send)) try(sendMaster(send, FALSE), silent = TRUE)
    .Call(C_mc_exit, as.integer(exit.code))
}

## used by mccollect, mclapply
children <- function(select)
{
    p <- .Call(C_mc_children)
    if (!missing(select)) p <- p[p %in% processID(select)]
    ## FIXME: this is not the meaning of this class as returned by mcfork
    lapply(p, function(x)
           structure(list(pid = x), class = c("childProcess", "process")))
}

## not used
childrenDescriptors <- function(index = 0L)
    .Call(C_mc_fds, as.integer(index))

## not used
masterDescriptor <- function() .Call(C_mc_master_fd)

## used by mclapply
isChild <- function() .Call(C_mc_is_child)

## used by mccollect, mclapply
closeStdout <- function(to.null=FALSE) .Call(C_mc_close_stdout, to.null)

## not used
closeStderr <- function(to.null=FALSE) .Call(C_mc_close_stderr, to.null)

## not used
closeFD <- function(fds) .Call(C_mc_close_fds, as.integer(fds))

## not used
closeAll <- function(includeStd = FALSE)
{
    if (!isChild()) {
        warning("closeAll() is a no-op in the master process", domain = NA)
        return(invisible(FALSE))
    }
    fds <- masterDescriptor()
    if (identical(fds, -1L)) fds <- integer(0)
    if (includeStd) fds <- c(1L, 2L, fds)
    mf <- max(fds) + 16L # take a few more ...
    ## close all but those that we actually use
    closeFD((1:mf)[-fds])
}

# used by mcparallel, mclapply, mcmapply
mcaffinity <- function(affinity = NULL) .Call(C_mc_affinity, affinity)

# used by mcparallel
mcinteractive <- function(interactive) .Call(C_mc_interactive, interactive)

# used by mclapply, pvec
prepareCleanup <- function() .Call(C_mc_prepare_cleanup)

# used by mclapply, pvec, mccollect
cleanup <- function(kill = TRUE, detach = TRUE, shutdown = FALSE)
    .Call(C_mc_cleanup, kill, detach, shutdown)
#  File src/library/parallel/R/unix/mclapply.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2019 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

### Derived from multicore version 0.1-6 by Simon Urbanek

mclapply <- function(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE,
                     mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
                     mc.cleanup = TRUE, mc.allow.recursive = TRUE,
                     affinity.list = NULL)
{
    cores <- as.integer(mc.cores)
    if((is.na(cores) || cores < 1L) && is.null(affinity.list))
        stop("'mc.cores' must be >= 1")
    .check_ncores(cores)

    if (isChild() && !isTRUE(mc.allow.recursive))
        return(lapply(X = X, FUN = FUN, ...))

    ## Follow lapply
    if(!is.vector(X) || is.object(X)) X <- as.list(X)
    if(!is.null(affinity.list) && length(affinity.list) < length(X))
        stop("affinity.list and X must have the same length")

    if(mc.set.seed) mc.reset.stream()
    if(length(X) < 2) {
        old.aff <- mcaffinity()
        mcaffinity(affinity.list[[1]])
        res <- lapply(X = X, FUN = FUN, ...)
        mcaffinity(old.aff)
        return(res)
    }

    if (length(X) < cores) cores <- length(X)
    if (cores < 2L && is.null(affinity.list))
	return(lapply(X = X, FUN = FUN, ...))

    jobs <- list()
    ## all processes created from now on will be terminated by cleanup
    prepareCleanup()
    on.exit(cleanup(mc.cleanup))
    if (!mc.preschedule) {              # sequential (non-scheduled)
        FUN <- match.fun(FUN)
        if (length(X) <= cores && is.null(affinity.list)) { # we can use one-shot parallel
            jobs <- lapply(seq_along(X),
                           function(i) mcparallel(FUN(X[[i]], ...),
                                                  name = names(X)[i],
                                                  mc.set.seed = mc.set.seed,
                                                  silent = mc.silent))
            res <- mccollect(jobs)
            if (length(res) == length(X)) names(res) <- names(X)
            has.errors <- sum(sapply(res, inherits, "try-error"))
        } else { # more complicated, we have to wait for jobs selectively
            sx <- seq_along(X)
            res <- vector("list", length(sx))
            names(res) <- names(X)
            fin <- rep(FALSE, length(X)) # values finished
            if (!is.null(affinity.list)) {
                ## build matrix for job mapping with affinity.list
                ## entry i,j is true if item i is allowed to run on core j
                cores <- max(unlist(x = affinity.list, recursive = TRUE))
                d0 <- logical(cores)
                cpu.map <- lapply(sx, function (i){
                    data <- d0
                    data[as.vector(affinity.list[[i]])] <- TRUE
                    data
                })
                ava <- do.call(rbind, cpu.map)
            } else {
                ## build matrix for job mapping without affinity.list
                ## all entries true
                ava <- matrix(TRUE, nrow = length(X), ncol = cores)
            }
            jobid <- integer(cores)
            ## choose first job for each core to start
            for (i in 1:cores) {
                jobid[i] <- match(TRUE, ava[,i]) # = which(ava[, i])[1]
                ava[jobid[i],] <- FALSE
            }
            ## remove unused cores from matrix
            if(anyNA(jobid)) {
                unused <- which(is.na(jobid))
                jobid <- jobid[-unused]
                ava   <- ava[, -unused, drop = FALSE]
            }
	    ## NOTE: we have to wrap the result in list() since readChild()
	    ## doesn't serialize raw vectors so if the result is a raw
	    ## it won't be serialized and we can't tell (see #17779)
	    ## This also allows us to distinguish try() in user code vs
	    ## our own (eventually)
            jobs <- lapply(jobid,
                           function(i) mcparallel(list(FUN(X[[i]], ...)),
                                                  mc.set.seed = mc.set.seed,
                                                  silent = mc.silent,
                                                  mc.affinity = affinity.list[[i]]))
            jobsp <- processID(jobs)
            has.errors <- 0L
            delivered.result <- 0L
            while (!all(fin)) {
                s <- selectChildren(jobs[!is.na(jobsp)], -1)
                if (is.null(s)) break   # no children -> no hope (should not happen)
                if (is.integer(s))
                    for (ch in s) {
                        ji <- match(TRUE, jobsp == ch)
                        ci <- jobid[ji]
                        r <- readChild(ch)
                        if (is.raw(r)) {
                            child.res <- unserialize(r)
                            if (inherits(child.res, "try-error"))
                                has.errors <- has.errors + 1L
			    ## unwrap the result
			    if (is.list(child.res)) child.res <- child.res[[1]]
                            ## we can't just assign it since a NULL
                            ## assignment would remove it from the list
                            if (!is.null(child.res)) res[[ci]] <- child.res
                            delivered.result <- delivered.result + 1L
                        } else {
                            fin[ci] <- TRUE
                            ## the job has finished, so we must not run
                            ## select on its fds again
                            jobsp[ji] <- jobid[ji] <- NA
                            if (any(ava)) { # still something to do,
                                ## look for first job which is allowed to
                                ## run on the now idling core and spawn it
                                nexti <- which.max(ava[, ji])
                                if(!is.na(nexti)) {
                                    jobid[ji] <- nexti
                                    jobs[[ji]] <- mcparallel(list(FUN(X[[nexti]], ...)),
                                                             mc.set.seed = mc.set.seed,
                                                             silent = mc.silent,
                                                             mc.affinity = affinity.list[[nexti]])
                                    jobsp[ji] <- processID(jobs[[ji]])
                                    ava[nexti,] <- FALSE
                                }
                            }
                        }
                    }
            }
            nores <- length(X) - delivered.result
            if (nores > 0)
                warning(sprintf(ngettext(nores,
                                         "%d parallel function call did not deliver a result",
                                         "%d parallel function calls did not deliver results"),
                                nores),
                        domain = NA)
        }
        if (has.errors)
            warning(gettextf("%d function calls resulted in an error",
                             has.errors), domain = NA)
        return(res)
    }

    ## mc.preschedule = TRUE from here on.
    if(!is.null(affinity.list))
        warning("'mc.preschedule' must be false if 'affinity.list' is used")
    sindex <- lapply(seq_len(cores),
                     function(i) seq(i, length(X), by = cores))
    schedule <- lapply(seq_len(cores),
                       function(i) X[seq(i, length(X), by = cores)])
    ch <- list()
    res <- vector("list", length(X))
    names(res) <- names(X)
    cp <- rep(0L, cores)
    fin <- rep(FALSE, cores)
    dr <- rep(FALSE, cores)
    inner.do <- function(core) {
        S <- schedule[[core]]
        f <- mcfork()
        if (isTRUE(mc.set.seed)) mc.advance.stream()
        if (inherits(f, "masterProcess")) { # this is the child process
            on.exit(mcexit(1L, structure("fatal error in wrapper code", class="try-error")))
            if (isTRUE(mc.set.seed)) mc.set.stream()
            if (isTRUE(mc.silent)) closeStdout(TRUE)
            sendMaster(try(lapply(X = S, FUN = FUN, ...), silent = TRUE))
            mcexit(0L)
        }
        jobs[[core]] <<- ch[[core]] <<- f
        cp[core] <<- processID(f)
        NULL
    }
    job.res <- lapply(seq_len(cores), inner.do)
    ac <- cp[cp > 0]
    has.errors <- integer(0)
    while (!all(fin)) {
        s <- selectChildren(ac[!fin], -1)
        if (is.null(s)) break # no children -> no hope we get anything (should not happen)
        if (is.integer(s))
            for (ch in s) {
                a <- readChild(ch)
                if (is.integer(a)) {
                    core <- which(cp == a)
                    fin[core] <- TRUE
                } else if (is.raw(a)) {
                    core <- which(cp == attr(a, "pid"))
                    job.res[[core]] <- ijr <- unserialize(a)
                    if (inherits(ijr, "try-error"))
                        has.errors <- c(has.errors, core)
                    dr[core] <- TRUE
                } else if (is.null(a)) {
                    # the child no longer exists (should not happen)
                    core <- which(cp == ch)
                    fin[core] <- TRUE
                }
            }
    }
    for (i in seq_len(cores)) {
        this <- job.res[[i]]
        if (inherits(this, "try-error")) { ## length-1 result
            for (j in sindex[[i]]) res[[j]] <- this
        } else
            ## we can't just assign it since a NULL
            ## assignment would remove it from the list
            if (!is.null(this)) res[sindex[[i]]] <- this
    }
    nores <- cores - sum(dr)
    if (nores > 0)
        warning(sprintf(ngettext(nores,
                                 "scheduled core %s did not deliver a result, all values of the job will be affected",
                                 "scheduled cores %s did not deliver results, all values of the jobs will be affected"),
                        paste(which(dr == FALSE), collapse = ", ")),
                domain = NA)
    if (length(has.errors)) {
        if (length(has.errors) == cores)
            warning("all scheduled cores encountered errors in user code")
        else
            warning(sprintf(ngettext(length(has.errors),
                                     "scheduled core %s encountered error in user code, all values of the job will be affected",
                                     "scheduled cores %s encountered errors in user code, all values of the jobs will be affected"),
                            paste(has.errors, collapse = ", ")),
                    domain = NA)
    }
    res
}
#  File src/library/parallel/R/unix/mcmapply.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2020 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

mcmapply <-
    function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE,
             mc.preschedule = TRUE, mc.set.seed = TRUE,
             mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
             mc.cleanup = TRUE, affinity.list = NULL)
{
    FUN <- match.fun(FUN)
    dots <- list(...)
    if(!length(dots)) return(list())
    lens <- lengths(dots)
    n <- max(lens)
    if(n && min(lens) == 0L)
        stop("Zero-length inputs cannot be mixed with those of non-zero length")
    answer <- if(n < 2L){ 
      ## ensure that it runs on the right core 
      if(!is.null(affinity.list)){
        save <- mcaffinity()
        mcaffinity(affinity.list[[1]])
      }
      answer <- .mapply(FUN, dots, MoreArgs)
        # .mapply will not catch errors thrown by FUN
      if(!is.null(affinity.list)) mcaffinity(save)
      answer
    } else {    
        ## recycle shorter vectors
        X <- if (!all(lens == n))
            lapply(dots, function(x) rep(x, length.out = n))
        else dots
        do_one <- function(indices, ...) {
            dots <- lapply(X, function(x) x[indices])
            .mapply(FUN, dots, MoreArgs)
        }
        answer <- mclapply(seq_len(n), do_one, mc.preschedule = mc.preschedule,
                           mc.set.seed = mc.set.seed, mc.silent = mc.silent,
                           mc.cores = mc.cores, mc.cleanup = mc.cleanup, 
                           affinity.list = affinity.list)
        answer <- lapply(answer, function(x) {
            if (inherits(x, "try-error")) {
                SIMPLIFY <<- FALSE # protect attributes from simplify2array()
                list(x) # protect attributes from c()
            } else
                x
        })
        do.call(c, answer)
    }
    if (USE.NAMES && length(dots)) {
        if (is.null(names1 <- names(dots[[1L]])) && is.character(dots[[1L]]))
            names(answer) <- dots[[1L]]
        else if (!is.null(names1))
            names(answer) <- names1
    }
    if (!identical(SIMPLIFY, FALSE))
        simplify2array(answer, higher = (SIMPLIFY == "array"))
    else answer
}

mcMap <- function (f, ...)
{
    f <- match.fun(f)
    mcmapply(f, ..., SIMPLIFY = FALSE, mc.silent = TRUE)
}
#  File src/library/parallel/R/unix/mcparallel.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2019 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

### Derived from multicore version 0.1-6 by Simon Urbanek

mcparallel <- function(expr, name, mc.set.seed = TRUE, silent = FALSE, mc.affinity = NULL, mc.interactive = FALSE, detached = FALSE)
{
    f <- mcfork(detached)
    env <- parent.frame()
    if (isTRUE(mc.set.seed)) mc.advance.stream()
    if (inherits(f, "masterProcess")) {
        on.exit(mcexit(1L, structure("fatal error in wrapper code",
                                  class = "try-error")))
        if (isTRUE(mc.set.seed)) mc.set.stream()
        mc.interactive <- as.logical(mc.interactive)
        if (isTRUE(mc.interactive)) mcinteractive(TRUE)
        if (isTRUE(!mc.interactive)) mcinteractive(FALSE)
        if (!is.null(mc.affinity)) mcaffinity(mc.affinity)
        if (isTRUE(silent)) closeStdout(TRUE)
	if (detached) {
	    on.exit(mcexit(1L))
	    eval(expr, env)
	    mcexit(0L)
	}
	sendMaster(try(eval(expr, env), silent = TRUE), FALSE)
        mcexit(0L)
    }
    if (!missing(name) && !is.null(name)) f$name <- as.character(name)[1L]
    class(f) <- c("parallelJob", class(f))
    f
}

mccollect <- function(jobs, wait = TRUE, timeout = 0, intermediate = FALSE)
{
    if (missing(jobs)) jobs <- children()
    if (!length(jobs)) return (NULL)
    if (isTRUE(intermediate)) intermediate <- utils::str
    pids <- if (inherits(jobs, "process") || is.list(jobs))
        processID(jobs) else jobs
    if (!length(pids)) return(NULL)
    if (!is.numeric(pids)) stop("invalid 'jobs' argument")
    pids <- as.integer(pids)
    pnames <- as.character(pids)
    if (!inherits(jobs, "process") && is.list(jobs))
        for(i in seq(jobs))
            if (!is.null(jobs[[i]]$name))
                pnames[i] <- as.character(jobs[[i]]$name)

    if (!wait) {
        s <- selectChildren(jobs, timeout)
        if (is.logical(s) || !length(s)) return(NULL) ## select error
        res <- lapply(s, function(x) NULL)
        delivered.result <- 0
        for (i in seq_along(s)) {
            x <- s[i]
            r <- readChild(x)
            if (is.raw(r)) {
                rmChild(x) ## avoid zombie process without waiting
                ## unserialize(r) might be null
                res[i] <- list(unserialize(r))
                delivered.result <- delivered.result + 1L
            }
        }
        names(res) <- pnames[match(s, pids)]
        expected.result <- length(s)
    } else {
        res <- lapply(pids, function(x) NULL)
        names(res) <- pnames
        fin <- rep(FALSE, length(pids))
        delivered.result <- 0
        while (!all(fin)) {
            s <- selectChildren(pids[!fin], -1)
            if (is.integer(s)) {
                for (pid in s) {
                    r <- readChild(pid)
                    if (is.raw(r)) {
                        ## unserialize(r) might be null
                        res[which(pid == pids)] <- list(unserialize(r))
                        delivered.result <- delivered.result + 1L
                    } else
                        ## child exiting or error
                        fin[pid == pids] <- TRUE 
                }
                if (is.function(intermediate)) intermediate(res)
            } else
                ## should not happen (select error)
                if (all(is.na(match(pids, processID(children()))))) break
        }
	expected.result <- length(pids)
        
    }
    nores <- expected.result - delivered.result
    if (nores > 0)
        warning(sprintf(ngettext(nores,
                                 "%d parallel job did not deliver a result",
                                 "%d parallel jobs did not deliver results"),
                        nores),
                domain = NA)
    cleanup(kill = FALSE, detach = FALSE) # compact children
    res
}
#  File src/library/parallel/R/unix/pvec.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2017 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

### Derived from multicore version 0.1-6 by Simon Urbanek

pvec <- function(v, FUN, ..., mc.set.seed = TRUE, mc.silent = FALSE,
                 mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE)
{
    if (!is.vector(v)) stop("'v' must be a vector")

    cores <- as.integer(mc.cores)
    if(cores < 1L) stop("'mc.cores' must be >= 1")
    if(cores == 1L) return(FUN(v, ...))
    .check_ncores(cores)

    if(mc.set.seed) mc.reset.stream()

    n <- length(v)
    l <- if (n <= cores) as.list(v) else {
        ## compute the scheduling, making it as fair as possible
        il <- as.integer(n / cores)
        xc <- n - il * cores
        sl <- rep(il, cores)
        if (xc) sl[1:xc] <- il + 1L
        si <- cumsum(c(1L, sl))
        se <- si + c(sl, 0L) - 1L
        lapply(seq_len(cores), function(ix) v[si[ix]:se[ix]])
    }
    jobs <- NULL
    ## all processes created from now on will be terminated by cleanup
    prepareCleanup()
    on.exit(cleanup())
    FUN <- match.fun(FUN)
    ## may have more cores than tasks ....
    jobs <- lapply(seq_len(min(n, cores)),
                   function(i) mcparallel(FUN(l[[i]], ...), name = i,
                                          mc.set.seed = mc.set.seed,
                                          silent = mc.silent))
    res <- mccollect(jobs)
    names(res) <- NULL
    res <- do.call(c, res)
    if (length(res) != n)
        warning("some results may be missing, folded or caused an error")
    res
}
#  File src/library/parallel/R/worker.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2020 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

## Derived from snow 0.3-6 by Luke Tierney

## NB: there is also workerCommand in snowSOCK.R
workCommand <- function(master)
{
    tryCatch({
        msg <- recvData(master)
        # cat(paste("Type:", msg$type, "\n"))

        if (msg$type == "DONE") {
            closeNode(master)
            FALSE
        } else if (msg$type == "EXEC") {
            success <- TRUE
            ## This uses the message rather than the exception since
            ## the exception class/methods may not be available on the
            ## master.
            handler <- function(e) {
                success <<- FALSE
                structure(conditionMessage(e),
                          class = c("snow-try-error","try-error"))
            }
            t1 <- proc.time()
            value <- tryCatch(do.call(msg$data$fun, msg$data$args, quote = TRUE),
                              error = handler)
            t2 <- proc.time()
            value <- list(type = "VALUE", value = value, success = success,
                          time = t2 - t1, tag = msg$data$tag)
            msg <- NULL ## release for GC
            sendData(master, value)
            value <- NULL ## release for GC
            TRUE
        } else {
            ## unknown command / shutdown
            TRUE
        }
    }, interrupt = function(e) TRUE)
}

workLoop <- function(master)
{
    if (!is.null(master))
        while(workCommand(master)) {}
}

## NB: this only sinks the connections, not C-level stdout/err.
sinkWorkerOutput <- function(outfile)
{
    if (nzchar(outfile)) {
        if (.Platform$OS.type == "windows" && outfile == "/dev/null")
            outfile <- "nul:"
        ## all the workers log to the same file.
        outcon <- file(outfile, open = "a")
        sink(outcon)
        sink(outcon, type = "message")
    }
}

#  File src/library/parallel/R/zzz.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2014 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

.noGenerics <- TRUE

if (.Platform$OS.type == "windows")
    utils::globalVariables(c("mc_pids", "clean_pids"), add = TRUE)

## dummy, just so we can register a finalizer
.fin.env <- new.env()

.onLoad <- function(libname, pkgname)
{
    initDefaultClusterOptions(libname)
    cores <- getOption("mc.cores", NULL)
    if(is.null(cores) && !is.na(nc <- as.integer(Sys.getenv("MC_CORES"))))
        options("mc.cores" = nc)
    if(.Platform$OS.type == "unix") reg.finalizer(.fin.env, clean_pids, TRUE)
}

.onUnload <-
function(libpath)
    library.dynam.unload("parallel", libpath)
