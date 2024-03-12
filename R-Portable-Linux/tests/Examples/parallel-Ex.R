pkgname <- "parallel"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('parallel')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("RngStream")
### * RngStream

flush(stderr()); flush(stdout())

### Name: RNGstreams
### Title: Implementation of Pierre L'Ecuyer's RngStreams
### Aliases: nextRNGStream nextRNGSubStream clusterSetRNGStream
###   mc.reset.stream
### Keywords: distribution sysdata

### ** Examples

RNGkind("L'Ecuyer-CMRG")
set.seed(123)
(s <- .Random.seed)
## do some work involving random numbers.
nextRNGStream(s)
nextRNGSubStream(s)



cleanEx()
nameEx("children")
### * children

flush(stderr()); flush(stdout())

### Name: mcchildren
### Title: Low-level Functions for Management of Forked Processes
### Aliases: children readChild readChildren selectChildren sendChildStdin
###   sendMaster mckill
### Keywords: interface

### ** Examples
## Not run: 
##D p  <- mcparallel(scan(n = 1, quiet = TRUE))
##D sendChildStdin(p, "17.4\n")
##D mccollect(p)[[1]]
## End(Not run)


cleanEx()
nameEx("clusterApply")
### * clusterApply

flush(stderr()); flush(stdout())

### Name: clusterApply
### Title: Apply Operations using Clusters
### Aliases: clusterApply clusterApplyLB clusterCall clusterEvalQ
###   clusterExport clusterMap clusterSplit parApply parCapply parLapply
###   parRapply parSapply parLapplyLB parSapplyLB

### ** Examples


cleanEx()
nameEx("detectCores")
### * detectCores

flush(stderr()); flush(stdout())

### Name: detectCores
### Title: Detect the Number of CPU Cores
### Aliases: detectCores

### ** Examples

detectCores()
detectCores(logical = FALSE)



cleanEx()
nameEx("mcfork")
### * mcfork

flush(stderr()); flush(stdout())

### Name: mcfork
### Title: Fork a Copy of the Current R Process
### Aliases: mcfork mcexit
### Keywords: interface

### ** Examples


cleanEx()
nameEx("mclapply")
### * mclapply

flush(stderr()); flush(stdout())

### Name: mclapply
### Title: Parallel Versions of 'lapply' and 'mapply' using Forking
### Aliases: mclapply mcmapply mcMap
### Keywords: interface

### ** Examples


cleanEx()
nameEx("mcparallel")
### * mcparallel

flush(stderr()); flush(stdout())

### Name: mcparallel
### Title: Evaluate an R Expression Asynchronously in a Separate Process
### Aliases: mccollect mcparallel
### Keywords: interface

### ** Examples

p <- mcparallel(1:10)
q <- mcparallel(1:20)
# wait for both jobs to finish and collect all results
res <- mccollect(list(p, q))

## IGNORE_RDIFF_BEGIN
## reports process ids, so not reproducible
p <- mcparallel(1:10)
mccollect(p, wait = FALSE, 10) # will retrieve the result (since it's fast)
mccollect(p, wait = FALSE)     # will signal the job as terminating
mccollect(p, wait = FALSE)     # there is no longer such a job
## IGNORE_RDIFF_END

## Don't show: 
set.seed(123, "L'Ecuyer"); mc.reset.stream()
## End(Don't show)
# a naive parallel lapply can be created using mcparallel alone:
jobs <- lapply(1:10, function(x) mcparallel(rnorm(x), name = x))
mccollect(jobs)



cleanEx()
nameEx("pvec")
### * pvec

flush(stderr()); flush(stdout())

### Name: pvec
### Title: Parallelize a Vector Map Function using Forking
### Aliases: pvec
### Keywords: interface

### ** Examples

x <- pvec(1:1000, sqrt)
stopifnot(all(x == sqrt(1:1000)))



cleanEx()
nameEx("splitIndices")
### * splitIndices

flush(stderr()); flush(stdout())

### Name: splitIndices
### Title: Divide Tasks for Distribution in a Cluster
### Aliases: splitIndices
### Keywords: utility

### ** Examples

splitIndices(20, 3)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
