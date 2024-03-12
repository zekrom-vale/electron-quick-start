pkgname <- "compiler"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('compiler')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("compile")
### * compile

flush(stderr()); flush(stdout())

### Name: compile
### Title: Byte Code Compiler
### Aliases: compile cmpfun cmpfile loadcmp disassemble enableJIT
###   compilePKGS getCompilerOption setCompilerOptions
### Keywords: programming

### ** Examples

oldJIT <- enableJIT(0)
# a simple example
f <- function(x) x+1
fc <- cmpfun(f)
fc(2)
disassemble(fc)

# old R version of lapply
la1 <- function(X, FUN, ...) {
    FUN <- match.fun(FUN)
    if (!is.list(X))
	X <- as.list(X)
    rval <- vector("list", length(X))
    for(i in seq_along(X))
	rval[i] <- list(FUN(X[[i]], ...))
    names(rval) <- names(X)		  # keep `names' !
    return(rval)
}
# a small variation
la2 <- function(X, FUN, ...) {
    FUN <- match.fun(FUN)
    if (!is.list(X))
	X <- as.list(X)
    rval <- vector("list", length(X))
    for(i in seq_along(X)) {
        v <- FUN(X[[i]], ...)
        if (is.null(v)) rval[i] <- list(v)
        else rval[[i]] <- v
    }
    names(rval) <- names(X)		  # keep `names' !
    return(rval)
}
# Compiled versions
la1c <- cmpfun(la1)
la2c <- cmpfun(la2)
# some timings
x <- 1:10
y <- 1:100
## Don't show: 
for (i in 1:10000) lapply(x, is.null)
for (i in 1:10000) la1(x, is.null)
for (i in 1:10000) la1c(x, is.null)
for (i in 1:10000) la2(x, is.null)
for (i in 1:10000) la2c(x, is.null)
for (i in 1:1000) lapply(y, is.null)
for (i in 1:1000) la1(y, is.null)
for (i in 1:1000) la1c(y, is.null)
for (i in 1:1000) la2(y, is.null)
for (i in 1:1000) la2c(y, is.null)
## End(Don't show)
enableJIT(oldJIT)


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
