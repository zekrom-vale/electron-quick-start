pkgname <- "splines"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('splines')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("asVector")
### * asVector

flush(stderr()); flush(stdout())

### Name: asVector
### Title: Coerce an Object to a Vector
### Aliases: asVector
### Keywords: models

### ** Examples

require(stats)
ispl <- interpSpline( weight ~ height,  women )
pred <- predict(ispl)
class(pred)
utils::str(pred)
asVector(pred)



cleanEx()
nameEx("backSpline")
### * backSpline

flush(stderr()); flush(stdout())

### Name: backSpline
### Title: Monotone Inverse Spline
### Aliases: backSpline
### Keywords: models

### ** Examples

require(graphics)
ispl <- interpSpline( women$height, women$weight )
bspl <- backSpline( ispl )
plot( bspl )                   # plots over the range of the knots
points( women$weight, women$height )



cleanEx()
nameEx("bs")
### * bs

flush(stderr()); flush(stdout())

### Name: bs
### Title: B-Spline Basis for Polynomial Splines
### Aliases: bs
### Keywords: smooth

### ** Examples

require(stats); require(graphics)
bs(women$height, df = 5)
summary(fm1 <- lm(weight ~ bs(height, df = 5), data = women))

## example of safe prediction
plot(women, xlab = "Height (in)", ylab = "Weight (lb)")
ht <- seq(57, 73, length.out = 200)
lines(ht, predict(fm1, data.frame(height = ht)))
## Don't show: 
## Consistency:
x <- c(1:3, 5:6)
stopifnot(identical(bs(x), bs(x, df = 3)),
          identical(bs(x, df = 4), bs(x, df = 4, knots = NULL)), # not true till 2.15.2
          !is.null(kk <- attr(bs(x), "knots")), # not true till 1.5.1
          length(kk) == 0)
## End(Don't show)


cleanEx()
nameEx("interpSpline")
### * interpSpline

flush(stderr()); flush(stdout())

### Name: interpSpline
### Title: Create an Interpolation Spline
### Aliases: interpSpline
### Keywords: models

### ** Examples

require(graphics); require(stats)
ispl <- interpSpline( women$height, women$weight )
ispl2 <- interpSpline( weight ~ height,  women )
# ispl and ispl2 should be the same
plot( predict( ispl, seq( 55, 75, length.out = 51 ) ), type = "l" )
points( women$height, women$weight )
plot( ispl )    # plots over the range of the knots
points( women$height, women$weight )
splineKnots( ispl )



cleanEx()
nameEx("ns")
### * ns

flush(stderr()); flush(stdout())

### Name: ns
### Title: Generate a Basis Matrix for Natural Cubic Splines
### Aliases: ns
### Keywords: smooth

### ** Examples

require(stats); require(graphics)
ns(women$height, df = 5)
summary(fm1 <- lm(weight ~ ns(height, df = 5), data = women))

## To see what knots were selected
attr(terms(fm1), "predvars")

## example of safe prediction
plot(women, xlab = "Height (in)", ylab = "Weight (lb)")
ht <- seq(57, 73, length.out = 200) ; nD <- data.frame(height = ht)
lines(ht, p1 <- predict(fm1, nD))
stopifnot(all.equal(p1, predict(update(fm1, . ~
                            splines::ns(height, df=5)), nD)))
          # not true in R < 3.5.0
## Don't show: 
## Consistency:
x <- c(1:3, 5:6)
stopifnot(identical(ns(x), ns(x, df = 1)),
          identical(ns(x, df = 2),
                    ns(x, df = 2, knots = NULL)), # not true till 2.15.2
          !is.null(kk <- attr(ns(x), "knots")), # not true till 1.5.1
          length(kk) == 0)
## End(Don't show)



cleanEx()
nameEx("periodicSpline")
### * periodicSpline

flush(stderr()); flush(stdout())

### Name: periodicSpline
### Title: Create a Periodic Interpolation Spline
### Aliases: periodicSpline
### Keywords: models

### ** Examples

require(graphics); require(stats)
xx <- seq( -pi, pi, length.out = 16 )[-1]
yy <- sin( xx )
frm <- data.frame( xx, yy )
pispl <- periodicSpline( xx, yy, period = 2 * pi)
pispl2 <- periodicSpline( yy ~ xx, frm, period = 2 * pi )
stopifnot(all.equal(pispl, pispl2))  # pispl and pispl2 are the same

plot( pispl )          # displays over one period
points( yy ~ xx, col = "brown")
plot( predict( pispl, seq(-3*pi, 3*pi, length.out = 101) ), type = "l" )



cleanEx()
nameEx("polySpline")
### * polySpline

flush(stderr()); flush(stdout())

### Name: polySpline
### Title: Piecewise Polynomial Spline Representation
### Aliases: polySpline as.polySpline
### Keywords: models

### ** Examples

require(graphics)
ispl <- polySpline(interpSpline( weight ~ height,  women, bSpline = TRUE))
## IGNORE_RDIFF_BEGIN
print( ispl )   # print the piecewise polynomial representation
## IGNORE_RDIFF_END
plot( ispl )    # plots over the range of the knots
points( women$height, women$weight )



cleanEx()
nameEx("predict.bSpline")
### * predict.bSpline

flush(stderr()); flush(stdout())

### Name: predict.bSpline
### Title: Evaluate a Spline at New Values of x
### Aliases: predict.bSpline predict.nbSpline predict.pbSpline
###   predict.npolySpline predict.ppolySpline
### Keywords: models

### ** Examples

require(graphics); require(stats)
ispl <- interpSpline( weight ~ height,  women )
opar <- par(mfrow = c(2, 2), las = 1)
plot(predict(ispl, nseg = 201),     # plots over the range of the knots
     main = "Original data with interpolating spline", type = "l",
     xlab = "height", ylab = "weight")
points(women$height, women$weight, col = 4)
plot(predict(ispl, nseg = 201, deriv = 1),
     main = "First derivative of interpolating spline", type = "l",
     xlab = "height", ylab = "weight")
plot(predict(ispl, nseg = 201, deriv = 2),
     main = "Second derivative of interpolating spline", type = "l",
     xlab = "height", ylab = "weight")
plot(predict(ispl, nseg = 401, deriv = 3),
     main = "Third derivative of interpolating spline", type = "l",
     xlab = "height", ylab = "weight")
par(opar)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("predict.bs")
### * predict.bs

flush(stderr()); flush(stdout())

### Name: predict.bs
### Title: Evaluate a Spline Basis
### Aliases: predict.bs predict.ns
### Keywords: smooth

### ** Examples

require(stats)
basis <- ns(women$height, df = 5)
newX <- seq(58, 72, length.out = 51)
# evaluate the basis at the new data
predict(basis, newX)



cleanEx()
nameEx("splineDesign")
### * splineDesign

flush(stderr()); flush(stdout())

### Name: splineDesign
### Title: Design Matrix for B-splines
### Aliases: splineDesign spline.des
### Keywords: models

### ** Examples

require(graphics)
splineDesign(knots = 1:10, x = 4:7)
splineDesign(knots = 1:10, x = 4:7, derivs = 1)
## visualize band structure

knots <- c(1,1.8,3:5,6.5,7,8.1,9.2,10)  # 10 => 10-4 = 6 Basis splines
x <- seq(min(knots)-1, max(knots)+1, length.out = 501)
bb <- splineDesign(knots, x = x, outer.ok = TRUE)

plot(range(x), c(0,1), type = "n", xlab = "x", ylab = "",
     main =  "B-splines - sum to 1 inside inner knots")
mtext(expression(B[j](x) *"  and "* sum(B[j](x), j == 1, 6)), adj = 0)
abline(v = knots, lty = 3, col = "light gray")
abline(v = knots[c(4,length(knots)-3)], lty = 3, col = "gray10")
lines(x, rowSums(bb), col = "gray", lwd = 2)
matlines(x, bb, ylim = c(0,1), lty = 1)



cleanEx()
nameEx("splineKnots")
### * splineKnots

flush(stderr()); flush(stdout())

### Name: splineKnots
### Title: Knot Vector from a Spline
### Aliases: splineKnots
### Keywords: models

### ** Examples

ispl <- interpSpline( weight ~ height, women )
splineKnots( ispl )



cleanEx()
nameEx("splineOrder")
### * splineOrder

flush(stderr()); flush(stdout())

### Name: splineOrder
### Title: Determine the Order of a Spline
### Aliases: splineOrder
### Keywords: models

### ** Examples

splineOrder( interpSpline( weight ~ height, women ) )



cleanEx()
nameEx("xyVector")
### * xyVector

flush(stderr()); flush(stdout())

### Name: xyVector
### Title: Construct an 'xyVector' Object
### Aliases: xyVector
### Keywords: models

### ** Examples

require(stats); require(graphics)
ispl <- interpSpline( weight ~ height, women )
weights <- predict( ispl, seq( 55, 75, length.out = 51 ))
class( weights )
plot( weights, type = "l", xlab = "height", ylab = "weight" )
points( women$height, women$weight )
weights



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
