#  File src/library/grDevices/R/Hershey.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
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

Hershey <-
    list(typeface =
         c("serif", "sans serif", "script",
           "gothic english", "gothic german", "gothic italian",
           "serif symbol", "sans serif symbol"),
         fontindex =
         c("plain", "italic", "bold", "bold italic",
           "cyrillic", "oblique cyrillic", "EUC"),
## List of valid combinations : ../man/Hershey.Rd
## *checking* of allowed combinations is done in
## (via max{#}) in    FixupVFont() ../../../main/plot.c
## The basic "table" really is in  ../../../modules/vfonts/g_fontdb.c

         allowed = rbind(cbind(1L, 1L:7L), cbind(2L, 1L:4L), cbind(3L,1L:3L),
                         cbind(4L:6L, 1L), cbind(7L, 1L:4L), cbind(8L,1L:2L))
         )
#  File src/library/grDevices/R/cairo.R
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

svg <- function(filename = if(onefile) "Rplots.svg" else "Rplot%03d.svg",
                width = 7, height = 7, pointsize = 12,
                onefile = FALSE, family = "sans", bg = "white",
                antialias = c("default", "none", "gray", "subpixel"),
                symbolfamily)
{
    if(!checkIntFormat(filename)) stop("invalid 'filename'")
    if(!capabilities("cairo"))
        stop("svg: Cairo-based devices are not available for this platform")
    antialiases <- eval(formals()$antialias)
    antialias <- match(match.arg(antialias, antialiases), antialiases)
    if (missing(symbolfamily)) symbolfamily <- symbolfamilyDefault(family)
    invisible(.External(C_devCairo, filename, 4L, 72*width, 72*height,
                        pointsize, bg, NA_integer_, antialias, onefile,
                        family, 300, checkSymbolFont(symbolfamily)))
}

cairo_pdf <- function(filename = if(onefile) "Rplots.pdf" else "Rplot%03d.pdf",
                      width = 7, height = 7, pointsize = 12,
                      onefile = FALSE, family = "sans", bg = "white",
                      antialias = c("default", "none", "gray", "subpixel"),
                      fallback_resolution = 300, symbolfamily)
{
    if(!checkIntFormat(filename)) stop("invalid 'filename'")
    if(!capabilities("cairo"))
        stop("cairo_pdf: Cairo-based devices are not available for this platform")
    antialiases <- eval(formals()$antialias)
    antialias <- match(match.arg(antialias, antialiases), antialiases)
    if (missing(symbolfamily)) symbolfamily <- symbolfamilyDefault(family)
    invisible(.External(C_devCairo, filename, 6L, 72*width, 72*height,
                        pointsize, bg, NA_integer_, antialias, onefile,
                        family, fallback_resolution,
                        checkSymbolFont(symbolfamily)))
}

cairo_ps <- function(filename = if(onefile) "Rplots.ps" else "Rplot%03d.ps",
                     width = 7, height = 7, pointsize = 12,
                     onefile = FALSE, family = "sans", bg = "white",
                     antialias = c("default", "none", "gray", "subpixel"),
                     fallback_resolution = 300, symbolfamily)
{
    if(!checkIntFormat(filename)) stop("invalid 'filename'")
    if(!capabilities("cairo"))
        stop("cairo_ps: Cairo-based devices are not available for this platform")
    antialiases <- eval(formals()$antialias)
    antialias <- match(match.arg(antialias, antialiases), antialiases)
    if (missing(symbolfamily)) symbolfamily <- symbolfamilyDefault(family)
    invisible(.External(C_devCairo, filename, 7L, 72*width, 72*height,
                        pointsize, bg, NA_integer_, antialias, onefile,
                        family, fallback_resolution,
                        checkSymbolFont(symbolfamily)))
}

cairoVersion <- function() .Call(C_cairoVersion)

pangoVersion <- function() .Call(C_pangoVersion)

cairoFT <- function() .Call(C_cairoFT)

## Modelled on utils::compareVersion()
comparePangoVersion <- function(a, b) {
    a <- as.integer(strsplit(a, "[.]")[[1L]])
    b <- as.integer(strsplit(b, "[.]")[[1L]])
    for (k in seq_along(a)) {
        if (k <= length(b)) {
            if (a[k] > b[k])
                return(1)
            else if (a[k] < b[k])
                return(-1L)
        } else {
            return(1L)
        }
    }
    if (length(b) > length(a))
        return(-1L)
    else
        return(0L)
}

symbolType1support <- function() {
    pangoVersion <- grSoftVersion()["pango"]
    pangoVersion == "" ||
        comparePangoVersion(pangoVersion, "1.44") < 0
}

cairoSymbolFont <- function(family, usePUA = TRUE) {
    font <- as.character(family)
    attr(font, "usePUA") <- as.logical(usePUA)
    class(font) <- "CairoSymbolFont"
    font
}

checkSymbolFont <- function(x) {
    if (inherits(x, "CairoSymbolFont")) {
        x
    } else {
        cairoSymbolFont(x)
    }
}

symbolfamilyDefault <- function(family) {
    if (symbolType1support()) {
        if (grSoftVersion()["cairoFT"] == "yes") {
            if (.Platform$OS.type == "windows") {
                "Standard Symbols L"
            } else {
                "Symbol"
            }
        } else {
            "Symbol"
        }
    } else {
        cairoSymbolFont(family, usePUA=FALSE)
    }
}

#  File src/library/grDevices/R/calc.R
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

#### Functions that calculate useful stuff for plotting
#### BUT which do not do any actual drawing
#### Useful for both graphics and grid to have access to

boxplot.stats <- function(x, coef = 1.5, do.conf = TRUE, do.out = TRUE)
{
    if(coef < 0) stop("'coef' must not be negative")
    nna <- !is.na(x)
    n <- sum(nna)                       # including +/- Inf
    stats <- stats::fivenum(x, na.rm = TRUE)
    iqr <- diff(stats[c(2, 4)])
    if(coef == 0)
	do.out <- FALSE
    else { ## coef > 0
	out <- if(!is.na(iqr)) { x < (stats[2L] - coef * iqr) |
				 x > (stats[4L] + coef * iqr)
			     } else !is.finite(x)
	if(any(out[nna], na.rm = TRUE))
	    stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
    }
    conf <- if(do.conf) stats[3L] + c(-1.58, 1.58) * iqr / sqrt(n)
    list(stats = stats, n = n, conf = conf,
	 out = if(do.out) x[out & nna] else numeric())
}

## Contour lines
contourLines <-
function (x = seq(0, 1, length.out = nrow(z)),
          y = seq(0, 1, length.out = ncol(z)),
	  z, nlevels = 10, levels = pretty(range(z, na.rm = TRUE), nlevels))
{
    ## FIXME: This "validation" code for the x, y, z values
    ## should be put in a function for contourLines, contour,
    ## image (and persp?) to share.  Unfortunately, an xyz.coords
    ## already exists which isn't really compatible with the
    ## desired behaviour here.
    if (missing(z)) {
	if (!missing(x)) {
	    if (is.list(x)) {
		z <- x$z; y <- x$y; x <- x$x
	    } else {
		z <- x
		x <- seq.int(0, 1, length.out = nrow(z))
	    }
	} else stop("no 'z' matrix specified")
    } else if (is.list(x)) {
	y <- x$y
	x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
	stop("increasing 'x' and 'y' values expected")
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1)
	stop("no proper 'z' matrix specified")
    if (1.0 * length(x) * length(y) != length(z))
        stop("dimensions of 'x', 'y' and 'z' do not match")
    .External2(C_contourLines, x, y, z, levels)
}

chull <- function(x, y = NULL)
{
    X <- xy.coords(x, y, recycle = TRUE, setLab = FALSE)
    x <- cbind(X$x, X$y)
    if(any(!is.finite(x))) stop("finite coordinates are needed")
    if(nrow(x) == 0) return(integer())
    if(nrow(x) == 1) return(1L)
    res <- .Call(C_chull, x)
    ## if this is called on multiple copies of a single point
    ## res is of length one.
    if (length(res) < 2L) return(res)
    ## fix up order: needed in rare cases: PR#15127
    xx <- sweep(x[res, ], 2L, colMeans(x[res, ]))
    angs <- atan2(xx[, 2L], -xx[, 1L])
    res[order(angs)]
}

nclass.Sturges <- function(x) ceiling(log2(length(x)) + 1)

nclass.scott <- function(x)
{
    h <- 3.5 * sqrt(stats::var(x)) * length(x)^(-1/3)
    if(h > 0) max(1, ceiling(diff(range(x))/h)) else 1L
}

nclass.FD <- function(x, digits = 5)
{
    h <- 2 * stats::IQR(x. <- signif(x, digits=digits))
    if (h == 0) {
	x. <- sort(x.)
	al <- 1/4; al.min <- 1/512 # try quantiles 1/8, 1/16, ... 1/512
	while(h == 0 && (al <- al/2) >= al.min)
	    h <- diff(stats::quantile(x., c(al, 1-al), names = FALSE)) / (1 - 2*al)
    }
    if (h == 0) ## revert to Scott's:
	h <- 3.5 * sqrt(stats::var(x))
    if (h > 0) ceiling(diff(range(x))/h * length(x)^(1/3)) else 1L
}


## Sunflower Plot computation:
## Used to be part of ../../graphics/R/sunflowerplot.R :
xyTable <- function(x, y = NULL, digits)
{
    ## Compute number := multiplicities of (x[i], y[i])

    x <- xy.coords(x, y, setLab = FALSE)

    ## get rid of rounding fuzz:
    y <- signif(x$y, digits=digits)
    x <- signif(x$x, digits=digits)
    n <- length(x)
    number <-
	if(n > 0) {
	    orderxy <- order(x, y)
	    x <- x[orderxy]
	    y <- y[orderxy]
	    first <- c(TRUE, (x[-1L] != x[-n]) | (y[-1L] != y[-n]))
	    x <- x[first]
	    y <- y[first]
	    diff(c((1L:n)[first], n + 1L))
	}
	else integer()

    list(x = x, y = y, number = number)
}

axisTicks <- function(usr, log, axp = NULL, nint = 5) {
    if(is.null(axp))
	axp <- unlist(.axisPars(usr, log=log, nintLog=nint), use.names=FALSE)
    .Call(C_R_CreateAtVector, axp, if(log) 10^usr else usr, nint, log)
}

.axisPars <- function(usr, log = FALSE, nintLog = 5) {
    .Call(C_R_GAxisPars, usr, log, nintLog)
}

#  File src/library/grDevices/R/gradients.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 2019      The R Foundation
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

## Create R objects defining clipping paths

#######################################################
## MUST match C API in
## ../../../main/clippath.c and
## ../../../include/R_ext/GraphicsEngine.h
#######################################################

.clipPath <- function(fun, rule) {
    attr(fun, "rule") <- .ruleIndex(rule)
    fun
}

.setClipPath <- function(path, index) {
    .External(C_setClipPath, path, index)
}
#  File src/library/grDevices/R/cm.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
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

cm <- function(x) 2.54*x

#  File src/library/grDevices/R/colorRamp.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
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


colorRampPalette <- function(colors,...)
{
    ramp <- colorRamp(colors,...)
    function(n) {
        x <- ramp(seq.int(0, 1, length.out = n))
        if (ncol(x) == 4L)
            rgb(x[, 1L], x[, 2L], x[, 3L], x[, 4L], maxColorValue = 255)
        else rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255)
    }

}

colorRamp <- function(colors, bias = 1, space = c("rgb","Lab"),
                    interpolate = c("linear","spline"), alpha = FALSE)
{
    if (bias <= 0) stop("'bias' must be positive")
    if (!missing(space) && alpha)
        stop("'alpha' must be false if 'space' is specified")

    colors <- t(col2rgb(colors, alpha = alpha)/255)
    space <- match.arg(space)
    interpolate <- match.arg(interpolate)

    if (space == "Lab")
        colors <- convertColor(colors, from = "sRGB", to = "Lab")


    interpolate <- switch(interpolate,
                          linear = stats::approxfun,
                          spline = stats::splinefun)

    if((nc <- nrow(colors)) == 1L) {
        colors <- colors[c(1L, 1L) ,]
        nc <- 2L
    }
    x <- seq.int(0, 1, length.out = nc)^bias
    palette <- c(interpolate(x, colors[, 1L]),
                 interpolate(x, colors[, 2L]),
                 interpolate(x, colors[, 3L]),
                 if(alpha) interpolate(x, colors[, 4L]))

    roundcolor <- function(rgb) ## careful to preserve matrix:
	pmax(pmin(rgb, 1), 0)

    if (space == "Lab")
        function(x)
            roundcolor(convertColor(cbind(palette[[1L]](x),
                                          palette[[2L]](x),
                                          palette[[3L]](x),
                                          if(alpha) palette[[4L]](x)),
                                    from = "Lab", to = "sRGB"))*255
    else
        function(x)
            roundcolor(cbind(palette[[1L]](x),
                             palette[[2L]](x),
                             palette[[3L]](x),
                             if(alpha) palette[[4L]](x)))*255
}
#  File src/library/grDevices/R/colorstuff.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2023 The R Core Team
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

colours <- colors <- function(distinct = FALSE)
{
    c <- .Call(C_colors)
    if(distinct) c[!duplicated(t(col2rgb(c)))] else c
}

col2rgb <- function(col, alpha = FALSE)
{
    ## R-level trap for now.
    if(any(as.character(col) %in% "0"))
        stop("numerical color values must be positive", domain = NA)
    if (is.factor(col)) col <- as.character(col)
    .Call(C_col2rgb, col, alpha)
}

gray <- function(level, alpha)
    .Call(C_gray, level, if (missing(alpha)) NULL else alpha)
grey <- gray

rgb <- function(red, green, blue, alpha, names = NULL, maxColorValue = 1)
{
    ## Only 'red' given
    if(missing(green) && missing(blue)) {
	if(is.matrix(red) || is.data.frame(red)) {
	    red <- data.matrix(red)
	    if(ncol(red) < 3L) stop("at least 3 columns needed")
	    green <- red[,2L]; blue <- red[,3L]; red <- red[,1L]
	}
    }

    .Call(C_rgb, red, green, blue, if (missing(alpha)) NULL else alpha,
          maxColorValue, names)
}

hsv <- function(h = 1, s = 1, v = 1, alpha)
    .Call(C_hsv, h, s, v, if(missing(alpha)) NULL else alpha)

hcl <- function (h = 0, c = 35, l = 85, alpha, fixup = TRUE)
    .Call(C_hcl, h, c, l, if(missing(alpha)) NULL else alpha, fixup)


rgb2hsv <- function(r, g = NULL, b = NULL, maxColorValue = 255)
{
    rgb <- if(is.null(g) && is.null(b)) as.matrix(r) else rbind(r, g, b)
    if(!is.numeric(rgb)) stop("rgb matrix must be numeric")
    d <- dim(rgb)
    if(d[1L] != 3L) stop("rgb matrix must have 3 rows")
    n <- d[2L]
    if(n == 0L) return(cbind(c(h = 1, s = 1, v = 1))[, 0L])
    ## else:
    rgb <- rgb/maxColorValue
    if(any(0 > rgb) || any(rgb > 1))
        stop("rgb values must be in [0, maxColorValue]")

    .Call(C_RGB2hsv, rgb)
}

## A quick little ''rainbow'' function -- improved by MM
## doc in	../man/palettes.Rd
rainbow <- function (n, s = 1, v = 1, start = 0, end = max(1,n - 1)/n,
                     alpha, rev = FALSE)
{
    if ((n <- as.integer(n[1L])) > 0) {
	if(start == end || any(c(start,end) < 0)|| any(c(start,end) > 1))
	    stop("'start' and 'end' must be distinct and in [0, 1].")
	alpha <- if(missing(alpha) || is.null(alpha)) NULL else rep_len(alpha, n)
	cols <- hsv(h = seq.int(start, (start > end)*1 + end,
				length.out = n) %% 1,
		    s, v, alpha)
        if(rev) rev(cols) else cols
    } else character()
}

topo.colors <- function (n, alpha, rev = FALSE)
{
    if ((n <- as.integer(n[1L])) > 0) {
	j <- n %/% 3
	k <- n %/% 3
	i <- n - j - k
	alpha <- if(missing(alpha) || is.null(alpha)) NULL else rep_len(alpha, n)
	cols <- c(if(i > 0) hsv(h = seq.int(from = 43/60, to = 31/60,
                                            length.out = i), alpha = alpha[seq_len(i)]),
                  if(j > 0) hsv(h = seq.int(from = 23/60, to = 11/60,
                                            length.out = j), alpha = alpha[i + seq_len(j)]),
                  if(k > 0) hsv(h = seq.int(from = 10/60, to =  6/60,
                                            length.out = k), alpha = alpha[i + j + seq_len(k)],
                                s = seq.int(from = 1, to = 0.3,
                                            length.out = k), v = 1))
        if(rev) rev(cols) else cols
    } else character()
}

terrain.colors <- function (n, alpha, rev = FALSE)
{
    if ((n <- as.integer(n[1L])) > 0) {
	k <- n%/%2
	h <- c(4/12, 2/12, 0/12)
	s <- c(1, 1, 0)
	v <- c(0.65, 0.9, 0.95)
	alpha <- if(missing(alpha) || is.null(alpha)) NULL else rep_len(alpha, n)
	cols <- c(hsv(h = seq.int(h[1L], h[2L], length.out = k),
                      s = seq.int(s[1L], s[2L], length.out = k),
                      v = seq.int(v[1L], v[2L], length.out = k), alpha = alpha[seq_len(k)]),
                  hsv(h = seq.int(h[2L], h[3L], length.out = n - k + 1)[-1L],
                      s = seq.int(s[2L], s[3L], length.out = n - k + 1)[-1L],
                      v = seq.int(v[2L], v[3L], length.out = n - k + 1)[-1L],
                      alpha = alpha[k + seq_len(n - k)]))
        if(rev) rev(cols) else cols
    } else character()
}

heat.colors <- function (n, alpha, rev = FALSE)
{
    if ((n <- as.integer(n[1L])) > 0) {
	j <- n %/% 4
	i <- n - j
	alpha <- if(missing(alpha) || is.null(alpha)) NULL else rep_len(alpha, n)
	cols <- c(rainbow(i, start = 0, end = 1/6, alpha = alpha[seq_len(i)]),
                  if (j > 0)
                      hsv(h = 1/6,
                          s = seq.int(from = 1-1/(2*j), to = 1/(2*j),
                                      length.out = j),
                          v = 1, alpha = alpha[i + seq_len(n - i)]))
        if(rev) rev(cols) else cols
    } else character()
}

cm.colors <- function (n, alpha, rev = FALSE)
{
    if ((n <- as.integer(n[1L])) > 0L) {
	even.n <- n %% 2L == 0L
	k <- n %/% 2L
	l1 <- k + 1L - even.n
	l2 <- n - k + even.n
	alpha <- if(missing(alpha) || is.null(alpha)) NULL else rep_len(alpha, n)
	cols <- c(if(l1 > 0L)
                      hsv(h =  6/12,
                          s = seq.int(.5, if(even.n) .5/k else 0,
                                      length.out = l1),
                          v = 1, alpha = alpha[seq_len(l1)]),
                  if(l2 > 1)
                      hsv(h = 10/12, s = seq.int(0, 0.5, length.out = l2)[-1L],
                          v = 1, alpha = alpha[seq_len(n)[-seq_len(l1)]]))
        if(rev) rev(cols) else cols
    } else character()
}

gray.colors <- function(n, start = 0.3, end = 0.9, gamma = 2.2, alpha,
                        rev = FALSE) {
    alpha <- if(missing(alpha) || is.null(alpha)) NULL else rep_len(alpha, n)
    cols <- gray(seq.int(from = start^gamma,
                         to = end^gamma, length.out = n)^(1/gamma),
                 alpha)
    if(rev) rev(cols) else cols
}

grey.colors <- gray.colors

##' match palette name (not exported)
palette.match <- function(pal) {
    fx <- function(x) tolower(gsub("[-, _, \\,, (, ), \\ , \\.]", "", x))
    charmatch(fx(pal), fx(names(.palette_colors_hex)))
}

palette <- function (value)
{
    ## if value missing return current palette (visibly)
    if (missing(value)) return(.Call(C_palette, character()))
    
    ## in case value is just a single string, select the corresponding set
    ## colors with "default" handled at C level
    if (length(value) == 1L && value != "default") {
        n <- palette.match(value)
	if (!is.na(n)) value <- .palette_colors_hex[[n]]
    }

    ## set new palette value, return old one invisibly
    ## if a .Device is open, record the .Call.graphics
    if (.Device == "null device") {
      invisible(.Call(C_palette, value))
    } else {
      invisible(.Call.graphics(C_palette, value))    
    }
}

## palette.colors() is a function for accessing the colors behind palette()
## directly. palette.pals() shows the available names (a la hcl.pals()).
palette.pals <- function() names(.palette_colors_hex)

palette.colors <- function(n = NULL, palette = "Okabe-Ito",
                           alpha, recycle = FALSE)
{
    ## number of colors
    if (!is.null(n)) {
        n <- as.integer(n[1L])
        if (n < 1L) return(character())
    }

    p <- palette.match(palette)
    if (is.na(p)) stop("'palette' does not match any given palette")
    if (p < 1L) stop("'palette' is ambiguous")
 
    ## select n colors from palette
    nc <- length(cols <- .palette_colors_hex[[p]])
    if (is.null(n))
        n <- nc
    else if(n > nc) {
        if(recycle) {
            cols <- cols[rep_len(seq_len(nc), n)]
        } else {
            warning(sprintf("'n' set to %s, the maximum available for %s palette",
                            nc, palette))
            n <- nc
        }
    }
    else if(n < nc)
        cols <- cols[seq_len(n)]

    ## add alpha if specified as number:
    if (!(missing(alpha) || is.null(alpha))) {
        alpha <- rep_len(alpha, n)
        alpha <- format(as.hexmode(round(alpha * 255 + 0.0001)),
                        width = 2L, upper.case = TRUE)
	cols <- paste0(cols, alpha)
    }
    cols
}

## underlying hex codes for palette color sets
.palette_colors_hex <- list(
    ## default in R <= 3.6.x
    "R3" = c("#000000", "#FF0000", "#00CD00", "#0000FF",
             "#00FFFF", "#FF00FF", "#FFFF00", "#BEBEBE"),
    ## rgb(
    ##     r = c(0, 255,   0,   0,   0, 255, 255, 190),
    ##     g = c(0,   0, 205,   0, 255,   0, 255, 190),
    ##     b = c(0,   0,   0, 255, 255, 255,   0, 190),
    ##     maxColorValue = 255
    ## ),

    ## new default in R >= 4.0.0		       
    "R4" = c("#000000", "#DF536B", "#61D04F", "#2297E6",
             "#28E2E5", "#CD0BBC", "#F5C710", "#9E9E9E"),
    ## hcl(h = c(0,   5, 125, 245, 195, 315,  65,   0),
    ##     c = c(0, 100,  90,  85,  63, 105,  94,   0),
    ##     l = c(0,  55,  75,  60,  82,  48,  82,  65)
    ## ),

    ## scales::hue_pal (Hadley Wickham)
    ## re-ordered for RGBCMY plus black/gray
    "ggplot2" = c("#000000", "#F8766D", "#00BA38", "#619CFF",
                  "#00BFC4", "#F564E3", "#B79F00", "#9E9E9E"),
    ## hcl(
    ##     h = c(0,  15, 135, 255, 195, 315,  75,   0),
    ##     c = c(0, 100, 100, 100, 100, 100, 100,   0),
    ##     l = c(0,  65,  65,  65,  65,  65,  65,  65)
    ## ),

    ## Masataka Okabe & Kei Ito
    ## http://jfly.iam.u-tokyo.ac.jp/color/
    "Okabe-Ito" = c(black = "#000000", orange = "#E69F00", skyblue = "#56B4E9", 
        bluishgreen = "#009E73", yellow = "#F0E442", blue = "#0072B2", 
        vermillion = "#D55E00", reddishpurple = "#CC79A7", gray = "#999999"),

    ## ColorBrewer.org (Mark A. Harrower & Cynthia A. Brewer)
    ## http://ColorBrewer2.org/
    "Accent" = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", 
        "#F0027F", "#BF5B17", "#666666"),
    "Dark 2" = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
        "#E6AB02", "#A6761D", "#666666"),
    "Paired" = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
        "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
        "#B15928"),
    "Pastel 1" = c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6",
        "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2"),
    "Pastel 2" = c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9",
        "#FFF2AE", "#F1E2CC", "#CCCCCC"),
    "Set 1" = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
        "#FFFF33", "#A65628", "#F781BF", "#999999"),
    "Set 2" = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
        "#FFD92F", "#E5C494", "#B3B3B3"),
    "Set 3" = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
        "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5",
        "#FFED6F"),

    ## Tableau 10 (Maureen Stone & Cristy Miller)
    ## https://www.tableau.com/about/blog/2016/7/colors-upgrade-tableau-10-56782
    "Tableau 10" = c(blue = "#4E79A7", orange = "#F28E2B", red = "#E15759",
        lightteal = "#76B7B2", green = "#59A14F", yellow = "#EDC948",
        purple = "#B07AA1", pink = "#FF9DA7", brown = "#9C755F",
        lightgray = "#BAB0AC"),
    "Classic Tableau" = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
        "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"),
	
    ## Polychrome
    ## (Kevin R. Coombes, Guy Brock, Zachary B. Abrams, Lynne V. Abruzzo)
    ## https://doi.org/10.18637/jss.v090.c01
    "Polychrome 36" = c(darkpurplishgray = "#5A5156", purplishwhite = "#E4E1E3",
        vividred = "#F6222E", vividpurple = "#FE00FA",
        vividyellowishgreen = "#16FF32",
	strongpurplishblue = "#3283FE", vividorangeyellow = "#FEAF16",
	vividpurplishred = "#B00068", brilliantgreen = "#1CFFCE",
	vividyellowgreen = "#90AD1C", vividblue = "#2ED9FF",
        brilliantpurple = "#DEA0FD",
	vividviolet = "#AA0DFE", strongpink = "#F8A19F", strongblue = "#325A9B",
	strongreddishorange = "#C4451C", vividgreen = "#1C8356",
	lightolivebrown = "#85660D", vividreddishpurple = "#B10DA1",
	vividgreenishyellow = "#FBE426", vividyellowishgreen = "#1CBE4F", 
        vividred = "#FA0087", vividpurplishred = "#FC1CBF",
        paleyellow = "#F7E1A0", 
        strongreddishpurple = "#C075A6", vividviolet = "#782AB6",
	vividyellowgreen = "#AAF400", verylightblue = "#BDCDFF",
	strongreddishbrown = "#822E1C", verylightyellowishgreen = "#B5EFB5", 
        verylightbluishgreen = "#7ED7D1", deepgreenishblue = "#1C7F93", 
        vividpurple = "#D85FF7", deeppurple = "#683B79",
        brilliantblue = "#66B0FF", 
        vividviolet = "#3B00FB"),
    "Alphabet" = c(amethyst = "#AA0DFE", blue = "#3283FE", caramel = "#85660D", 
        damson = "#782AB6", ebony = "#565656", forest = "#1C8356",
        green = "#16FF32", honey = "#F7E1A0", iron = "#E2E2E2",
        jade = "#1CBE4F", kingcrab = "#C4451C",  lavender = "#DEA0FD",
        magenta = "#FE00FA", navy = "#325A9B", orange = "#FEAF16",
        pink = "#F8A19F", quagmire = "#90AD1C", red = "#F6222E", 
        sea = "#1CFFCE", turquoise = "#2ED9FF", ultraviolet = "#B10DA1", 
        violet = "#C075A6", wine = "#FC1CBF", xanthin = "#B00068",
        yellow = "#FBE426", zinnia = "#FA0087")
)

## An unexported version that works with internal representation as 'rcolor'
## We could avoid this if we knew at R level whether the display list was
## enabled or inhibited: but we do need to record a call to C_palette2.
recordPalette <- function()
    .Call.graphics(C_palette2, .Call(C_palette2, NULL))

#  File src/library/grDevices/R/convertColor.R
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

## Simplified and faster version of `ifelse` with constraints:
##
## * test, yes, and no must be the same length
## * test must be logical
## * yes and no must be numeric
## * if test is NA, then `no` is returned, which is particularly okay if it is
##   the case that when test is NA so is NO, as tends to be the case here.

.ifelse <- function(test, yes, no) {
  test.w <- which(test)
  no[test.w] <- yes[test.w]
  no
}
## Benchmarks show x ^ 3 is much slower than x * x * x

pow3 <- function(x) x * x * x

## easyRGB scales Y=100 for white
## brucelindbloom uses XYZ in [0,1], so multiply by 100 to convert

## white points in xyY format (Y=1 omitted)
white.points <- cbind(A = c(x = 0.44757, y = 0.40745),
                      B = c(x = 0.34842, y = 0.35161),
                      C = c(x = 0.31006, y = 0.31616),
                      D50 = c(x = 0.34574, y = 0.35867),
                      D55 = c(x = 0.33250, y = 0.34761),
                      D65 = c(x = 0.3137, y = 0.3291),
                      E = c(x = 1/3, y = 1/3))
## converting these:
c2to3 <- function(col) c(col[1L]/col[2L], 1, (1 - sum(col[1L:2L]))/col[2L])

## http://www.brucelindbloom.com/index.html?Equations.html


make.rgb <-
    function(red, green, blue, name = NULL, white = "D65", gamma = 2.2)
{
    whitexyz <- c2to3(white.points[, white])
    rgb <- rbind(c2to3(red),
                 c2to3(green),
                 c2to3(blue))
    S <- drop(whitexyz %*% solve(rgb))
    M <- S * rgb

    if (is.numeric(gamma) && length(gamma) == 1) {
        dogamma <- function(x) x %^% gamma
        ungamma <- function(x) x %^% (1/gamma)
    } else if (gamma == "sRGB") {
        dogamma <- function(x) .ifelse(x < 0.04045,
                                       x/12.92,
                                       ((x+0.055)/1.055)^2.4)
        ungamma <- function(x) .ifelse(x <= 0.0031308,
                                       12.92*x,
                                       1.055*x %^% (1/2.4)-0.055)
    } else stop("'gamma' must be a scalar or 'sRGB'")

    toXYZ <- function(rgb,...) { dogamma(rgb) %*% M }
    toRGB <- function(xyz,...) {
      res <- ungamma(xyz %*% solve(M))
      # for backward compatibily, return vector if input is vector
      if(nrow(res) == 1L) res[1L, ,drop=TRUE] else res
    }
    if (is.null(name)) name <- deparse(sys.call())[1L]
    RGBcolorConverter(toXYZ = toXYZ, fromXYZ = toRGB, gamma = gamma,
                      white = white, name = name, vectorized = TRUE)
}

print.colorConverter <- function(x,...) {
    cat(gettextf("Color space converter: %s", x$name), "\n", sep = "")
    if (!is.null(x$reference.white))
        cat(gettextf("Reference white: %s", x$reference.white), "\n", sep = "")
    invisible(x)
}

print.RGBcolorConverter <- function(x,...) {
    print.colorConverter(x, ...)
    if (!is.null(x$gamma))
        cat(gettextf("display gamma = %s", format(x$gamma)), "\n", sep = "")
    invisible(x)
}

chromaticAdaptation <- function(xyz, from, to) {
    ## Von Kries scaling algorithm
    Ma <- matrix(c( 0.40024, -0.22630, 0.,
                    0.70760,  1.16532, 0.,
                   -0.08081,  0.04570, 0.91822), nrow = 3L, byrow = TRUE)
    from.cone <- drop(from %*% Ma)
    to.cone   <- drop(to %*% Ma)
    ## M <- Ma %*% diag(to.cone/from.cone) %*% solve(Ma)
    M <- (Ma * rep(to.cone/from.cone, each=3)) %*% solve(Ma)
    xyz %*% M
}

vectorizeConverter <- function(converter) {
    function(color, white) {
        nr <- if(is.null(nrow(color))) length(color) else nrow(color)
        res <- apply(color, 1L, converter, white)
        if (is.null(nrow(res)))
            matrix(res, nrow=nr)
        else
            t(res)
    }
}

colorConverter <-
    function(toXYZ, fromXYZ, name, white = NULL, vectorized=FALSE)
{
    stopifnot(isTRUE(as.logical(vectorized) %in% c(TRUE, FALSE)))

    toXYZv <- if (!vectorized)
                  vectorizeConverter(toXYZ)
              else toXYZ
    fromXYZv <- if (!vectorized)
                    vectorizeConverter(fromXYZ)
              else fromXYZ

    ## redundant `white` / `reference.white` for backwards-compatibility
    ## see https://bugs.r-project.org/show_bug.cgi?id=17473
    rval <- list(toXYZ = toXYZv, fromXYZ = fromXYZv,
                 name = name, white = white, reference.white = white)
    class(rval) <- "colorConverter"
    rval
}
RGBcolorConverter <- function(..., gamma) {
    rval <- colorConverter(...)
    rval[['gamma']] <- gamma
    class(rval) <- c("RGBcolorConverter", class(rval))
    rval
}
colorspaces <-
    list("XYZ" =
         colorConverter(toXYZ = function(x,w) x,
                        fromXYZ = function(x,w) x,
                        white = NULL,name = "XYZ",vectorized = TRUE),

         "Apple RGB" =
         make.rgb(red = c(0.6250,0.3400),
                  green = c(0.2800,0.5950),
                  blue = c(0.1550,0.0700),gamma = 1.8,
                  white = "D65", name = "Apple RGB"),

         "sRGB" =
         make.rgb(red = c(0.6400, 0.3300),
                  green = c(0.3000,0.6000),
                  blue = c(0.1500,0.0600), gamma = "sRGB",
                  white = "D65", name = "sRGB"),

         "CIE RGB" =
         make.rgb(red = c(0.7350,0.2650),
                  green = c(0.2740,0.7170),
                  blue = c(0.1670,0.0090), gamma = 2.2,
                  white = "E", name = "CIE RGB"),

         "Lab" =
         colorConverter(fromXYZ = function(XYZ, white) {
             stopifnot(length(XYZ) == 3 || ncol(XYZ) == 3L)
             white <- rep(white, length.out=3L)
             if (is.null(nrow(XYZ))) XYZ <- matrix(XYZ, nrow = 1L)

             epsilon <- 216/24389
             kappa <- 24389/27

             xyzr <- cbind(XYZ[,1L] / white[1L],
                           XYZ[,2L] / white[2L],
                           XYZ[,3L] / white[3L])
             fxyz <- .ifelse(xyzr <= epsilon, (kappa*xyzr+16)/116, xyzr^(1/3))

             res <- cbind(L = 116*fxyz[,2L]-16,
                          a = 500*(fxyz[,1L]-fxyz[,2L]),
                          b = 200*(fxyz[,2L]-fxyz[,3L]))
             if(nrow(res) == 1L) res[1L, ,drop=TRUE] else res
         },
         toXYZ = function(Lab, white) {
             stopifnot(ncol(Lab) == 3L || length(Lab)==3)
             white <- rep(white, length.out=3L)
             if (is.null(nrow(Lab))) Lab <- matrix(Lab, nrow = 1L)

             epsilon <- 216/24389
             kappa <- 24389/27

             L <- Lab[,1L]
             yr <- .ifelse(L < kappa*epsilon,
                           L/kappa, pow3((L+16)/116))
             fy <- (.ifelse(yr <= epsilon, kappa*yr, L)+16)/116
             fx <- Lab[,2L]/500+fy
             fz <- fy-Lab[,3L]/200

             fz3 <- pow3(fz)
             fx3 <- pow3(fx)
             zr <- .ifelse(fz3 <= epsilon, (116*fz-16)/kappa, fz3)
             xr <- .ifelse(fx3 <= epsilon, (116*fx-16)/kappa, fx3)

             res <- cbind(X = xr*white[1], Y = yr*white[2], Z = zr*white[3])

             if(nrow(res) == 1L) res[1L, ,drop=TRUE] else res
         }, name = "Lab", white = NULL, vectorized = TRUE),

         "Luv" =
         colorConverter(fromXYZ = function(XYZ, white) {
             epsilon <- 216/24389
             kappa <- 24389/27

             yr <- XYZ[,2L]/white[2L]

             denom  <- rowSums(cbind(XYZ[,1L], XYZ[,2L]*15, XYZ[,3L]*3))
             wdenom <- sum(white*c(1,15,3))

             one <- rep_len(1, length(denom))
             u1 <- .ifelse(denom == 0, one, 4*XYZ[,1L]/denom)
             v1 <- .ifelse(denom == 0, one, 9*XYZ[,2L]/denom)
             ur <- 4*white[1L]/wdenom
             vr <- 9*white[2L]/wdenom

             L <- .ifelse(yr <= epsilon, kappa*yr, 116*(yr^(1/3))-16)
             res <- cbind(L = L, u = 13*L*(u1-ur), v = 13*L*(v1-vr))
             if(nrow(res) == 1L) res[1L, ,drop=TRUE] else res
         }, toXYZ = function(Luv,white) {
             epsilon <- 216/24389
             kappa <- 24389/27

             u0 <- 4*white[1L]/(white[1L]+15*white[2L]+3*white[3L])
             v0 <- 9*white[2L]/(white[1L]+15*white[2L]+3*white[3L])

             L <- Luv[,1L]
             Y <- .ifelse(L <= kappa*epsilon,
                          L/kappa, pow3((L+16)/116))
             a <- (52*L/(Luv[,2L]+13*L*u0)-1)/3
             b <- -5*Y
             c <- -1/3
             d <- Y*(39*L/(Luv[,3L]+13*L*v0)-5)

             X <- (d-b)/(a-c)
             Z <- X*a+b

             res <- cbind(X = X,Y = Y,Z = Z)

             res[which(L == 0L),] <- c(0,0,0)
             if(nrow(res) == 1L) res[1L, ,drop=TRUE] else res
         }, name = "Luv", white = NULL, vectorized = TRUE)

         ) # colorspaces


`%^%` <- function(a,b) sign(a) * (abs(a) ^ b)

convertColor <-
    function(color, from, to,
             from.ref.white = NULL, to.ref.white = NULL,
             scale.in = 1, scale.out = 1, clip = TRUE)
{
  if (is.data.frame(color)) {
      color <- as.matrix(color)
  }

  if (is.character(from))
      from <- colorspaces[[match.arg(from, names(colorspaces))]]
  if (!inherits(from,"colorConverter"))
      stop("'from' must be a \"colorConverter\" object or a string")
  if (is.character(to))
      to <- colorspaces[[match.arg(to, names(colorspaces))]]
  if (!inherits(to,"colorConverter"))
      stop("'to' must be a \"colorConverter\" object or a string")

  ## Need a reference white. If both the definition and the argument
  ## specify one they must agree.

  if (is.null(from.ref.white))
      from.ref.white <- from$reference.white
  else if (!is.null(from$reference.white) &&
           from.ref.white != from$reference.white)
      stop(gettextf("'from.ref.white' disagrees with definition of %s",
                    from$name), domain = NA)

  if (is.null(to.ref.white))
      to.ref.white <- to$reference.white
  else if (!is.null(to$reference.white) &&
           to.ref.white != to$reference.white)
      stop(gettextf("'to.ref.white' disagrees with definition of %s",
                    to$name), domain = NA)

  if (is.null(to.ref.white) && is.null(from.ref.white))
      to.ref.white <- from.ref.white <- "D65"

  if (is.null(to.ref.white))
      to.ref.white <- from.ref.white
  if (is.null(from.ref.white))
      from.ref.white <- to.ref.white

  from.ref.white <- c2to3(white.points[, from.ref.white])
  to.ref.white   <- c2to3(white.points[, to.ref.white])

  if (is.null(nrow(color)))
    color <- matrix(color, nrow = 1L)

  if (!is.null(scale.in))
      color <- color/scale.in

  trim <- function(rgb) {
      rgb <- round(rgb,5)
      if (is.na(clip))
          rgb[rgb < 0 | rgb > 1] <- NaN
      else if(clip) {
          rgb[rgb < 0] <- 0
          rgb[rgb > 1] <- 1
      }
      rgb
  }

  xyz <- from$toXYZ(color, from.ref.white)

  if (is.null(nrow(xyz)))
    xyz <- matrix(xyz, nrow = 1L)

  if (!isTRUE(all.equal(from.ref.white, to.ref.white))) {
      mc <- match.call()
      if (is.null(mc$from.ref.white) || is.null(mc$to.ref.white)) {
          warning("color spaces use different reference whites")
      }
      xyz <- chromaticAdaptation(xyz, from.ref.white, to.ref.white)
  }

  rval <- to$fromXYZ(xyz, to.ref.white)

  if(is.null(nrow(rval)))
      rval <- t(rval)

  if (inherits(to,"RGBcolorConverter"))
      rval <- trim(rval)

  if (is.null(scale.out))
      rval
  else
      rval*scale.out
}

##' @title Modify a vector of colors by "screwing" any of (r,g,b,alpha)
##'   by multification by a factor
##' @param col vector of colors, in any format that col2rgb() accepts
##' @param alpha.f factor modifying the opacity alpha; typically in [0,1]
##' @param red.f   factor modifying "red"ness
##' @param green.f factor modifying "green"ness
##' @param blue.f  factor modifying "blue"ness
##' @return From rgb(), a color vector of the same length as 'col'.
##' @author Thomas Lumley, Luke Tierney, Martin Maechler, Duncan Murdoch...
adjustcolor <- function(col, alpha.f = 1, red.f = 1, green.f = 1,
                        blue.f = 1, offset = c(0,0,0,0),
                        transform = diag(c(red.f, green.f, blue.f, alpha.f)))
{
    stopifnot(exprs = {
        length(offset) %% 4L == 0L
        !is.null(d <- dim(transform))
        d == c(4L, 4L)
    })
    x <- col2rgb(col, alpha = TRUE)/255
    x[] <- pmax(0, pmin(1,
                        transform %*% x +
                        matrix(offset, nrow = 4L, ncol = ncol(x))))
    rgb(x[1L,], x[2L,], x[3L,], x[4L,])
}
#  File src/library/grDevices/R/device.R
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


dev.interactive <- function(orNone = FALSE)
{
    if(!interactive()) return(FALSE)
    .known_interactive_devices <- deviceIsInteractive()
    if(.Device %in% .known_interactive_devices) return(TRUE)
    if(!(orNone && .Device == "null device")) return(FALSE)
    ## at this point we have mo active device.
    newdev <- getOption("device")
    if(is.character(newdev)) newdev %in% .known_interactive_devices
    else { # a function
        if(.Platform$OS.type == "windows") identical(newdev, windows)
        else identical(newdev, X11) || identical(newdev, quartz)
    }
}

deviceIsInteractive <- local({
    .known_interactive_devices <- 
        c("X11", "X11cairo", "quartz", "windows", "JavaGD", "CairoWin",
          "CairoX11")
    function(name = NULL) {
        if(length(name)) {
            if(!is.character(name))
                stop("'name' must be a character vector")
            .known_interactive_devices <<-
                unique(c(.known_interactive_devices, name))
            invisible(.known_interactive_devices)
        } else
            .known_interactive_devices
    }
})


`%||%` <- function(L,R) if(is.null(L)) R else L

dev.list <- function()
{
    n <- get0(".Devices") %||% list("null device")
    n <- unlist(n)
    i <- seq_along(n)[n != ""]
    names(i) <- n[i]
    i <- i[-1L]
    if(length(i) == 0L) NULL else i
}

dev.cur <- function()
{
    if(!exists(".Devices"))
	.Devices <- list("null device")
    num.device <- .External(C_devcur)
    names(num.device) <- .Devices[[num.device]]
    num.device
}

dev.set <-
    function(which = dev.next())
{
    which <- .External(C_devset, as.integer(which))
    names(which) <- .Devices[[which]]
    which
}

dev.next <-
    function(which = dev.cur())
{
    if(!exists(".Devices"))
	.Devices <- list("null.device")
    num.device <- .External(C_devnext, as.integer(which))
    names(num.device) <- .Devices[[num.device]]
    num.device
}

dev.prev <-
    function(which = dev.cur())
{
    if(!exists(".Devices"))
	.Devices <- list("null device")
    num.device <- .External(C_devprev, as.integer(which))
    names(num.device) <- .Devices[[num.device]]
    num.device
}

dev.off <-
    function(which = dev.cur())
{
    if(which == 1)
	stop("cannot shut down device 1 (the null device)")
    .External(C_devoff, as.integer(which))
    dev.cur()
}

dev.copy <- function(device, ..., which = dev.next())
{
    if(!missing(which) && !missing(device))
	stop("cannot supply 'which' and 'device' at the same time")
    old.device <- dev.cur()
    if(old.device == 1)
	stop("cannot copy from the null device")
    if(missing(device)) {
	if(which == 1)
	    stop("cannot copy to the null device")
	else if(which == dev.cur())
	    stop("cannot copy device to itself")
	dev.set(which)
    }
    else {
	if(!is.function(device))
	    stop("'device' should be a function")
	else device(...)
    }
    ## protect against failure
    on.exit(dev.set(old.device))
    .External(C_devcopy, old.device)
    on.exit()
    dev.cur()
}

dev.print <- function(device = postscript, ...)
{
    current.device <- dev.cur()
    nm <- names(current.device)[1L]
    if(nm == "null device") stop("no device to print from")
    if(!dev.displaylist())
        stop("can only print from a screen device")
    oc <- match.call()
    oc[[1L]] <- quote(grDevices::dev.copy)
    oc$device <- device
    din <- graphics::par("din"); w <- din[1L]; h <- din[2L]
    if(missing(device)) { ## safe way to recognize postscript
        if(is.null(oc$file)) oc$file <- ""
        hz0 <- oc$horizontal
        hz <- if(is.null(hz0)) ps.options()$horizontal else eval.parent(hz0)
        paper <- oc$paper
        if(is.null(paper)) paper <- ps.options()$paper
        if(paper == "default") paper <- getOption("papersize")
        paper <- tolower(paper)
        switch(paper,
               a4 = 	 {wp <- 8.27; hp <- 11.69},
               legal =	 {wp <- 8.5;  hp <- 14.0},
               executive={wp <- 7.25; hp <- 10.5},
               { wp <- 8.5; hp <- 11}) ## default is "letter"

        wp <- wp - 0.5; hp <- hp - 0.5  # allow 0.25" margin on each side.
        if(!hz && is.null(hz0) && h < wp && wp < w && w < hp) {
            ## fits landscape but not portrait
            hz <- TRUE
        } else if (hz && is.null(hz0) && w < wp && wp < h && h < hp) {
            ## fits portrait but not landscape
            hz <- FALSE
        } else {
            h0 <- if(hz) wp else hp
            if(h > h0) { w <- w * h0/h; h <- h0 }
            w0 <- if(hz) hp else wp
            if(w > w0) { h <- h * w0/w; w <- w0 }
        }
        if(is.null(oc$pointsize)) {
            pt <- ps.options()$pointsize
            oc$pointsize <- pt * w/din[1L]
        }
        if(is.null(hz0)) oc$horizontal <- hz
        if(is.null(oc$width)) oc$width <- w
        if(is.null(oc$height)) oc$height <- h
    } else {
        devname <- deparse(substitute(device))
        if(devname %in% c("png", "jpeg", "bmp") &&
           is.null(oc$width) && is.null(oc$height))
            warning("need to specify one of 'width' and 'height'")
        if(is.null(oc$width))
            oc$width <- if(!is.null(oc$height)) w/h * eval.parent(oc$height) else w
        if(is.null(oc$height))
            oc$height <- if(!is.null(oc$width)) h/w * eval.parent(oc$width) else h
    }
    ## protect against failure (PR#9801)
    on.exit(dev.set(current.device))
    dev.off(eval.parent(oc))
}

dev.copy2eps <- function(...)
{
    current.device <- dev.cur()
    nm <- names(current.device)[1L]
    if(nm == "null device") stop("no device to print from")
    if(!dev.displaylist())
        stop("can only print from a screen device")
    oc <- match.call()
    oc[[1L]] <- quote(grDevices::dev.copy)
    oc$device <- postscript
    oc$onefile <- FALSE
    oc$horizontal <- FALSE
    if(is.null(oc$paper))
        oc$paper <- "special"
    din <- dev.size("in"); w <- din[1L]; h <- din[2L]
    if(is.null(oc$width))
        oc$width <- if(!is.null(oc$height)) w/h * eval.parent(oc$height) else w
    if(is.null(oc$height))
        oc$height <- if(!is.null(oc$width)) h/w * eval.parent(oc$width) else h
    if(is.null(oc$file)) oc$file <- "Rplot.eps"
    ## protect against failure (PR#9801)
    on.exit(dev.set(current.device))
    dev.off(eval.parent(oc))
}

dev.copy2pdf <- function(..., out.type = "pdf")
{
    out.type <- match.arg(out.type, c("pdf", "quartz", "cairo"))
    current.device <- dev.cur()
    nm <- names(current.device)[1L]
    if(nm == "null device") stop("no device to print from")
    if(!dev.displaylist())
        stop("can only print from a screen device")
    oc <- match.call()
    oc[[1L]] <- quote(grDevices::dev.copy)
    if(out.type == "quartz" && capabilities("aqua")) {
        oc$device <- quartz
        oc$type <- "pdf"
    } else if(out.type == "cairo" && capabilities("cairo")) {
        oc$device <- cairo_pdf
        oc$onefile <- FALSE # future-proofing
    } else {
        oc$device <- pdf
        ## the defaults in pdf() are all customizable, so we override
        ## even those which are the ultimate defaults.
        oc$onefile <- FALSE
        if(is.null(oc$paper)) oc$paper <- "special"
    }
    oc$out.type <- NULL
    din <- dev.size("in"); w <- din[1L]; h <- din[2L]
    if(is.null(oc$width))
        oc$width <- if(!is.null(oc$height)) w/h * eval.parent(oc$height) else w
    if(is.null(oc$height))
        oc$height <- if(!is.null(oc$width)) h/w * eval.parent(oc$width) else h
    if(is.null(oc$file)) oc$file <- "Rplot.pdf"
    ## protect against failure (PR#9801)
    on.exit(dev.set(current.device))
    dev.off(eval.parent(oc))
}

dev.control <- function(displaylist = c("inhibit", "enable"))
{
    if(dev.cur() <= 1)
        stop("dev.control() called without an open graphics device")
    if(!missing(displaylist)) {
        displaylist <- match.arg(displaylist)
	.External(C_devcontrol, displaylist == "enable")
    } else stop("argument is missing with no default")
    invisible()
}

dev.displaylist <- function()
{
    if(dev.cur() <= 1)
        stop("dev.displaylist() called without an open graphics device")
    .External(C_devdisplaylist)
}

## This records graphics ops and manipulates visibility, so needs to stay .Internal
recordGraphics <- function(expr, list, env)
  .Internal(recordGraphics(substitute(expr), list, env))

graphics.off <- function ()
{
    while ((which <- dev.cur()) != 1) dev.off(which)
    invisible()
}

dev.new <- function(..., noRStudioGD = FALSE)
{
    dev <- getOption("device")
    if(!is.character(dev) && !is.function(dev))
        stop("invalid setting for 'getOption(\"device\")'")
    if(noRStudioGD && is.character(dev) && dev == "RStudioGD")
        dev <- .select_device()
    if(is.character(dev)) {
        ## this is documented to be searched for from workspace,
        ## then in the grDevices namespace.
        ## We could restrict the search to functions, but the C
        ## code in devices.c does not.
        dev <- get0(dev, .GlobalEnv) %||%
               get0(dev, asNamespace("grDevices")) %||%
            stop(gettextf("device '%s' not found", dev), domain=NA)
    }
    ## only include named args in the devices's arglist
    a <- list(...)
    a2 <- names(formals(dev))
    a <- a[names(a) %in% a2]
    if(identical(dev, pdf)) {
        ## Take care not to open device on top of another.
        if(is.null(a[["file"]]) && file.exists("Rplots.pdf")) {
            fe <- file.exists(tmp <- paste0("Rplots", 1L:999, ".pdf"))
            if(all(fe)) stop("no suitable unused file name for pdf()")
            message(gettextf("dev.new(): using pdf(file=\"%s\")", tmp[!fe][1L]),
                    domain=NA)
            a$file <- tmp[!fe][1L]
        }
    } else if(identical(dev, postscript)) {
        ## Take care not to open device on top of another.
        if(is.null(a[["file"]]) && file.exists("Rplots.ps")) {
            fe <- file.exists(tmp <- paste0("Rplots", 1L:999, ".ps"))
            if(all(fe)) stop("no suitable unused file name for postscript()")
            message(gettextf("dev.new(): using postscript(file=\"%s\")",
                             tmp[!fe][1L]), domain=NA)
            a$file <- tmp[!fe][1L]
        }
    } else if (!is.null(a[["width"]]) && !is.null(a[["height"]]) &&
               (identical(dev, png) || identical(dev, jpeg) ||
                identical(dev, bmp) || identical(dev, tiff))) {
        ## some people want dev.new(width=12, height=7) to be portable
        if(is.null(a[["units"]]) && is.null(a[["res"]])) {
            a$units <- "in"
            a$res <- 72
        }
    }
    do.call(dev, a)
}

## Check string to contain at most one valid integer format (starting with '%')
## 's' is typically the file(name) or title argument of a graphical device.
checkIntFormat <- function(s)
{
    ## OK if no unescaped %, so first remove those
    s <- gsub("%%", "", s, fixed=TRUE)
    if(length(grep("%", s)) == 0L) return(TRUE)
    ## now remove at most one valid(ish) integer format
    s <- sub("%[#0 ,+-]*[0-9.]*[diouxX]", "", s)
    length(grep("%", s, fixed=TRUE)) == 0L
}

devAskNewPage <- function(ask=NULL) .External2(C_devAskNewPage, ask)

dev.size <- function(units = c("in", "cm", "px"))
{
    units <- match.arg(units)
    size <- .External(C_devsize)
    if(units == "px") size else size * graphics::par("cin")/graphics::par("cra") *
        if(units == "cm") 2.54 else 1
}

dev.hold <- function(level = 1L) .External(C_devholdflush, max(0L, level))
dev.flush <- function(level = 1L) .External(C_devholdflush, -max(0L, level))

dev.capture <- function(native = FALSE) .External(C_devcapture, native)

dev.capabilities <- function(what = NULL)
{
    ncap <- 13
    template <- vector("list", ncap)
    capabilities <- .External(C_devcap, template)
    ## The device may have filled in some capabilities so check it is still
    ## the right sort of structure
    if (!(is.list(capabilities) &&
          length(capabilities) == ncap &&
          all(sapply(capabilities, class) == "integer")))
        stop("Invalid capabilities - alert the device maintainer")

    z <- vector("list", ncap)
    names(z) <- c("semiTransparency",
                  "transparentBackground",
                  "rasterImage",
                  "capture",
                  "locator",
                  "events",
                  "patterns",
                  "clippingPaths",
                  "masks",
                  "compositing",
                  "transformations",
                  "paths",
                  "glyphs")
    z[[1L]] <- c(NA, FALSE, TRUE)[capabilities[[1L]] + 1L]
    z[[2L]] <- c(NA, "no", "fully", "semi")[capabilities[[2L]] + 1L]
    z[[3L]] <- c(NA, "no", "yes", "non-missing")[capabilities[[3L]] + 1L]
    z[[4L]] <- c(NA, FALSE, TRUE)[capabilities[[4L]] + 1L]
    z[[5L]] <- c(NA, FALSE, TRUE)[capabilities[[5L]] + 1L]
    z[[6L]] <- c( "", 
                  if (capabilities[[6L]][1L]) "MouseDown",
                  if (capabilities[[6L]][2L]) "MouseMove",
                  if (capabilities[[6L]][3L]) "MouseUp",
                  if (capabilities[[6L]][4L]) "Keybd" )[-1L]
    ## Patterns
    if (length(capabilities[[7]]) == 1 && is.na(capabilities[[7]])) {
        z[[7]] <- NA
    } else if (length(capabilities[[7]]) == 1 && capabilities[[7]] < 1) {
        z[[7]] <- FALSE
    } else if (length(capabilities[[7]]) > 0) {
        z[[7]] <- patternTypes[capabilities[[7]]]
    } else {
        z[[7]] <- NA
    }
    ## Clipping paths
    if (is.na(capabilities[[8]]))
        z[[8]] <- NA
    else 
        z[[8]] <- as.logical(capabilities[[8]])
    ## Masks
    if (length(capabilities[[9]]) == 1 && is.na(capabilities[[9]])) {
        z[[9]] <- NA
    } else if (length(capabilities[[9]]) == 1 && capabilities[[9]] < 1) {
        z[[9]] <- FALSE
    } else if (length(capabilities[[9]]) > 0) {
        z[[9]] <- maskTypes[capabilities[[9]]]
    } else {
        z[[9]] <- NA
    }
    ## Compositing operators
    if (length(capabilities[[10]]) == 1 && is.na(capabilities[[10]])) {
        z[[10]] <- NA
    } else if (length(capabilities[[10]]) == 1 && capabilities[[10]] < 1) {
        z[[10]] <- FALSE
    } else if (length(capabilities[[10]]) > 0) {
        z[[10]] <- compositingOperators[capabilities[[10]]]
    } else {
        z[[10]] <- NA
    }
    ## Transforms
    if (is.na(capabilities[[11]]))
        z[[11]] <- NA
    else 
        z[[11]] <- as.logical(capabilities[[11]])
    ## Paths
    if (is.na(capabilities[[12]]))
        z[[12]] <- NA
    else 
        z[[12]] <- as.logical(capabilities[[12]])
    ## Glyphs
    if (is.na(capabilities[[13]]))
        z[[13]] <- NA
    else 
        z[[13]] <- as.logical(capabilities[[13]])

    if (!is.null(what)) z[charmatch(what, names(z), 0L)] else z
}

## for use in dev.new and .onLoad
.select_device <- function() {
    ## Use device functions rather than names to make it harder to get masked.
    if(!nzchar(defdev <- Sys.getenv("R_DEFAULT_DEVICE"))) defdev <- pdf
    if(interactive()) {
        if(nzchar(intdev <- Sys.getenv("R_INTERACTIVE_DEVICE"))) intdev
        else {
            if(.Platform$OS.type == "windows") windows
            else {
                ## This detects if quartz() was built and if we are
                ## running at the macOS console (both of which have to
                ## be true under R.app).
                if(.Platform$GUI == "AQUA" ||.Call(C_makeQuartzDefault)) quartz
                else if(nzchar(Sys.getenv("DISPLAY"))
                        && .Platform$GUI %in% c("X11", "Tk")) X11
                else defdev
            }
        }
    } else defdev
}
#  File src/library/grDevices/R/gevents.R
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

setGraphicsEventHandlers <- function(which=dev.cur(),
				     ...)
    setGraphicsEventEnv(which, as.environment(list(...)))

setGraphicsEventEnv <- function(which=dev.cur(), env) {
    which <- as.integer(which)
    stopifnot(length(which) == 1)
    result <- getGraphicsEventEnv(which)
    env$which <- which
    .External2(C_setGraphicsEventEnv, which, env)
    invisible(result)
}

getGraphicsEventEnv <- function(which=dev.cur()) {
    which <- as.integer(which)
    stopifnot(length(which) == 1)
    .External2(C_getGraphicsEventEnv, which)
}

getGraphicsEvent <- function(prompt = "Waiting for input",
                 onMouseDown = NULL, onMouseMove = NULL, onMouseUp = NULL,
                 onKeybd = NULL, onIdle = NULL, consolePrompt = prompt) {
    if (!interactive()) return(NULL)
    if (!missing(prompt) || !missing(onMouseDown) || !missing(onMouseMove)
     || !missing(onMouseUp) || !missing(onKeybd) || !missing(onIdle)) {
        setGraphicsEventHandlers(prompt=prompt, onMouseDown=onMouseDown,
          onMouseMove=onMouseMove, onMouseUp=onMouseUp, onKeybd=onKeybd,
          onIdle=onIdle)
    }
    .External2(C_getGraphicsEvent, consolePrompt)
}

## Generate object containing information about typeset glyphs
## All x/y and width/height and anchors should be in "big" pts (1/72 inch)

mapCharWeight <- function(x) {
    if (is.na(x))
        x
    else 
        switch(as.character(x),
               normal=400,
               bold=700,
               as.numeric(x))
}

mapWeight <- function(x) {
    if (is.numeric(x)) {
        if (min(x, na.rm=TRUE) < 0 || max(x, na.rm=TRUE) > 1000)
            stop("Invalid span weight")
        x
    } else {
        sapply(x, mapCharWeight, USE.NAMES=FALSE)
    }   
}

fontStyles <- c("normal", "italic", "oblique")

mapStyle <- function(x) {
    ## NA passes through
    match(x, fontStyles)
}

invertStyle <- function(x) {
    fontStyles[x]
}

################################################################################
## glyph dimensions, anchors, and justification
glyphWidth <- function(w, label="width", left="left") {
    if (!length(w) ||
        length(w) != length(label) ||
        length(label) != length(left))
        stop("length of arguments must match (and be greater than 0)")
    w <- as.numeric(w)
    if (any(!is.finite(w))) stop("Invalid glyph width(s)")
    names(w) <- as.character(label)
    attr(w, "anchor") <- as.character(left)
    class(w) <- "GlyphWidth"
    w
}

glyphWidthLeft <- function(w, label) {
    if (!inherits(w, "GlyphWidth")) stop("Invalid glyph width")
    if (!label %in% names(w)) {
        warning("Unknown width; using left anchor")
        return("left")
    }
    which <- match(label, names(w))
    attr(w, "anchor")[which]
}

glyphHeight <- function(h, label="height", bottom="bottom") {
    if (!length(h) ||
        length(h) != length(label) ||
        length(label) != length(bottom))
        stop("length of arguments must match (and be greater than 0)")
    h <- as.numeric(h)
    if (any(!is.finite(h))) stop("Invalid glyph height(s)")
    names(h) <- as.character(label)
    attr(h, "anchor") <- as.character(bottom)
    class(h) <- "GlyphHeight"
    h
}

glyphHeightBottom <- function(h, label) {
    if (!inherits(h, "GlyphHeight")) stop("Invalid glyph height")
    if (!label %in% names(h)) {
        warning("Unknown height; using bottom anchor")
        return("bottom")
    }
    which <- match(label, names(h))
    attr(h, "anchor")[which]
}

glyphAnchor <- function(value, label) {
    if (!length(value) ||
        length(value) != length(label))
        stop("length of arguments must match (and be greater than 0)")
    value <- as.numeric(value)
    if (any(!is.finite(value))) stop("Invalid glyph anchor")
    names(value) <- as.character(label)
    class(value) <- "GlyphAnchor"
    value
}

glyphJust <- function(just, ...) {
    UseMethod("glyphJust")
}
glyphJust.GlyphJust <- function(just, ...) {
    just
}
glyphJust.character <- function(just, ...) {
    class(just) <- "GlyphJust"
    just
}
glyphJust.numeric <- function(just, which=NULL, ...) {
    if (is.null(which)) {
        which <- names(just)
    }
    names(just) <- which
    class(just) <- "GlyphJust"
    just
}

################################################################################
## glyph font
glyphFont <- function(file, index,
                      family, weight, style,
                      PSname=NA) {
    file <- as.character(file)
    nafile <- is.na(file)
    if (any(nchar(file[!nafile], "bytes") > 500))
        warning("Font file longer than 500 will be truncated")
    index <- as.integer(index)
    family <- as.character(family)
    nafamily <- is.na(family)
    if (any(nchar(family[!nafamily], "bytes") > 200))
        warning("Font family longer than 200 will be truncated")
    weight <- mapWeight(weight)
    style <- mapStyle(style)
    PSname <- as.character(PSname)
    ## Missing PSname values are "estimated"
    naPS <- is.na(PSname)
    if (any(naPS)) {
        PSbold <- ifelse(weight >= 700, "Bold", "")
        PSstyle <- ifelse(style > 1,
                          ifelse(style > 2, "Oblique", "Italic"),
                          "")
        face <- paste0(PSbold, PSstyle)
        PSname[naPS] <-
            ifelse(nchar(file[naPS]),
                   sub("([^.]+)\\.[[:alnum:]]+$", "\\1",
                       basename(file[naPS])),
                   paste0(family[naPS],
                          ifelse(nchar(PSstyle), paste0("-", PSstyle), "")))
    }
    if (any(nchar(PSname, "bytes") > 200))
        warning("PostScript font name longer than 200 will be truncated")
    ## Check that family-weight-style and file and PSname all line up
    families <- rle(paste0(family, weight, style))$lengths
    files <- rle(file)$lengths
    names <- rle(PSname)$lengths
    if (!(all(families == files) && all(files == names)))
        stop("Font information is inconsistent")
    
    font <- list(file=file, index=index,
                 family=family, weight=weight, style=style,
                 PSname=PSname)
    class(font) <- "RGlyphFont"
    font
}

print.RGlyphFont <- function(x, ...) {
    cat(paste0(x$family, " wgt: ", x$weight, " style: ", invertStyle(x$style),
               "\n  (", x$file, " [", x$index, "])\n"))
}

glyphFontList <- function(...) {
    fonts <- list(...)
    if (!length(fonts))
        stop("List must include at least one font")
    if (!all(sapply(fonts, function(x) inherits(x, "RGlyphFont"))))
        stop("Invalid glyph font")
    class(fonts) <- "RGlyphFontList"
    fonts
}

################################################################################
## glyph information
glyphInfo <- function(id, x, y, font, size,
                      fontList,
                      width, height,
                      hAnchor, vAnchor,
                      col=NA) {
    id <- as.integer(id)
    x <- as.numeric(x)
    y <- as.numeric(y)
    ## Check font
    font <- as.integer(font)
    if (!inherits(fontList, "RGlyphFontList"))
        stop("Invalid font list")
    if (any(is.na(font)) || !all(font %in% seq_along(fontList)))
        stop("Unknown font")
    size <- as.numeric(size)
    ## Check colour (allow any R colour spec)
    nacol <- is.na(col)
    if (any(!nacol)) {
        rgb <- col2rgb(col[!nacol], alpha=TRUE)
        col[!nacol] <- rgb(rgb[1,], rgb[2,], rgb[3,], rgb[4,],
                           maxColorValue=255)
    }
    ## Check width/height
    if (!inherits(width, "GlyphWidth"))
        width <- glyphWidth(width)
    if (!inherits(height, "GlyphHeight"))
        height <- glyphHeight(height)
    ## Check anchors (and provide defaults if missing)
    if (missing(hAnchor))
        hAnchor <- glyphAnchor(c(min(x), min(x) + width[1],
                                 min(x) + width[1]/2),
                               label=c("left", "right", "centre"))
    if (missing(vAnchor))
        vAnchor <- glyphAnchor(c(min(y), min(y) + height[1],
                                 min(y) + height[1]/2),
                               label=c("bottom", "top", "centre"))
    if (!inherits(hAnchor, "GlyphAnchor"))
        hAnchor <- glyphAnchor(hAnchor, names(hAnchor))
    if (!inherits(vAnchor, "GlyphAnchor"))
        vAnchor <- glyphAnchor(vAnchor, names(vAnchor))
    hNames <- names(hAnchor)
    vNames <- names(vAnchor)
    if (!("left" %in% hNames && "bottom" %in% vNames))
        stop('There must be anchors named "left" and "bottom"')
    if (!"right" %in% hNames)
        hAnchor <- c(hAnchor,
                     right=unname(hAnchor["left"]) + unname(width[1]))
    if (!"top" %in% vNames)
        vAnchor <- c(vAnchor,
                     top=unname(vAnchor["bottom"]) + unname(height[1]))
    if (!"centre" %in% hNames)
        hAnchor <- c(hAnchor,
                     centre=unname(hAnchor["left"]) + unname(width[1]/2))
    if (!"centre" %in% vNames)
        vAnchor <- c(vAnchor,
                     centre=unname(vAnchor["bottom"]) + unname(height[1]/2))
    if (!"center" %in% hNames)
        hAnchor <- c(hAnchor, center=unname(hAnchor["centre"]))
    if (!"center" %in% vNames)
        vAnchor <- c(vAnchor, center=unname(vAnchor["centre"]))
    ## Build glyph info
    dropNA <- !(is.na(id) | is.na(x) | is.na(y) |
                ## is.na(font) already checked
                is.na(size))
    glyphs <- data.frame(id, x, y, font, size)[dropNA, ]
    if (nrow(glyphs) < 1)
        stop("Invalid glyph info")
    ## Colour can be NA
    if (inherits(glyphs, "omit")) {
        glyphs$colour <- col[-attr(glyphs, "na.action")]
    } else {
        glyphs$colour <- col
    }
    ## Construct final structure
    info <- list(glyphs=glyphs, fonts=fontList,
                 width=width, height=height,
                 hAnchor=hAnchor, vAnchor=vAnchor)
    class(info) <- c("RGlyphInfo")
    info
}



## Must match order in src/include/R_ext/GraphicsEngine.h
compositingOperators <- c(## Porter-Duff
                          "clear", "source", "over",
                          "in", "out", "atop",
                          "dest", "dest.over", "dest.in",
                          "dest.out", "dest.atop", "xor",
                          "add", "saturate",
                          ## PDF "blend modes"
                          "multiply", "screen", "overlay",
                          "darken", "lighten", "color.dodge", "color.burn",
                          "hard.light", "soft.light", "difference", "exclusion"
                          )

.opIndex <- function(x) {
    op <- match(x, compositingOperators)
    if (is.na(op))
        stop("Invalid compositing operator")
    as.integer(op)
}

.defineGroup <- function(source, op, destination) {
    .External(C_defineGroup, source, .opIndex(op), destination)
}

.useGroup <- function(ref, trans) {
    .External(C_useGroup, ref, trans)
}

.devUp <- function() {
    .External(C_devUp)
}
.hcl_colors_parameters <- as.data.frame(rbind(structure(numeric(0L), dim = c(0L, 16L), dimnames = list(NULL,
                      c("type", "h1", "h2", "h3", "c1", "c2", "c3", "l1", "l2", "l3", "p1", "p2", "p3", "p4", "cmax1", "cmax2"))),
    ##  1 : qualitative
    "Pastel 1"      = c(     1,    0,   NA,   NA,   35,   NA,   NA,   85,   NA,   NA,   NA,   NA,   NA,   NA,      NA,      NA),
    "Dark 2"        = c(     1,    0,   NA,   NA,   50,   NA,   NA,   60,   NA,   NA,   NA,   NA,   NA,   NA,      NA,      NA),
    "Dark 3"        = c(     1,    0,   NA,   NA,   80,   NA,   NA,   60,   NA,   NA,   NA,   NA,   NA,   NA,      NA,      NA),
    "Set 2"         = c(     1,    0,   NA,   NA,   60,   NA,   NA,   70,   NA,   NA,   NA,   NA,   NA,   NA,      NA,      NA),
    "Set 3"         = c(     1,   10,   NA,   NA,   50,   NA,   NA,   80,   NA,   NA,   NA,   NA,   NA,   NA,      NA,      NA),
    "Warm"          = c(     1,   90,  -30,   NA,   50,   NA,   NA,   70,   NA,   NA,   NA,   NA,   NA,   NA,      NA,      NA),
    "Cold"          = c(     1,  270,  150,   NA,   50,   NA,   NA,   70,   NA,   NA,   NA,   NA,   NA,   NA,      NA,      NA),
    "Harmonic"      = c(     1,   60,  240,   NA,   50,   NA,   NA,   70,   NA,   NA,   NA,   NA,   NA,   NA,      NA,      NA),
    "Dynamic"       = c(     1,   30,   NA,   NA,   50,   NA,   NA,   70,   NA,   NA,   NA,   NA,   NA,   NA,      NA,      NA),
    ##  2 : sequential
    "Grays"         = c(     2,    0,   NA,   NA,    0,   NA,   NA,   10,   98,   NA,  1.3,   NA,   NA,   NA,      NA,      NA),
    "Light Grays"   = c(     2,    0,   NA,   NA,    0,   NA,   NA,   30,   90,   NA,  1.5,   NA,   NA,   NA,      NA,      NA),
    "Blues 2"       = c(     2,  260,   NA,   NA,   80,   NA,   NA,   30,   90,   NA,  1.5,   NA,   NA,   NA,      NA,      NA),
    "Blues 3"       = c(     2,  245,   NA,   NA,   50,   NA,   NA,   20,   98,   NA,  0.8,  1.4,   NA,   NA,      75,      NA),
    "Purples 2"     = c(     2,  270,   NA,   NA,   70,   NA,   NA,   25,   95,   NA,  1.2,   NA,   NA,   NA,      NA,      NA),
    "Purples 3"     = c(     2,  270,   NA,   NA,   50,   NA,   NA,   20,   98,   NA,  0.9,  1.4,   NA,   NA,      75,      NA),
    "Reds 2"        = c(     2,   10,   NA,   NA,   85,   NA,   NA,   25,   95,   NA,  1.3,   NA,   NA,   NA,      NA,      NA),
    "Reds 3"        = c(     2,   10,   NA,   NA,   65,   NA,   NA,   20,   97,   NA,  1.1,  1.3,   NA,   NA,     150,      NA),
    "Greens 2"      = c(     2,  135,   NA,   NA,   45,   NA,   NA,   35,   95,   NA,  1.3,   NA,   NA,   NA,      NA,      NA),
    "Greens 3"      = c(     2,  135,   NA,   NA,   35,   NA,   NA,   25,   98,   NA,    1,  1.5,   NA,   NA,      70,      NA),
    "Oslo"          = c(     2,  250,   NA,   NA,    0,    0,   NA,   99,    1,   NA,    1,   NA,   NA,   NA,      70,      NA),
    "Purple-Blue"   = c(     2,  300,  200,   NA,   60,    0,   NA,   25,   95,   NA,  0.7,  1.3,   NA,   NA,      NA,      NA),
    "Red-Purple"    = c(     2,   10,  -80,   NA,   80,    5,   NA,   25,   95,   NA,  0.7,  1.3,   NA,   NA,      NA,      NA),
    "Red-Blue"      = c(     2,    0, -100,   NA,   80,   40,   NA,   40,   75,   NA,    1,    1,   NA,   NA,      NA,      NA),
    "Purple-Orange" = c(     2,  -83,   20,   NA,   65,   18,   NA,   32,   90,   NA,  0.5,    1,   NA,   NA,      NA,      NA),
    "Purple-Yellow" = c(     2,  320,   80,   NA,   60,   20,   NA,   30,   95,   NA,  0.7,  1.3,   NA,   NA,      65,      NA),
    "Blue-Yellow"   = c(     2,  265,   80,   NA,   60,   10,   NA,   25,   95,   NA,  0.7,    2,   NA,   NA,      NA,      NA),
    "Green-Yellow"  = c(     2,  140,   80,   NA,   50,   10,   NA,   40,   97,   NA,  0.7,  1.8,   NA,   NA,      NA,      NA),
    "Red-Yellow"    = c(     2,   10,   85,   NA,   80,   10,   NA,   25,   95,   NA,  0.4,  1.3,   NA,   NA,      NA,      NA),
    "Heat"          = c(     2,    0,   90,   NA,   80,   30,   NA,   30,   90,   NA,  0.2,    2,   NA,   NA,      NA,      NA),
    "Heat 2"        = c(     2,    0,   90,   NA,  100,   30,   NA,   50,   90,   NA,  0.2,    1,   NA,   NA,      NA,      NA),
    "Terrain"       = c(     2,  130,    0,   NA,   80,    0,   NA,   60,   95,   NA,  0.1,    1,   NA,   NA,      NA,      NA),
    "Terrain 2"     = c(     2,  130,   30,   NA,   65,    0,   NA,   45,   90,   NA,  0.5,  1.5,   NA,   NA,      NA,      NA),
    "Viridis"       = c(     2,  300,   75,   NA,   40,   95,   NA,   15,   90,   NA,    1,  1.1,   NA,   NA,      NA,      NA),
    "Plasma"        = c(     2, -100,  100,   NA,   60,  100,   NA,   15,   95,   NA,    2,  0.9,   NA,   NA,      NA,      NA),
    "Inferno"       = c(     2, -100,   85,   NA,    0,   65,   NA,    1,   98,   NA,  1.1,  0.9,   NA,   NA,     120,      NA),
    "Rocket"        = c(     2,  -70,   60,   NA,    0,   10,   NA,    2,   97,   NA,  0.8,  0.8,   NA,   NA,     130,      NA),
    "Mako"          = c(     2,  325,  130,   NA,    0,   18,   NA,    2,   95,   NA,  1.0,  1.0,   NA,   NA,      70,      NA),
    "Dark Mint"     = c(     2,  240,  130,   NA,   30,   33,   NA,   25,   95,   NA,    1,   NA,   NA,   NA,      NA,      NA),
    "Mint"          = c(     2,  205,  140,   NA,   40,   12,   NA,   34,   94,   NA,  0.5,    1,   NA,   NA,      NA,      NA),
    "BluGrn"        = c(     2,  215,  120,   NA,   25,   30,   NA,   31,   88,   NA,  0.7,  1.1,   NA,   NA,      45,      NA),
    "Teal"          = c(     2,  240,  180,   NA,   35,   15,   NA,   35,   92,   NA,  0.6,  1.1,   NA,   NA,      40,      NA),
    "TealGrn"       = c(     2,  220,  125,   NA,   44,   50,   NA,   49,   90,   NA,  0.8,  1.2,   NA,   NA,      60,      NA),
    "Emrld"         = c(     2,  224,  105,   NA,   23,   55,   NA,   25,   92,   NA,  1.5,    1,   NA,   NA,      NA,      NA),
    "BluYl"         = c(     2,  250,   90,   NA,   40,   55,   NA,   33,   98,   NA,  0.5,    1,   NA,   NA,      NA,      NA),
    "ag_GrnYl"      = c(     2,  225,   87,   NA,   27,   86,   NA,   34,   92,   NA,  0.9,   NA,   NA,   NA,      NA,      NA),
    "Peach"         = c(     2,   15,   50,   NA,  128,   30,   NA,   55,   90,   NA,  1.1,   NA,   NA,   NA,      NA,      NA),
    "PinkYl"        = c(     2,   -4,   80,   NA,  100,   47,   NA,   55,   96,   NA,    1,   NA,   NA,   NA,      NA,      NA),
    "Burg"          = c(     2,  -10,   10,   NA,   40,   40,   NA,   25,   85,   NA,  1.2,    1,   NA,   NA,      75,      NA),
    "BurgYl"        = c(     2,  -10,   55,   NA,   45,   30,   NA,   30,   90,   NA,  0.7,    1,   NA,   NA,      80,      NA),
    "RedOr"         = c(     2,   -3,   53,   NA,   75,   42,   NA,   44,   86,   NA,  0.8,    1,   NA,   NA,      90,      NA),
    "OrYel"         = c(     2,    5,   72,   NA,  120,   49,   NA,   56,   87,   NA,    1,   NA,   NA,   NA,     125,      NA),
    "Purp"          = c(     2,  270,  300,   NA,   55,   20,   NA,   42,   92,   NA,  0.6,    1,   NA,   NA,      60,      NA),
    "PurpOr"        = c(     2,  -83,   20,   NA,   55,   18,   NA,   32,   90,   NA,  0.6,    1,   NA,   NA,      65,      NA),
    "Sunset"        = c(     2,  -80,   78,   NA,   60,   55,   NA,   40,   91,   NA,  0.8,    1,   NA,   NA,      75,      NA),
    "Magenta"       = c(     2,  312,  358,   NA,   50,   24,   NA,   27,   85,   NA,  0.6,  1.1,   NA,   NA,      65,      NA),
    "SunsetDark"    = c(     2,  -35,   50,   NA,   55,   60,   NA,   30,   90,   NA,  1.2,    1,   NA,   NA,     120,      NA),
    "ag_Sunset"     = c(     2,  -85,   70,   NA,   70,   45,   NA,   25,   85,   NA,  0.6,    1,   NA,   NA,     105,      NA),
    "BrwnYl"        = c(     2,  -20,   70,   NA,   30,   20,   NA,   20,   90,   NA,    1,  1.1,   NA,   NA,      60,      NA),
    "YlOrRd"        = c(     2,    5,   85,   NA,   75,   40,   NA,   25,   99,   NA,  1.6,  1.3,   NA,   NA,     180,      NA),
    "YlOrBr"        = c(     2,   20,   85,   NA,   50,   20,   NA,   25,   99,   NA,  1.3,  1.5,   NA,   NA,     150,      NA),
    "OrRd"          = c(     2,    0,   60,   NA,   90,   10,   NA,   25,   97,   NA,    1,  1.5,   NA,   NA,     135,      NA),
    "Oranges"       = c(     2,   20,   55,   NA,   70,   10,   NA,   30,   97,   NA,  1.2,  1.3,   NA,   NA,     150,      NA),
    "YlGn"          = c(     2,  160,   85,   NA,   25,   20,   NA,   25,   99,   NA,  1.2,  1.6,   NA,   NA,      70,      NA),
    "YlGnBu"        = c(     2,  270,   90,   NA,   40,   25,   NA,   15,   99,   NA,    2,  1.5,   NA,   NA,      90,      NA),
    "Reds"          = c(     2,    0,   35,   NA,   65,    5,   NA,   20,   97,   NA,  1.1,  1.3,   NA,   NA,     150,      NA),
    "RdPu"          = c(     2,  -70,   40,   NA,   45,    5,   NA,   15,   97,   NA,    1,  1.3,   NA,   NA,     100,      NA),
    "PuRd"          = c(     2,   20,  -95,   NA,   60,    5,   NA,   20,   97,   NA,  1.6,  1.1,   NA,   NA,     140,      NA),
    "Purples"       = c(     2,  275,  270,   NA,   55,    5,   NA,   20,   99,   NA,  1.3,  1.3,   NA,   NA,      70,      NA),
    "PuBuGn"        = c(     2,  160,  320,   NA,   25,    5,   NA,   25,   98,   NA,  1.4,  1.2,   NA,   NA,      70,      NA),
    "PuBu"          = c(     2,  240,  260,   NA,   30,    5,   NA,   25,   98,   NA,  1.5,  1.2,   NA,   NA,      70,      NA),
    "Greens"        = c(     2,  135,  115,   NA,   35,    5,   NA,   25,   98,   NA,    1,  1.5,   NA,   NA,      70,      NA),
    "BuGn"          = c(     2,  125,  200,   NA,   30,    5,   NA,   25,   98,   NA,  1.4,  1.6,   NA,   NA,      65,      NA),
    "GnBu"          = c(     2,  265,   95,   NA,   55,   10,   NA,   25,   97,   NA,  1.3,  1.7,   NA,   NA,      80,      NA),
    "BuPu"          = c(     2,  320,  200,   NA,   40,    5,   NA,   15,   98,   NA,  1.2,  1.3,   NA,   NA,      65,      NA),
    "Blues"         = c(     2,  260,  220,   NA,   45,    5,   NA,   25,   98,   NA,  1.2,  1.3,   NA,   NA,      70,      NA),
    "Lajolla"       = c(     2,   90,  -20,   NA,   40,    5,   NA,   99,    5,   NA,  0.7,  0.8,   NA,   NA,     100,      NA),
    "Turku"         = c(     2,   10,  120,   NA,   20,    0,   NA,   95,    1,   NA,  1.7,  0.8,   NA,   NA,      55,      NA),
    "Hawaii"        = c(     2,  -30,  200,   NA,   70,   35,   NA,   30,   92,   NA,  0.3,    1,   NA,   NA,      75,      NA),
    "Batlow"        = c(     2,  270,  -40,   NA,   35,   35,   NA,   12,   88,   NA,  0.6,  1.1,   NA,   NA,      75,      NA),
    ##  3 : diverging
    "Blue-Red"      = c(     3,  260,    0,   NA,   80,   NA,   NA,   30,   90,   NA,  1.5,   NA,   NA,   NA,      NA,      NA),
    "Blue-Red 2"    = c(     3,  260,    0,   NA,  100,   NA,   NA,   50,   90,   NA,    1,   NA,   NA,   NA,      NA,      NA),
    "Blue-Red 3"    = c(     3,  255,   12,   NA,   50,   NA,   NA,   20,   97,   NA,    1,  1.3,   NA,   NA,      80,      NA),
    "Red-Green"     = c(     3,  340,  128,   NA,   60,   NA,   NA,   30,   97,   NA,  0.8,  1.5,   NA,   NA,      80,      NA),
    "Purple-Green"  = c(     3,  300,  128,   NA,   30,   NA,   NA,   20,   95,   NA,    1,  1.4,   NA,   NA,      65,      NA),
    "Purple-Brown"  = c(     3,  270,   40,   NA,   30,   NA,   NA,   20,   98,   NA,  0.8,  1.2,   NA,   NA,      70,      NA),
    "Green-Brown"   = c(     3,  180,   55,   NA,   40,   NA,   NA,   25,   97,   NA,  0.8,  1.4,   NA,   NA,      65,      NA),
    "Blue-Yellow 2" = c(     3,  265,   80,   NA,   80,   NA,   NA,   40,   95,   NA,  1.2,   NA,   NA,   NA,      NA,      NA),
    "Blue-Yellow 3" = c(     3,  265,   80,   NA,   80,   NA,   NA,   70,   95,   NA,  0.5,    2,   NA,   NA,      NA,      NA),
    "Green-Orange"  = c(     3,  130,   43,   NA,  100,   NA,   NA,   70,   90,   NA,    1,   NA,   NA,   NA,      NA,      NA),
    "Cyan-Magenta"  = c(     3,  180,  330,   NA,   59,   NA,   NA,   75,   95,   NA,  1.5,   NA,   NA,   NA,      NA,      NA),
    "Tropic"        = c(     3,  195,  325,   NA,   70,   NA,   NA,   55,   95,   NA,    1,   NA,   NA,   NA,      NA,      NA),
    "Broc"          = c(     3,  240,   85,   NA,   30,   NA,   NA,   15,   98,   NA,  0.9,   NA,   NA,   NA,      45,      NA),
    "Cork"          = c(     3,  245,  125,   NA,   30,   NA,   NA,   15,   95,   NA,  0.9,  1.1,   NA,   NA,      55,      NA),
    "Vik"           = c(     3,  240,   55,   NA,   45,   NA,   NA,   15,   95,   NA,  0.8,  1.1,   NA,   NA,      65,      NA),
    "Berlin"        = c(     3,  240,   15,   NA,   60,   NA,   NA,   75,    5,   NA,  1.2,  1.5,   NA,   NA,      80,      NA),
    "Lisbon"        = c(     3,  240,   85,   NA,   30,   NA,   NA,   98,    8,   NA,    1,   NA,   NA,   NA,      45,      NA),
    "Tofino"        = c(     3,  260,  120,   NA,   45,   NA,   NA,   90,    5,   NA,  0.8,    1,   NA,   NA,      55,      NA),
    ##  4 : divergingx
    "ArmyRose"      = c(     4,    0,   NA,   93,   73,   18,   47,   58,   98,   52,  1.5,  0.8,  0.8,    1,      NA,      55),
    "Earth"         = c(     4,   43,   82,  221,   61,   30,   45,   50,   92,   52,    1,    1,  0.8,    1,      NA,      10),
    "Fall"          = c(     4,  133,   77,   21,   20,   35,  100,   35,   95,   50,    1,   NA,  1.5,   NA,      NA,      NA),
    "Geyser"        = c(     4,  192,   77,   21,   40,   35,  100,   50,   95,   50,    1,    1,  1.2,    1,      20,      NA),
    "TealRose"      = c(     4,  190,   77,    0,   50,   25,   80,   55,   92,   55,  1.5,    1,  1.8,    1,      15,      NA),
    "Temps"         = c(     4,  191,   80,   -4,   43,   50,   78,   55,   89,   54,  1.6,    1,    1,    1,      57,      85),
    "PuOr"          = c(     4,   40,   NA,  270,   70,    0,   30,   30,   98,   10,  0.6,  1.4,  1.5,  1.3,     100,      65),
    "RdBu"          = c(     4,   20,   NA,  230,   60,    0,   50,   20,   98,   15,  1.4,  1.2,  1.5,  1.5,     125,      90),
    "RdGy"          = c(     4,    5,   50,   50,   60,    0,    0,   20,   98,   20,  1.2,  1.2,    1,  1.2,     125,      NA),
    "PiYG"          = c(     4,  340,   NA,  115,   75,    0,   50,   30,   98,   35,  1.3,  1.4,  0.8,  1.5,     100,      80),
    "PRGn"          = c(     4,  300,   NA,  128,   30,    0,   30,   15,   97,   25,  1.3,  1.2,  0.9,  1.5,      65,      65),
    "BrBG"          = c(     4,   55,   NA,  180,   40,    0,   30,   25,   97,   20,  0.8,  1.4,  0.8,  1.4,      75,      45),
    "RdYlBu"        = c(     4,   10,   85,  260,  105,   45,   70,   35,   98,   35,  1.5,  1.2,  0.6,  1.2,     150,      10),
    "RdYlGn"        = c(     4,   10,   85,  140,  105,   45,   50,   35,   98,   35,  1.5,  1.2,  0.8,  1.2,     150,      75),
    "Spectral"      = c(     4,    0,   85,  270,   90,   45,   65,   37,   98,   37,    1,  1.2,    1,  1.2,     120,      NA),
    "Zissou 1"      = c(     4,  218,   71,   12,   46,   88,  165,   59,   82,   52,  0.2,    1,    3,    1,      33,      NA),
    "Cividis"       = c(     4,  255,   NA,   75,   30,    0,   95,   13,   52,   92,  1.1,    1,    1,   NA,      47,      NA),
    "Roma"          = c(     4,   10,  120,  265,   80,   25,   80,   25,   92,   25,  0.4,  1.5,  1.0,  1.2,      NA,      NA)
))
.hcl_colors_parameters$type <- factor(.hcl_colors_parameters$type,
  labels = c("qualitative", "sequential", "diverging", "divergingx"))

hcl.pals <- function(type = NULL) {
    if (is.null(type)) return(rownames(.hcl_colors_parameters))
    type <- match.arg(tolower(type), levels(.hcl_colors_parameters$type))
    rownames(.hcl_colors_parameters)[.hcl_colors_parameters$type == type]
}

## palette function a la rainbow(n, ...), heat.colors(n) etc.
hcl.colors <- function(n, palette = "viridis",
                       alpha = NULL, rev = FALSE, fixup = TRUE)
{
    ## empty palette
    n <- as.integer(n[1L])
    if(n < 1L) return(character())

    ## expand alpha to length n (if any)
    if(!is.null(alpha)) alpha <- rep_len(alpha, n)

    ## match palette (ignoring case, space, -, _)
    fx <- function(x) tolower(gsub("[-, _, \\,, (, ), \\ , \\.]", "", x))
    p <- charmatch(fx(palette), fx(rownames(.hcl_colors_parameters)))
    if(is.na(p)) stop("'palette' does not match any given palette")
    if(p < 1L) stop("'palette' is ambiguous")
    p <- .hcl_colors_parameters[p, ]
    p.type <- as.character(p$type[1L]) # one of the factor levels
    p <- unlist(p[1L, -1L, drop = TRUE])

    ## trajectories
    lintrj <- function(i, p1, p2) p2 - (p2 - p1) * i
    tritrj <- function(i, j, p1, p2, pm) ifelse(i <= j,
        p2 - (p2 - pm) * i/j,
        pm - (pm - p1) * abs((i - j)/(1 - j)))
    seqhcl <- function(i, h1, h2, c1, c2, l1, l2, p1, p2, cmax, alpha) {
        j <- 1/(1 + abs(cmax - c1) / abs(cmax - c2))
        if (!is.na(j) && (j <= 0 || 1 <= j)) j <- NA
        hcl(h = lintrj(i, h1, h2),
            c = if(is.na(j)) lintrj(i^p1, c1, c2) else tritrj(i^p1, j, c1, c2, cmax),
	    l = lintrj(i^p2, l1, l2),
	    alpha = alpha,
	    fixup = fixup)
    }

    ## adapt defaults and set up HCL colors
    switch(p.type,
           "qualitative" = {
	if(is.na(p[["h2"]])) p[["h2"]] <- p[["h1"]] + 360 * (n - 1L)/n
        ## h/c/l trajectories
        i <- seq.int(1, 0, length.out = n)
	col <- hcl(h = lintrj(i, p[["h1"]], p[["h2"]]),
                   c = p[["c1"]], l = p[["l1"]],
		   alpha = alpha, fixup = fixup)
    }, "sequential" = {
	if(is.na(p[["h2"]])) p[["h2"]] <- p[["h1"]]
	if(is.na(p[["c2"]])) p[["c2"]] <- 0
	if(is.na(p[["p2"]])) p[["p2"]] <- p[["p1"]]
        ## h/c/l trajectories
        i <- seq.int(1, 0, length.out = n)
	col <- seqhcl(i, p[["h1"]], p[["h2"]], p[["c1"]], p[["c2"]],
			 p[["l1"]], p[["l2"]], p[["p1"]], p[["p2"]], p[["cmax1"]], alpha)
    }, "diverging" = {
	if(is.na(p[["p2"]])) p[["p2"]] <- p[["p1"]]
        ## h/c/l trajectories
        n2 <- ceiling(n/2)
	i <- if(n == 1L) 0 else seq.int(1, by = -2/(n - 1), length.out = n2)
	col <- c(seqhcl(i, p[["h1"]], p[["h1"]], p[["c1"]], 0, p[["l1"]], p[["l2"]], p[["p1"]], p[["p2"]], p[["cmax1"]], alpha[seq_len(n2)]),
	     rev(seqhcl(i, p[["h2"]], p[["h2"]], p[["c1"]], 0, p[["l1"]], p[["l2"]], p[["p1"]], p[["p2"]], p[["cmax1"]], alpha[rev(seq_len(n))[seq_len(n2)]])))
        if(n%/%2 < n2) col <- col[-n2]
    }, "divergingx" = {
	if(is.na(p[["p2"]])) p[["p2"]] <- p[["p1"]]
	if(is.na(p[["p4"]])) p[["p4"]] <- p[["p2"]]
        ## h/c/l trajectories
        n2 <- ceiling(n/2)
	i <- if(n == 1L) 0 else seq.int(1, by = -2/(n - 1), length.out = n2)
	col <- c(seqhcl(i, p[["h1"]], if(is.na(p[["h2"]])) p[["h1"]] else p[["h2"]], p[["c1"]], p[["c2"]],
			   p[["l1"]], p[["l2"]], p[["p1"]], p[["p2"]], p[["cmax1"]], alpha[seq_len(n2)]),
		 rev(seqhcl(i, p[["h3"]], if(is.na(p[["h2"]])) p[["h3"]] else p[["h2"]], p[["c3"]], p[["c2"]],
			       p[["l3"]], p[["l2"]], p[["p3"]], p[["p4"]], p[["cmax2"]], alpha[rev(seq_len(n))[seq_len(n2)]])))
        if(n%/%2 < n2) col <- col[-n2]
    },
    stop("wrong 'type'; should never happen, please report!"))

    if(rev) rev(col) else col
}

maskTypes <- c("alpha", "luminance")

.typeIndex <- function(x) {
    type <- match(x, maskTypes)
    if (is.na(type))
        stop("Invalid mask type")
    as.integer(type)
}

.mask <- function(fun, type) {
    attr(fun, "type") <- .typeIndex(type)
    fun
}

.setMask <- function(mask, ref) {    
    .External(C_setMask, mask, ref)
}

#  File src/library/grDevices/R/gradients.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 2019      The R Foundation
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

## Create R objects defining paths

fillRules <- c("winding", "evenodd")

.ruleIndex <- function(x) {
    rule <- match(x, fillRules)
    if (is.na(rule))
        stop("Invalid fill rule")
    as.integer(rule)
}


#  File src/library/grDevices/R/gradients.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 2019      The R Foundation
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

## Create R objects defining patterns

#######################################################
## MUST match C #defines in
## ../../../main/patterns.c and
## ../../../include/R_ext/GraphicsEngine.h
#######################################################

patternTypes <- c("LinearGradient", "RadialGradient", "TilingPattern")

extendOptions <- c("pad", "repeat", "reflect", "none")

pattern <- function(type, ...) {
    type <- match(type, patternTypes)
    if (is.na(type))
        stop("Invalid pattern type")
    
    pattern <- c(list(type=as.integer(type)), list(...))
    class(pattern) <- c(type, "Pattern")
    pattern
}

.linearGradientPattern <- function(colours = c("black", "white"),
                                   stops = seq(0, 1,
                                               length.out = length(colours)),
                                   x1 = 0, y1 = 0, x2 = 1, y2 = 1,
                                   extend = "pad") {
    ## Vectorising colours & stops
    nstops <- max(length(colours), length(stops))
    colours <- rep(colours, length.out = nstops)
    stops <- rep(stops, length.out = nstops)

    extend <- match(extend, extendOptions)
    if (is.na(extend))
        stop("Invalid 'extend' value")
    
    ## Make sure we really have colours in some form
    ## AND convert to "#RRGGBB" format
    RGBA <- col2rgb(colours, alpha=TRUE)
    colours <- rgb(t(RGBA[1:3,]), alpha=RGBA[4,], maxColorValue=255)
    
    pattern("LinearGradient",
            x1 = as.numeric(x1), y1 = as.numeric(y1),
            x2 = as.numeric(x2), y2 = as.numeric(y2),
            stops = as.numeric(stops), colours = colours,
            extend = as.integer(extend))
}

.radialGradientPattern <- function(colours = c("black", "white"),
                                   stops = seq(0, 1,
                                               length.out = length(colours)),
                                   cx1 = 0, cy1 = 0, r1 = 0,
                                   cx2 = 1, cy2 = 1, r2 = .5,
                                   extend = "pad") {
    ## Vectorising colours & stops
    nstops <- max(length(colours), length(stops))
    colours <- rep(colours, length.out = nstops)
    stops <- rep(stops, length.out = nstops)

    extend <- match(extend, extendOptions)
    if (is.na(extend))
        stop("Invalid 'extend' value")
    
    ## Make sure we really have colours in some form
    ## AND convert to "#RRGGBB" format
    RGBA <- col2rgb(colours, alpha=TRUE)
    colours <- rgb(t(RGBA[1:3,]), alpha=RGBA[4,], maxColorValue=255)
    
    pattern("RadialGradient",
            cx1 = as.numeric(cx1), cy1 = as.numeric(cy1),
            r1 = as.numeric(r1),
            cx2 = as.numeric(cx2), cy2 = as.numeric(cy2),
            r2 = as.numeric(r2),
            stops = as.numeric(stops), colours = colours,
            extend = as.integer(extend))
}

## (x, y) is (left, bottom)
.tilingPattern <- function(fun, x, y, width, height, extend) {
    extend <- match(extend, extendOptions)
    if (is.na(extend))
        stop("Invalid 'extend' value")
    
    pattern("TilingPattern",
            f = fun,
            x = as.numeric(x), y = as.numeric(y),
            width = as.numeric(width), height = as.numeric(height),
            extend = as.integer(extend))
}

.setPattern <- function(pattern) {
    .External(C_setPattern, pattern)
}

#  File src/library/grDevices/R/pictex.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
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

pictex <-
    function(file = "Rplots.tex", width = 5, height = 4, debug = FALSE,
	     bg = "white", fg = "black")
{
    .External(C_PicTeX, file, bg, fg, width, height, as.logical(debug))

    graphics::par(mar = c(5,4,2,4)+0.1)
}
#  File src/library/grDevices/R/postscript.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2023 The R Core Team
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

## An environment not exported from namespace:graphics used to pass
## .PostScript.Options and .PDF.options to the windows() device for
## use in its menus, and also to hide the variables.
.PSenv <- new.env()

check.options <-
    function(new, name.opt, reset = FALSE, assign.opt = FALSE,
	     envir = .GlobalEnv, check.attributes = c("mode", "length"),
	     override.check = FALSE)
{
    lnew <- length(new)
    if(lnew != length(newnames <- names(new)))
	stop(gettextf("invalid arguments in '%s' (need named args)",
                      deparse(sys.call(sys.parent()))), domain = NA)
    if(!is.character(name.opt))
	stop("'name.opt' must be character, name of an existing list")
    if(reset) {
	if(exists(name.opt, envir=envir, inherits=FALSE)) {
	    if(length(utils::find(name.opt)) > 1)
		rm(list=name.opt, envir=envir)

	} else stop(gettextf("cannot reset non-existent '%s'", name.opt),
                    domain = NA)
    }
    old <- get(name.opt, envir=envir, inherits=FALSE)
    if(!is.list(old))
	stop(gettextf("invalid options in '%s'", name.opt), domain = NA)
    oldnames <- names(old)
    if(lnew > 0) {
	matches <- pmatch(newnames, oldnames)
	if(any(is.na(matches)))
	    stop(sprintf(ngettext(as.integer(sum(is.na(matches))),
                                 "invalid argument name %s in '%s'",
                                 "invalid argument names %s in '%s'"),
                         paste(sQuote(newnames[is.na(matches)]),
                               collapse=", "),
                         deparse(sys.call(sys.parent()))),
                 domain = NA)
	else { #- match(es) found:  substitute if appropriate
	    i.match <- oldnames[matches]
	    prev <- old[i.match]
	    doubt <- rep.int(FALSE, length(prev))
	    for(fn in check.attributes)
		if(any(ii <- sapply(prev, fn) != sapply(new, fn))) {
                    ## skip 'fonts';
                    ii <- ii & (names(prev) != "fonts")
                    if(!any(ii)) next
		    doubt <- doubt | ii
		    do.keep <- ii & !override.check
		    warning(paste(sQuote(paste0(fn, "(", names(prev[ii]), ")" )),
                                  collapse = " and "), " ",
                            ngettext(as.integer(sum(ii)),
                                     "differs between new and previous",
                                     "differ between new and previous"),
                            if(any(do.keep)) {
                                paste0("\n\t ==> ",
                                       gettextf("NOT changing %s",
                                                paste(sQuote(names(prev[do.keep])),
                                                      collapse=" & ")))
                            } else "",
                            domain = NA, call. = FALSE)
		}
	    names(new) <- NULL
	    if(any(doubt)) {
		ii <- !doubt | override.check
		old[i.match[ii]] <- new[ii]
	    } else old[i.match] <- new

	}
	if(assign.opt) assign(name.opt, old, envir=envir)
    }
    old
}

ps.options <- function(..., reset = FALSE, override.check = FALSE)
{
    ## do initialization if needed
    initPSandPDFfonts()
    old <- get(".PostScript.Options", envir = .PSenv)
    if(reset) {
        assign(".PostScript.Options",
               get(".PostScript.Options.default", envir = .PSenv),
               envir = .PSenv)
    }
    new <- list(...)
    if(m <- match("append", names(new), 0L)) {
        warning("argument 'append' is for back-compatibility and will be ignored",
                immediate. = TRUE)
        new <- new[-m]
    }
    assign <- length(new) > 0
    check.options(new, name.opt = ".PostScript.Options", envir = .PSenv,
                  assign.opt = assign, override.check = override.check)
    if(reset || assign) invisible(old) else old
}

setEPS <- function(...)
{
    dots <- list(...)
    args <- list(width = 7, height = 7)
    args[names(dots)] <- dots
    force <- list(onefile = FALSE, horizontal = FALSE, paper = "special")
    args[names(force)] <- force
    do.call("ps.options", args)
}

setPS <- function(...)
{
    dots <- list(...)
    args <- list(width = 0, height = 0)
    args[names(dots)] <- dots
    force <- list(onefile = TRUE, horizontal = TRUE, paper = "default")
    args[names(force)] <- force
    do.call("ps.options", args)
}

pdf.options <- function(..., reset=FALSE)
{
    ## do initialization if needed
    initPSandPDFfonts()
    old <- get(".PDF.Options", envir = .PSenv)
    if(reset) {
        assign(".PDF.Options",
               get(".PDF.Options.default", envir = .PSenv),
               envir = .PSenv)
    }
    l... <- length(new <- list(...))
    check.options(new, name.opt = ".PDF.Options", envir = .PSenv,
                  assign.opt = l... > 0)
    if(reset || l... > 0) invisible(old) else old
}

guessEncoding <- function(family)
{
    # Three special families have special encodings, regardless of locale
    if (!missing(family) &&
        family %in% c("ComputerModern", "ComputerModernItalic")) {
        switch(family,
               "ComputerModern" = "TeXtext.enc",
               "ComputerModernItalic" = "TeXtext.enc")
    }  else {
        switch(.Platform$OS.type,
               "windows" = {
                   switch(utils::localeToCharset()[1L],
                          "ISO8859-2" = "CP1250.enc",
                          "ISO8859-7" = "CP1253.enc", # Greek
                          "ISO8859-13" = "CP1257.enc",
                          "CP1251" = "CP1251.enc", # Cyrillic
                          "WinAnsi.enc")
               },
               { lc <- utils::localeToCharset()
                 if(length(lc) == 1L)
                     switch(lc,
                            "ISO8859-1" = "ISOLatin1.enc",
                            "ISO8859-2" = "ISOLatin2.enc",
                            "ISO8859-5" = "Cyrillic.enc",
                            "ISO8859-7" = "Greek.enc",
                            "ISO8859-13" = "ISOLatin7.enc",
                            "ISO8859-15" = "ISOLatin9.enc",
                            "KOI8-R" = "KOI8-R.enc",
                            "KOI8-U" = "KOI8-U.enc",
                            "ISOLatin1.enc")
                 else if(lc[1L] == "UTF-8")
                     switch(lc[2L],
                            "ISO8859-1" = "ISOLatin1.enc", # what about Euro?
                            "ISO8859-2" = "ISOLatin2.enc",
                            "ISO8859-5" = "Cyrillic.enc",
                            "ISO8859-7" = "Greek.enc",
                            "ISO8859-13" = "ISOLatin7.enc",
                            "ISOLatin1.enc")
                 else "ISOLatin1.enc"})
    }
}

##--> source in devPS.c :

postscript <- function(file = if(onefile) "Rplots.ps" else "Rplot%03d.ps",
                       onefile, family, title , fonts, encoding, bg, fg,
                       width, height, horizontal, pointsize,
                       paper, pagecentre, print.it, command, colormodel,
                       useKerning, fillOddEven)
{
    ## do initialization if needed
    initPSandPDFfonts()

    new <- list()
    if(!missing(onefile)) new$onefile <- onefile
    ## 'family' is handled separately
    if(!missing(title)) new$title <- title
    if(!missing(fonts)) new$fonts <- fonts
    if(!missing(encoding)) new$encoding <- encoding
    if(!missing(bg)) new$bg <- bg
    if(!missing(fg)) new$fg <- fg
    if(!missing(width)) new$width <- width
    if(!missing(height)) new$height <- height
    if(!missing(horizontal)) new$horizontal <- horizontal
    if(!missing(pointsize)) new$pointsize <- pointsize
    if(!missing(paper)) new$paper <- paper
    if(!missing(pagecentre)) new$pagecentre <- pagecentre
    if(!missing(print.it)) new$print.it <- print.it
    if(!missing(command)) new$command <- command
    if(!missing(colormodel)) new$colormodel <- colormodel
    if(!missing(useKerning)) new$useKerning <- useKerning
    if(!missing(fillOddEven)) new$fillOddEven <- fillOddEven

    old <- check.options(new, name.opt = ".PostScript.Options", envir = .PSenv)

    if(is.null(old$command) || old$command == "default")
        old$command <- if(!is.null(cmd <- getOption("printcmd"))) cmd else ""

    # need to handle this case before encoding
    if(!missing(family) &&
       (inherits(family, "Type1Font") || inherits(family, "CIDFont"))) {
        enc <- family$encoding
        if(inherits(family, "Type1Font") &&!is.null(enc) && enc != "default"
           && (is.null(old$encoding) || old$encoding  == "default"))
            old$encoding <- enc
        family <- family$metrics
    }
    if(is.null(old$encoding) || old$encoding  == "default")
        old$encoding <- guessEncoding(family)

    ## handle family separately as length can be 1, 4, or 5
    if(!missing(family)) {
        # Case where family is a set of AFMs
        if(length(family) == 4L) {
            family <- c(family, "Symbol.afm")
        } else if (length(family) == 5L) {
            ## nothing to do
        } else if (length(family) == 1L) {
            ## If family has been specified, match with a font in the
            ## font database (see postscriptFonts())
            ## and pass in a device-independent font name.
            ## NOTE that in order to match, we need both family name
            ## and encoding to match.
            pf <- postscriptFonts(family)[[1L]]
            if(is.null(pf))
              stop(gettextf("unknown family '%s'", family), domain = NA)
            matchFont(pf, old$encoding)
        } else
            stop("invalid 'family' argument")
        old$family <- family
    }
    if(grepl("[\n\r\f\177]", old$title))
        ## title with these characters generates corrupt postscript file
        stop(gettextf("'title' argument \"%s\" contains invalid characters",
                      old$title), domain = NA)

    onefile <- old$onefile # for 'file'
    if(!checkIntFormat(file))
        stop(gettextf("invalid 'file' argument '%s'", file), domain = NA)
    if (!is.na(pmatch("MacRoman", old$encoding)))
        warning("encoding MacRoman is deprecated")
    .External(C_PostScript,
              file, old$paper, old$family, old$encoding, old$bg, old$fg,
              old$width, old$height, old$horizontal, old$pointsize,
              onefile, old$pagecentre, old$print.it, old$command,
              old$title, old$fonts, old$colormodel, old$useKerning,
              old$fillOddEven)
    # if .ps.prolog is searched for and fails, NULL got returned.
    invisible()
}

xfig <- function (file = if(onefile) "Rplots.fig" else "Rplot%03d.fig",
                  onefile = FALSE, encoding = "none",
                  paper = "default", horizontal = TRUE,
                  width = 0, height = 0, family = "Helvetica",
                  pointsize = 12, bg = "transparent", fg = "black",
                  pagecentre = TRUE,
                  defaultfont = FALSE, textspecial = FALSE)
{
    ## do initialization if needed
    initPSandPDFfonts()

    if(!checkIntFormat(file))
        stop(gettextf("invalid 'file' argument '%s'", file), domain = NA)
    .External(C_XFig, file, paper, family, bg, fg,
              width, height, horizontal, pointsize,
              onefile, pagecentre, defaultfont, textspecial, encoding)
    invisible()
}

pdf <- function(file = if(onefile) "Rplots.pdf" else "Rplot%03d.pdf",
                width, height, onefile, family, title, fonts, version,
                paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
                useDingbats, useKerning, fillOddEven, compress)
{
    ## do initialization if needed
    initPSandPDFfonts()

    new <- list()
    if(!missing(width)) new$width <- width
    if(!missing(height)) new$height <- height
    if(!missing(onefile)) new$onefile <- onefile
    ## 'family' is handled separately
    if(!missing(title)) new$title <- title
    if(!missing(fonts)) new$fonts <- fonts
    if(!missing(version)) new$version <- version
    if(!missing(paper)) new$paper <- paper
    if(!missing(encoding)) new$encoding <- encoding
    if(!missing(bg)) new$bg <- bg
    if(!missing(fg)) new$fg <- fg
    if(!missing(pointsize)) new$pointsize <- pointsize
    if(!missing(pagecentre)) new$pagecentre <- pagecentre
    if(!missing(colormodel)) new$colormodel <- colormodel
    if(!missing(useDingbats)) new$useDingbats <- useDingbats
    if(!missing(useKerning)) new$useKerning <- useKerning
    if(!missing(fillOddEven)) new$fillOddEven <- fillOddEven
    if(!missing(compress)) new$compress <- compress

    old <- check.options(new, name.opt = ".PDF.Options", envir = .PSenv)

    ## need to handle this before encoding
    if(!missing(family) &&
       (inherits(family, "Type1Font") || inherits(family, "CIDFont"))) {
        enc <- family$encoding
        if(inherits(family, "Type1Font") &&!is.null(enc) && enc != "default"
           && (is.null(old$encoding) || old$encoding  == "default"))
            old$encoding <- enc
        family <- family$metrics
    }
    if(is.null(old$encoding) || old$encoding  == "default")
        old$encoding <- guessEncoding()
    ## handle family separately as length can be 1, 4, or 5
    if(!missing(family)) {
        # Case where family is a set of AFMs
        if(length(family) == 4L) {
            family <- c(family, "Symbol.afm")
        } else if (length(family) == 5L) {
            ## nothing to do
        } else if (length(family) == 1L) {
            ## If family has been specified, match with a font in the
            ## font database (see postscriptFonts())
            ## and pass in a device-independent font name.
            ## NOTE that in order to match, we need both family name
            ## and encoding to match.
            pf <- pdfFonts(family)[[1L]]
            if(is.null(pf))
              stop(gettextf("unknown family '%s'", family), domain = NA)
            matchFont(pf, old$encoding)
        } else
            stop("invalid 'family' argument")
        old$family <- family
    }
    # Extract version
    version <- old$version
    versions <- c("1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "2.0")
    if (version %in% versions)
        version <- as.integer(strsplit(version, "[.]")[[1L]])
    else
        stop("invalid PDF version")

    onefile <- old$onefile # needed to set 'file'
    if(!checkIntFormat(file))
        stop(gettextf("invalid 'file' argument '%s'", file), domain = NA)
    if (!is.na(pmatch("MacRoman", old$encoding)))
        warning("encoding MacRoman is deprecated")
    .External(C_PDF,
              file, old$paper, old$family, old$encoding, old$bg, old$fg,
              old$width, old$height, old$pointsize, onefile, old$pagecentre,
              old$title, old$fonts, version[1L], version[2L],
              old$colormodel, old$useDingbats, old$useKerning,
              old$fillOddEven, old$compress)
    invisible()
}

.ps.prolog <- c(
"/gs  { gsave } bind def",
"/gr  { grestore } bind def",
"/ep  { showpage gr gr } bind def",
"/m   { moveto } bind def",
"/l  { rlineto } bind def",
"/np  { newpath } bind def",
"/cp  { closepath } bind def",
"/f   { fill } bind def",
"/o   { stroke } bind def",
"/c   { newpath 0 360 arc } bind def",
"/r   { 4 2 roll moveto 1 copy 3 -1 roll exch 0 exch rlineto 0 rlineto -1 mul 0 exch rlineto closepath } bind def",
"/p1  { stroke } bind def",
"/p2  { gsave bg fill grestore newpath } bind def",
"/p3  { gsave bg fill grestore stroke } bind def",
"/p6  { gsave bg eofill grestore newpath } bind def",
"/p7  { gsave bg eofill grestore stroke } bind def",
"/t   { 5 -2 roll moveto gsave rotate",
"       1 index stringwidth pop",
"       mul neg 0 rmoveto show grestore } bind def",
"/ta  { 4 -2 roll moveto gsave rotate show } bind def",
"/tb  { 2 -1 roll 0 rmoveto show } bind def",
"/cl  { grestore gsave newpath 3 index 3 index moveto 1 index",
"       4 -1 roll lineto  exch 1 index lineto lineto",
"       closepath clip newpath } bind def",
"/rgb { setrgbcolor } bind def",
"/s   { scalefont setfont } bind def")

.ps.prolog.srgb <- c(## From PLRM 3rd Ed pg 225
"/sRGB { [ /CIEBasedABC",
"          << /DecodeLMN",
"               [ { dup 0.03928 le",
"                        {12.92321 div}",
"                        {0.055 add 1.055 div 2.4 exp }",
"                     ifelse",
"                 } bind dup dup",
"               ]",
"             /MatrixLMN [0.412457 0.212673 0.019334",
"                         0.357576 0.715152 0.119192",
"                         0.180437 0.072175 0.950301]",
"             /WhitePoint [0.9505 1.0 1.0890]",
"           >>",
"         ] setcolorspace } bind def"
)

####################
# PostScript font database
#
# PostScript fonts may be either Type1 or CID-keyed fonts
# (the latter provides support for CJK fonts)
####################

assign(".PostScript.Fonts", list(), envir = .PSenv)

checkFont <- function(font) UseMethod("checkFont")

checkFont.default <- function(font) stop("Invalid font type")

# A Type1 font family has a name, plus a vector of 4 or 5 directories
# for font metric afm files, plus an encoding file

# Check that the font has the correct structure and information
# Already checked that it had a name
checkFont.Type1Font <- function(font) {
    if (is.null(font$family) || !is.character(font$family))
        stop("invalid family name in font specification")
    if (is.null(font$metrics) || !is.character(font$metrics) ||
        length(font$metrics) < 4L)
        stop("invalid metric information in font specification")
        ## Add default symbol font metric if none provided
    if (length(font$metrics) == 4L)
        font$metrics <- c(font$metrics, "Symbol.afm")
    if (is.null(font$encoding) || !is.character(font$encoding))
        stop("invalid encoding in font specification")
    font
}

# A CID-keyed font family has a name, four afm files,
# a CMap name, a CMap encoding, and (for now at least) a
# PDF chunk
# (I really hope we can dispense with the latter!)
checkFont.CIDFont <- function(font) {
    if (!inherits(font, "CIDFont"))
        stop("Not a CID font")
    if (is.null(font$family) || !is.character(font$family))
        stop("invalid family name in font specification")
    if (is.null(font$metrics) || !is.character(font$metrics) ||
        length(font$metrics) < 4L)
        stop("invalid metric information in font specification")
        ## Add default symbol font metric if none provided
    if (length(font$metrics) == 4L)
        font$metrics <- c(font$metrics, "Symbol.afm")
    if (is.null(font$cmap) || !is.character(font$cmap))
        stop("invalid CMap name in font specification")
    if (is.null(font$cmapEncoding) || !is.character(font$cmapEncoding))
        stop("invalid 'cmapEncoding' in font specification")
    if (is.null(font$pdfresource) || !is.character(font$pdfresource))
        stop("invalid PDF resource in font specification")
    font
}

isPDF <- function(fontDBname) {
  switch(fontDBname,
         .PostScript.Fonts=FALSE,
         .PDF.Fonts=TRUE,
         stop("Invalid font database name"))
}

checkFontInUse <- function(names, fontDBname) {
    for (i in names)
        if (.Call(C_Type1FontInUse, i, isPDF(fontDBname))
            || .Call(C_CIDFontInUse, i, isPDF(fontDBname)))
            stop(gettextf("font %s already in use", i), domain = NA)
    invisible()
}

setFonts <- function(fonts, fontNames, fontDBname) {
    fonts <- lapply(fonts, checkFont)
    fontDB <- get(fontDBname, envir=.PSenv)
    existingFonts <- fontNames %in% names(fontDB)
    if (sum(existingFonts) > 0) {
        checkFontInUse(fontNames[existingFonts], fontDBname)
        fontDB[fontNames[existingFonts]] <- fonts[existingFonts]
    }
    if (sum(existingFonts) < length(fontNames))
        fontDB <- c(fontDB, fonts[!existingFonts])
    assign(fontDBname, fontDB, envir=.PSenv)
}

printFont <- function(font) UseMethod("printFont")

printFont.Type1Font <- function(font)
    paste0(font$family, "\n    (", paste(font$metrics, collapse = " "),
           "\n    ", font$encoding, "\n")

printFont.CIDFont <- function(font)
    paste0(font$family, "\n    (", paste(font$metrics, collapse = " "),
           ")\n    ", font$CMap, "\n    ", font$encoding, "\n")

printFonts <- function(fonts)
    cat(paste0(names(fonts), ": ", unlist(lapply(fonts, printFont)),
               collapse = ""))

# If no arguments specified, return entire font database
# If no named arguments specified, all args should be font names
# on which to get info from the database
# Else, must specify new fonts to enter into database (all
# of which must be valid PostScript font descriptions and
# all of which must be named args)
postscriptFonts <- function(...)
{
    ## do initialization if needed: not recursive
    initPSandPDFfonts()
    ndots <- length(fonts <- list(...))
    if (ndots == 0L)
        get(".PostScript.Fonts", envir=.PSenv)
    else {
        fontNames <- names(fonts)
        nnames <- length(fontNames)
        if (nnames == 0L) {
            if (!all(sapply(fonts, is.character)))
                stop(gettextf("invalid arguments in '%s' (must be font names)",
                              "postscriptFonts"), domain = NA)
            else
                get(".PostScript.Fonts", envir=.PSenv)[unlist(fonts)]
        } else {
            if (ndots != nnames)
                stop(gettextf("invalid arguments in '%s' (need named args)",
                              "postscriptFonts"), domain = NA)
            setFonts(fonts, fontNames, ".PostScript.Fonts")
        }
    }
}

# Create a valid postscript font description
Type1Font <- function(family, metrics, encoding="default")
{
    font <- list(family=family, metrics=metrics, encoding=encoding)
    class(font) <- "Type1Font"
    checkFont(font)
}

CIDFont <- function(family, cmap, cmapEncoding, pdfresource="")
{
    font <- list(family=family, metrics=c("", "", "", ""), cmap=cmap,
                 cmapEncoding=cmapEncoding, pdfresource=pdfresource)
    class(font) <- "CIDFont"
    checkFont(font)
}


####################
# PDF font database
#
# PDF fonts may be either Type1 or CID-keyed fonts
# (the latter provides support for CJK fonts)
#
# PDF font database separate from PostScript one because
# some standard CID fonts are different
####################

assign(".PDF.Fonts", list(), envir = .PSenv)

pdfFonts <- function(...)
{
    ## do initialization if needed: not recursive
    initPSandPDFfonts()
    ndots <- length(fonts <- list(...))
    if (ndots == 0L)
        get(".PDF.Fonts", envir=.PSenv)
    else {
        fontNames <- names(fonts)
        nnames <- length(fontNames)
        if (nnames == 0L) {
            if (!all(sapply(fonts, is.character)))
                stop(gettextf("invalid arguments in '%s' (must be font names)",
                              "pdfFonts"), domain = NA)
            else
                get(".PDF.Fonts", envir=.PSenv)[unlist(fonts)]
        } else {
            if (ndots != nnames)
                stop(gettextf("invalid arguments in '%s' (need named args)",
                              "pdfFonts"), domain = NA)
            setFonts(fonts, fontNames, ".PDF.Fonts")
        }
    }
}

# Match an encoding
# NOTE that if encoding in font database is "default", that is a match
matchEncoding <- function(font, encoding) UseMethod("matchEncoding")

matchEncoding.Type1Font <- function(font, encoding) {
    ## the trailing .enc is optional
    font$encoding %in% c("default", encoding, paste0(encoding, ".enc"))
}

# Users should not be specifying a CID font AND an encoding
# when starting a new device
matchEncoding.CIDFont <- function(font, encoding) TRUE

# Match a font name (and an encoding)
matchFont <- function(font, encoding) {
    if (is.null(font))
        stop("unknown font")
    if (!matchEncoding(font, encoding))
        stop(gettextf("font encoding mismatch '%s'/'%s'",
                      font$encoding, encoding), domain=NA)
}

# Function to initialise default PostScript and PDF fonts
# Called at first use
#   a) because that's a sensible place to do initialisation of package globals
#   b) because it does not work to do it BEFORE then.  In particular,
#      if the body of this function is evaluated when the R code of the
#      package is sourced, then the method dispatch on checkFont() does
#      not work because when the R code is sourced, the S3 methods in
#      this package have not yet been registered.
#      Also, we want the run-time locale not the install-time locale.

initPSandPDFfonts <- function() {
    if(exists(".PostScript.Options", envir = .PSenv, inherits=FALSE)) return()

assign(".PostScript.Options",
    list(onefile = TRUE,
         family = "Helvetica",
         title = "R Graphics Output",
         fonts = NULL,
	 encoding = "default",
	 bg	= "transparent",
	 fg	= "black",
	 width	= 0,
	 height = 0,
         horizontal = TRUE,
	 pointsize  = 12,
         paper	= "default",
         pagecentre = TRUE,
	 print.it   = FALSE,
	 command    = "default",
         colormodel = "srgb",
         useKerning = TRUE,
         fillOddEven = FALSE), envir = .PSenv)
assign(".PostScript.Options.default",
       get(".PostScript.Options", envir = .PSenv),
       envir = .PSenv)

assign(".PDF.Options",
    list(width	= 7,
	 height = 7,
         onefile = TRUE,
         family = "Helvetica",
         title = "R Graphics Output",
         fonts = NULL,
         version = "1.4",
         paper = "special",
         encoding = "default",
	 bg	= "transparent",
	 fg	= "black",
	 pointsize  = 12,
	 pagecentre = TRUE,
         colormodel = "srgb",
         useDingbats = FALSE,
         useKerning = TRUE,
         fillOddEven = FALSE,
         compress = TRUE), envir = .PSenv)
assign(".PDF.Options.default",
       get(".PDF.Options", envir = .PSenv),
       envir = .PSenv)


postscriptFonts(# Default Serif font is Times
                serif = Type1Font("Times",
                  c("Times-Roman.afm", "Times-Bold.afm",
                    "Times-Italic.afm", "Times-BoldItalic.afm",
                    "Symbol.afm")),
                # Default Sans Serif font is Helvetica
                sans = Type1Font("Helvetica",
                  c("Helvetica.afm", "Helvetica-Bold.afm",
                    "Helvetica-Oblique.afm", "Helvetica-BoldOblique.afm",
                    "Symbol.afm")),
                # Default Monospace font is Courier
                mono = Type1Font("Courier",
                  c("Courier.afm", "Courier-Bold.afm",
                    "Courier-Oblique.afm", "Courier-BoldOblique.afm",
                    "Symbol.afm")),
                # Remainder are standard Adobe fonts that
                # should be present on PostScript devices
                AvantGarde = Type1Font("AvantGarde",
                  c("agw_____.afm", "agd_____.afm",
                    "agwo____.afm", "agdo____.afm",
                    "Symbol.afm")),
                Bookman = Type1Font("Bookman",
                  c("bkl_____.afm", "bkd_____.afm",
                    "bkli____.afm", "bkdi____.afm",
                    "Symbol.afm")),
                Courier = Type1Font("Courier",
                  c("Courier.afm", "Courier-Bold.afm",
                    "Courier-Oblique.afm", "Courier-BoldOblique.afm",
                    "Symbol.afm")),
                Helvetica = Type1Font("Helvetica",
                  c("Helvetica.afm", "Helvetica-Bold.afm",
                    "Helvetica-Oblique.afm", "Helvetica-BoldOblique.afm",
                    "Symbol.afm")),
                "Helvetica-Narrow" = Type1Font("Helvetica-Narrow",
                  c("hvn_____.afm", "hvnb____.afm",
                    "hvno____.afm", "hvnbo___.afm",
                    "Symbol.afm")),
                NewCenturySchoolbook = Type1Font("NewCenturySchoolbook",
                  c("ncr_____.afm", "ncb_____.afm",
                    "nci_____.afm", "ncbi____.afm",
                    "Symbol.afm")),
                Palatino = Type1Font("Palatino",
                  c("por_____.afm", "pob_____.afm",
                    "poi_____.afm", "pobi____.afm",
                    "Symbol.afm")),
                Times = Type1Font("Times",
                  c("Times-Roman.afm", "Times-Bold.afm",
                    "Times-Italic.afm", "Times-BoldItalic.afm",
                    "Symbol.afm")),
                # URW equivalents
                URWGothic = Type1Font("URWGothic",
                  c("a010013l.afm", "a010015l.afm",
                    "a010033l.afm", "a010035l.afm",
                    "s050000l.afm")),
                URWBookman = Type1Font("URWBookman",
                  c("b018012l.afm", "b018015l.afm",
                    "b018032l.afm", "b018035l.afm",
                    "s050000l.afm")),
                NimbusMon = Type1Font("NimbusMon",
                  c("n022003l.afm", "n022004l.afm",
                    "n022023l.afm", "n022024l.afm",
                    "s050000l.afm")),
                NimbusSan = Type1Font("NimbusSan",
                  c("n019003l.afm", "n019004l.afm",
                    "n019023l.afm", "n019024l.afm",
                    "s050000l.afm")),
                URWHelvetica = Type1Font("URWHelvetica",
                  c("n019003l.afm", "n019004l.afm",
                    "n019023l.afm", "n019024l.afm",
                    "s050000l.afm")),
                NimbusSanCond = Type1Font("NimbusSanCond",
                  c("n019043l.afm", "n019044l.afm",
                    "n019063l.afm", "n019064l.afm",
                    "s050000l.afm")),
                CenturySch = Type1Font("CenturySch",
                  c("c059013l.afm", "c059016l.afm",
                    "c059033l.afm", "c059036l.afm",
                    "s050000l.afm")),
                URWPalladio = Type1Font("URWPalladio",
                  c("p052003l.afm", "p052004l.afm",
                    "p052023l.afm", "p052024l.afm",
                    "s050000l.afm")),
                NimbusRom = Type1Font("NimbusRom",
                  c("n021003l.afm", "n021004l.afm",
                    "n021023l.afm", "n021024l.afm",
                    "s050000l.afm")),
                URWTimes = Type1Font("URWTimes",
                  c("n021003l.afm", "n021004l.afm",
                    "n021023l.afm", "n021024l.afm",
                    "s050000l.afm")),
                ## And Monotype Arial
                ArialMT = Type1Font("ArialMT",
                  c("ArialMT.afm", "ArialMT-Bold.afm",
                    "ArialMT-Italic.afm", "ArialMT-BoldItalic.afm",
                    "Symbol.afm"))
                )

## All of the above Type1 fonts are the same for PostScript and PDF
do.call("pdfFonts", postscriptFonts())

## add ComputerModern to postscript only
postscriptFonts(# Computer Modern as recoded by Brian D'Urso
                ComputerModern = Type1Font("ComputerModern",
                  c("CM_regular_10.afm", "CM_boldx_10.afm",
                    "CM_italic_10.afm", "CM_boldx_italic_10.afm",
                    "CM_symbol_10.afm"), encoding = "TeXtext.enc"),
                 ComputerModernItalic = Type1Font("ComputerModernItalic",
                  c("CM_regular_10.afm", "CM_boldx_10.afm", "cmti10.afm",
                    "cmbxti10.afm", "CM_symbol_10.afm"),
                 encoding = "TeXtext.enc")
                )


# CJK fonts
postscriptFonts(Japan1 = CIDFont("HeiseiKakuGo-W5", "EUC-H", "EUC-JP"),
                Japan1HeiMin = CIDFont("HeiseiMin-W3", "EUC-H", "EUC-JP"),
                Japan1GothicBBB =
                CIDFont("GothicBBB-Medium", "EUC-H", "EUC-JP"),
                Japan1Ryumin = CIDFont("Ryumin-Light", "EUC-H", "EUC-JP"),
                Korea1 = CIDFont("Baekmuk-Batang", "KSCms-UHC-H", "CP949"),
                Korea1deb = CIDFont("Batang-Regular", "KSCms-UHC-H", "CP949"),
                CNS1 = CIDFont("MOESung-Regular", "B5pc-H", "CP950"),
                GB1 = CIDFont("BousungEG-Light-GB", "GBK-EUC-H", "GBK"))

pdfFonts(Japan1 = CIDFont("KozMinPro-Regular-Acro", "EUC-H", "EUC-JP",
           paste("/FontDescriptor",
                 "<<",
                 "  /Type /FontDescriptor",
                 "  /CapHeight 740 /Ascent 1075 /Descent -272 /StemV 72",
                 "  /FontBBox [-195 -272 1110 1075]",
                 "  /ItalicAngle 0 /Flags 6 /XHeight 502",
                 "  /Style << /Panose <000001000500000000000000> >>",
                 ">>",
                 "/CIDSystemInfo << /Registry(Adobe) /Ordering(Japan1) /Supplement  2 >>",
                 "/DW 1000",
                 "/W [",
                 "   1 632 500 ",
                 "   8718 [500 500] ",
                 "]\n",
                 sep = "\n      ")),
         Japan1HeiMin = CIDFont("HeiseiMin-W3-Acro", "EUC-H", "EUC-JP",
           paste("/FontDescriptor",
                 "<<",
                 "  /Type /FontDescriptor",
                 "  /CapHeight 709 /Ascent 723 /Descent -241 /StemV 69",
                 "  /FontBBox [-123 -257 1001 910]",
                 "  /ItalicAngle 0 /Flags 6 /XHeight 450",
                 "  /Style << /Panose <000002020500000000000000> >>",
                 ">>",
                 "/CIDSystemInfo << /Registry(Adobe) /Ordering(Japan1) /Supplement  2 >>",
                 "/DW 1000",
                 "/W [",
                 "   1 632 500 ",
                 "   8718 [500 500] ",
                 "]\n",
                 sep = "\n      ")),
         Japan1GothicBBB = CIDFont("GothicBBB-Medium", "EUC-H", "EUC-JP",
           paste("/FontDescriptor",
                 "<<",
                 "  /Type /FontDescriptor",
                 "  /CapHeight 737 /Ascent 752 /Descent -271 /StemV 99",
                 "  /FontBBox [-22 -252 1000 892]",
                 "  /ItalicAngle 0 /Flags 4",
                 "  /Style << /Panose <0801020b0500000000000000> >>",
                 ">>",
                 "/CIDSystemInfo << /Registry(Adobe) /Ordering(Japan1) /Supplement  2 >>",
                 "/DW 1000",
                 "/W [",
                 "   1 632 500",
                 "   8718 [500 500]",
                 "]\n",
                 sep = "\n      ")),
         Japan1Ryumin = CIDFont("Ryumin-Light", "EUC-H", "EUC-JP",
           paste("/FontDescriptor",
                 "<<",
                 "  /Type /FontDescriptor",
                 "  /CapHeight 709 /Ascent 723 /Descent -241 /StemV 69",
                 "  /FontBBox [-54 -305 1000 903]",
                 "  /ItalicAngle 0 /Flags 6",
                 "  /Style << /Panose <010502020300000000000000> >>",
                 ">>",
                 "/CIDSystemInfo << /Registry(Adobe) /Ordering(Japan1) /Supplement  2 >>",
                 "/DW 1000",
                 "/W [",
                 "   1 632 500",
                 "   8718 [500 500]",
                 "]\n",
                 sep = "\n      ")),
         Korea1 = CIDFont("HYSMyeongJoStd-Medium-Acro", "KSCms-UHC-H", "CP949",
           paste("/FontDescriptor",
                 "<<",
                 "  /Type /FontDescriptor",
                 "  /CapHeight 720 /Ascent 880 /Descent -148 /StemV 59",
                 "  /FontBBox [-28 -148 1001 880]",
                 "  /ItalicAngle 0 /Flags 6 /XHeight 468",
                 "  /Style << /Panose <000001000600000000000000> >>",
                 ">>",
                 "/CIDSystemInfo << /Registry(Adobe) /Ordering(Korea1) /Supplement 1 >>",
                 "/DW 1000",
                 "/W [",
                 "   1 94 500",
                 "   97 [500] ",
                 "   8094 8190 500",
                 "]\n",
                 sep = "\n      ")),
         Korea1deb = CIDFont("HYGothic-Medium-Acro", "KSCms-UHC-H", "CP949",
           paste("/FontDescriptor",
                 "<<",
                 "  /Type /FontDescriptor",
                 "  /CapHeight 737 /Ascent 752 /Descent -271 /StemV 58",
                 "  /FontBBox [-6 -145 1003 880]",
                 "  /ItalicAngle 0 /Flags 4 /XHeight 553",
                 "  /Style << /Panose <000001000600000000000000> >>",
                 ">>",
                 "/CIDSystemInfo << /Registry(Adobe) /Ordering(Korea1) /Supplement 1 >>",
                 "/DW 1000",
                 "/W [",
                 "   1 94 500",
                 "   97 [500] ",
                 "   8094 8190 500",
                 "]\n",
                 sep = "\n      ")),
         CNS1 = CIDFont("MSungStd-Light-Acro", "B5pc-H", "CP950",
           paste("/FontDescriptor",
                 "<<",
                 "  /Type /FontDescriptor",
                 "  /CapHeight 662 /Ascent 1071 /Descent -249 /StemV 66",
                 "  /FontBBox [-160 -249 1015 1071]",
                 "  /ItalicAngle 0 /Flags 6 /XHeight 400",
                 "  /Style << /Panose <000001000600000000000000> >>",
                 ">>",
                 "/CIDSystemInfo << /Registry(Adobe) /Ordering(CNS1) /Supplement  0 >>",
                 "/DW 1000",
                 "/W [",
                 "     1 33 500",
                 "     34 [749 673 679 679 685 671 738 736 333 494 739 696 902 720 750 674 746 672 627 769 707 777 887 709 616]",
                 "     60 65 500",
                 "     66 [500 511 502 549 494 356 516 550 321 321 510 317 738 533 535 545 533 376 443 261 529 742 534 576 439]",
                 "     92 95 500",
                 "     13648 13742 500",
                 "     17603 [500]",
                 "]\n",
                 sep = "\n      ")),
         GB1 = CIDFont("STSong-Light-Acro", "GBK-EUC-H", "GBK",
           paste("/FontDescriptor",
                 "<<",
                 "  /Type /FontDescriptor",
                 "  /CapHeight 626 /Ascent 905 /Descent -254 /StemV 48",
                 "  /FontBBox [-134 -254 1001 905]",
                 "  /ItalicAngle 0 /Flags 6 /XHeight 416",
                 "  /Style << /Panose <000000000400000000000000> >>",
                 ">>",
                 "/CIDSystemInfo << /Registry(Adobe) /Ordering(GB1) /Supplement  2 >>",
                 "/DW 1000",
                 "/W [",
                 "     1 95 500",
                 "     814 939 500",
                 "     7712 7716 500",
                 "     22127 22357 500",
                 "]\n",
                 sep = "\n      ")))
}

# Call ghostscript to process postscript or pdf file to embed fonts
# (could also be used to convert ps or pdf to any supported format)
embedFonts <- function(file, # The ps or pdf file to convert
                       format, # Default guessed from file suffix
                       outfile = file, # By default overwrite file
                       fontpaths = character(),
                       options = character() # Additional options to ghostscript
                       )
{
    if(!is.character(file) || length(file) != 1L || !nzchar(file))
        stop("'file' must be a non-empty character string")
    gsexe <- tools::find_gs_cmd()
    if(!nzchar(gsexe)) stop("GhostScript was not found")
    if(.Platform$OS.type == "windows") gsexe <- shortPathName(gsexe)
    suffix <- gsub(".+[.]", "", file)
    if (missing(format))
        format <- switch(suffix,
                         ps = , eps = "ps2write",
                         pdf = "pdfwrite")
    if (!is.character(format)) stop("invalid output format")
    check_gs_type(gsexe, format)
    tmpfile <- tempfile("Rembed")
    if (length(fontpaths))
        fontpaths <-
            paste0("-sFONTPATH=",
                   shQuote(paste(fontpaths, collapse = .Platform$path.sep)))
    args <- c(paste0("-dNOPAUSE -dBATCH -q -dAutoRotatePages=/None -sDEVICE=", format),
              paste0(" -sOutputFile=", shQuote(tmpfile)),
              fontpaths, options, shQuote(file))
    ret <- system2(gsexe, args)
    if(ret != 0)
        stop(gettextf("status %d in running command '%s'", ret, cmd),
             domain = NA)
    if(outfile != file) args[2] <- paste0(" -sOutputFile=", shQuote(outfile))
    cmd <- paste(c(shQuote(gsexe), args), collapse = " ")
    file.copy(tmpfile, outfile, overwrite = TRUE)
    invisible(cmd)
}

## 'file' is the pdf file to convert
## 'glyphInfo' is RGlyphInfo
## 'outfile' is the new pdf file
## 'options' are additional options to ghostscript
embedGlyphs <- function(file, glyphInfo, outfile = file,
                        options = character()) {
    if (!is.character(file) || length(file) != 1L || !nzchar(file))
        stop("'file' must be a non-empty character string")
    if (!inherits(glyphInfo, "RGlyphInfo"))
        stop("Invalid 'glyphInfo'")
    gsexe <- tools::find_gs_cmd()
    if(!nzchar(gsexe)) stop("GhostScript was not found")
    if(.Platform$OS.type == "windows") gsexe <- shortPathName(gsexe)
    format <- "pdfwrite"
    check_gs_type(gsexe, format)
    tmpfile <- tempfile("Rembed")
    ## Generate cidfmap to relate font names to font files
    cidfmap <- file.path(tempdir(), "cidfmap")
    fontfile <- unique(sapply(glyphInfo$fonts, function(x) x$file))
    fontname <- unique(sapply(glyphInfo$fonts, function(x) x$PSname))
    writeLines(paste0("/", fontname,
                      " << /FileType /TrueType /Path (", fontfile,
                      ") /SubfontID 0 /CSI [(Identity) 0] >>;"),
               cidfmap)
    args <- c(paste0("-dNOPAUSE -dBATCH -q -dAutoRotatePages=/None ",
                     "-sDEVICE=", format),
              paste0("-sOutputFile=", shQuote(tmpfile)),
              ## Make sure ghostscript can see the cidfmap
              paste0("-I", shQuote(tempdir())),
              options, shQuote(file))
    ret <- system2(gsexe, args)
    if (ret != 0)
        stop(gettextf("status %d in running command '%s'", ret, cmd),
             domain = NA)
    if (outfile != file)
        args[2] <- paste0(" -sOutputFile=", shQuote(outfile))
    cmd <- paste(c(shQuote(gsexe), args), collapse = " ")
    file.copy(tmpfile, outfile, overwrite = TRUE)
    invisible(cmd)
}
#  File src/library/grDevices/R/prettyDate.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2023 The R Core Team
#
# Original code Copyright (C) 2010 Felix Andrews
# Modifications Copyright (C) 2010 The R Core Team
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

##' S3 method =:  pretty.Date() and pretty.POSIXt [in ../NAMESPACE]
prettyDate <- function(x, n = 5, min.n = n %/% 2, sep = " ", ...)
{
    stopifnot(min.n <= n)
    isDate <- inherits(x, "Date")
    x <- as.POSIXct(x)
    if (isDate) # the timezone *does* matter
	attr(x, "tzone") <- "GMT"
    zz <- rx <- range(x, na.rm = TRUE)
    D <- diff(nzz <- as.numeric(zz))
    MIN <- 60
    HOUR <- MIN * 60
    DAY <- HOUR * 24
    YEAR <- DAY * 365.25
    MONTH <- YEAR / 12
    makeOutput <- function(at, s, round = TRUE, do) {
	structure(if(isDate)
		      if(round) as.Date(round(at, units = "days")) else at
		  else as.POSIXct(at),
		  labels = format(at, s$format),
          format = s$format)
    }
    if(isDate && D <= n * DAY) { # D <= 'n days' & Date  ==> use days
	zz <- as.Date(zz)
	r <- round(n - D/DAY)
	m <- max(0, r %/% 2)
        m2 <- m + (r %% 2)
	while(length(dd <- seq.Date(zz[1] - m, zz[2] + m2, by = "1 day")) < min.n + 1)
	    if(m < m2) m <- m+1 else m2 <- m2+1
	return(makeOutput(dd, round = FALSE, ## "1 DSTday" from steps:
			  list(format = paste("%b", "%d", sep = sep))))
    }
    else if(D < 1) { # unique values / sub-second ranges: [? or use "1 ms" steps below?]
	m <- min(30, max(D == 0, n/2))
	zz <- structure(c(floor(nzz[1] - m), ceiling(nzz[2] + m)),
			class = class(x), tzone = attr(x, "tzone"))
    }
    xspan <- as.numeric(diff(zz), units = "secs")
    ## specify the set of pretty timesteps
    steps <-
        list("1 sec" = list(1, format = "%S", start = "mins"),
             "2 secs" = list(2),
             "5 secs" = list(5),
             "10 secs" = list(10),
             "15 secs" = list(15),
             "30 secs" = list(30, format = "%H:%M:%S"),
             "1 min" = list(1*MIN, format = "%H:%M"),
             "2 mins" = list(2*MIN, start = "hours"),
             "5 mins" = list(5*MIN),
             "10 mins" = list(10*MIN),
             "15 mins" = list(15*MIN),
             "30 mins" = list(30*MIN),
             ## "1 hour" = list(1*HOUR),
	     "1 hour" = list(1*HOUR, format = if (xspan <= DAY) "%H:%M"
					      else paste("%b %d", "%H:%M", sep = sep)),
             "3 hours" = list(3*HOUR, start = "days"),
             "6 hours" = list(6*HOUR, format = paste("%b %d", "%H:%M", sep = sep)),
             "12 hours" = list(12*HOUR),
             "1 DSTday" = list(1*DAY, format = paste("%b", "%d", sep = sep)),
             "2 DSTdays" = list(2*DAY),
             "1 week" = list(7*DAY, start = "weeks"),
             "halfmonth" = list(MONTH/2, start = "months"),
             ## "1 month" = list(1*MONTH, format = "%b"),
	     "1 month" = list(1*MONTH, format = if (xspan < YEAR) "%b"
						else paste("%b", "%Y", sep = sep)),
             "3 months" = list(3*MONTH, start = "years"),
             "6 months" = list(6*MONTH, format = "%Y-%m"),
             "1 year" = list(1*YEAR, format = "%Y"),
             "2 years" = list(2*YEAR, start = "decades"),
             "5 years" = list(5*YEAR),
             "10 years" = list(10*YEAR),
             "20 years" = list(20*YEAR, start = "centuries"),
             "50 years" = list(50*YEAR),
             "100 years" = list(100*YEAR),
             "200 years" = list(200*YEAR),
             "500 years" = list(500*YEAR),
             "1000 years" = list(1000*YEAR))
    ## carry forward 'format' and 'start' to following steps
    for (i in seq_along(steps)) {
        if (is.null(steps[[i]]$format))
            steps[[i]]$format <- steps[[i-1]]$format
        if (is.null(steps[[i]]$start))
            steps[[i]]$start <- steps[[i-1]]$start
        steps[[i]]$spec <- names(steps)[i]
    }
    ## crudely work out number of steps in the given interval
    nsteps <- xspan / vapply(steps, `[[`, numeric(1), 1L, USE.NAMES=FALSE)
    init.i <- init.i0 <- which.min(abs(nsteps - n))
    ## calculate actual number of ticks in the given interval
    calcSteps <- function(s, lim = range(zz)) {
        startTime <- trunc_POSIXt(lim[1], units = s$start) ## FIXME: should be trunc() eventually
        at <- seqDtime(startTime, end = lim[2], by = s$spec)
	if(anyNA(at)) { at <- at[!is.na(at)]; if(!length(at)) return(at) }
	r1 <- sum(at <= lim[1])
	r2 <- length(at) + 1 - sum(at >= lim[2])
	if(r2 == length(at) + 1) { # not covering at right -- add point at right
	    nat <- seqDtime(at[length(at)], by = s$spec, length=2)[2]
	    if(is.na(nat) || !(nat > at[length(at)])) # failed
		r2 <- length(at)
	    else
		at[r2] <- nat
	}
	## Now we could see if we are *smaller* than 'n+1' and add even more at[] on both sides
	at[r1:r2]
    }
    init.at <- calcSteps(st.i <- steps[[init.i]])
    ## bump it up if below acceptable threshold
    R <- TRUE # R := TRUE iff "right"
    L.fail <- R.fail <- FALSE
    while ((init.n <- length(init.at) - 1L) < min.n) {
	if(init.i == 1L) { ## keep steps[[1]]
	    ## add new interval right or left
            if(R) {
                nat <- seqDtime(init.at[length(init.at)], by = st.i$spec, length=2)[2]
                R.fail <- is.na(nat) || !(nat > init.at[length(init.at)])
                if(!R.fail)
                    init.at[length(init.at) + 1] <- nat
            } else { # left
                nat <- seqDtime(init.at[1], by = paste0("-",st.i$spec), length=2)[2]
                L.fail <- is.na(nat) || !(nat < init.at[1])
                if(!L.fail) {
                    init.at[seq_along(init.at) + 1] <- init.at
                    init.at[1] <- nat
                }
            }
            if(R.fail && L.fail)
                stop("failed to add more ticks; 'min.n' too large?")
	    R <- !R # alternating right <-> left
	} else { # smaller step sizes
	    init.i <- init.i - 1L
	    init.at <- calcSteps(st.i <- steps[[init.i]])
	}
    }
    if (init.n == n) ## perfect
        return(makeOutput(init.at, st.i))
    ## else : have a difference dn :
    dn <- init.n - n
    if(dn > 0L) {  ## too many ticks
	## ticks "outside", on left and right, keep at least one on each side
	nl <- sum(init.at <= rx[1]) - 1L
	nr <- sum(init.at >= rx[2]) - 1L
	if(nl > 0L || nr > 0L) {
	    n.c <- nl+nr # number of removable ticks
	    if(dn < n.c) { # remove dn, not all
		nl <- round(dn * nl/n.c)
		nr <- dn - nl
	    }
	    ## remove nl on left,  nr on right:
	    init.at <- init.at[-c(seq_len(nl), length(init.at)+1L-seq_len(nr))]
	}
    } else { ## too few ticks
        ## warning("trying to add more points -- not yet implemented")
        ## but after all, 'n' is approximate
	## init.at <- calcSteps(st.i, "more ticks")
    }
    if ((dn <- length(init.at) - 1L - n) == 0L  ## perfect
	|| (dn > 0L && init.i < init.i0) # too many, but we tried init.i + 1 already
        || (dn < 0L && init.i == 1)) # too few, but init.i = 1
	return(makeOutput(init.at, st.i))

    new.i <- if (dn > 0L) ## too many ticks
		 min(init.i + 1L, length(steps))
	     else ## too few ticks (and init.i > 1):
		 init.i - 1L
    new.at <- calcSteps(steps[[new.i]])
    new.n <- length(new.at) - 1L
    ## work out whether new.at or init.at is better
    if (new.n < min.n)
        new.n <- -Inf
    if (abs(new.n - n) < abs(dn))
	makeOutput(new.at, steps[[new.i]])
    else
	makeOutput(init.at, st.i)
}


## Utility, a generalization/special case of seq.POSIXct() / seq.Date()
seqDtime <- function(beg, end, by, length=NULL) {
    if(missing(by) || !identical(by, "halfmonth"))
        return( seq(beg, end, by = by, length.out=length) )
    ## else  by == "halfmonth" => can only go forward (!)
    if(is.null(length)) {
        l2 <- NULL; i <- TRUE
    } else {
        l2 <- ceiling(length/2); i <- seq_len(length)
    }
    at <- seq(beg, end, by = "months", length.out = l2)
    at2 <- as.POSIXlt(at)
    stopifnot(length(md <- unique(at2$mday)) == 1)
    at <- as.POSIXct(at)
    ## intersperse at and at2 := 15-day-shifted( at ), via rbind():
    if(md == 1) {
        at2$mday <- 15L
    } else if(md >= 15) { # (md == 16 may happen; not seen yet)
        at2$mday <- 1L
        at2$mon <- at2$mon + 1L
        ## at2 now has wrong 'yday','wday',.. and we rely on as.POSIXct():
    } else if(md < 15) { ## e.g., southern hemisphere, seen 14
        at2$mday <- md + 14L # consistent w (1 -> 15) in 1st case; ok even in Feb.
    }
    at2$isdst <- -1L
    at2 <- rbind(at, as.POSIXct(at2), deparse.level = 0L)
    structure(at2[i], class = class(at), tzone = attr(at, "tzone"))
}


## utility function, extending the base function trunc.POSIXt.
## Ideally this should replace the original, but that should be done
## with a little more thought (what about round.POSIXt etc.?)

trunc_POSIXt <-
    function(x, units = c("secs", "mins", "hours", "days",
                "weeks", "months", "years", "decades", "centuries"),
             start.on.monday = TRUE)
{
    x <- as.POSIXlt(x)
    if (units %in% c("secs", "mins", "hours", "days"))
	return(trunc.POSIXt(x, units))
    x <- trunc.POSIXt(x, "days")
    if (length(x$sec))
        switch(units,
               weeks = {
                   x$mday <- x$mday - x$wday
                   if (start.on.monday)
                       x$mday <- x$mday + ifelse(x$wday > 0L, 1L, -6L)
               },
               months = {
                   x$mday <- 1
               },
               years = {
                   x$mday <- 1
                   x$mon <- 0
               },
               decades = {
                   x$mday <- 1
                   x$mon <- 0
                   x$year <- (x$year %/% 10) * 10
               },
               centuries = {
                   x$mday <- 1
                   x$mon <- 0
                   x$year <- (x$year %/% 100) * 100
               })
    x
}
#  File src/library/grDevices/R/raster.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2015 The R Core Team
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


# A raster object is a character vector
# of colour strings
# plus a number of rows and columns
# The vector gives colours in ROW ORDER,
# starting from the TOP row
#
# due to the inherent inefficiency of
# raster implementation the graphics
# routines also support "nativeRaster"
# which is the native R representation
# (integer matrix) of colors in the same
# order as raster, suitable for practical
# use with images

is.raster <- function(x)
    inherits(x, "raster")

as.raster <- function(x, ...)
    UseMethod("as.raster")

as.raster.raster <- function(x, ...)  x

as.raster.logical <- function(x, max = 1, ...)
    as.raster(matrix(x, ...), max)

as.raster.raw <- function(x, max = 255L, ...)
    as.raster(matrix(x, ...), max=max)

as.raster.numeric <- as.raster.logical

as.raster.character <- as.raster.logical


as.raster.matrix <- function(x, max = 1, ...)
{
    if (is.character(x)) {
        ## Assume to be color names
        r <- t(x)
    } else if (is.numeric(x) || is.logical(x)) {
        ## Assume greyscale or b&w values
        ## We have to use rgb() indirectly as it
        ## doesn't hande NAs correctly
        tx <- t(x)
        tx.na <- which(is.na(tx))
        if (length(tx.na)) tx[tx.na] <- 0
        r <- rgb(tx, tx, tx, maxColorValue = max)
        if (length(tx.na)) r[tx.na] <- NA
    } else if (is.raw(x)) { ## non NA's here
	storage.mode(x) <- "integer"
	tx <- t(x)
	r <- rgb(tx, tx, tx, maxColorValue = 255L)
    } else
        stop("a raster matrix must be character, or numeric, or logical")
    ## Transposed, but retain original dimensions
    dim(r) <- dim(x)
    class(r) <- "raster"
    r
}

as.raster.array <- function(x, max = 1, ...)
{
    if (!is.numeric(x)) {
	if (is.raw(x)) {
	    storage.mode(x) <- "integer" # memory x 4 (!)
	    max <- 255L
	} else
	    stop("a raster array must be numeric")
    }
    if (length(d <- dim(x)) != 3L)
        stop("a raster array must have exactly 3 dimensions")
    r <- array(if (d[3L] == 3L)
        rgb(t(x[,,1L]), t(x[,,2L]), t(x[,,3L]), maxColorValue = max)
    else if (d[3L] == 4L)
        rgb(t(x[,,1L]), t(x[,,2L]), t(x[,,3L]), t(x[,,4L]), maxColorValue = max)
    else
        stop("a raster array must have exactly 3 or 4 planes"),
    dim = d[1:2])
    class(r) <- "raster"
    r
}

# Conversion to (character) matrix
as.matrix.raster <- function(x, ...)
{
    dim <- dim(x)
    matrix(x, nrow = dim[1L], ncol = dim[2L], byrow = TRUE)
}

is.na.raster <- function(x) is.na(as.matrix(x))
anyNA.raster <- function(x, recursive = FALSE) anyNA(as.matrix(x))

# FIXME:
# It would be useful to have conversion to array (rgb[a])
# so that people could play with numeric machinations
# with raster images

print.raster <- function(x, ...) print(as.matrix(x), ...)


# Subsetting methods
# Non-standard because raster is ROW-wise
# Try to piggy-back on existing methods as much as possible
# IGNORE 'drop' -- i.e. use "drop = FALSE" -- in all cases, but  m[i]
`[.raster` <- function(x, i, j, drop, ...)
{
    mdrop <- missing(drop)
    nA <- nargs() - (!mdrop)
    if(!mdrop && !isFALSE(drop))
        warning("'drop' is always implicitly FALSE in '[.raster'")
    m <- as.matrix(x)
    m <-
	if (missing(i)) {
	    if(missing(j)) m[ , drop = FALSE] else m[, j, drop = FALSE]
	} else if (missing(j)) {
	    if (nA == 2) ## is.matrix(i) || is.logical(i))
		return(m[i]) # behave as a matrix and directly return character vector
	    else if(nA == 3) m[i, , drop = FALSE]
	    else stop("invalid raster subsetting")
	} else m[i, j, drop = FALSE]
    as.raster(m)
}

`[<-.raster` <- function(x, i, j, value)
{
    nA <- nargs()
    m <- as.matrix(x)
    if (missing(i)) {
	if(missing(j)) m[] <- value else m[, j] <- value
    } else if (missing(j)) {
	if (nA == 3) ## typically is.matrix(i) || is.logical(i))
	    m[i] <- value
	else if(nA == 4) m[i, ] <- value
	else stop("invalid raster subassignment")
    } else m[i, j] <- value
    as.raster(m)
}

Ops.raster <- function(e1, e2)
{
    if (.Generic %in% c("==", "!=")) {
        ## Allows comparison of rasters with each other or with colour names
        if (is.raster(e1)) e1 <- as.matrix(e1)
        if (is.raster(e2)) e2 <- as.matrix(e2)
        ## The result is a logical matrix
        get(.Generic)(e1, e2)
    } else {
        stop("operator not meaningful for raster objects")
    }
}

#  File src/library/grDevices/R/recordplot.R
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

# 'load' and 'attach' should be character vectors of package names
recordPlot <- function(load=NULL, attach=NULL)
{
    if(dev.cur() == 1L)
        stop("no current device to record from")
    res <- .External2(C_getSnapshot)
    attr(res, "pid") <- Sys.getpid()
    attr(res, "Rversion") <- getRversion()
    attr(res, "load") <- as.character(load)
    attr(res, "attach") <- as.character(attach)
    class(res) <- "recordedplot"
    res
}

replayPlot <- function(x, reloadPkgs=FALSE)
{
    if(!inherits(x, "recordedplot"))
        stop(gettextf("argument is not of class %s", dQuote("recordedplot")),
             domain = NA)
    pid <- attr(x, "pid") ## added after R 3.0.2
    if (doRestore <- (is.null(pid) || pid != Sys.getpid())) {
        # This is a "recordedplot" loaded from another session
        x <- restoreRecordedPlot(x, reloadPkgs)
    }
    restore <- function(e) {
        if(doRestore)
            stop("invalid \"recordedplot\": ", conditionMessage(e))
        ## else: typically deserialized recordedplot from this session:
        .External2(C_playSnapshot,
                   restoreRecordedPlot(x, reloadPkgs))
    }
    r <- tryCatch(.External2(C_playSnapshot, x), error = restore)
    invisible(r)
}

print.recordedplot <- function(x, ...)
{
    replayPlot(x)
    invisible(x)
}

# If this is a "recordedplot" that has been saved and reloaded
# (possibly across sessions) then we need to ...
# - warn if have R version mismatch
# - restore NativeSymbolInfo on each element of the snapshot display list
# - bail out gracefully if something is not right
restoreRecordedPlot <- function(plot, reloadPkgs) {
    snapshotRversion <- attr(plot, "Rversion")
    if (is.null(snapshotRversion)) {
        warning("snapshot recorded in different R version (pre 3.3.0)")
    } else if (snapshotRversion != getRversion()) {
        warning(gettextf("snapshot recorded in different R version (%s)",
                         snapshotRversion))
    }
    # Ensure that all graphics systems in the snapshot are available
    # (snapshots only started recording pkgName in R 3.3.0)
    # Similar for any 'pkgs' saved with the snapshot
    if (!is.null(snapshotRversion) &&
        snapshotRversion >= R_system_version("3.3.0")) {
        for (i in seq_along(plot)[-1]) { # " 2:length(plot) "
            library(attr(plot[[i]], "pkgName"), character.only=TRUE)
        }
        if (reloadPkgs) {
            load <- attr(plot, "load")
            for (i in load) {
                loadNamespace(i)
            }
            attach <- attr(plot, "attach")
            for (i in attach) {
                library(i, character.only=TRUE)
            }
        }
    }
    # The display list is the first component of the snapshot
    for (i in seq_along(plot[[1]])) {
        # get the symbol then test if it's a native symbol
        symbol <- plot[[c(1L, i, 2L, 1L)]]
        if (inherits(symbol, "NativeSymbolInfo")) {
            # determine the dll that the symbol lives in
            name <- (if(!is.null(symbol$package))
                         symbol$package else symbol$dll)[["name"]]
            pkgDLL <- getLoadedDLLs()[[name]]
            # reconstruct the native symbol and assign it into the plot
            # This will error out if it fails to find the symbol, which
            # is some protection against running "recordedplot" in
            # R session where the recorded function does not exist!
            nativeSymbol <- getNativeSymbolInfo(name = symbol$name,
                                                PACKAGE = pkgDLL,
                                                withRegistrationInfo = TRUE)
            # Check that the 'numParameters' matches.
            # If it does not, we would also receive a redundant WARNING
            # about R version or graphics engine version mismatch,
            # but this mismatch is serious enough to put a STOP to things.
            if (nativeSymbol$numParameters != symbol$numParameters) {
                stop("snapshot contains invalid graphics call")
            }
            plot[[c(1L, i, 2L, 1L)]] <- nativeSymbol
        }
    }
    plot
}
#  File src/library/grDevices/R/smooth2d.R
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


## need some standard blues to plot ; output of brewer.pal(9, "Blues"):
blues9 <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
	    "#4292C6", "#2171B5", "#08519C", "#08306B")

.smoothScatterCalcDensity <- function(x, nbin, bandwidth, range.x)
{
    if (length(nbin) == 1)
	nbin <- c(nbin, nbin)
    if (!is.numeric(nbin) || length(nbin) != 2)
	stop("'nbin' must be numeric of length 1 or 2")

    if (missing(bandwidth)) { ## cheap
	bandwidth <- diff(apply(x, 2, stats::quantile,
				probs = c(0.05, 0.95),
                                na.rm = TRUE, names = FALSE)) / 25
	bandwidth[bandwidth==0] <- 1
    }
    else {
	if(!is.numeric(bandwidth)) stop("'bandwidth' must be numeric")
	if(any(bandwidth <= 0)) stop("'bandwidth' must be positive")
    }
    if(requireNamespace("KernSmooth", quietly=TRUE))
    ## create density map
    rv <- KernSmooth::bkde2D(x, bandwidth=bandwidth, gridsize=nbin,
			     range.x=range.x)
    else stop("Must have the ('Recommended') package \"KernSmooth\" installed")
    rv$bandwidth <- bandwidth
    rv
}

densCols <- function(x, y = NULL, nbin = 128, bandwidth,
		     colramp = colorRampPalette(blues9[-(1:3)]))
{
    ## similar as in plot.default
    xy <- xy.coords(x, y, setLab = FALSE)

    ## deal with NA, etc
    select <- is.finite(xy$x) & is.finite(xy$y)
    x <- cbind(xy$x, xy$y)[select, ]

    ## create density map
    map <- .smoothScatterCalcDensity(x, nbin, bandwidth)

    ## bin  x- and y- values
    mkBreaks <- function(u) u - diff(range(u))/(length(u)-1)/2
    xbin <- cut(x[,1], mkBreaks(map$x1), labels = FALSE)
    ybin <- cut(x[,2], mkBreaks(map$x2), labels = FALSE)

    dens <- map$fhat[cbind(xbin, ybin)]
    dens[is.na(dens)] <- 0

    ## transform densities to colors
    colpal <- cut(dens, length(dens), labels = FALSE)
    cols   <- rep(NA_character_, length(select))
    cols[select] <- colramp(length(dens))[colpal]

    cols
}
#  File src/library/grDevices/R/unix/dev2bitmap.R
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

dev2bitmap <- function(file, type = "png16m", height = 7, width = 7, res = 72,
                       units = "in", pointsize, ...,
                       method = c("postscript", "pdf"), taa = NA, gaa = NA)
{
    if(missing(file)) stop("'file' is missing with no default")
    if(!is.character(file) || length(file) != 1L || !nzchar(file))
        stop("'file' must be a non-empty character string")
    method <- match.arg(method)
    units <- match.arg(units, c("in", "px", "cm", "mm"))
    height <- switch(units, "in"=1, "cm"=1/2.54, "mm"=1/25.4, "px"=1/res) * height
    width <- switch(units, "in"=1, "cm"=1/2.54, "mm"=1/25.4, "px"=1/res) * width
    gsexe <- tools::find_gs_cmd()
    if(!nzchar(gsexe)) stop("GhostScript was not found")
    check_gs_type(gsexe, type)
    if(missing(pointsize)) pointsize <- 1.5*min(width, height)
    tmp <- tempfile("Rbit")
    on.exit(unlink(tmp))
    din <- graphics::par("din"); w <- din[1L]; h <- din[2L]
    if(missing(width) && !missing(height)) width <- w/h * height
    if(missing(height) && !missing(width)) height <- h/w * width

    current.device <- dev.cur()
    if(method == "pdf")
        dev.off(dev.copy(device = pdf, file = tmp, width = width,
                         height = height,
                         pointsize = pointsize, paper = "special", ...))
    else
        dev.off(dev.copy(device = postscript, file = tmp, width = width,
                         height = height,
                         pointsize = pointsize, paper = "special",
                         horizontal = FALSE, ...))
    dev.set(current.device)
    extra <- ""
    if (!is.na(taa)) extra <- paste0(" -dTextAlphaBits=", taa)
    if (!is.na(gaa)) extra <- paste0(extra, " -dGraphicsAlphaBits=", gaa)
    cmd <- paste0(shQuote(gsexe), " -dNOPAUSE -dBATCH -q -sDEVICE=", type,
                  " -r", res,
                  " -dAutoRotatePages=/None",
                  " -g", ceiling(res*width), "x", ceiling(res*height),
                  extra,
                  " -sOutputFile=", shQuote(file), " ", tmp)
    system(cmd)
    invisible()
}

bitmap <- function(file, type = "png16m", height = 7, width = 7, res = 72,
                   units = "in", pointsize, taa = NA, gaa = NA, ...)
{
    if(missing(file)) stop("'file' is missing with no default")
    if(!is.character(file) || length(file) != 1L || !nzchar(file))
        stop("'file' must be a non-empty character string")
    units <- match.arg(units, c("in", "px", "cm", "mm"))
    height <- switch(units, "in"=1, "cm"=1/2.54, "mm"=1/25.4, "px"=1/res) * height
    width <- switch(units, "in"=1, "cm"=1/2.54, "mm"=1/25.4, "px"=1/res) * width
    gsexe <- tools::find_gs_cmd()
    if(!nzchar(gsexe)) stop("GhostScript was not found")
    check_gs_type(gsexe, type)
    if(missing(pointsize)) pointsize <- 1.5*min(width, height)
    extra <- ""
    if (!is.na(taa)) extra <- paste0(" -dTextAlphaBits=", taa)
    if (!is.na(gaa)) extra <- paste0(extra, " -dGraphicsAlphaBits=", gaa)
    cmd <- paste0("|", shQuote(gsexe),
                  " -dNOPAUSE -dBATCH -q -sDEVICE=", type,
                  " -r", res,
                  " -dAutoRotatePages=/None",
                  " -g", ceiling(res*width), "x", ceiling(res*height),
                  extra,
                  " -sOutputFile=", shQuote(file), " -")
    postscript(file = cmd, width = width, height = height,
               pointsize = pointsize, paper = "special", horizontal = FALSE, ...)
    invisible()
}


## unexported
check_gs_type <- function(gsexe, type)
{
    gshelp <- system(paste(gsexe, "-help"), intern = TRUE)
    st <- grep("^Available", gshelp)
    en <- grep("^Search", gshelp)
    if(!length(st) || !length(en))
        warning("unrecognized format of gs -help")
    else {
        gsdevs <- gshelp[(st+1L):(en-1L)]
        devs <- c(strsplit(gsdevs, " "), recursive = TRUE)
        if(match(type, devs, 0L) == 0L) {
            op <- options(warning.length = 8000L)
            on.exit(options(op))
            stop(gettextf("device '%s' is not available\n", type),
                 gettextf("Available devices are:\n%s",
                          paste(gsdevs, collapse = "\n")),
                 domain = NA)
        }
    }
}
#  File src/library/grDevices/R/unix/png.R
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

.geometry <- function(width, height, units, res)
{
    width <- as.numeric(width)
    if (length(width) == 0 || !is.finite(width[1]))
        stop("invalid width")
    else width <- width[1]
    height <- as.numeric(height)
    if (length(height) == 0 || !is.finite(height[1]))
        stop("invalid height")
    else height <- height[1]

    units <- match.arg(units, c("in", "px", "cm", "mm"))
    if(units != "px" && is.na(res))
        stop("'res' must be specified unless 'units = \"px\"'")
    width <- switch(units,
                    "in" = res,
                    "cm" = res/2.54,
                    "mm" = res/25.4,
                    "px" = 1) * width
    height <- switch(units,
                     "in" = res,
                     "cm" = res/2.54,
                     "mm" = res/25.4,
                     "px" = 1) * height
    list(width = width, height = height)
}

png <- function(filename = "Rplot%03d.png",
                width = 480, height = 480, units = "px",
                pointsize = 12, bg = "white", res = NA, ...,
                type = c("cairo", "cairo-png", "Xlib", "quartz"), antialias)
{
    if(!is.character(filename) || length(filename) != 1L || !nzchar(filename))
        stop("'filename' must be a non-empty character string")
    if(!checkIntFormat(filename)) stop("invalid 'filename'")
    g <- .geometry(width, height, units, res)
    new <- list(...)
    if(missing(type)) type <- getOption("bitmapType")
    type <- match.arg(type)
    if(!missing(antialias)) new$antialias <- match.arg(antialias, aa.cairo)
    d <- check.options(new, name.opt = ".X11.Options", envir = .X11env)
    antialias <- match(d$antialias, aa.cairo)
    if(type == "quartz") {
        if(capabilities("aqua")) {
            width <- g$width/ifelse(is.na(res), 72, res);
            height <- g$height/ifelse(is.na(res), 72, res);
            invisible(.External(C_Quartz, "png", path.expand(filename),
                                width, height, pointsize, d$family,
                                d$antialias != "none", "", bg,
                                "white", if(is.na(res)) NULL else res))
        } else warning('type = "quartz" is unavailable. trying "Xlib"')
    } else if (type == "cairo") {
        if(capabilities("cairo")) {
        invisible(.External(C_devCairo, filename, 2L, g$width, g$height,
                            pointsize, bg, res, antialias, 100L, d$family, 300,
                            optionSymbolFont(d$symbolfamily)))
        } else warning('type = "cairo" is unavailable. trying "Xlib"')
    } else if (type == "cairo-png") {
        if(capabilities("cairo")) {
            invisible(.External(C_devCairo, filename, 5L, g$width, g$height,
                                pointsize, bg, res, antialias, 100L, d$family,
                                300, optionSymbolFont(d$symbolfamily)))
        } else warning('type = "cairo-png" is unavailable. trying "Xlib"')
    } else
        invisible(.External2(C_X11,
                             paste0("png::", filename),
                             g$width, g$height, pointsize, d$gamma,
                             d$colortype, d$maxcubesize, bg, bg, d$fonts, res,
                             0L, 0L, "", 0, 0, d$family,
                             optionSymbolFont(d$symbolfamily)))
}

jpeg <- function(filename = "Rplot%03d.jpeg",
                 width = 480, height = 480, units = "px",
                 pointsize = 12, quality = 75,
                 bg = "white", res = NA, ...,
                 type = c("cairo", "Xlib", "quartz"), antialias)
{
    if(!is.character(filename) || length(filename) != 1L || !nzchar(filename))
        stop("'filename' must be a non-empty character string")
    if(!checkIntFormat(filename)) stop("invalid 'filename'")
    g <- .geometry(width, height, units, res)
    new <- list(...)
    type <- if(!missing(type)) match.arg(type) else getOption("bitmapType")
    if(!missing(antialias)) new$antialias <- match.arg(antialias, aa.cairo)
    d <- check.options(new, name.opt = ".X11.Options", envir = .X11env)
    antialias <- match(d$antialias, aa.cairo)
    if(type == "quartz") {
        if(capabilities("aqua")) {
            width <- g$width/ifelse(is.na(res), 72, res);
            height <- g$height/ifelse(is.na(res), 72, res);
            invisible(.External(C_Quartz, "jpeg", path.expand(filename),
                                width, height, pointsize, d$family,
                                d$antialias != "none", "", bg,
                                "white", if(is.na(res)) NULL else res))
        } else warning('type = "quartz" is unavailable. trying "Xlib"')
    } else if (type == "cairo") {
        if(capabilities("cairo")) {
            invisible(.External(C_devCairo, filename, 3L, g$width, g$height,
                                pointsize, bg, res, antialias, quality, d$family,
                                300, optionSymbolFont(d$symbolfamily)))
        } else warning('type = "cairo" is unavailable. trying "Xlib"')
    } else
        invisible(.External2(C_X11,
                            paste0("jpeg::", quality, ":", filename),
                            g$width, g$height, pointsize, d$gamma,
                            d$colortype, d$maxcubesize, bg, bg, d$fonts, res,
                            0L, 0L, "", 0, 0, d$family,
                            optionSymbolFont(d$symbolfamily)))
}

tiff <- function(filename = "Rplot%03d.tiff",
                 width = 480, height = 480, units = "px", pointsize = 12,
                 compression = c("none", "rle", "lzw", "jpeg", "zip",
                                 "lzw+p", "zip+p"),
                 bg = "white", res = NA, ...,
                 type = c("cairo", "Xlib", "quartz"), antialias)
{
    if(!is.character(filename) || length(filename) != 1L || !nzchar(filename))
        stop("'filename' must be a non-empty character string")
    if(!checkIntFormat(filename)) stop("invalid 'filename'")
    g <- .geometry(width, height, units, res)
    new <- list(...)
    type <- if(!missing(type)) match.arg(type) else getOption("bitmapType")
    if(!missing(antialias)) new$antialias <- match.arg(antialias, aa.cairo)
    d <- check.options(new, name.opt = ".X11.Options", envir = .X11env)
    antialias <- match(d$antialias, aa.cairo)
    comp <- switch( match.arg(compression),
                   "none" = 1L, "rle" = 2L, "lzw" = 5L, "jpeg" = 7L, "zip" = 8L,
                   "lzw+p" = 15L, "zip+p" = 18L)
    if(type == "quartz") {
        if(capabilities("aqua")) {
            width <- g$width/ifelse(is.na(res), 72, res);
            height <- g$height/ifelse(is.na(res), 72, res);
            if (comp != 1L)
                warning('compression is not supported for type = "quartz"')
            invisible(.External(C_Quartz, "tiff", path.expand(filename),
                                width, height, pointsize, d$family,
                                d$antialias != "none", "", bg,
                                "white", if(is.na(res)) NULL else res))
        } else warning('type = "quartz" is unavailable. trying "Xlib"')
    } else if (type == "cairo") {
        if(capabilities("cairo")) {
            invisible(.External(C_devCairo, filename, 8L, g$width, g$height,
                                pointsize, bg, res, antialias, comp, d$family,
                                300, optionSymbolFont(d$symbolfamily)))
        } else warning('type = "cairo" is unavailable. trying "Xlib"')
    } else
        invisible(.External2(C_X11,
                             paste0("tiff::", comp, ":", filename),
                             g$width, g$height, pointsize, d$gamma,
                             d$colortype, d$maxcubesize, bg, bg, d$fonts, res,
                             0L, 0L, "", 0, 0, d$family,
                             optionSymbolFont(d$symbolfamily)))
}

bmp <- function(filename = "Rplot%03d.bmp",
                width = 480, height = 480, units = "px", pointsize = 12,
                bg = "white", res = NA, ...,
                type = c("cairo", "Xlib", "quartz"), antialias)
{
    if(!is.character(filename) || length(filename) != 1L || !nzchar(filename))
        stop("'filename' must be a non-empty character string")
    if(!checkIntFormat(filename)) stop("invalid 'filename'")
    g <- .geometry(width, height, units, res)
    new <- list(...)
    type <- if(!missing(type)) match.arg(type) else getOption("bitmapType")
    if(!missing(antialias)) new$antialias <- match.arg(antialias, aa.cairo)
    d <- check.options(new, name.opt = ".X11.Options", envir = .X11env)
    antialias <- match(d$antialias, aa.cairo)
    if(type == "quartz") {
        if(capabilities("aqua")) {
        width <- g$width/ifelse(is.na(res), 72, res);
        height <- g$height/ifelse(is.na(res), 72, res);
        invisible(.External(C_Quartz, "bmp", path.expand(filename),
                            width, height, pointsize, d$family,
                            d$antialias != "none", "", bg,
                            "white", if(is.na(res)) NULL else res))
        } else warning('type = "quartz" is unavailable. trying "Xlib"')
    } else if (type == "cairo") {
        if(capabilities("cairo")) {
            invisible(.External(C_devCairo, filename, 9L, g$width, g$height,
                                pointsize, bg, res, antialias, 100L, d$family,
                                300, optionSymbolFont(d$symbolfamily)))
        } else warning('type = "cairo" is unavailable. trying "Xlib"')
    } else
        invisible(.External2(C_X11, paste0("bmp::", filename),
                             g$width, g$height, pointsize, d$gamma,
                             d$colortype, d$maxcubesize, bg, bg, d$fonts, res,
                             0L, 0L, "", 0, 0, d$family,
                             optionSymbolFont(d$symbolfamily)))
}

grSoftVersion <- function() {
    bm <- .Call(C_bmVersion)
    if(nzchar(bm[3L])) bm[3L] <- strsplit(bm[3L], "\n")[[1L]][1L]
    c(cairo = cairoVersion(), cairoFT = cairoFT(), pango = pangoVersion(), bm)
}
#  File src/library/grDevices/R/unix/quartz.R
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

.Quartzenv <- new.env()

assign(".quartz.Options",
       list(title = "Quartz %d",
            width = 7, height = 7, pointsize = 12,
            family = "Helvetica",
            antialias = TRUE,
            type = "native",
            bg = "transparent", canvas = "white",
            dpi = NA_real_),
       envir = .Quartzenv)

assign(".quartz.Options.default",
       get(".quartz.Options", envir = .Quartzenv),
       envir = .Quartzenv)

quartz.options <- function(..., reset = FALSE)
{
    old <- get(".quartz.Options", envir = .Quartzenv)
    if(reset) {
        assign(".quartz.Options",
               get(".quartz.Options.default", envir = .Quartzenv),
               envir = .Quartzenv)
    }
    l... <- length(new <- list(...))
    check.options(new, name.opt = ".quartz.Options", envir = .Quartzenv,
                  assign.opt = l... > 0L)
    if(reset || l... > 0L) invisible(old) else old
}

quartz <- function(title, width, height, pointsize, family, antialias,
                   type, file = NULL, bg, canvas, dpi)
{
    if (missing(type) || type %in% c("", "native", "Cocoa")) {
        check <- Sys.getenv("_R_CHECK_SCREEN_DEVICE_", "")
        msg <- "screen devices should not be used in examples etc"
        if (identical(check, "stop"))
            stop(msg, domain = NA)
        else if (identical(check, "warn"))
            warning(msg, immediate. = TRUE, noBreaks. = TRUE, domain = NA)
    }

    new <- list()
    if(!missing(title)) new$title <- title
    if(!missing(width)) new$width <- width
    if(!missing(height)) new$height <- height
    if(!missing(pointsize)) new$pointsize <- pointsize
    if(!missing(family)) new$family <- family
    if(!missing(antialias)) new$antialias <- antialias
    if(!missing(bg)) new$bg <- bg
    if(!missing(canvas)) new$canvas <- canvas
    if(!missing(type)) new$type <- type
    if(!missing(dpi)) new$dpi <- dpi
    if(!checkIntFormat(new$title)) stop("invalid 'title'")
    if(!is.null(file) && !checkIntFormat(file)) stop("invalid 'file'")
    d <- check.options(new, name.opt = ".quartz.Options", envir = .Quartzenv)
    .External(C_Quartz, d$type, file, d$width, d$height, d$pointsize, d$family,
              d$antialias, d$title, d$bg, d$canvas,
              if(is.na(d$dpi)) NULL else d$dpi)
    invisible()
}

#########
# QUARTZ font database
# To map device-independent font to device-specific font
#########

# Each font family has only a name
assign(".Quartz.Fonts", list(), envir = .Quartzenv)

# Check that the font has the correct structure and information
checkQuartzFont <- function(font) {
    if (!is.character(font) || length(font) != 4)
        stop("invalid Quartz font:  must be 4 strings")
    font
}

setQuartzFonts <- function(fonts, fontNames) {
    fonts <- lapply(fonts, checkQuartzFont)
    fontDB <- get(".Quartz.Fonts", envir=.Quartzenv)
    existingFonts <- fontNames %in% names(fontDB)
    if (sum(existingFonts) > 0L)
        fontDB[fontNames[existingFonts]] <- fonts[existingFonts]
    if (sum(existingFonts) < length(fontNames))
        fontDB <- c(fontDB, fonts[!existingFonts])
    assign(".Quartz.Fonts", fontDB, envir=.Quartzenv)
}

printFont <- function(font) {
    paste0(font, "\n")
}

printFonts <- function(fonts) {
    cat(paste0(names(fonts), ": ", unlist(lapply(fonts, printFont)),
               collapse=""))
}

# If no arguments spec'ed, return entire font database
# If no named arguments spec'ed, all args should be font names
# to get info on from the database
# Else, must specify new fonts to enter into database (all
# of which must be valid PostScript font descriptions and
# all of which must be named args)
quartzFonts <- function(...) {
    ndots <- length(fonts <- list(...))
    if (ndots == 0L)
        get(".Quartz.Fonts", envir=.Quartzenv)
    else {
        fontNames <- names(fonts)
        nnames <- length(fontNames)
        if (nnames == 0L) {
            if (!all(sapply(fonts, is.character)))
                stop("invalid arguments in 'quartzFonts' (must be font names)")
            else
                get(".Quartz.Fonts", envir=.Quartzenv)[unlist(fonts)]
        } else {
            if (ndots != nnames)
                stop("invalid arguments in 'quartzFonts' (need named args)")
            setQuartzFonts(fonts, fontNames)
        }
    }
}

# Create a valid quartz font description
quartzFont <- function(family) {
    checkQuartzFont(family)
}

quartzFonts(# Default Serif font is Times
            serif = quartzFont(c("Times-Roman", "Times-Bold",
            "Times-Italic", "Times-BoldItalic")),
            ## Default Sans Serif font is Helvetica,
            ## even the device default is Arial
            sans = quartzFont(c("Helvetica", "Helvetica-Bold",
            "Helvetica-Oblique", "Helvetica-BoldOblique")),
            ## Default Monospace font is Courier
            mono = quartzFont(c("Courier", "Courier-Bold",
            "Courier-Oblique", "Courier-BoldOblique")))

## Formerly for R.app only
quartz.save <- function(file, type = 'png', device = dev.cur(), dpi = 100, ...)
{
    ## modified version of dev.copy2pdf
    dev.set(device)
    current.device <- dev.cur()
    nm <- names(current.device)[1L]
    if (nm == "null device") stop("no device to print from")
    if (!dev.displaylist()) stop("can only print from a screen device")
    oc <- match.call()
    oc[[1L]] <- as.name("dev.copy")
    oc$file <- NULL
    oc$device <- quartz
    oc$type <- type
    if(missing(file)) file <- paste0("Rplot.", type)
    oc$file <- file
    oc$dpi <- dpi
    din <- dev.size("in")
    w <- din[1L]
    h <- din[2L]
    if (is.null(oc$width))
        oc$width <- if (!is.null(oc$height)) w/h * eval.parent(oc$height) else w
    if (is.null(oc$height))
        oc$height <- if (!is.null(oc$width)) h/w * eval.parent(oc$width) else h
    on.exit(dev.set(current.device))
    dev.off(eval.parent(oc))
}
#  File src/library/grDevices/R/unix/x11.R
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

## An environment not exported from namespace:grDevices used to
## pass .X11.Fonts to the X11 device.
.X11env <- new.env()

assign(".X11.Options",
       list(display = "",
            width = NA_real_, height = NA_real_, pointsize = 12,
            bg = "transparent", canvas = "white",
            gamma = 1,
            colortype = "true", maxcubesize = 256,
            fonts = c("-adobe-helvetica-%s-%s-*-*-%d-*-*-*-*-*-*-*",
            "-adobe-symbol-medium-r-*-*-%d-*-*-*-*-*-*-*"),
            family = "sans",
            symbolfamily = "default",
            xpos = NA_integer_, ypos = NA_integer_,
	    title = "", type = "cairo", antialias = "default"),
       envir = .X11env)

assign(".X11.Options.default",
       get(".X11.Options", envir = .X11env),
       envir = .X11env)

aa.cairo  <- c("default", "none", "gray", "subpixel")

X11.options <- function(..., reset = FALSE)
{
    old <- get(".X11.Options", envir = .X11env)
    if(reset) {
        assign(".X11.Options",
               get(".X11.Options.default", envir = .X11env),
               envir = .X11env)
    }
    l... <- length(new <- list(...))
    check.options(new, name.opt = ".X11.Options", envir = .X11env,
                  assign.opt = l... > 0)
    if(reset || l... > 0) invisible(old) else old
}

x11 <-
X11 <- function(display = "", width, height, pointsize, gamma,
                bg, canvas, fonts, family,
                xpos, ypos, title, type, antialias, symbolfamily)
{
    if(display != "XImage") { # used by tkrplot
        check <- Sys.getenv("_R_CHECK_SCREEN_DEVICE_", "")
        msg <- "screen devices should not be used in examples etc"
        if (identical(check, "stop"))
            stop(msg, domain = NA)
        else if (identical(check, "warn"))
            warning(msg, immediate. = TRUE, noBreaks. = TRUE, domain = NA)
    }

    if(display == "" && .Platform$GUI == "AQUA" &&
       is.na(Sys.getenv("DISPLAY", NA))) Sys.setenv(DISPLAY = ":0")

    new <- list()
    if(!missing(display)) new$display <- display
    if(!missing(width)) new$width <- width
    if(!missing(height)) new$height <- height
    if(!missing(gamma)) new$gamma <- gamma
    if(!missing(pointsize)) new$pointsize <- pointsize
    if(!missing(bg)) new$bg <- bg
    if(!missing(canvas)) new$canvas <- canvas
    if(!missing(xpos)) new$xpos <- xpos
    if(!missing(ypos)) new$ypos <- ypos
    if(!missing(title)) new$title <- title
    if(!checkIntFormat(new$title)) stop("invalid 'title'")
    if(!missing(type)) {
        new$type <- match.arg(type, c("Xlib", "cairo", "nbcairo", "dbcairo"))
        if(!capabilities("cairo") && type != "Xlib")
            warning("cairo-based types are not supported on this build - using \"Xlib\"")
    }
    if(!missing(family)) new$family <- family
    if(!missing(fonts)) new$fonts <- fonts
    if(!missing(antialias) && type != "Xlib")
        new$antialias <- match.arg(antialias, aa.cairo)
    if(!missing(symbolfamily)) new$symbolfamily <- symbolfamily
    d <- check.options(new, name.opt = ".X11.Options", envir = .X11env)
    if(d$type == "Xlib" && !missing(family)) {
        fns <- X11Fonts()
        if (! family %in% names(fns))
            stop('unknown family for X11(type = "XLib")')
        d$fonts[1] <- fns[[family]]
    }
    type <-
	if(capabilities("cairo"))
            switch(d$type, "cairo" = 1L, "nbcairo" = 2L, "dbcairo" = 3L, 0L)
	else 0L
    ## Aargh -- trkplot has a trapdoor and does not set type.
    if (display == "XImage") type <- 0L
    antialias <- match(d$antialias, aa.cairo)
    if (grepl("darwin", R.version$os))
        utils:::check_for_XQuartz(file.path(R.home("modules"), "R_X11.so"))
    .External2(C_X11, d$display, d$width, d$height, d$pointsize, d$gamma,
               d$colortype, d$maxcubesize, d$bg, d$canvas, d$fonts,
               NA_integer_, d$xpos, d$ypos, d$title,
               type, antialias, d$family, optionSymbolFont(d$symbolfamily))
    invisible()
}


####################
# X11 font database
####################

assign(".X11.Fonts", list(), envir = .X11env)

X11FontError <- function(errDesc)
    stop("invalid X11 font specification: ", errDesc)


# Check that the font has the correct structure and information
# Already checked that it had a name
checkX11Font <- function(font)
{
    if (!is.character(font))
        X11FontError("must be a string")
    ## Check it has the right format
    if (length(grep("(-[^-]+){14}", font)) > 0) {
        ## Force the %s and %d substitution formats into the right spots
        font <- sub("((-[^-]+){2})(-[^-]+){2}((-[^-]+){2})(-[^-]+)((-[^-]+){7})",
                    "\\1-%s-%s\\4-%d\\7", font, perl = TRUE)
    } else {
        X11FontError("incorrect format")
    }
    font
}

setX11Fonts <- function(fonts, fontNames)
{
    fonts <- lapply(fonts, checkX11Font)
    fontDB <- get(".X11.Fonts", envir=.X11env)
    existingFonts <- fontNames %in% names(fontDB)
    if (sum(existingFonts) > 0)
        fontDB[fontNames[existingFonts]] <- fonts[existingFonts]
    if (sum(existingFonts) < length(fontNames))
        fontDB <- c(fontDB, fonts[!existingFonts])
    assign(".X11.Fonts", fontDB, envir=.X11env)
}

printFont <- function(font) paste0(font, "\n")


printFonts <- function(fonts)
    cat(paste0(names(fonts), ": ", unlist(lapply(fonts, printFont)),
               collapse=""))

# If no arguments spec'ed, return entire font database
# If no named arguments spec'ed, all args should be font names
# to get info on from the database
# Else, must specify new fonts to enter into database (all
# of which must be valid X11 font descriptions and
# all of which must be named args)
X11Fonts <- function(...)
{
    ndots <- length(fonts <- list(...))
    if (ndots == 0)
        get(".X11.Fonts", envir=.X11env)
    else {
        fontNames <- names(fonts)
        nnames <- length(fontNames)
        if (nnames == 0) {
            if (!all(sapply(fonts, is.character)))
                stop("invalid arguments in 'X11Fonts' (must be font names)")
            else
                get(".X11.Fonts", envir=.X11env)[unlist(fonts)]
        } else {
            if (ndots != nnames)
                stop("invalid arguments in 'X11Fonts' (need named args)")
            setX11Fonts(fonts, fontNames)
        }
    }
}

# Create a valid X11 font description
X11Font <- function(font) checkX11Font(font)

X11Fonts(# Default Serif font is Times
         serif = X11Font("-*-times-%s-%s-*-*-%d-*-*-*-*-*-*-*"),
         # Default Sans Serif font is Helvetica
         sans = X11Font("-*-helvetica-%s-%s-*-*-%d-*-*-*-*-*-*-*"),
         # Default Monospace font is Courier
         mono = X11Font("-*-courier-%s-%s-*-*-%d-*-*-*-*-*-*-*"),
         Times = X11Font("-adobe-times-%s-%s-*-*-%d-*-*-*-*-*-*-*"),
         Helvetica = X11Font("-adobe-helvetica-%s-%s-*-*-%d-*-*-*-*-*-*-*"),
         CyrTimes = X11Font("-cronyx-times-%s-%s-*-*-%d-*-*-*-*-*-*-*"),
         CyrHelvetica = X11Font("-cronyx-helvetica-%s-%s-*-*-%d-*-*-*-*-*-*-*"),
         Arial = X11Font("-monotype-arial-%s-%s-*-*-%d-*-*-*-*-*-*-*"),
         Mincho = X11Font("-*-mincho-%s-%s-*-*-%d-*-*-*-*-*-*-*")
         )

savePlot <- function(filename = paste0("Rplot.", type),
                     type = c("png", "jpeg", "tiff", "bmp"),
                     device = dev.cur())
{
    type <- match.arg(type)
    devlist <- dev.list()
    devcur <- match(device, devlist, NA)
    if(is.na(devcur)) stop("no such device")
    devname <- names(devlist)[devcur]
    if(devname != "X11cairo")
        stop("can only copy from 'X11(type=\"*cairo\")' devices")
    invisible(.External2(C_savePlot, filename, type, device))
}


optionSymbolFont <- function(family) {
    if (family == "default") {
        if (symbolType1support()) {
            cairoSymbolFont("symbol")
        } else {
            cairoSymbolFont("sans", usePUA = FALSE)
        }
    } else {
        checkSymbolFont(family)
    }
}
#  File src/library/grDevices/R/utils.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2023 The R Core Team
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

n2mfrow <- function(nr.plots, asp = 1)
{
  nr.plots <- as.integer(nr.plots)
  if(asp == 1 && nr.plots <= 12L) { # if only for back compatibility
      if   (nr.plots <=  3L)
          c(nr.plots, 1L)        # nr.p = 1, 2, 3
      else if (nr.plots <=  6L)
           c((nr.plots+1L) %/% 2L, 2L)# = 4, 5, 6
      else c((nr.plots+2L) %/% 3L, 3L)# = 7, 8,.., 12
  }
  else {
      rc <- c(nrow <- ceiling(sqrt(nr.plots / asp)),
              ceiling(nr.plots / nrow))
      ## improve (if easy):
      do <- TRUE
      while(do && prod(rc) > nr.plots) {
          if(do <- prod(n <- rc - 0:1) >= nr.plots)
              rc <- n
          else if(do <- prod(n <- rc - 1:0) >= nr.plots)
              rc <- n
      }
      rc
  }
}

extendrange <- function(x, r = range(x, na.rm = TRUE), f = 0.05)
{
    ## Purpose: extend a range by a factor 'f' - on each side
    if(!missing(r) && length(r) != 2)
        stop("'r' must be a \"range\", hence of length 2")
    f <- if(length(f) == 1L) c(-f,f) else c(-f[1L], f[2L])
    r + f * diff(r)
}

trans3d <- function(x,y,z, pmat, continuous = FALSE, verbose = TRUE) {
    tr <- cbind(x,y,z,1, deparse.level=0L) %*% pmat
    if(continuous && (n <- nrow(tr)) >= 2) {
        st4 <- sign(tr[,4])
        if((s1 <- st4[1]) != st4[n]) { # have a sign change ==> cut off at sign switch
            if((last <- (which.min(st4 == s1) - 1L)) >= 1L) { # needed? -- safe programming!
                if(verbose) message(sprintf("points cut off after point[%d]", last))
                tr <- tr[seq_len(last), , drop=FALSE]
            }
        }
    }
    list(x = tr[,1]/tr[,4],
	 y = tr[,2]/tr[,4])
}
#  File src/library/grDevices/R/xyz.coords.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2022 The R Core Team
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

## Both xy.coords() and xyz.coords()  --- should be kept in sync!

warnLogCoords <- function(coord, n) {
    warning(warningCondition(sprintf(
        ngettext(n,
                 "%d %s value <= 0 omitted from logarithmic plot",
                 "%d %s values <= 0 omitted from logarithmic plot"),
        n, coord),
        call = sys.call(-1L), coord = coord,
        class = c("log_le_0", "plot_coords")))
}

xy.coords <-
    function(x, y = NULL, xlab = NULL, ylab = NULL, log = NULL, recycle = FALSE,
             setLab = TRUE)
{
    if(is.null(y)) {
	if(is.null(ylab)) ylab <- xlab
	if(is.language(x)) {
	    if (inherits(x, "formula") && length(x) == 3) {
                if(setLab) {
                    ylab <- deparse(x[[2L]])
                    xlab <- deparse(x[[3L]])
                }
		y <- eval(x[[2L]], environment(x))
		x <- eval(x[[3L]], environment(x))
	    }
	    else stop("invalid first argument")
	}
	else if(inherits(x, "ts")) {
	    y <- if(is.matrix(x)) x[,1] else x
	    x <- stats::time(x)
            if(setLab) xlab <- "Time"
	}
	else if(is.complex(x)) {
	    y <- Im(x)
	    x <- Re(x)
            if(setLab) {
                xlab <- paste0("Re(", ylab, ")")
                ylab <- paste0("Im(", ylab, ")")
            }
	}
	else if(is.matrix(x) || is.data.frame(x)) {
	    x <- data.matrix(x)
	    if(ncol(x) == 1) {
		if(setLab) xlab <- "Index"
		y <- x[,1]
		x <- seq_along(y)
	    }
	    else {
		colnames <- dimnames(x)[[2L]]
                if(setLab) {
                    if(is.null(colnames)) {
                        xlab <- paste0(ylab, "[,1]")
                        ylab <- paste0(ylab, "[,2]")
                    }
                    else {
                        xlab <- colnames[1L]
                        ylab <- colnames[2L]
                    }
                }
		y <- x[,2]
		x <- x[,1]
	    }
	}
	else if(is.list(x)) {
            if (all(c("x", "y") %in% names(x))) {
                if(setLab) {
                    xlab <- paste0(ylab, "$x")
                    ylab <- paste0(ylab, "$y")
                }
                y <- x[["y"]]
                x <- x[["x"]]
            } else
                stop("'x' is a list, but does not have components 'x' and 'y'")
	}
	else {
	    if(is.factor(x)) x <- as.numeric(x)
	    if(setLab) xlab <- "Index"
	    y <- x
	    x <- seq_along(x)
	}
    }
    ## to allow e.g. lines, points, identify to be used with plot.POSIXlt
    if(inherits(x, "POSIXt")) x <- as.POSIXct(x)

    if(length(x) != length(y)) {
	if(recycle) {
	    if((nx <- length(x)) < (ny <- length(y)))
		x <- rep_len(x, ny)
	    else
		y <- rep_len(y, nx)
	}
	else
	    stop("'x' and 'y' lengths differ")
    }

    if(length(log) && log != "") {
	log <- strsplit(log, NULL)[[1L]]
	if("x" %in% log && any(ii <- x <= 0 & !is.na(x))) {
	    warnLogCoords("x", n = as.integer(sum(ii)))
	    x[ii] <- NA
	}
	if("y" %in% log && any(ii <- y <= 0 & !is.na(y))) {
	    warnLogCoords("y", n = as.integer(sum(ii)))
	    y[ii] <- NA
	}
    }
    list(x=as.double(x), y=as.double(y), xlab=xlab, ylab=ylab)
}

xyz.coords <- function(x, y=NULL, z=NULL, xlab=NULL, ylab=NULL, zlab=NULL,
		       log = NULL, recycle = FALSE, setLab = TRUE)
{
    ## Only x
    if(is.null(y)) {
	if (is.language(x)) {
	    if (inherits(x, "formula") && length(x) == 3
		&& length(rhs <- x[[3L]]) == 3) {
                if(setLab) {
                    zlab <- deparse(x[[2L]])
                    ylab <- deparse(rhs[[3L]])
                    xlab <- deparse(rhs[[2L]])
                }
		pf <- parent.frame()
		z <- eval(x[[2L]],   environment(x), pf)
		y <- eval(rhs[[3L]], environment(x), pf)
		x <- eval(rhs[[2L]], environment(x), pf)
	    }
	    else stop("invalid first argument [bad language object]")
	}
	else if(is.matrix(x) || is.data.frame(x)) {
	    x <- data.matrix(x)
	    if(ncol(x) < 2) stop("at least 2 columns needed")
	    if(ncol(x) == 2) {
		if(setLab) xlab <- "Index"
		y <- x[,1]
		z <- x[,2]
		x <- seq_along(y)
	    }
	    else { ## >= 3 columns
		colnames <- dimnames(x)[[2L]]
                if(setLab) {
                    if(is.null(colnames)) {
                        zlab <- paste0(xlab,"[,3]")
                        ylab <- paste0(xlab,"[,2]")
                        xlab <- paste0(xlab,"[,1]")
                    }
                    else {
                        xlab <- colnames[1L]
                        ylab <- colnames[2L]
                        zlab <- colnames[3L]
                    }
                }
		y <- x[,2]
		z <- x[,3]
		x <- x[,1]
	    }
	}
	else if(is.list(x)) {
            if (all(c("x", "y", "z") %in% names(x))) {
                if(setLab) {
                    zlab <- paste0(xlab,"$z")
                    ylab <- paste0(xlab,"$y")
                    xlab <- paste0(xlab,"$x")
                }
                y <- x[["y"]]
                z <- x[["z"]]
                x <- x[["x"]]
            } else
                stop("'x' is a list, but does not have components 'x', 'y'  and 'z'")
        }
    }

    ## Only x, y
    if(!is.null(y) && is.null(z)) {
	if(is.complex(x)) {
	    z <- y
	    y <- Im(x)
	    x <- Re(x)
            if(setLab) {
                zlab <- ylab
                ylab <- paste0("Im(", xlab, ")")
                xlab <- paste0("Re(", xlab, ")")
            }
	}
	else if(is.complex(y)) {
	    z <- x
	    x <- Re(y)
	    y <- Im(y)
            if(setLab) {
                zlab <- xlab
                xlab <- paste0("Re(", ylab, ")")
                ylab <- paste0("Im(", ylab, ")")
            }
	}
	else {
	    if(is.factor(x)) x <- as.numeric(x)
	    if(is.factor(y)) y <- as.numeric(y)
            if(setLab) xlab <- "Index"
	    z <- y
	    y <- x
	    x <- seq_along(x)
	}
    }

    ## Lengths and recycle
    if(((xl <- length(x)) != length(y)) || (xl != length(z))) {
	if(recycle) {
	    ml <- max(xl, (yl <- length(y)), (zl <- length(z)))
	    if(xl < ml && !is.null(x)) x <- rep_len(x, ml)
	    if(yl < ml && !is.null(y)) y <- rep_len(y, ml)
	    if(zl < ml && !is.null(z)) z <- rep_len(z, ml)
	}
	else stop("'x', 'y' and 'z' lengths differ")
    }

    ## log
    if(length(log) && log != "") {
	log <- strsplit(log, NULL)[[1L]]
	if("x" %in% log && any(ii <- x <= 0 & !is.na(x))) {
	    warnLogCoords("x", n = sum(ii))
	    x[ii] <- NA
	}
	if("y" %in% log && any(ii <- y <= 0 & !is.na(y))) {
	    warnLogCoords("y", n = sum(ii))
	    y[ii] <- NA
	}
	if("z" %in% log && any(ii <- z <= 0 & !is.na(z))) {
	    warnLogCoords("z", n = sum(ii))
	    z[ii] <- NA
	}
    }
    list(x=as.double(x), y=as.double(y), z=as.double(z),
	 xlab=xlab, ylab=ylab, zlab=zlab)
}
#  File src/library/grDevices/R/zzz.R
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

if (.Platform$OS.type == "windows") {
    utils::globalVariables(c("C_cairoProps", "C_makeQuartzDefault"))
    utils::suppressForeignCheck(c("C_cairoProps", "C_makeQuartzDefault"))
}

.onLoad <- function(libname, pkgname)
{
    if (.Platform$OS.type != "windows" && !.Call(C_cairoProps, 2L))
        X11.options(type = "Xlib")

    extras <- if(.Platform$OS.type == "windows")
        list(windowsTimeouts = c(100L,500L)) else
        list(bitmapType = if(capabilities("aqua")) "quartz"
        else if(.Call(C_cairoProps, 2L)) "cairo" else "Xlib")
    op.grDevices <- c(list(locatorBell = TRUE, device.ask.default = FALSE),
                      extras, list(device = .select_device()))
    toset <- !(names(op.grDevices) %in% names(.Options))
    if(any(toset)) options(op.grDevices[toset])
}

.onUnload <- function(libpath)
    library.dynam.unload("grDevices", libpath)


### Used by text, mtext, strwidth, strheight, title, axis,
### L_text and L_textBounds, all of which
### coerce SYMSXPs and LANGSXPs to EXPRSXPs
### We don't want to use as.expression here as that is generic
### even though is.language no longer is

### Possibly later have
### if (is.language(x)) x
### else if(isS4(x)) methods::as(x, "character")
### else if(is.object(x)) as.character(x)
### else x

as.graphicsAnnot <- function(x)
    if(is.language(x) || !is.object(x)) x else as.character(x)
