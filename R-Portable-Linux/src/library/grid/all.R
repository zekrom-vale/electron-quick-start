
createClipPath <- function(clip) {
    force(clip)
    pathFun <- function() {
        grid.draw(clip$grob, recording=FALSE)
    }
    path <- list(f=.clipPath(pathFun, clip$rule), ref=NULL)
    class(path) <- "GridClipPath"
    path
}   

isClipPath <- function(x) {
    inherits(x, "GridClipPath")
}

## "resolve" clipping paths
resolveClipPath <- function(path) {
    ref <- .setClipPath(path$f, path$ref)
    resolvedClipPath(path, ref)
}

resolvedClipPath <- function(path, ref) {
    UseMethod("resolvedClipPath")
}

resolvedClipPath.GridClipPath <- function(path, ref) {
    path$ref <- ref
    class(path) <- c("GridResolvedClipPath", class(path))
    path
}    

resolvedClipPath.GridResolvedClipPath <- function(path, ref) {
    path$ref <- ref
    path
}    

unresolveClipPath <- function(path) {
    UseMethod("unresolveClipPath")
}
    
## Unresolved clipPaths just pass through
unresolveClipPath.GridClipPath <- function(path) {
    path
}

unresolveClipPath.GridResolvedClipPath <- function(path) {
    path <- list(f=path$f, ref=NULL)
    class(path) <- "GridClipPath"
    path
}

#  File src/library/grid/R/components.R
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

grid.collection <- function(..., gp=gpar(), draw=TRUE, vp=NULL) {
    .Defunct("gTree")
}

######################################
# AXES
######################################

# Axes are extended from the "gTree" class
# This means that the standard (e.g., draw.details)
# methods for gTrees will apply

# The children of an axis are fixed to be:

# NOTE that the `at' parameter is numeric (i.e., NOT a unit) for
# grid.xaxis and grid.yaxis.  These functions assume a unit for the `at'
# values rather than letting the user specify a unit.

validDetails.axis <- function(x) {
  if (!is.null(x$at)) {
    x$at <- as.numeric(x$at)
    if (length(x$at) < 1 || !all(is.finite(x$at)))
      stop("invalid 'at' location in 'axis'")
  }
  if (!is.logical(x$label)) {
    # labels specified
    # Can only spec labels if at is not NULL
    if (is.null(x$at))
      stop("invalid to specify axis labels when 'at' is NULL")
    # Must be either language object or string
    x$label <- as.graphicsAnnot(x$label)
    # Must be same number of labels as "at" locations
    if (length(x$label) != length(x$at))
      stop("'labels' and 'at' locations must have same length")
  }
  x$main <- as.logical(x$main)
  x
}

makeContent.xaxis <- function(x) {
    # If x$at is NULL, then we must calculate the
    # tick marks on-the-fly
    if (is.null(x$at)) {
        x$at <- grid.pretty(current.viewport()$xscale)
        # Add the new output as children
        x <- addGrob(x, make.xaxis.major(x$at, x$main))
        x <- addGrob(x, make.xaxis.ticks(x$at, x$main))
        x <- updateXlabels(x)
        # Apply any edits relevant to children
        x <- applyEdits(x, x$edits)
    }
    x
}

# NOTE that this can't be for all axes because it needs to
# call make.XAXIS.ticks and make.XAXIS.labels
editDetails.xaxis <- function(x, specs) {
  slot.names <- names(specs)
  if ("at" %in% slot.names) {
    # NOTE that grid.edit has already set x$at to the new value
    # We might set at to NULL to get ticks recalculated at redraw
    if (is.null(x$at)) {
      x <- removeGrob(x, "major", warn=FALSE)
      x <- removeGrob(x, "ticks", warn=FALSE)
      x <- removeGrob(x, "labels", warn=FALSE)
    } else {
      x <- addGrob(x, make.xaxis.major(x$at, x$main))
      x <- addGrob(x, make.xaxis.ticks(x$at, x$main))
      x <- updateXlabels(x)
    }
  }
  if ("label" %in% slot.names) {
    if (!is.null(x$at))
      x <- updateXlabels(x)
  }
  if ("main" %in% slot.names)
    if (!is.null(x$at)) {
      x <- addGrob(x, make.xaxis.major(x$at, x$main))
      x <- addGrob(x, make.xaxis.ticks(x$at, x$main))
      x <- updateXlabels(x)
    }
  x
}

make.xaxis.major <- function(at, main) {
  if (main)
    y <- c(0, 0)
  else
    y <- c(1, 1)
  linesGrob(unit(c(min(at), max(at)), "native"),
            unit(y, "npc"), name="major")
}

make.xaxis.ticks <- function(at, main) {
  if (main) {
    tick.y0 <- unit(0, "npc")
    tick.y1 <- unit(-.5, "lines")
  }
  else {
    tick.y0 <- unit(1, "npc")
    tick.y1 <- unit(1, "npc") + unit(.5, "lines")
  }
  segmentsGrob(unit(at, "native"), tick.y0,
               unit(at, "native"), tick.y1,
               name="ticks")
}

make.xaxis.labels <- function(at, label, main) {
  # FIXME:  labels only character versions of "at"
  if (main)
    label.y <- unit(-1.5, "lines")
  else
    label.y <- unit(1, "npc") + unit(1.5, "lines")
  if (is.logical(label))
    labels <- as.character(at)
  else
    labels <- label
  textGrob(labels, unit(at, "native"), label.y,
           just="centre", rot=0,
           check.overlap=TRUE, name="labels")
}

updateXlabels <- function(x) {
  if (is.logical(x$label) && !x$label)
    removeGrob(x, "labels", warn=FALSE)
  else
    addGrob(x, make.xaxis.labels(x$at, x$label, x$main))
}

xaxisGrob <- function(at=NULL, label=TRUE, main=TRUE,
                      edits=NULL,
                      name=NULL, gp=gpar(), vp=NULL) {
  grid.xaxis(at=at, label=label, main=main,
             edits=edits,
             name=name, gp=gp, draw=FALSE, vp=vp)
}

# The "main" x-axis is on the bottom when vp$origin is "bottom.*"
# and on the top when vp$origin is "top.*"
grid.xaxis <- function(at=NULL, label=TRUE, main=TRUE,
                       edits=NULL, name=NULL, gp=gpar(),
                       draw=TRUE, vp=NULL) {
  if (is.null(at)) {
    # We do not have enough information to make the ticks and labels
    major <- NULL
    ticks <- NULL
    labels <- NULL
  } else {
    major <- make.xaxis.major(at, main)
    ticks <- make.xaxis.ticks(at, main)
    if (is.logical(label) && length(label) == 0)
	stop("logical 'label' supplied of length 0")
    if (is.logical(label) && !label)
      labels <- NULL
    else
      labels <- make.xaxis.labels(at, label, main)
  }
  xg <- applyEdits(gTree(at=at, label=label, main=main,
                         children=gList(major, ticks, labels),
                         edits=edits,
                         name=name, gp=gp, vp=vp,
                         cl=c("xaxis", "axis")),
                   edits)
  if (draw)
    grid.draw(xg)
  invisible(xg)
}

makeContent.yaxis <- function(x) {
    # If x$at is NULL, then we must calculate the
    # tick marks on-the-fly
    if (is.null(x$at)) {
        x$at <- grid.pretty(current.viewport()$yscale)
        # Add the new output as children
        x <- addGrob(x, make.yaxis.major(x$at, x$main))
        x <- addGrob(x, make.yaxis.ticks(x$at, x$main))
        x <- updateYlabels(x)
        # Apply any edits relevant to children
        x <- applyEdits(x, x$edits)
    }
    x
}

editDetails.yaxis <- function(x, specs) {
  slot.names <- names(specs)
  if ("at" %in% slot.names) {
    if (is.null(x$at)) {
      x <- removeGrob(x, "major", warn=FALSE)
      x <- removeGrob(x, "ticks", warn=FALSE)
      x <- removeGrob(x, "labels", warn=FALSE)
    } else {
      x <- addGrob(x, make.yaxis.major(x$at, x$main))
      x <- addGrob(x, make.yaxis.ticks(x$at, x$main))
      x <- updateYlabels(x)
    }
  }
  if ("label" %in% slot.names) {
    if (!is.null(x$at))
      x <- updateYlabels(x)
  }
  if ("main" %in% slot.names)
    if (!is.null(x$at)) {
      x <- addGrob(x, make.yaxis.major(x$at, x$main))
      x <- addGrob(x, make.yaxis.ticks(x$at, x$main))
      x <- updateYlabels(x)
    }
  x
}

make.yaxis.major <- function(at, main) {
  if (main)
    x <- c(0, 0)
  else
    x <- c(1, 1)
  linesGrob(unit(x, "npc"), unit(c(min(at), max(at)), "native"),
            name="major")
}

make.yaxis.ticks <- function(at, main) {
  if (main) {
    tick.x0 <- unit(0, "npc")
    tick.x1 <- unit(-.5, "lines")
  }
  else {
    tick.x0 <- unit(1, "npc")
    tick.x1 <- unit(1, "npc") + unit(.5, "lines")
  }
  segmentsGrob(tick.x0, unit(at, "native"),
               tick.x1, unit(at, "native"),
               name="ticks")
}

make.yaxis.labels <- function(at, label, main) {
  if (main) {
    hjust <- "right"
    label.x <- unit(-1, "lines")
  }
  else {
    hjust <- "left"
    label.x <- unit(1, "npc") + unit(1, "lines")
  }
  just <- c(hjust, "centre")
  if (is.logical(label))
    labels <- as.character(at)
  else
    labels <- label
  textGrob(labels, label.x, unit(at, "native"),
           just=just, rot=0, check.overlap=TRUE, name="labels")
}

updateYlabels <- function(x) {
  if (is.logical(x$label) && !x$label)
    removeGrob(x, "labels", warn=FALSE)
  else
    addGrob(x, make.yaxis.labels(x$at, x$label, x$main))
}

yaxisGrob <- function(at=NULL, label=TRUE, main=TRUE,
                      edits=NULL,
                      name=NULL, gp=gpar(), vp=NULL) {
  grid.yaxis(at=at, label=label, main=main, edits=edits,
             name=name, gp=gp, draw=FALSE, vp=vp)
}

# The "main" y-axis is on the left when vp$origin is "*.left"
# and on the right when vp$origin is "*.right"
grid.yaxis <- function(at=NULL, label=TRUE, main=TRUE,
                       edits=NULL,
                       name=NULL, gp=gpar(),
                       draw=TRUE, vp=NULL) {
  if (is.null(at)) {
    # We do not have enough information to make the ticks and labels
    major <- NULL
    ticks <- NULL
    labels <- NULL
  } else {
    major <- make.yaxis.major(at, main)
    ticks <- make.yaxis.ticks(at, main)
    if (is.logical(label) && length(label) == 0)
	stop("logical 'label' supplied of length 0")
    if (is.logical(label) && !label)
      labels <- NULL
    else
      labels <- make.yaxis.labels(at, label, main)
  }
  yg <- applyEdits(gTree(at=at, label=label, main=main,
                         children=gList(major, ticks, labels),
                         edits=edits,
                         name=name, gp=gp, vp=vp,
                         cl=c("yaxis", "axis")),
                   edits)
  if (draw)
    grid.draw(yg)
  invisible(yg)
}

######################################
# Simple "side-effect" plotting functions
######################################

grid.grill <- function(h=unit(seq(0.25, 0.75, 0.25), "npc"),
                       v=unit(seq(0.25, 0.75, 0.25), "npc"),
                       default.units="npc",
                       gp=gpar(col="grey"), vp=NULL) {
  if (!is.unit(h))
    h <- unit(h, default.units)
  if (!is.unit(v))
    v <- unit(v, default.units)
  # FIXME:  Should replace for loop and call to grid.lines with call to grid.segments
  # once the latter exists
  if (!is.null(vp))
    pushViewport(vp)
  grid.segments(v, unit(0, "npc"), v, unit(1, "npc"), gp=gp)
  grid.segments(unit(0, "npc"), h, unit(1, "npc"), h, gp=gp)
  if (!is.null(vp))
    popViewport()
}


## Functions to calculate a set of points around the perimeter
## (or along the length) of a grob

## grobCoords() is a user-level function that emulates drawing
## set up behaviour (pushing viewports and setting graphical parameters)

## grobPoints() does not perform any set up and is for use
## within other 'grid' functions when set up has already been done,
## e.g., within resolveFill()

################################################################################
## Functions for creating coords data structures

validCoords <- function(x) {
    is.list(x) && length(x) > 0 &&
        is.numeric(x$x) && is.numeric(x$y) && length(x$x) == length(x$y)
        
}

validGrobCoords <- function(x) {
    is.list(x) && length(x) > 0 &&
        all(sapply(x, inherits, "GridCoords"))
}

validGTreeCoords <- function(x) {
    is.list(x) && length(x) > 0 &&
        all(sapply(x, inherits, "GridGrobCoords") |
            sapply(x, inherits, "GridGTreeCoords"))
}

coordPrintIndent <- "  "

## Public function for creating valid 'grid' points/coords result
gridCoords <- function(x, y) {
    coords <- list(x=as.numeric(x), y=as.numeric(y))
    if (validCoords(coords)) {
        class(coords) <- "GridCoords"
        coords
    } else
        stop("Invalid coordinates")
}

print.GridCoords <- function(x, indent="", ...) {
    if (length(x$x) > 3) {
        dots <- "..."
    } else {
        dots <- ""
    }
    cat(paste0(indent, "x:"),
        head(x$x, 3), dots, paste0("[", length(x$x), " values]\n"))
    cat(paste0(indent, "y:"),
        head(x$y, 3), dots, paste0("[", length(x$y), " values]\n"))    
}

gridGrobCoords <- function(x, name, rule = NULL) {
    if (validGrobCoords(x)) {
        class(x) <- "GridGrobCoords"
        attr(x, "name") <- name
        attr(x, "rule") <- rule
        x
    } else
        stop("Invalid grob coordinates")
}

print.GridGrobCoords <- function(x, indent="", ...) {
    if (is.null(names(x))) {
        names <- 1:length(x)
    } else {
        names <- names(x)
    }
    rule <- attr(x, "rule")
    if (is.null(rule)) {
        fillrule <- ""
    } else {
        fillrule <- paste0(" (fill: ", rule, ")")
    }
    cat(paste0(indent, "grob"), attr(x, "name"), fillrule, "\n")
    for (i in seq_along(x)) {
        cat(paste0(indent, coordPrintIndent, "shape"), names[i], "\n")
        print(x[[i]], indent=paste0(indent, coordPrintIndent, coordPrintIndent))
    }
}

gridGTreeCoords <- function(x, name) {
    if (validGTreeCoords(x)) {
        class(x) <- "GridGTreeCoords"
        attr(x, "name") <- name
        x
    } else {
        stop("Invalid gTree coordinates")
    }
}

print.GridGTreeCoords <- function(x, indent="", ...) {
    cat(paste0(indent, "gTree"), attr(x, "name"), "\n")
    for (i in seq_along(x)) {
        print(x[[i]], indent=paste0(indent, coordPrintIndent))
    }
}

################################################################################
## Calculate bounding box of coordinates

getX <- function(x, ...) {
    UseMethod("getX")
}

getY <- function(x, ...) {
    UseMethod("getY")
}

getX.GridCoords <- function(x, ...) {
    x$x
}

getY.GridCoords <- function(x, ...) {
    x$y
}

getX.GridGrobCoords <- function(x, subset = NULL, ...) {
    if (is.null(subset)) {
        x <- unlist(lapply(x, getX, ...))
    } else {
        x <- unlist(lapply(x[subset], getX, ...))
    }
}

getY.GridGrobCoords <- function(x, subset = NULL, ...) {
    if (is.null(subset)) {
        x <- unlist(lapply(x, getY, ...))
    } else {
        x <- unlist(lapply(x[subset], getY, ...))
    }
}

getX.GridGTreeCoords <- function(x, ...) {
    x <- unlist(lapply(x, getX, ...))
}

getY.GridGTreeCoords <- function(x, ...) {
    y <- unlist(lapply(x, getY, ...))
}

coordsBBox <- function(x, subset = NULL) {
    xx <- getX(x, subset)
    yy <- getY(x, subset)
    list(left = min(xx),
         bottom = min(yy),
         width = diff(range(xx)),
         height = diff(range(yy)))
}

################################################################################
## Support functions for calculating gTree coords

toDevice <- function(x) {
    if (isEmptyCoords(x)) return(x)
    UseMethod("toDevice")
}

toDevice.GridCoords <- function(x) {
    pts <- deviceLoc(unit(x$x, "in"), unit(x$y, "in"),
                        valueOnly=TRUE)
    gridCoords(pts$x, pts$y)
}

toDevice.GridGrobCoords <- function(x) {
    pts <- lapply(x, toDevice)
    gridGrobCoords(pts, attr(x, "name"), attr(x, "rule"))
}

toDevice.GridGTreeCoords <- function(x) {
    pts <- lapply(x, toDevice)
    gridGTreeCoords(pts, attr(x, "name"))
}

fromDevice <- function(x, trans) {
    UseMethod("fromDevice")
}

fromDevice.GridCoords <- function(x, trans) {
    ptsMatrix <- cbind(x$x, x$y, 1) %*% solve(trans)
    gridCoords(x=ptsMatrix[,1], y=ptsMatrix[,2])
}

fromDevice.GridGrobCoords <- function(x, trans) {
    pts <- lapply(x, fromDevice, trans)
    gridGrobCoords(pts, attr(x, "name"), attr(x, "rule"))
}

fromDevice.GridGTreeCoords <- function(x, trans) {
    pts <- lapply(x, fromDevice, trans)
    gridGTreeCoords(pts, attr(x, "name"))
}

################################################################################
## Empty coordinates

emptyCoords <- gridCoords(x = 0, y = 0)

emptyGrobCoords <- function(name) {
    gridGrobCoords(list("0"=emptyCoords), name)
}

emptyGTreeCoords <- function(name) {
    gridGTreeCoords(list(emptyGrobCoords("0")), name)
}

isEmptyCoords <- function(coords) {
    UseMethod("isEmptyCoords")
}

isEmptyCoords.GridCoords <- function(coords) {
    identical(coords, emptyCoords) 
}

isEmptyCoords.GridGrobCoords <- function(coords) {
    all(sapply(coords, identical, emptyCoords))
}

isEmptyCoords.GridGTreeCoords <- function(coords) {
    all(sapply(coords, isEmptyCoords))    
}

################################################################################
## Determine default 'closed' value
isClosed <- function(x, ...) {
    UseMethod("isClosed")
}

isClosedTRUE <- function(x, ...) {
    TRUE
}

isClosedFALSE <- function(x, ...) {
    FALSE
}

isClosed.default <- isClosedTRUE

isClosed.move.to <- isClosedFALSE
isClosed.line.to <- isClosedFALSE
isClosed.lines <- isClosedFALSE
isClosed.polyline <- isClosedFALSE
isClosed.segments <- isClosedFALSE
isClosed.beziergrob <- isClosedFALSE

isClosed.xspline <- function(x, ...) {
    if (x$open)
        FALSE
    else
        TRUE
}

isClosed.points <- function(x, ...) {
    switch(as.character(x$pch),
           "3"=, ## plus
           "4"=, ## times
           "8"=FALSE, ## plus-times
           TRUE)
}

################################################################################
## grobCoords()
##   Do drawing set up then calculate points

grobCoords <- function(x, closed, ...) {
    UseMethod("grobCoords")
}

grobCoords.grob <- function(x, closed=isClosed(x), ...) {
    vp <- x$vp
    trans <- current.transform()
    # Same set up as drawGrob()
    dlon <- grid.Call(C_setDLon, FALSE)
    on.exit(grid.Call(C_setDLon, dlon))
    tempgpar <- grid.Call(C_getGPar)
    on.exit(grid.Call(C_setGPar, tempgpar), add=TRUE)
    # Same drawing context set up as drawGrob()
    # including enforcing the drawing context
    x <- preDraw(x)
    # Same drawing content set up as drawGrob() ...
    x <- makeContent(x)
    ## Does this grob change the viewport ?
    ## (including has preDraw() changed the viewport)
    vpgrob <- !is.null(x$vp) || !identical(vp, x$vp)
    # BUT NO DRAWING
    ## Polygon outline in inches
    pts <- grobPoints(x, closed, ...)
    if (vpgrob && !isEmptyCoords(pts)) {
        ## Calc locations on device
        pts <- gridGrobCoords(lapply(pts, toDevice), x$name,
                              attr(pts, "rule"))
    }
    # Same context clean up as drawGrob()
    postDraw(x)
    if (vpgrob && !isEmptyCoords(pts)) {
        ## Transform back to locations
        pts <- gridGrobCoords(lapply(pts, fromDevice, trans), x$name,
                              attr(pts, "rule"))
    }
    pts
}

## "gTree"s
grobCoords.gList <- function(x, closed=isClosed(x), ...) {
    gridGTreeCoords(lapply(x, grobCoords, closed, ...),
                    grobName())
}

grobCoords.gTree <- function(x, closed=isClosed(x), ...) {
    vp <- x$vp
    trans <- current.transform()
    # Same set up as drawGTree()
    dlon <- grid.Call(C_setDLon, FALSE)
    on.exit(grid.Call(C_setDLon, dlon))
    tempgrob <- grid.Call(C_getCurrentGrob)
    tempgpar <- grid.Call(C_getGPar)
    on.exit({ grid.Call(C_setGPar, tempgpar)
              grid.Call(C_setCurrentGrob, tempgrob)
            }, add=TRUE)
    # Same drawing context set up as drawGTree(),
    # including enforcing the drawing context
    x <- preDraw(x)
    # Same drawing content set up as drawGTree() ...
    x <- makeContent(x)
    ## Does this grob change the viewport ?
    ## (including has preDraw() changed the viewport)
    vpgrob <- !is.null(x$vp) || !identical(vp, x$vp)
    ## Polygon outline in inches
    pts <- gridGTreeCoords(unname(lapply(x$children[x$childrenOrder],
                                         grobCoords, closed, ...)),
                           x$name)
    if (vpgrob && !isEmptyCoords(pts)) {
        ## Calc locations on device
        pts <- gridGTreeCoords(lapply(pts, toDevice), x$name)
    }
    # Same context clean up as drawGTree()
    postDraw(x)
    if (vpgrob && !isEmptyCoords(pts)) {
        ## Transform back to locations
        pts <- gridGTreeCoords(lapply(pts, fromDevice, trans), x$name)
    }
    pts
}

################################################################################
## grobPoints()
##   No drawing set up

grobPoints <- function(x, closed, ...) {
    UseMethod("grobPoints")
}

grobPoints.grob <- function(x, closed, ...) {
    emptyGrobCoords(x$name)
}

grobPoints.move.to <- function(x, closed, ...) {
    emptyGrobCoords(x$name)
}

grobPoints.line.to <- function(x, closed, ...) {
    emptyGrobCoords(x$name)
}

grobPoints.circle <- function(x, closed=TRUE, ..., n=100) {
    if (closed) {
        cx <- convertX(x$x, "in", valueOnly=TRUE)
        cy <- convertY(x$y, "in", valueOnly=TRUE)
        r <- pmin(convertWidth(x$r, "in", valueOnly=TRUE),
                  convertHeight(x$r, "in", valueOnly=TRUE))
        t <- seq(0, 2*pi, length.out=n+1)[-(n+1)]
        ## Recycle via cbind()
        circs <- cbind(cx, cy, r)
        n <- nrow(circs)
        pts <- lapply(1:n,
                      function(i) {
                          gridCoords(x=circs[i, 1] + circs[i, 3]*cos(t),
                                     y=circs[i, 2] + circs[i, 3]*sin(t))
                      })
        names(pts) <- 1:n
        gridGrobCoords(pts, x$name)
    } else {
        emptyGrobCoords(x$name)
    }
}

grobPoints.lines <- function(x, closed=FALSE, ..., n=100) {
    if (closed) {
        emptyGrobCoords(x$name)
    } else {
        xx <- convertX(x$x, "in", valueOnly=TRUE)
        yy <- convertY(x$y, "in", valueOnly=TRUE)
        ## Recycle via cbind()
        lines <- cbind(xx, yy)
        gridGrobCoords(list("1"=gridCoords(x=lines[,1], y=lines[,2])), x$name)
    }
}

grobPoints.polyline <- function(x, closed=FALSE, ...) {
    if (closed) {
        emptyGrobCoords(x$name)
    } else {
        ## polylineGrob() ensures that x/y same length
        xx <- convertX(x$x, "in", valueOnly=TRUE)
        yy <- convertY(x$y, "in", valueOnly=TRUE)
        pts <- list(x=xx, y=yy)
        if (is.null(x$id) && is.null(x$id.lengths)) {
            gridGrobCoords(list("1"=do.call(gridCoords, pts)), x$name)
        } else {
            if (is.null(x$id)) {
                n <- length(x$id.lengths)
                id <- rep(1L:n, x$id.lengths)
            } else {
                n <- length(unique(x$id))
                id <- x$id
            }
            if (n > 1) {
                gridGrobCoords(lapply(split(as.data.frame(pts), id),
                                  function(z) do.call(gridCoords, z)),
                               x$name)
            } else {
                gridGrobCoords(list("1"=do.call(gridCoords, pts)), x$name)
            }
        }
    }    
}

## NOTE that grid.polygon() does not provide ability to set fill rule
## (and neither does dev->polygon());  some devices allow a global
## device fill rule (!), e.g., pdf(), postscript(), windows(),
## but we can't do anything about that.
## If you want proper control, use grid.path() instead
## (which does have a fill rule arg).
grobPoints.polygon <- function(x, closed=TRUE, ...) {
    if (closed) {
        ## polygonGrob() ensures that x/y same length
        xx <- convertX(x$x, "in", valueOnly=TRUE)
        yy <- convertY(x$y, "in", valueOnly=TRUE)
        pts <- list(x=xx, y=yy)
        if (is.null(x$id) && is.null(x$id.lengths)) {
            gridGrobCoords(list("1"=do.call(gridCoords, pts)), x$name)
        } else {
            if (is.null(x$id)) {
                n <- length(x$id.lengths)
                id <- rep(1L:n, x$id.lengths)
            } else {
                n <- length(unique(x$id))
                id <- x$id
            }
            if (n > 1) {
                gridGrobCoords(lapply(split(as.data.frame(pts), id),
                                  function(z) do.call(gridCoords, z)),
                               x$name)
            } else {
                gridGrobCoords(list("1"=do.call(gridCoords, pts)), x$name)
            }
        }
    } else {
        emptyGrobCoords(x$name)
    }
}

xyListFromMatrix <- function(m, xcol, ycol) {
    n <- nrow(m)
    lapply(1:n,
           function(i) {
               gridCoords(x=m[i, xcol], y=m[i, ycol])
           })
}

grobPoints.pathgrob <- function(x, closed=TRUE, ...) {
    if (closed) {
        ## pathGrob() ensures that x/y same length
        xx <- convertX(x$x, "in", valueOnly=TRUE)
        yy <- convertY(x$y, "in", valueOnly=TRUE)
        pts <- list(x=xx, y=yy)
        hasMultiple <- !(is.null(x$pathId) && is.null(x$pathId.lengths))
        if (hasMultiple) {
            if (is.null(x$pathId)) {
                n <- length(x$pathId.lengths)
                pathId <- rep(1L:n, x$pathId.lengths)
            } else {
                pathId <- x$pathId
            }
        }
        if (is.null(x$id) && is.null(x$id.lengths)) {
            if (hasMultiple) {
                gridGrobCoords(lapply(split(as.data.frame(pts), pathId),
                                  function(z) do.call(gridCoords, z)),
                               x$name, x$rule)
            } else {
                gridGrobCoords(list("1"=do.call(gridCoords, pts)),
                               x$name, x$rule)
            }
        } else {
            if (is.null(x$id)) {
                n <- length(x$id.lengths)
                id <- rep(1L:n, x$id.lengths)
            } else {
                n <- length(unique(x$id))
                id <- x$id
            }
            if (hasMultiple) {
                pts <- unlist(mapply(split,
                                     split(as.data.frame(pts), pathId),
                                     split(id, pathId),
                                     SIMPLIFY=FALSE),
                              recursive=FALSE)
                names(pts) <- gsub("[.][0-9]+$", "", names(pts))
                gridGrobCoords(lapply(pts,
                                  function(z) do.call(gridCoords, z)),
                               x$name, x$rule)
            } else {
                pts <- split(as.data.frame(pts), id)
                names(pts) <- rep(1, length(pts))
                gridGrobCoords(lapply(pts,
                                  function(z) do.call(gridCoords, z)),
                               x$name, x$rule)
            }
        }
    } else {
        emptyGrobCoords(x$name)
    }
}

grobPoints.rect <- function(x, closed=TRUE, ...) {
    if (closed) {
        hjust <- resolveHJust(x$just, x$hjust)
        vjust <- resolveVJust(x$just, x$vjust)
        w <- convertWidth(x$width, "in", valueOnly=TRUE)
        h <- convertHeight(x$height, "in", valueOnly=TRUE)
        left <- convertX(x$x, "in", valueOnly=TRUE) - hjust*w
        bottom <- convertY(x$y, "in", valueOnly=TRUE) - vjust*h
        right <- left + w
        top <- bottom + h
        ## Recycle via cbind()
        rects <- cbind(left, right, bottom, top)
        pts <- xyListFromMatrix(rects, c(1, 1, 2, 2), c(3, 4, 4, 3))
        names(pts) <- seq_along(pts)
        gridGrobCoords(pts, x$name)
    } else {
        emptyGrobCoords(x$name)
    }
}

grobPoints.segments <- function(x, closed=FALSE, ...) {
    if (closed) {
        emptyGrobCoords(x$name)
    } else {
        x0 <- convertX(x$x0, "in", valueOnly=TRUE)
        x1 <- convertX(x$x1, "in", valueOnly=TRUE)
        y0 <- convertY(x$y0, "in", valueOnly=TRUE)
        y1 <- convertY(x$y1, "in", valueOnly=TRUE)
        ## Recycle via cbind()        
        xy <- cbind(x0, x1, y0, y1)
        pts <- xyListFromMatrix(xy, 1:2, 3:4)
        names(pts) <- seq_along(pts)
        gridGrobCoords(pts, x$name)
    }
}

grobPoints.xspline <- function(x, closed=!x$open, ...) {
    if ((closed && !x$open) ||
        (!closed && x$open)) {
        ## xsplinePoints() takes care of multiple X-splines
        trace <- xsplinePoints(x)
        if ("x" %in% names(trace)) {
            ## Single X-spline
            gridGrobCoords(list("1"=gridCoords(x=as.numeric(trace$x),
                                           y=as.numeric(trace$y))),
                           x$name)
        } else {
            pts <- lapply(trace,
                          function(t) {
                              gridCoords(x=as.numeric(t$x), y=as.numeric(t$y))
                          })
            names(pts) <- seq_along(pts)
            gridGrobCoords(pts, x$name)
        }
    } else {
        emptyGrobCoords(x$name)
    }
}

## beziergrob covered by splinegrob (via makeContent)

## Just return a bounding box for the text (if closed=TRUE)
grobPoints.text <- function(x, closed=TRUE, ...) {
    if (closed) {
        bounds <- grid.Call(C_textBounds, as.graphicsAnnot(x$label),
                            x$x, x$y,
                            resolveHJust(x$just, x$hjust),
                            resolveVJust(x$just, x$vjust),
                            x$rot, as.numeric(NA))
        if (is.null(bounds))
            emptyGrobCoords(x$name)
        else {
            left <- bounds[5]
            bottom <- bounds[6]
            right <- left + bounds[3]
            top <- bottom + bounds[4]
            gridGrobCoords(list("1"=gridCoords(x=c(left, left, right, right),
                                           y=c(bottom, top, top, bottom))),
                           x$name)
        }
    } else {
        emptyGrobCoords(x$name)
    }
}

grobPoints.points <- function(x, closed=TRUE, ...) {
    closed <- as.logical(closed)
    if (is.na(closed)) 
        stop("Closed must not be a missing value")
    pts <- grid.Call(C_pointsPoints, x$x, x$y, x$pch, x$size, closed)
    if (is.null(pts) ||
        all(sapply(pts, is.null))) {
        emptyGrobCoords(x$name)
    } else {
        names <- attr(pts, "coordNames")
        pts <- lapply(pts,
                      function(x) {
                          if (is.null(x))
                              emptyCoords
                          else {
                              names(x) <- c("x", "y")
                              do.call(gridCoords, x)
                          }
                      })
        names(pts) <- names
        gridGrobCoords(pts, x$name)
    }
}

## Do not treat these as open or closed shapes (for now at least)
grobPoints.rastergrob <- function(x, closed, ...) {
    emptyGrobCoords(x$name)
}

grobPoints.clip <- function(x, closed, ...) {
    emptyGrobCoords(x$name)
}

grobPoints.null <- function(x, closed, ...) {
    emptyGrobCoords(x$name)
}

## Collections of grobs

## NOTE that these generate coordinates from their children
## and they must call grobCoords() rather than grobPoints()
## on those children so that the children can perform any
## relevant set up

grobPoints.gList <- function(x, closed=TRUE, ...) {
    if (length(x) > 0) {
        gridGTreeCoords(lapply(x, grobCoords, closed, ...), grobName())
    } else {
        emptyGTreeCoords(grobName())
    }
}

grobPoints.gTree <- function(x, closed=TRUE, ...) {
    if (length(x$children) > 0) {
        pts <- lapply(x$children[x$childrenOrder], grobCoords, closed, ...)
        gridGTreeCoords(unname(pts), x$name)
    } else {
        emptyGTreeCoords(x$name)
    }
}

#  File src/library/grid/R/curve.R
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


###############################
# CURVE primitive
###############################

calcOrigin <- function(x1, y1, x2, y2, origin, hand) {
    # Positive origin means origin to the "right"
    # Negative origin means origin to the "left"
    xm <- (x1 + x2)/2
    ym <- (y1 + y2)/2
    dx <- x2 - x1
    dy <- y2 - y1
    slope <- dy/dx
    oslope <- -1/slope
    # The origin is a point somewhere along the line between
    # the end points, rotated by 90 (or -90) degrees
    # Two special cases:
    # If slope is non-finite then the end points lie on a vertical line, so
    # the origin lies along a horizontal line (oslope = 0)
    # If oslope is non-finite then the end points lie on a horizontal line,
    # so the origin lies along a vertical line (oslope = Inf)
    tmpox <- ifelse(!is.finite(slope),
                    xm,
                    ifelse(!is.finite(oslope),
                           xm + origin*(x2 - x1)/2,
                           xm + origin*(x2 - x1)/2))
    tmpoy <- ifelse(!is.finite(slope),
                    ym + origin*(y2 - y1)/2,
                    ifelse(!is.finite(oslope),
                           ym,
                           ym + origin*(y2 - y1)/2))
    # ALWAYS rotate by -90 about midpoint between end points
    # Actually no need for "hand" because "origin" also
    # encodes direction
    # sintheta <- switch(hand, left=-1, right=1)
    sintheta <- -1
    ox <- xm - (tmpoy - ym)*sintheta
    oy <- ym + (tmpox - xm)*sintheta

    list(x=ox, y=oy)
}

# Given ncp*ncurve vector of values, ncurve vector of start values,
# ncurve vector of end values, ncurve vector of end logicals,
# combine start or end values with original values based on logicals
interleave <- function(ncp, ncurve, val, sval, eval, e) {
    sval <- rep(sval, length.out=ncurve)
    eval <- rep(eval, length.out=ncurve)
    result <- matrix(NA, ncol=ncurve, nrow=ncp+1)
    m <- matrix(val, ncol=ncurve)
    for (i in 1L:ncurve) {
        if (e[i])
            result[,i] <- c(m[,i], eval[i])
        else
            result[,i] <- c(sval[i], m[,i])
    }
    as.numeric(result)
}

# Calculate a "square" set of end points to calculate control points from
# NOTE: end points may be vector
calcSquareControlPoints <- function(x1, y1, x2, y2,
                                    curvature, angle, ncp,
                                    debug=FALSE) {
    dx <- x2 - x1
    dy <- y2 - y1
    slope <- dy/dx

    # FIXME:  There MUST be a more compact way of calculating the
    # new end point!
    end <- (slope > 1 |
            (slope < 0 & slope > -1))
    if (curvature < 0)
        end <- !end
    startx <- ifelse(end,
                     x1,
                     ifelse(abs(slope) > 1, x2 - dx, x2 - sign(slope)*dy))
    starty <- ifelse(end,
                     y1,
                     ifelse(abs(slope) > 1, y2 - sign(slope)*dx, y2 - dy))
    endx <- ifelse(end,
                   ifelse(abs(slope) > 1, x1 + dx, x1 + sign(slope)*dy),
                   x2)
    endy <- ifelse(end,
                   ifelse(abs(slope) > 1, y1 + sign(slope)*dx, y1 + dy),
                   y2)

    cps <- calcControlPoints(startx, starty, endx, endy,
                             curvature, angle, ncp,
                             debug)

    # Intereave control points and extra "square" control points
    ncurve <- length(x1)
    cps$x <- interleave(ncp, ncurve, cps$x, startx, endx, end)
    cps$y <- interleave(ncp, ncurve, cps$y, starty, endy, end)

    list(x=cps$x, y=cps$y, end=end)
}

# Find origin of rotation
# Rotate around that origin
calcControlPoints <- function(x1, y1, x2, y2, curvature, angle, ncp,
                              debug=FALSE) {
    # Negative curvature means curve to the left
    # Positive curvature means curve to the right
    # Special case curvature = 0 (straight line) has been handled
    xm <- (x1 + x2)/2
    ym <- (y1 + y2)/2
    dx <- x2 - x1
    dy <- y2 - y1
    slope <- dy/dx

    # Calculate "corner" of region to produce control points in
    # (depends on 'angle', which MUST lie between 0 and 180)
    # Find by rotating start point by angle around mid point
    if (is.null(angle)) {
        # Calculate angle automatically
        angle <- ifelse(slope < 0,
                        2*atan(abs(slope)),
                        2*atan(1/slope))
    } else {
        angle <- angle/180*pi
    }
    sina <- sin(angle)
    cosa <- cos(angle)
    # FIXME:  special case of vertical or horizontal line ?
    cornerx <- xm + (x1 - xm)*cosa - (y1 - ym)*sina
    cornery <- ym + (y1 - ym)*cosa + (x1 - xm)*sina

    # Debugging
    if (debug) {
        grid.points(cornerx, cornery, default.units="inches",
                    pch=16, size=unit(3, "mm"),
                    gp=gpar(col="grey"))
    }

    # Calculate angle to rotate region by to align it with x/y axes
    beta <- -atan((cornery - y1)/(cornerx - x1))
    sinb <- sin(beta)
    cosb <- cos(beta)
    # Rotate end point about start point to align region with x/y axes
    newx2 <- x1 + dx*cosb - dy*sinb
    newy2 <- y1 + dy*cosb + dx*sinb

    # Calculate x-scale factor to make region "square"
    # FIXME:  special case of vertical or horizontal line ?
    scalex <- (newy2 - y1)/(newx2 - x1)
    # Scale end points to make region "square"
    newx1 <- x1*scalex
    newx2 <- newx2*scalex

    # Calculate the origin in the "square" region
    # (for rotating start point to produce control points)
    # (depends on 'curvature')
    # 'origin' calculated from 'curvature'
    ratio <- 2*(sin(atan(curvature))^2)
    origin <- curvature - curvature/ratio
    # 'hand' also calculated from 'curvature'
    if (curvature > 0)
        hand <- "right"
    else
        hand <- "left"
    oxy <- calcOrigin(newx1, y1, newx2, newy2, origin, hand)
    ox <- oxy$x
    oy <- oxy$y

    # Calculate control points
    # Direction of rotation depends on 'hand'
    dir <- switch(hand,
                  left=-1,
                  right=1)
    # Angle of rotation depends on location of origin
    maxtheta <- pi + sign(origin*dir)*2*atan(abs(origin))
    theta <- seq(0, dir*maxtheta,
                 dir*maxtheta/(ncp + 1))[c(-1, -(ncp + 2))]
    costheta <- cos(theta)
    sintheta <- sin(theta)
    # May have BOTH multiple end points AND multiple
    # control points to generate (per set of end points)
    # Generate consecutive sets of control points by performing
    # matrix multiplication
    cpx <- ox + ((newx1 - ox) %*% t(costheta)) -
        ((y1 - oy) %*% t(sintheta))
    cpy <- oy + ((y1 - oy) %*% t(costheta)) +
        ((newx1 - ox) %*% t(sintheta))

    # Reverse transformations (scaling and rotation) to
    # produce control points in the original space
    cpx <- cpx/scalex
    sinnb <- sin(-beta)
    cosnb <- cos(-beta)
    finalcpx <- x1 + (cpx - x1)*cosnb - (cpy - y1)*sinnb
    finalcpy <- y1 + (cpy - y1)*cosnb + (cpx - x1)*sinnb

    # Debugging
    if (debug) {
        ox <- ox/scalex
        fox <- x1 + (ox - x1)*cosnb - (oy - y1)*sinnb
        foy <- y1 + (oy - y1)*cosnb + (ox - x1)*sinnb
        grid.points(fox, foy, default.units="inches",
                    pch=16, size=unit(1, "mm"),
                    gp=gpar(col="grey"))
        grid.circle(fox, foy, sqrt((ox - x1)^2 + (oy - y1)^2),
                    default.units="inches",
                    gp=gpar(col="grey"))
    }

    list(x=as.numeric(t(finalcpx)), y=as.numeric(t(finalcpy)))
}

# Debugging
cbDiagram <- function(x1, y1, x2, y2, cps) {
    grid.segments(x1, y1, x2, y2,
                gp=gpar(col="grey"),
                default.units="inches")
    grid.points(x1, y1, pch=16, size=unit(1, "mm"),
                gp=gpar(col="green"),
                default.units="inches")
    grid.points(x2, y2, pch=16, size=unit(1, "mm"),
                gp=gpar(col="red"),
                default.units="inches")
    grid.points(cps$x, cps$y, pch=16, size=unit(1, "mm"),
                default.units="inches",
                gp=gpar(col="blue"))
}

straightCurve <- function(x1, y1, x2, y2, arrow, debug) {
    if (debug) {
        xm <- (x1 + x2)/2
        ym <- (y1 + y2)/2
        cbDiagram(x1, y1, x2, y2, list(x=xm, y=ym))
    }

    segmentsGrob(x1, y1, x2, y2,
                 default.units="inches",
                 arrow=arrow, name="segment")
}

# Return a gTree (even if it only has one grob as a child)
# because that is the only way to get more than one child
# to draw
calcCurveGrob <- function(x, debug) {
    x1 <- x$x1
    x2 <- x$x2
    y1 <- x$y1
    y2 <- x$y2
    curvature <- x$curvature
    angle <- x$angle
    ncp <- x$ncp
    shape <- x$shape
    square <- x$square
    squareShape <- x$squareShape
    inflect <- x$inflect
    arrow <- x$arrow
    open <- x$open

    # Calculate a set of control points based on:
    # 'curvature', ' angle', and 'ncp',
    # and the start and end point locations.

    # The origin is a point along the perpendicular bisector
    # of the line between the end points.

    # The control points are found by rotating the end points
    # about the origin.

    # Do everything in inches to make things easier.
    # Because this is within a makeContent() method,
    # the conversions will not be an
    # issue (in terms of device resizes).
    x1 <- convertX(x1, "inches", valueOnly=TRUE)
    y1 <- convertY(y1, "inches", valueOnly=TRUE)
    x2 <- convertX(x2, "inches", valueOnly=TRUE)
    y2 <- convertY(y2, "inches", valueOnly=TRUE)

    # Outlaw identical end points
    if (any(x1 == x2 & y1 == y2))
        stop("end points must not be identical")

    # Rep locations to allow multiple curves from single call
    maxn <- max(length(x1),
                length(y1),
                length(x2),
                length(y2))
    x1 <- rep(x1, length.out=maxn)
    y1 <- rep(y1, length.out=maxn)
    x2 <- rep(x2, length.out=maxn)
    y2 <- rep(y2, length.out=maxn)
    if (!is.null(arrow))
        arrow <- rep(arrow, length.out=maxn)

    if (curvature == 0) {
        children <- gList(straightCurve(x1, y1, x2, y2, arrow, debug))
    } else {
        # Treat any angle less than 1 or greater than 179 degrees
        # as a straight line
        # Takes care of some nasty limit effects as well as simplifying
        # things
        if (angle < 1 || angle > 179) {
            children <- gList(straightCurve(x1, y1, x2, y2, arrow, debug))
        } else {
            # Handle 'square' vertical and horizontal lines
            # separately
            if (square && any(x1 == x2 | y1 == y2)) {
                subset <- x1 == x2 | y1 == y2
                straightGrob <- straightCurve(x1[subset], y1[subset],
                                               x2[subset], y2[subset],
                                               arrow, debug)
                # Remove these from the curves to draw
                x1 <- x1[!subset]
                x2 <- x2[!subset]
                y1 <- y1[!subset]
                y2 <- y2[!subset]
                if (!is.null(arrow))
                    arrow <- arrow[!subset]
            } else {
                straightGrob <- NULL
            }
            ncurve <- length(x1)
            # If nothing to draw, we're done
            if (ncurve == 0) {
                children <- gList(straightGrob)
            } else {
                if (inflect) {
                    xm <- (x1 + x2)/2
                    ym <- (y1 + y2)/2
                    shape1 <- rep(rep(shape, length.out=ncp), ncurve)
                    shape2 <- rev(shape1)
                    if (square) {
                      # If 'square' then add an extra control point
                        cps1 <- calcSquareControlPoints(x1, y1, xm, ym,
                                                        curvature, angle,
                                                        ncp,
                                                        debug=debug)
                        cps2 <- calcSquareControlPoints(xm, ym, x2, y2,
                                                        -curvature, angle,
                                                        ncp,
                                                        debug=debug)
                        shape1 <- interleave(ncp, ncurve, shape1,
                                             squareShape, squareShape,
                                             cps1$end)
                        shape2 <- interleave(ncp, ncurve, shape2,
                                             squareShape, squareShape,
                                             cps2$end)
                        ncp <- ncp + 1
                    } else {
                        cps1 <- calcControlPoints(x1, y1, xm, ym,
                                                  curvature, angle, ncp,
                                                  debug=debug)
                        cps2 <- calcControlPoints(xm, ym, x2, y2,
                                                  -curvature, angle, ncp,
                                                  debug=debug)
                    }

                    if (debug) {
                        cbDiagram(x1, y1, xm, ym, cps1)
                        cbDiagram(xm, ym, x2, y2, cps2)
                    }

                    idset <- 1L:ncurve
                    splineGrob <-
                        xsplineGrob(c(x1, cps1$x, xm, cps2$x, x2),
                                    c(y1, cps1$y, ym, cps2$y, y2),
                                    id=c(idset, rep(idset, each=ncp),
                                      idset, rep(idset, each=ncp),
                                      idset),
                                    default.units="inches",
                                    shape=c(rep(0, ncurve), shape1,
                                      rep(0, ncurve), shape2,
                                      rep(0, ncurve)),
                                    arrow=arrow, open=open,
                                    name="xspline")
                    if (is.null(straightGrob)) {
                        children <- gList(splineGrob)
                    } else {
                        children <- gList(straightGrob, splineGrob)
                    }
                } else {
                    shape <- rep(rep(shape, length.out=ncp), ncurve)
                    if (square) {
                      # If 'square' then add an extra control point
                        cps <- calcSquareControlPoints(x1, y1, x2, y2,
                                                       curvature, angle,
                                                       ncp,
                                                       debug=debug)
                        shape <- interleave(ncp, ncurve, shape,
                                            squareShape, squareShape,
                                            cps$end)
                        ncp <- ncp + 1
                    } else {
                        cps <- calcControlPoints(x1, y1, x2, y2,
                                                 curvature, angle, ncp,
                                                 debug=debug)
                    }
                    if (debug) {
                        cbDiagram(x1, y1, x2, y2, cps)
                    }

                    idset <- 1L:ncurve
                    splineGrob <- xsplineGrob(c(x1, cps$x, x2),
                                              c(y1, cps$y, y2),
                                              id=c(idset,
                                                rep(idset, each=ncp), idset),
                                              default.units="inches",
                                              shape=c(rep(0, ncurve), shape,
                                                rep(0, ncurve)),
                                              arrow=arrow, open=open,
                                              name="xspline")
                    if (is.null(straightGrob)) {
                        children <- gList(splineGrob)
                    } else {
                        children <- gList(straightGrob, splineGrob)
                    }
                }
            }
        }
    }
    gTree(children=children,
          name=x$name, gp=x$gp, vp=x$vp)
}

validDetails.curve <- function(x) {
    if ((!is.unit(x$x1) || !is.unit(x$y1)) ||
        (!is.unit(x$x2) || !is.unit(x$y2)))
        stop("'x1', 'y1', 'x2', and 'y2' must be units")
    x$curvature <- as.numeric(x$curvature)
    x$angle <- x$angle %% 180
    x$ncp <- as.integer(x$ncp)
    if (x$shape < -1 || x$shape > 1)
        stop("'shape' must be between -1 and 1")
    x$square <- as.logical(x$square)
    if (x$squareShape < -1 || x$squareShape > 1)
        stop("'squareShape' must be between -1 and 1")
    x$inflect <- as.logical(x$inflect)
    if (!is.null(x$arrow) && !inherits(x$arrow, "arrow"))
        stop("'arrow' must be an arrow object or NULL")
    x$open <- as.logical(x$open)
    x
}

makeContent.curve <- function(x) {
    calcCurveGrob(x, x$debug)
}

xDetails.curve <- function(x, theta) {
    cg <- calcCurveGrob(x, FALSE)
    # Could do better here
    # (result for more than 1 child is basically to give up)
    if (length(cg$children) == 1)
        xDetails(cg$children[[1]], theta)
    else
        xDetails(cg, theta)
}

yDetails.curve <- function(x, theta) {
    cg <- calcCurveGrob(x, FALSE)
    if (length(cg$children) == 1)
        yDetails(cg$children[[1]], theta)
    else
        yDetails(cg, theta)
}

widthDetails.curve <- function(x) {
    cg <- calcCurveGrob(x, FALSE)
    if (length(cg$children) == 1)
        widthDetails(cg$children[[1]])
    else
        widthDetails(cg)
}

heightDetails.curve <- function(x) {
    cg <- calcCurveGrob(x, FALSE)
    if (length(cg$children) == 1)
        heightDetails(cg$children[[1]])
    else
        heightDetails(cg)
}

curveGrob <- function(x1, y1, x2, y2, default.units="npc",
                      curvature=1, angle=90, ncp=1,
                      shape=0.5, square=TRUE, squareShape=1,
                      inflect=FALSE, arrow=NULL, open=TRUE,
                      debug=FALSE,
                      name=NULL, gp=gpar(), vp=NULL) {
    # FIXME:  add arg checking
    # FIXME:  angle MUST be between 0 and 180
    if (!is.unit(x1))
        x1 <- unit(x1, default.units)
    if (!is.unit(y1))
        y1 <- unit(y1, default.units)
    if (!is.unit(x2))
        x2 <- unit(x2, default.units)
    if (!is.unit(y2))
        y2 <- unit(y2, default.units)
    gTree(x1=x1, y1=y1, x2=x2, y2=y2,
          curvature=curvature, angle=angle, ncp=ncp,
          shape=shape, square=square, squareShape=squareShape,
          inflect=inflect, arrow=arrow, open=open, debug=debug,
          name=name, gp=gp, vp=vp,
          cl="curve")
}

grid.curve <- function(...) {
    grid.draw(curveGrob(...))
}

# Calculate the curvature to use if you want to produce control
# points lying along the arc of a circle that spans theta degrees
# (Use ncp=8 and shape=-1 to actually produce such an arc)
arcCurvature <- function(theta) {
    # Avoid limiting cases (just draw a straight line)
    if (theta < 1 || theta > 359)
        return(0)
    angle <- 0.5*theta/180*pi
    1/sin(angle) - 1/tan(angle)
}

#  File src/library/grid/R/debug.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2016 The R Core Team
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

### Label grobs in a scene


labelGrob <- function(grob, recurse, curdepth, depth, labelfun, ...) {
    UseMethod("labelGrob")
}

# The default grob label needs to do some calculations
# on sizes so need a drawDetails method to get the
# calculations right
drawDetails.groblabel <- function(x, ...) {
    gw <- convertWidth(grobWidth(x$grob), "inches", valueOnly=TRUE)
    gh <- convertHeight(grobHeight(x$grob), "inches", valueOnly=TRUE)
    grid.rect(grobX(x$grob, "west"), grobY(x$grob, "south"),
              unit(gw, "inches"), unit(gh, "inches"),
              just=c("left", "bottom"), gp=x$gp)
    tw <- convertWidth(stringWidth(x$grob$name), "inches", valueOnly=TRUE)
    th <- convertHeight(stringHeight(x$grob$name), "inches", valueOnly=TRUE)
    eps <- .01
    # If grob is REALLY short, draw horiz at normal cex
    if (gh < eps) {
        rot <- 0
        cex <- 1
    # If grob is REALLY thin, draw vertical at normal cex
    } else if (gw < eps) {
        rot <- 90
        cex <- 1
    } else {
        gratio <- gh/gw
        if (gratio > 1 && tw > gw) {
            rot <- 90
            wratio <- th/gw
            hratio <- tw/gh
        } else {
            rot <- 0
            wratio <- tw/gw
            hratio <- th/gh
        }
        if (wratio > 1 || hratio > 1) {
            cex <- 1/max(wratio, hratio)
        } else {
            cex <- 1
        }
    }
    if (is.null(x$gp)) {
        x$gp <- gpar(cex=cex)
    } else {
        if (is.null(x$gp$cex))
            x$gp$cex <- cex
    }
    if (is.null(x$otherArgs$rot))
        x$otherArgs$rot <- rot
    do.call("grid.text", c(list(label=x$grob$name,
                                x=grobX(x$grob, "north"),
                                y=grobY(x$grob, "west"),
                                gp=x$gp),
                           x$otherArgs))
}

grobLabel <- function(grob,
                      gp=gpar(col=rgb(1, 0, 0, .5),
                        fill=rgb(1, 0, 0, .2)),
                      ...) {
    grob(grob=grob, gp=gp, otherArgs=list(...),
         cl="groblabel")
}

labelGrob.grob <- function(grob, recurse, curdepth, depth, labelfun, ...) {
    if (is.null(depth) || curdepth %in% depth) {
        gTree(children=gList(grob,
                labelfun(grob, ...)),
              # Name new gTree same as old grob so that
              # setGrob() approach works below
              # (when 'gPath' is specified)
              name=grob$name)
    } else {
        grob
    }
}

labelGrob.gTree <- function(grob, recurse, curdepth, depth, labelfun, ...) {
    if (recurse) {
        newChildren <- do.call("gList",
                               lapply(grob$children,
                                      labelGrob,
                                      recurse, curdepth + 1, depth,
                                      labelfun, ...))
        grob <- setChildren(grob, newChildren)
    }
    if (is.null(depth) || curdepth %in% depth) {
        gTree(children=gList(grob,
                labelfun(grob, ...)),
              name=grob$name)
    } else {
        grob
    }
}

showGrob <- function(x=NULL,
                     gPath=NULL, strict=FALSE, grep=FALSE,
                     recurse=TRUE, depth=NULL,
                     labelfun=grobLabel, ...) {
    if (is.null(x)) {
        # Label all or part of current scene
        # The grid display list is NOT affected
        # To remove labels use grid.redraw()
        if (is.null(gPath)) {
            # Show the current scene
            dl <- grid.Call(C_getDisplayList)[1L : grid.Call(C_getDLindex)]
            grid.newpage(recording=FALSE)
            # -1 because first element on DL is ROOT viewport
            lapply(dl[-1],
                   function(y) {
                       # Modify the grob to add a label
                       if (is.grob(y))
                           y <- labelGrob(y, recurse, 1, depth, labelfun, ...)
                       # Draw either the original object or the modified grob
                       grid.draw(y, recording=FALSE)
                   })
        } else {
            # Only label the bit of the current scene specified by gPath
            grobToLabel <- grid.get(gPath, strict=strict, grep=grep)
            # NOTE: have to 'wrap' because otherwise the grobs in the
            # captured scene have been altered
            scene <- grid.grab(wrap=TRUE)
            modScene <- setGrob(scene, gPath,
                                labelGrob(grobToLabel, recurse, 1, depth,
                                          labelfun, ...),
                                strict=strict, grep=grep)
            grid.newpage(recording=FALSE)
            grid.draw(modScene, recording=FALSE)
        }
    } else {
        # Assume grob is not current scene so start a new page
        grid.newpage()
        grid.draw(x)
        showGrob(NULL, gPath, strict, grep, recurse, depth, labelfun, ...)
    }
    invisible()
}

#############
# Labelling viewports in a scene
#############

# FIXME:  some of this code for vpLists and vpStacks and vpTrees
# assumes that the components of a vpList or vpStack or the
# vpTree parent can ONLY be a viewport (when in fact they can
# also be a vpList, vpStack, or vpTree!)

# Label a viewport
# Get physical aspect ratio of vp to determine whether to rotate
# Shrink text to fit in vp
# (Assumes that we are currently occupying 'vp'
#  so that conversions are correct)
labelVP <- function(vp, col) {
    vw <- convertWidth(unit(1, "npc"), "inches", valueOnly=TRUE)
    vh <- convertHeight(unit(1, "npc"), "inches", valueOnly=TRUE)
    tw <- convertWidth(stringWidth(vp$name), "inches", valueOnly=TRUE)
    th <- convertHeight(stringHeight(vp$name), "inches", valueOnly=TRUE)
    eps <- .01
    # If viewport is REALLY short, draw horiz at normal cex
    if (vh < eps) {
        rot <- 0
        cex <- 1
    # If viewport is REALLY thin, draw vertical at normal cex
    } else if (vw < eps) {
        rot <- 90
        cex <- 1
    } else {
        vratio <- vh/vw
        if (vratio > 1 && tw > vw) {
            rot <- 90
            wratio <- th/vw
            hratio <- tw/vh
        } else {
            rot <- 0
            wratio <- tw/vw
            hratio <- th/vh
        }
        if (wratio > 1 || hratio > 1) {
            cex <- 1/max(wratio, hratio)
        } else {
            cex <- 1
        }
    }
    # Violate any clipping that is in effect
    pushViewport(viewport(clip="off"))
    grid.text(vp$name, rot=rot, gp=gpar(col=col, cex=cex))
    upViewport()
}

# Draw a "viewport"
drawVP <- function(vp, curDepth, depth, col, fill, label) {
    UseMethod("drawVP")
}

drawVP.viewport <- function(vp, curDepth, depth, col, fill, label) {
    if (vp$name != "ROOT" &&
        (is.null(depth) || curDepth %in% depth)) {
        pushViewport(vp)
        colIndex <- (curDepth - 1) %% length(col) + 1
        fillIndex <- (curDepth - 1) %% length(fill) + 1
        grid.rect(gp=gpar(col=col[colIndex], fill=fill[fillIndex]))
        if (label)
            labelVP(vp, col[colIndex])
        upViewport()
    }
}

drawVP.vpPath <- function(vp, curDepth, depth, col, fill, label) {
    if (is.null(depth) || curDepth %in% depth) {
        downViewport(vp)
        colIndex <- (curDepth - 1) %% length(col) + 1
        fillIndex <- (curDepth - 1) %% length(fill) + 1
        grid.rect(gp=gpar(col=col[colIndex], fill=fill[fillIndex]))
        if (label)
            labelVP(vp, col[colIndex])
        upViewport(depth(vp))
    }
}

drawVP.vpList <- function(vp, curDepth, depth, col, fill, label) {
    lapply(vp, drawVP, curDepth, depth, col, fill, label)
}

drawVP.vpStack <- function(vp, curDepth, depth, col, fill, label) {
    d <- depth(vp)
    for (i in 1:length(vp)) {
        this <- vp[[i]]
        drawVP(this, curDepth, depth, col, fill, label)
        curDepth <- curDepth + depth(this)
        pushViewport(this)
    }
    upViewport(d)
}

drawVP.vpTree <- function(vp, curDepth, depth, col, fill, label) {
    if (vp$parent$name == "ROOT") {
        lapply(vp$children, drawVP, curDepth, depth, col, fill, label)
    } else {
        pushViewport(vp$parent)
        if (is.null(depth) || curDepth %in% depth) {
            colIndex <- (curDepth - 1) %% length(col) + 1
            fillIndex <- (curDepth - 1) %% length(fill) + 1
            grid.rect(gp=gpar(col=col[colIndex], fill=fill[fillIndex]))
            if (label) {
                drawLabel <- is.null(vp$children) ||
                             (!is.null(depth) &&
                              curDepth == max(depth))
                if (drawLabel)
                    labelVP(vp$parent, col[colIndex])
            }
        }
        lapply(vp$children, drawVP, curDepth + 1, depth, col, fill, label)
        upViewport()
    }
}

# Draw all viewports in same viewport
showVP <- function(vp, newpage, cvpt, depth, col, fill,
                   label) {
    # If we've started a new page, we'll need the old
    # viewport tree to navigate within
    if (newpage) {
        pushViewport(cvpt)
        # "-1" for "ROOT"
        upViewport(depth(cvpt) - 1)
    }
    # Work off a vpTree, so convert vp if it's a vpPath
    showingPath <- inherits(vp, "vpPath")
    if (showingPath) {
        path <- vp
        downViewport(path)
        vp <- current.vpTree(all=FALSE)
        upViewport(1)
    }
    drawVP(vp, 1, depth, col, fill, label)
    if (showingPath)
        # "-1" because we went down the path then back up 1 originally
        upViewport(depth(path) - 1)
    invisible()
}

# Convert a "viewport" to a set of vpPaths
leafPaths <- function(vp) {
    UseMethod("leafPaths")
}

leafPaths.viewport <- function(vp) {
    if (vp$name == "ROOT")
        NULL
    else
        vp$name
}

leafPaths.vpList <- function(vp) {
    unlist(lapply(vp, leafPaths))
}

leafPaths.vpStack <- function(vp) {
    pathList <- lapply(vp, leafPaths)
    for (i in 1:length(pathList)) {
        if (i > 1) {
            pathList[[i]] <- paste(pathList[[i - 1]],
                                   pathList[[i]],
                                   sep=.grid.pathSep)
        }
    }
    unlist(pathList)
}

leafPaths.vpTree <- function(vp) {
    if (is.null(vp$children)) {
        if (vp$parent$name == "ROOT")
            NULL
        else
            vp$parent$name
    } else {
        pathList <- lapply(vp$children, leafPaths)
        if (vp$parent$name == "ROOT") {
            unlist(pathList)
        } else {
            paste(vp$parent$name,
                  unlist(pathList),
                  sep=.grid.pathSep)
        }
    }
}

leafPaths.vpPath <- function(vp) {
    as.character(vp)
}

# Draw a vpPath
drawPath <- function(path, depth, col, fill, label) {
    n <- depth(path)
    for (i in 1:n) {
        downViewport(path[i])
        if (is.null(depth) || i %in% depth) {
            colIndex <- (i - 1) %% length(col) + 1
            fillIndex <- (i - 1) %% length(fill) + 1
            grid.rect(gp=gpar(col=col[colIndex], fill=fill[fillIndex]))
            if (label) {
                if (is.null(depth))
                    drawLabel <- i == n
                else
                    drawLabel <- i == min(n, max(depth))
                if (drawLabel)
                    labelVP(current.viewport(), col[colIndex])
            }
        }
    }
    upViewport(n)
}

# Draw each leaf in separate viewports
# FIXME: allow control over number of rows and cols
# NOTE: this does NOT leave its viewports hanging around after
showVPmatrix <- function(vp, cvpt, depth, col, fill,
                         label, # Only the leaf viewports are labelled
                         nrow, ncol) {
    # Work off a vpPath, so convert vp if it's a "viewport"
    if (is.viewport(vp)) {
        paths <- leafPaths(vp)
    } else {
        # Should not happen
        stop("how did we get here?")
    }
    firstPath <- 0
    while (length(paths) - firstPath > 0) {
        if (firstPath > 0)
            grid.newpage()
        pushViewport(viewport(layout=grid.layout(nrow, ncol)))
        for (i in 1:nrow) {
            for (j in 1:ncol) {
                theLeaf <- firstPath + (i - 1)*nrow + j
                if (theLeaf <= length(paths)) {
                    thePath <- vpPath(paths[theLeaf])
                    pushViewport(viewport(layout.pos.row=i,
                                          layout.pos.col=j))
                    grid.rect(gp=gpar(col="grey80"))
                    # We may need the old vpTree to navigate within
                    # if 'vp' is a vpStack, or something similar, that
                    # contains a vpPath
                    if (!is.null(cvpt$children)) {
                        pushViewport(cvpt$children)
                        upViewport(depth(cvpt) - 1)
                    }
                    # Now push the viewport we are showing
                    pushViewport(vp)
                    upViewport(depth(vp))
                    # Now go to the particular viewport we
                    # are going to show
                    drawPath(thePath, depth, col, fill, label)
                    # Pop our placement within the layout
                    popViewport()
                }
            }
        }
        popViewport()
        firstPath <- firstPath + nrow*ncol
    }
}

showViewport <- function(vp=NULL, recurse=TRUE, depth=NULL,
                         newpage=FALSE, leaves=FALSE,
                         col=rgb(0, 0, 1, .2), fill=rgb(0, 0, 1, .1),
                         label=TRUE, nrow=3, ncol=nrow) {
    cvpt <- current.vpTree()
    if (is.null(vp))
        vp <- cvpt
    if (newpage == FALSE && leaves == TRUE)
        stop("must start new page if showing leaves separately")
    if (newpage) {
        grid.newpage()
    }
    if (!recurse)
        depth <- 1
    if (leaves) {
        # Special case of showing vpPath (i.e., only one viewport)
        # Ignores nrow & ncol
        if (inherits(vp, "vpPath"))
            showVP(vp, TRUE, cvpt, depth, col, fill, label)
        else
            showVPmatrix(vp, cvpt, depth, col, fill, label, nrow, ncol)
    } else {
        showVP(vp, newpage, cvpt, depth, col, fill, label)
    }
    invisible()
}
#  File src/library/grid/R/edit.R
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

# All args just used as args to editGrob
gEdit <- function(...) {
  edit <- list(...)
  class(edit) <- "gEdit"
  edit
}

is.gEdit <- function(x) {
  inherits(x, "gEdit")
}

applyEdit <- function(x, edit) {
  if (is.null(edit)) {
    x
  } else {
    if (!is.gEdit(edit))
      stop("invalid 'edit' information")
    # Intended to handle whether edit has gPath spec or not
    newx <- do.call("editGrob", c(list(x), edit))
    # If edit was specified for non-existent child, newx will be NULL
    if (is.null(newx))
      x
    else
      newx
  }
}

# A list of gEdit's to apply to the same grob
gEditList <- function(...) {
  edits <- list(...)
  if (!all(sapply(edits, is.gEdit)))
    stop("'gEditList' can only contain 'gEdit' objects")
  class(edits) <- "gEditList"
  edits
}

is.gEditList <- function(x) {
  inherits(x, "gEditList")
}

applyEdits <- function(x, edits) {
  if (is.null(edits)) {
    x
  } else {
    if (is.gEdit(edits))
      applyEdit(x, edits)
    else {
      if (!inherits(edits, "gEditList"))
        stop("invalid 'edit' information")
      for (i in edits)
        x <- applyEdits(x, i)
      x
    }
  }
}

#  File src/library/grid/R/frames.R
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

################
# frame class
################
# NOTE: make framevp separate slot (rather than combining with
# normal vp slot) so that it can be edited (e.g., by grid.pack)
frameGrob <- function(layout=NULL, name=NULL, gp=gpar(), vp=NULL) {
  framevp <- if(!is.null(layout)) viewport(layout=layout) # else NULL
  gTree(framevp=framevp, name=name, gp=gp, vp=vp, cl="frame")
}

# draw=TRUE will not draw anything, but will mean that
# additions to the frame are drawn
grid.frame <- function(layout=NULL, name=NULL, gp=gpar(), vp=NULL,
                       draw=TRUE) {
  fg <- frameGrob(layout=layout, name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(fg)
  invisible(fg)
}

makeContext.frame <- function(x) {
    if (!is.null(x$framevp)) {
        if (!is.null(x$vp)) {
            x$vp <- vpStack(x$vp, x$framevp)
        } else {
            x$vp <- x$framevp
        }
    }
    x
}

widthDetails.frame <- function(x) {
  if (is.null(x$framevp))
    unit(1, "null")
  else
    sum(layout.widths(viewport.layout(x$framevp)))
}

heightDetails.frame <- function(x) {
  if (is.null(x$framevp))
    unit(1, "null")
  else
    sum(layout.heights(viewport.layout(x$framevp)))
}

frameDim <- function(frame) {
  if (is.null(frame$framevp))
    rep(0, 2)
  else
    c(layout.nrow(viewport.layout(frame$framevp)),
      layout.ncol(viewport.layout(frame$framevp)))
}

################
# cellGrob class
################
cellViewport <- function(col, row, border) {
  vp <- viewport(layout.pos.col=col, layout.pos.row=row)
  if (!is.null(border))
    vpStack(vp,
            viewport(x=border[2L],
                     y=border[1L],
                     width =unit(1, "npc") - sum(border[c(2,4)]),
                     height=unit(1, "npc") - sum(border[c(1,3)]),
                     just=c("left", "bottom")))
  else vp
}

cellGrob <- function(col, row, border, grob, dynamic, vp) {
  gTree(col=col, row=row, border=border, dynamic=dynamic,
        children=gList(grob), cellvp=vp, cl="cellGrob")
}

makeContext.cellGrob <- function(x) {
    if (!is.null(x$cellvp)) {
        if (!is.null(x$vp)) {
            x$vp <- vpStack(x$vp, x$cellvp)
        } else {
            x$vp <- x$cellvp
        }
    }
    x
}

# For dynamically packed grobs, need to be able to
# recalculate cell sizes
widthDetails.cellGrob <- function(x) {
  if (x$dynamic)
    unit(1, "grobwidth", gPath(x$children[[1L]]$name))
  else
    unit(1, "grobwidth", x$children[[1L]])
}

heightDetails.cellGrob <- function(x) {
  if (x$dynamic)
    unit(1, "grobheight", gPath(x$children[[1L]]$name))
  else
    unit(1, "grobheight", x$children[[1L]])
}

################
# grid.place
################
# Place an object into an already existing cell of a frame ...
# ... for a grob on the display list
grid.place <- function(gPath, grob,
                       row=1, col=1,
                       redraw=TRUE) {
  grid.set(gPath,
           placeGrob(grid.get(gPath), grob, row, col),
           redraw)
}

# ... for a grob description
placeGrob <- function(frame, grob,
                      row=NULL, col=NULL) {
  if (!inherits(frame, "frame"))
    stop("invalid 'frame'")
  if (!is.grob(grob))
    stop("invalid 'grob'")
  dim <- frameDim(frame)
  if (is.null(row))
    row <- c(1, dim[1L])
  if (is.null(col))
    col <- c(1, dim[2L])
  if (min(row) < 1 || max(row) > dim[1L] ||
      min(col) < 1 || max(col) > dim[2L])
    stop("invalid 'row' and/or 'col' (no such cell in frame layout)")
  cgrob <- cellGrob(col, row, NULL, grob, FALSE,
                    cellViewport(col, row, NULL))
  addGrob(frame, cgrob)
}

################
# grid.pack
################
num.col.specs <- function(side, col, col.before, col.after) {
  4 - sum(is.null(side) || any(c("top", "bottom") %in% side),
          is.null(col), is.null(col.before), is.null(col.after))
}

# We are assuming that checking has been done so that only one
# of these specifications has been given
col.spec <- function(side, col, col.before, col.after, ncol) {
  if (!is.null(side)) {
    if (side == "left")
      col <- 1
    else if (side == "right")
      col <- ncol + 1
  }
  else if (!is.null(col.before))
    col <- col.before
  else if (!is.null(col.after))
    col <- col.after + 1
  col
}

# We are assuming that checking has been done so that only one
# of these specifications has been given
new.col <- function(side, col, col.before, col.after, ncol) {
  # Special case ncol==0 for first grob added to frame
  result <- TRUE
  if (!is.null(col)) {
    # It is an error to specify a range for col which is outside 1..ncol
    if (length(col) == 2)
      if (col[1L] < 1 || col[2L] > ncol)
        stop("'col' can only be a range of existing columns")
      else
        result <- FALSE
    # It is also an error to specify a single col outside 1..ncol+1
    else
      if (col < 1 || col > ncol + 1)
        stop("invalid 'col' specification")
      else
        result <- col == ncol+1
  }
  result
}

num.row.specs <- function(side, row, row.before, row.after) {
  4 - sum(is.null(side) || any(c("left", "right") %in% side),
          is.null(row), is.null(row.before), is.null(row.after))
}

# We are assuming that checking has been done so that only one
# of these specifications has been given
row.spec <- function(side, row, row.before, row.after, nrow) {
  if (!is.null(side)) {
    if (side == "top")
      row <- 1
    else if (side == "bottom")
      row <- nrow + 1
  }
  else if (!is.null(row.before))
    row <- row.before
  else if (!is.null(row.after))
    row <- row.after + 1
  row
}

# We are assuming that checking has been done so that only one
# of these specifications has been given
new.row <- function(side, row, row.before, row.after, nrow) {
  # Special case nrow==0 for first grob added to frame
  result <- TRUE
  if (!is.null(row)) {
    # It is an error to specify a range for row which is outside 1..nrow
    if (length(row) == 2)
      if (row[1L] < 1 || row[2L] > nrow)
        stop("'row' can only be a range of existing rows")
      else
        result <- FALSE
    # It is also an error to specify a single row outside 1..nrow+1
    else
      if (row < 1 || row > nrow + 1)
        stop("invalid 'row' specification")
      else
        result <- row == nrow+1
  }
  result
}

mod.dims <- function(dim, dims, index, new.index, nindex, force) {
  # If adding a new row/col, add the new width/height to the list
  if (new.index)
    if (index == 1)
      dims <- unit.c(dim, dims)
    else if (index == nindex)
      dims <- unit.c(dims, dim)
    else
      dims <- unit.c(dims[1L:(index-1)], dim, dims[index:nindex])
  # Otherwise, if force=TRUE, we override previous width/heights for the
  # row/col, otherotherwise, the width/height of the existing row/col
  # is the maximum of the previous width/height and the new width/height
  else {
    if (!force)
      dim <- max(dim, dims[index])
    if (index==1)
      if (nindex == 1)
        dims <- dim
      else
        dims <- unit.c(dim, dims[2:nindex])
    else if (index==nindex)
      dims <- unit.c(dims[1L:(nindex-1)], dim)
    else
      dims <- unit.c(dims[1L:(index-1)], dim, dims[(index+1):nindex])
  }
  dims
}

updateCol <- function(col, added.col) {
  old.col <- col
  # If grob$col is a range ...
  if (length(old.col) == 2) {
    if (added.col <= old.col[2L])
      col <- c(old.col[1L], old.col[2L] + 1)
  }
  else
    if (added.col <= old.col)
      col <- old.col + 1
  col
}

updateRow <- function(row, added.row) {
  old.row <- row
  # If grob$row is a range ...
  if (length(old.row) == 2) {
    if (added.row <= old.row[2L])
      row <- c(old.row[1L], old.row[2L] + 1)
  }
  else
    if (added.row <= old.row)
      row <- old.row + 1
  row
}

# FIXME:  Allow specification of respect for new row/col
# Pack a child grob within a frame grob ...
# (a special sort of editing just for frame grobs)
# ... for a grob on the display list
grid.pack <- function(gPath, grob, redraw=TRUE,
                      side=NULL,
                      row=NULL, row.before=NULL, row.after=NULL,
                      col=NULL, col.before=NULL, col.after=NULL,
                      width=NULL, height=NULL,
                      force.width=FALSE, force.height=FALSE,
                      border=NULL, dynamic=FALSE) {
  grid.set(gPath,
           packGrob(grid.get(gPath), grob, side,
                    row, row.before, row.after,
                    col, col.before, col.after,
                    width, height, force.width, force.height,
                    border),
           redraw)
}

packGrob <- function(frame, grob,
                     side=NULL,
                     row=NULL, row.before=NULL, row.after=NULL,
                     col=NULL, col.before=NULL, col.after=NULL,
                     width=NULL, height=NULL,
                     force.width=FALSE, force.height=FALSE,
                     border=NULL, dynamic=FALSE) {
  if (!inherits(frame, "frame"))
    stop("invalid 'frame'")
  if (!is.grob(grob))
    stop("invalid 'grob'")
  # col/row can be given as a range, but I only want to know
  # about the min and max
  if (!is.null(col) && length(col) > 1) {
    col <- range(col)
    col.range <- TRUE
  }
  else
    col.range <- FALSE
  if (!is.null(row) && length(row) > 1) {
    row <- range(row)
    row.range <- TRUE
  }
  else
    row.range <- FALSE

  frame.vp <- frame$framevp
  if (is.null(frame.vp))
    frame.vp <- viewport()
  lay <- viewport.layout(frame.vp)
  if (is.null(lay)) {
    ncol <- 0
    nrow <- 0
  } else {
    ncol <- layout.ncol(lay)
    nrow <- layout.nrow(lay)
  }

  # (i) Check that the specifications of the location of the grob
  # give a unique location
  ncs <- num.col.specs(side, col, col.before, col.after)
  # If user does not specify a col, assume it is all cols
  if (ncs == 0) {
    # Allow for fact that this might be first grob packed
    if (ncol > 0) {
      col <- c(1, ncol)
      col.range <- TRUE
    }
    else
      col <- 1
    ncs <- 1
  }
  if (ncs != 1)
    stop("cannot specify more than one of 'side=[\"left\", \"right\"]', 'col', 'col.before', or 'col.after'")
  nrs <- num.row.specs(side, row, row.before, row.after)
  # If user does not specify a row, assume it is all rows
  if (nrs == 0) {
    # Allow for fact that this might be first grob packed
    if (nrow > 0) {
      row <- c(1, nrow)
      row.range <- TRUE
    }
    else
      row <- 1
    nrs <- 1
  }
  if (nrs != 1)
    stop("must specify exactly one of 'side=[\"top\", \"bottom\"]', 'row', 'row.before', or 'row.after'")

  # (ii) Determine that location and check that it is valid
  new.col <- new.col(side, col, col.before, col.after, ncol)
  col <- col.spec(side, col, col.before, col.after, ncol)
  new.row <- new.row(side, row, row.before, row.after, nrow)
  row <- row.spec(side, row, row.before, row.after, nrow)

  # Wrap the child in a "cellGrob" to maintain additional info
  # (like row and col occupied in frame)
  # Need to do this here so can create widths/heights based on this cell grob
  if (!is.null(grob))
    cgrob <- cellGrob(col, row, border, grob, dynamic,
                      cellViewport(col, row, border))

  # (iii) If width and height are not given, take them from the child
  #       NOTE:  if dynamic is TRUE then use a gPath to the child
  if (is.null(width))
    if (is.null(grob))
      width <- unit(1, "null")
    else
      if (dynamic)
        width <- unit(1, "grobwidth", gPath(cgrob$name))
      else
        width <- unit(1, "grobwidth", cgrob)
  if (is.null(height))
    if (is.null(grob))
      height <- unit(1, "null")
    else
      if (dynamic)
        height <- unit(1, "grobheight", gPath(cgrob$name))
      else
        height <- unit(1, "grobheight", cgrob)
  # If there is a border, include it in the width/height
  if (!is.null(border)) {
    width <- sum(border[2L], width, border[4L])
    height <- sum(border[1L], height, border[3L])
  }

  # (iv) Update the frame.vp of the frame (possibly add new row/col,
  # possibly update existing widths/heights and respect)
  if (new.col) ncol <- ncol + 1
  if (new.row) nrow <- nrow + 1
  # If we are creating the frame.vp$layout for the first time then
  # we have to initialise the layout widths and heights
  if (is.null(lay)) {
    widths <- width
    heights <- height
  } else {
    # DO NOT modify widths/heights if the grob is being added to
    # multiple columns/rows
    if (col.range)
      widths <- layout.widths(lay)
    else
      widths <- mod.dims(width, layout.widths(lay), col, new.col, ncol,
                         force.width)
    if (row.range)
      heights <- layout.heights(lay)
    else
      heights <- mod.dims(height, layout.heights(lay), row, new.row, nrow,
                          force.height)
  }
  frame.vp$layout <- grid.layout(ncol=ncol, nrow=nrow,
                                 widths=widths, heights=heights)

  # Modify the locations (row, col) of existing children in the frame
  if (new.col || new.row) {
    for (i in childNames(frame)) {
      child <- getGrob(frame, i)
      if (new.col) {
        newcol <- updateCol(child$col, col)
        child <- editGrob(child, col=newcol,
                          cellvp=cellViewport(newcol, child$row, child$border))
      }
      if (new.row) {
        newrow <- updateRow(child$row, row)
        child <- editGrob(child, row=newrow,
                          cellvp=cellViewport(child$col, newrow, child$border))
      }
      frame <- addGrob(frame, child)
    }
  }

  # Add the new grob to the frame
  if (!is.null(grob)) {
    frame <- addGrob(frame, cgrob)
  }

  editGrob(frame, framevp=frame.vp)
}

#  File src/library/grid/R/function.R
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


validDetails.functiongrob <- function(x, ...) {
    if (x$n < 1)
        stop("invalid 'n'")
    if (!(is.character(x$range) && x$range %in% c("x", "y")))
        x$range <- as.numeric(x$range)
    if (!is.function(x$f))
        stop("invalid 'f'")
    x
}

genXY <- function(x) {
    if (is.numeric(x$range)) {
        range <- range(x$range)
    } else {
        if (x$range == "x")
            range <- current.viewport()$xscale
        else
            range <- current.viewport()$yscale
    }
    input <- seq(range[1], range[2], length.out=x$n)
    x$f(input)
}

makeContent.functiongrob <- function(x) {
    xy <- genXY(x)
    linesGrob(xy$x, xy$y, default.units=x$units,
              name=x$name, gp=x$gp, vp=x$vp)
}

xDetails.functiongrob <- function(x, theta) {
    xy <- genXY(x)
    xDetails(linesGrob(xy$x, xy$y, default.units=x$units), theta)
}

yDetails.functiongrob <- function(x, theta) {
    xy <- genXY(x)
    yDetails(linesGrob(xy$x, xy$y, default.units=x$units), theta)
}

widthDetails.functiongrob <- function(x) {
    xy <- genXY(x)
    widthDetails(linesGrob(xy$x, xy$y, default.units=x$units))
}

heightDetails.functiongrob <- function(x) {
    xy <- genXY(x)
    heightDetails(linesGrob(xy$x, xy$y, default.units=x$units))
}

functionGrob <- function(f, n=101, range="x", units="native",
                         name=NULL, gp=gpar(), vp=NULL) {
    grob(f=f, range=range, units=units, n=n,
         gp=gp, vp=vp, name=name, cl="functiongrob")
}

grid.function <- function(...) {
    grid.draw(functionGrob(...))
}

# Convenience wrappers
grid.abline <- function(intercept=0, slope=1, ...) {
    grid.function(function(x) list(x=x, y=intercept + slope*x), ...)
}

##############
# Tests
tests <- function() {

    # editing
    grid.newpage()
    pushViewport(viewport(xscale=c(0, 2*pi), yscale=c(-1, 1)))
    grid.function(function(x) list(x=x, y=sin(x)), name="fg")
    grid.edit("fg", n=10)
    grid.edit("fg", f=function(x) list(x=x, y=x))

    # x/y/width/height calculations
    grid.newpage()
    pushViewport(viewport(xscale=c(0, 2*pi), yscale=c(-2, 2)))
    grid.function(function(x) list(x=x, y=sin(x)), name="fg")
    grid.segments(0, 1,
                  grobX("fg", 135), grobY("fg", 135),
                  arrow=arrow())

}
#  File src/library/grid/R/gpar.R
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


# A "gpar" object is a list of graphics parameters
# A graphics parameter is a name-value pair

gpar <- function(...) {
  gp <- validGP(list(...))
  class(gp) <- "gpar"
  gp
}

is.gpar <- function(x) {
  inherits(x, "gpar")
}

print.gpar <- function(x, ...) {
  print(unclass(x), ...)
  invisible(x)
}

validGP <- function(gpars) {
  # Check a (non-NULL) gpar is not of length 0
  check.length <- function(gparname) {
    if (length(gpars[[gparname]]) == 0)
      stop(gettextf("'gpar' element '%s' must not be length 0", gparname),
           domain = NA)
  }
  # Check a gpar is numeric and not NULL
  numnotnull <- function(gparname) {
    if (!is.na(match(gparname, names(gpars)))) {
      if (is.null(gpars[[gparname]]))
        gpars[[gparname]] <<- NULL
      else {
        check.length(gparname)
        gpars[[gparname]] <<- as.numeric(gpars[[gparname]])
      }
    }
  }
  checkNA <- function(gparname) {
      if (!is.na(match(gparname, names(gpars)))) {
          if (any(is.na(gpars[[gparname]]))) {
              # ALL NA gets removed (ignored)
              if (all(is.na(gpars[[gparname]]))) {
                  gpars[[gparname]] <<- NULL
              } else {
                  stop(gettextf("mixture of missing and non-missing values for %s",
                                gparname),
                       domain=NA)
              }
          }
      }
  }
  # fontsize, lineheight, cex, lwd should be numeric and not NULL
  numnotnull("fontsize")
  checkNA("fontsize")
  numnotnull("lineheight")
  checkNA("lineheight")
  numnotnull("cex")
  checkNA("cex")
  numnotnull("lwd")
  numnotnull("lex")
  # gamma defunct in 2.7.0
  if ("gamma" %in% names(gpars)) {
    warning("'gamma' 'gpar' element is defunct")
    gpars$gamma <- NULL
  }
  numnotnull("alpha")
  # col and fill are converted in C code
  # BUT still want to check length > 0
  if (!is.na(match("col", names(gpars)))) {
      if (is.null(gpars$col))
          gpars$col <- NULL
      else
          check.length("col")
  }
  if (!is.na(match("fill", names(gpars)))) {
      if (is.null(gpars$fill)) {
          gpars$fill <- NULL
      } else {
          ## fill can be a simple colour (NA, integer, string)
          ## OR a "GridPattern" OR a list of "GridPattern"s
          if (!is.pattern(gpars$fill)) {
              if (is.list(gpars$fill)) {
                  if (!all(sapply(gpars$fill, is.pattern)))
                      stop("'fill' gpar list components must all be patterns")
                  class(gpars$fill) <- "GridPatternList"
              }
              check.length("fill")
          }
      }
  }
  # lty converted in C code
  # BUT still want to check for NULL and check length > 0
  if (!is.na(match("lty", names(gpars)))) {
    if (is.null(gpars$lty))
      gpars$lty <- NULL
    else
      check.length("lty")
  }
  if (!is.na(match("lineend", names(gpars)))) {
    if (is.null(gpars$lineend))
      gpars$lineend <- NULL
    else
      check.length("lineend")
  }
  if (!is.na(match("linejoin", names(gpars)))) {
    if (is.null(gpars$linejoin))
      gpars$linejoin <- NULL
    else
      check.length("linejoin")
  }
  # linemitre should be larger than 1
  numnotnull("linemitre")
  if (!is.na(match("linemitre", names(gpars)))) {
    if (any(gpars$linemitre < 1))
      stop("invalid 'linemitre' value")
  }
  # alpha should be 0 to 1
  if (!is.na(match("alpha", names(gpars)))) {
    if (any(gpars$alpha < 0) || any(gpars$alpha > 1))
      stop("invalid 'alpha' value")
  }
  # font should be integer and not NULL
  if (!is.na(match("font", names(gpars)))) {
    if (is.null(gpars$font))
      gpars$font <- NULL
    else {
      check.length("font")
      gpars$font <- as.integer(gpars$font)
    }
  }
  # fontfamily should be character
  if (!is.na(match("fontfamily", names(gpars)))) {
    if (is.null(gpars$fontfamily))
      gpars$fontfamily <- NULL
    else {
      check.length("fontfamily")
      gpars$fontfamily <- as.character(gpars$fontfamily)
      checkNA("fontfamily")
    }
  }
  # fontface can be character or integer;  map character to integer
  # store value in font
  # Illegal to specify both font and fontface
  if (!is.na(match("fontface", names(gpars)))) {
    if (!is.na(match("font", names(gpars))))
      stop("must specify only one of 'font' and 'fontface'")
    gpars$font <-
	if (is.null(gpars$fontface)) NULL # remove it
	else {
	    check.length("fontface")
	    if (is.numeric(gpars$fontface))
		as.integer(gpars$fontface)
	    else
		vapply(as.character(gpars$fontface),
		       function(ch) # returns integer
		       switch(ch,
			      plain = 1L,
			      bold  = 2L,
			      italic=, oblique = 3L,
			      bold.italic = 4L,
			      symbol= 5L,
					# These are Hershey variants
			      cyrillic=5L,
			      cyrillic.oblique=6L,
			      EUC   = 7L,
			      stop("invalid fontface ", ch)), 0L)
	}
    ## Remove fontface
    gpars$fontface <- NULL
  }
  gpars
}

# Method for subsetting "gpar" objects
`[.gpar` <- function(x, index, ...) {
    if (length(x) == 0)
        return(gpar())
    maxn <- do.call("max", lapply(x, length))
    newgp <- lapply(x, rep, length.out=maxn)
    newgp <- lapply(X = newgp, FUN = `[`, index, ...)
    class(newgp) <- "gpar"
    newgp
}

# possible gpar names
# The order must match the GP_* values in grid.h
.grid.gpar.names <- c("fill", "col", "gamma", "lty", "lwd", "cex",
                      "fontsize", "lineheight", "font", "fontfamily",
                      "alpha", "lineend", "linejoin", "linemitre",
                      "lex", "gradientFill",
                      # Keep fontface at the end because it is never
                      # used in C code (it gets mapped to font)
                      "fontface")

set.gpar <- function(gp, grob=NULL) {
  if (!is.gpar(gp))
    stop("argument must be a 'gpar' object")
  gp <- validGP(gp)
  temp <- grid.Call(C_getGPar)
  # gamma defunct in 2.7.0
  if ("gamma" %in% names(gp)) {
      warning("'gamma' 'gpar' element is defunct")
      gp$gamma <- NULL
  }
  # Special case "cex" (make it cumulative)
  if (match("cex", names(gp), nomatch=0L))
    tempcex <- temp$cex * gp$cex
  else
    tempcex <- temp$cex
  # Special case "alpha" (make it cumulative)
  if (match("alpha", names(gp), nomatch=0L))
    tempalpha <- temp$alpha * gp$alpha
  else
    tempalpha <- temp$alpha
  # Special case "lex" (make it cumulative)
  if (match("lex", names(gp), nomatch=0L))
    templex <- temp$lex * gp$lex
  else
    templex <- temp$lex
  ## resolve fill - could be a simple colour OR a "GridPattern"
  if (is.pattern(gp$fill)) {
      if (is.null(grob)) {
          class(gp$fill) <- c("GridViewportPattern", class(gp$fill))
      } else {
          if (inherits(grob, "gTree")) {
              class(gp$fill) <- c("GridGTreePattern", class(gp$fill))
          } else {
              class(gp$fill) <- c("GridGrobPattern", class(gp$fill))
          }
      }
  } else if (is.list(gp$fill)) {
      if (is.null(grob)) {
          ## Silently use just first pattern
          gp$fill <- gp$fill[[1]]
          class(gp$fill) <- c("GridViewportPatternList", class(gp$fill))
      } else {
          if (inherits(grob, "gTree")) {
              class(gp$fill) <- c("GridGTreePatternList", class(gp$fill))
          } else {
              class(gp$fill) <- c("GridGrobPatternList", class(gp$fill))
          }
      }      
  }
  # All other gpars
  temp[names(gp)] <- gp
  temp$cex <- tempcex
  temp$alpha <- tempalpha
  temp$lex <- templex
  if (is.null(grob)) {
      ## Do this as a .Call.graphics to get it onto the base display list
      grid.Call.graphics(C_setGPar, temp)
  } else {
      grid.Call(C_setGPar, temp)
  }
}

get.gpar <- function(names=NULL) {
  if (is.null(names)) {
    result <- grid.Call(C_getGPar)
    # drop gamma
    result$gamma <- NULL
  } else {
    if (!is.character(names) ||
        !all(names %in% .grid.gpar.names))
      stop("must specify only valid 'gpar' names")
    # gamma deprecated
    if ("gamma" %in% names) {
      warning("'gamma' 'gpar' element is defunct")
      names <- names[-match("gamma", names)]
    }
    result <- unclass(grid.Call(C_getGPar))[names]
  }
  class(result) <- "gpar"
  result
}

# When editing a gp slot, only update the specified gpars
# Assume gp is NULL or a gpar
# assume newgp is a gpar (and not NULL)
mod.gpar <- function(gp, newgp) {
  if (is.null(gp))
    gp <- newgp
  else
    gp[names(newgp)] <- newgp
  gp
}

#  File src/library/grid/R/grab.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2016 The R Core Team
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

#########
# Generate a gTree from the current display list
#
# Or from an expression
# (recording on to a null graphics device)
#########
rootVP <- function(pvp) {
  match(pvp$name, "ROOT", nomatch=FALSE)
}

# List the children of the current vp (as a vpList)
current.vpList <- function() {
  cpvp <- grid.Call(C_currentViewport)
  if (no.children(cpvp$children))
    NULL
  else
    vpListFromNode(cpvp)
}

current.vpNames <- function() {
  ls(grid.Call(C_currentViewport)$children)
}

# vp might be a viewport, or a vpList, or a vpStack, or a vpTree
vpExists <- function(vp) {
  UseMethod("vpExists")
}

vpExists.viewport <- function(vp) {
  exists(vp$name, .Call(C_currentViewport)$children)
}

vpExists.vpStack <- function(vp) {
  vpExists(vp[[1L]])
}

vpExists.vpList <- function(vp) {
  any(vapply(vp, vpExists, logical(1L), simplify=TRUE))
}

vpExists.vpTree <- function(vp) {
  vpExists(vp$parent)
}

# Handle vpPaths in a vpStack or vpTree
# Not a problem to downViewport() to a viewport that already exists
vpExists.vpPath <- function(vp) {
    FALSE
}

wrap <- function(x, ...) {
  UseMethod("wrap")
}

wrap.default <- function(x, ...) {
  if (!is.null(x))
    stop("invalid display list element")
  NULL
}

wrap.grob <- function(x, grobs=FALSE, ...) {
    if (grobs) {
        recordGrob(grid.draw(x), list(x=x))
    } else {
        x
    }
}

wrap.viewport <- function(x, ...) {
  recordGrob(pushViewport(vp), list(vp=x))
}

wrap.pop <- function(x, ...) {
  recordGrob(popViewport(n), list(n=x))
}

wrap.up <- function(x, ...) {
  recordGrob(upViewport(n), list(n=x))
}

wrap.vpPath <- function(x, ...) {
  recordGrob(downViewport(path), list(path=x))
}

# Grab the display list on the current device
# ... are passed to gTree
# If warn is 0, issue no warnings
# If warn is 1, issue warnings about situations that are definitely
#   NOT captured correctly (e.g., reuse of top-level grob name)
# If warn is 2, issue warnings about situations that
#   MAY not get captured correctly (e.g., top-level downViewport)
# If wrap is TRUE, grab will wrap all pushes and grobs
#   in a gTree
grabDL <- function(warn, wrap, wrap.grobs=FALSE, ...) {
  gList <- NULL
  dl.index <- grid.Call(C_getDLindex)
  if (dl.index > 1) {
      if (warn > 0 && !wrap.grobs) {
          names <- getNames()
          ## Check for overwriting existing grob
          if (length(unique(names)) != length(names))
              warning("one or more grobs overwritten (grab WILL not be faithful; try 'wrap.grobs = TRUE')")
      }
      if (!wrap) {
          grid.newpage(recording=FALSE)
      }
      ## Start at 2 because first element is viewport[ROOT]
      for (i in 2:dl.index) {
          ## Do all of this as a big ifelse rather than
          ## dispatching to a function call per element because
          ## we need to work with whole DL at times, not
          ## just individual elements
          elt <- grid.Call(C_getDLelt, as.integer(i - 1))
          if (wrap)
              gList <- addToGList(wrap(elt, grobs=wrap.grobs), gList)
          else {
              ## ####################
              ## grabGrob
              ## ####################
              if (inherits(elt, "grob")) {
                  ## Enforce grob$vp now and set grob$vp to NULL
                  ## Will be replaced later with full vpPath
                  tempvp <- elt$vp
                  if (warn > 1) {
                      ## Check to see if about to push a viewport
                      ## with existing viewport name
                      if (inherits(tempvp, "viewport") &&
                          vpExists(tempvp))
                          warning("viewport overwritten (grab MAY not be faithful)")
                  }
                  if (!is.null(tempvp))
                      tempdepth <- depth(tempvp)
                  grid.draw(tempvp, recording=FALSE)
                  ## vpPath after grob$vp slot has been pushed
                  ## Has to be recorded here in case grob drawing
                  ## pushes (and does not pop) more viewports
                  drawPath <- current.vpPath()
                  elt$vp <- NULL
                  grid.draw(elt, recording=FALSE)
                  if (warn > 1) {
                      ## Compare new vpPath
                      ## If not same, the grob has pushed some viewports
                      ## and not popped or upped them
                      pathSame <- TRUE
                      if (!(is.null(drawPath) && is.null(current.vpPath()))) {
                          if (is.null(drawPath))
                              pathSame <- FALSE
                          else if (is.null(current.vpPath()))
                              pathSame <- FALSE
                          else if (as.character(drawPath) !=
                                   as.character(current.vpPath()))
                              pathSame <- FALSE
                      }
                      if (!pathSame)
                          warning("grob pushed viewports and did not pop/up them (grab MAY not be faithful)")
                  }
                  elt$vp <- drawPath
                  if (!is.null(tempvp))
                      upViewport(tempdepth, recording=FALSE)
                  gList <- addToGList(elt, gList)
                  ## ####################
                  ## grabViewport
                  ## ####################
              } else if (inherits(elt, "viewport")) {
                  ## Includes viewports, vpLists, vpTrees, and vpStacks
                  ## Check to see if about to push a viewport
                  ## with existing viewport name
                  if (warn > 1) {
                      if (vpExists(elt))
                          warning("viewport overwritten (grab MAY not be faithful)")
                  }
                  grid.draw(elt, recording=FALSE)
                  ## ####################
                  ## grabPop
                  ## ####################
              } else if (inherits(elt, "pop")) {
                  ## Replace pop with up
                  upViewport(elt, recording=FALSE)

                  ## ####################
                  ## grabDefault
                  ## ####################
              } else {
                  grid.draw(elt, recording=FALSE)
              }
          } ## matches if (wrap)
      }
      ## Go to top level
      upViewport(0, recording=FALSE)
      gTree(children=gList, childrenvp=current.vpList(), ...)
  } else {
      NULL
  }
}

# expr is ignored if dev is NULL
# otherwise, it should be an expression, like postscript("myfile.ps")
grid.grab <- function(warn=2, wrap=wrap.grobs, wrap.grobs=FALSE, ...) {
  grabDL(warn, wrap, wrap.grobs, ...)
}

offscreen <- function(width, height) {
    pdf(file=NULL, width=width, height=height)
}

grid.grabExpr <- function(expr, warn=2, wrap=wrap.grobs, wrap.grobs=FALSE,
                          width=7, height=7, device=offscreen, ...) {
    ## Start an "offline" PDF device for this function
    cd <- dev.cur()
    device(width, height)
    grabd <- dev.cur()
    on.exit({ dev.set(grabd); dev.off(); dev.set(cd) })
    ## Run the graphics code in expr
    ## Rely on lazy evaluation for correct "timing"
    eval(expr)
    ## Grab the DL on the new device
    grabDL(warn, wrap, wrap.grobs, ...)
}

#########################
# A different sort of capture ...
# Just grab the screen raster image
#########################

grid.cap <- function() {
    # This does not need recording on the display list
    grid.Call(C_cap)
}


#  File src/library/grid/R/grid.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2016 The R Core Team
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


# FIXME:  all grid functions should check that .grid.started is TRUE
.grid.loaded <- FALSE

push.vp <- function(vp, recording) {
  UseMethod("push.vp")
}

push.vp.default <- function(vp, recording) {
  stop("only valid to push viewports")
}

push.vp.viewport <- function(vp, recording) {
  # Record on the display list
  if (recording)
    record(vp)
  # Store the entire set of gpar settings JUST PRIOR to push
  # We refer to this when calculating the viewport transform
  # We cannot simply rely on parent's gpar because we may be
  # being pushed from within a gTree which has enforced gpar
  # settings (i.e., the gTree$gp is enforced between this viewport
  # and the this viewport's parent$gp)
  vp$parentgpar <- grid.Call(C_getGPar)
  # Enforce gpar settings
  set.gpar(vp$gp)
  # Store the entire set of gpar settings for this viewport
  vp$gpar <- grid.Call(C_getGPar)
  # Pass in the pushedvp structure which will be used to store
  # things like the viewport transformation, parent-child links, ...
  # In C code, a pushedvp object is created, with a call to pushedvp(),
  # for the system to keep track of
  # (it happens in C code so that a "normal" vp gets recorded on the
  #  display list rather than a "pushedvp")
  grid.Call.graphics(C_setviewport, vp, TRUE)
}

# For all but the last viewport, push the
# viewport then pop it
# For the last viewport, just push
push.vp.vpList <- function(vp, recording) {
  push.vp.parallel <- function(vp, recording) {
    push.vp(vp, recording)
    upViewport(depth(vp), recording)
  }
  if (length(vp) == 1)
    push.vp(vp[[1L]], recording)
  else {
    lapply(vp[1L:(length(vp) - 1)], push.vp.parallel, recording)
    push.vp(vp[[length(vp)]], recording)
  }
}

# Push viewports in series
push.vp.vpStack <- function(vp, recording) {
  lapply(vp, push.vp, recording)
}

# Push parent
# Children are a vpList
push.vp.vpTree <- function(vp, recording) {
  # Special case if user has saved the entire vpTree
  # parent will be the ROOT viewport, which we don't want to
  # push (grid ensures it is ALWAYS there)
  if (!(vp$parent$name %in% "ROOT"))
    push.vp(vp$parent, recording)
  push.vp(vp$children, recording)
}

# "push"ing a vpPath is just a downViewport(..., strict=TRUE)
push.vp.vpPath <- function(vp, recording) {
    downViewport(vp, strict=TRUE, recording)
}

push.viewport <- function(..., recording=TRUE) {
    .Defunct("pushViewport")
}

pushViewport <- function(..., recording=TRUE) {
  if (missing(...))
    stop("must specify at least one viewport")
  else {
    vps <- list(...)
    lapply(vps, push.vp, recording)
  }
  invisible()
}

# Helper functions called from C
no.children <- function(children) {
  length(names(children)) == 0
}

child.exists <- function(name, children) {
  exists(name, envir=children, inherits=FALSE)
}

child.list <- function(children) {
  ls(children, all.names=TRUE) # sorted (needed ?)
}

pathMatch <- function(path, pathsofar, strict) {
  if (is.null(pathsofar))
    is.null(path)
  else {
    pattern <- paste0(if(strict) "^", path, "$")
    grepl(pattern, pathsofar)
  }
}

growPath <- function(pathsofar, name) {
  paste(pathsofar, name, sep=.grid.pathSep)
}

# Rather than pushing a new viewport, navigate down to one that has
# already been pushed
downViewport <- function(name, strict=FALSE, recording=TRUE) {
  UseMethod("downViewport")
}

# For interactive use, allow user to specify
# vpPath directly (i.e., w/o calling vpPath)
downViewport.default <- function(name, strict=FALSE, recording=TRUE) {
  name <- as.character(name)
  downViewport(vpPath(name), strict, recording=recording)
}

# Build vpPath from one (pushed) viewport up to another (pushed) viewport
# 'anc' is assumed to be an ancestor of 'desc'
# 'depth' is the depth that the final depth should have
buildPath <- function(desc, anc, depth) {
    path <- desc$name
    while (!identical(desc$parent, anc)) {
        if (is.null(desc$parent))
            stop("Down viewport failed to record on display list")
        desc <- desc$parent
        path <- c(desc$name, path)
    }
    result <- vpPath(path)
    if (depth(result) != depth)
        warning("Down viewport incorrectly recorded on display list")
    result
}

downViewport.vpPath <- function(name, strict=FALSE, recording=TRUE) {
    start <- grid.Call(C_currentViewport)
    if (name$n == 1)
        result <- grid.Call.graphics(C_downviewport, name$name, strict)
    else
        result <- grid.Call.graphics(C_downvppath,
                                     name$path, name$name, strict)
    # If the downViewport() fails, there is an error in C code
    # so none of the following code will be run

    # Enforce the gpar settings for the viewport
    pvp <- grid.Call(C_currentViewport)
    # Do not call set.gpar because set.gpar accumulates cex
    grid.Call.graphics(C_setGPar, pvp$gpar)
    # Record the viewport operation
    # ... including the depth navigated down
    if (recording) {
        attr(name, "depth") <- result
        # Record the strict path down
        path <- buildPath(pvp, start, result)
        record(path)
    }
    invisible(result)
}

# Similar to down.viewport() except it starts searching from the
# top-level viewport, so the result may be "up" or even "across"
# the current viewport tree
seekViewport <- function(name, recording=TRUE) {
  # up to the top-level
  upViewport(0, recording=recording)
  downViewport(name, recording=recording)
}

# Depth of the current viewport
vpDepth <- function() {
  pvp <- grid.Call(C_currentViewport)
  count <- 0
  while (!is.null(pvp$parent)) {
    pvp <- pvp$parent
    count <- count + 1
  }
  count
}

pop.viewport <- function(n=1, recording=TRUE) {
    .Defunct("popViewport")
}

popViewport <- function(n=1, recording=TRUE) {
  if (n < 0)
    stop("must pop at least one viewport")
  if (n == 0)
    n <- vpDepth()
  if (n > 0) {
    grid.Call.graphics(C_unsetviewport, as.integer(n))
    # Record on the display list
    if (recording) {
      class(n) <- "pop"
      record(n)
    }
  }
  invisible()
}

# Rather than removing the viewport from the viewport stack (tree),
# simply navigate up, leaving pushed viewports in place.
upViewport <- function(n=1, recording=TRUE) {
  if (n < 0)
    stop("must navigate up at least one viewport")
  if (n == 0) {
    n <- vpDepth()
    upPath <- current.vpPath()
  }
  if (n > 0) {
    path <- current.vpPath()
    upPath <- path[(depth(path) - n + 1):depth(path)]
    grid.Call.graphics(C_upviewport, as.integer(n))
    # Record on the display list
    if (recording) {
      class(n) <- "up"
      record(n)
    }
  }
  invisible(upPath)
}

# Return the full vpPath to the current viewport
current.vpPath <- function() {
  names <- NULL
  pvp <- grid.Call(C_currentViewport)
  while (!rootVP(pvp)) {
    names <- c(names, pvp$name)
    pvp <- pvp$parent
  }
  if (!is.null(names))
    vpPathFromVector(rev(names))
  else
    names
}

# Function to obtain the current viewport
current.viewport <- function() {
    # The system stores a pushedvp;  the user should only
    # ever see normal viewports, so convert.
    vpFromPushedvp(grid.Call(C_currentViewport))
}

# Return the parent of the current viewport
# (could be NULL)
current.parent <- function(n=1) {
    if (n < 1)
        stop("Invalid number of generations")
    vp <- grid.Call(C_currentViewport)
    generation <- 1
    while (generation <= n) {
        if (is.null(vp))
            stop("Invalid number of generations")
        vp <- vp$parent
        generation <- generation + 1
    }
    if (!is.null(vp))
        vpFromPushedvp(vp)
    else
        vp
}

vpListFromNode <- function(node) {
  vpListFromList(eapply(node$children, vpTreeFromNode, all.names=TRUE))
}

vpTreeFromNode <- function(node) {
  # If no children then just return viewport
  if (no.children(node$children))
    vpFromPushedvp(node)
  # Otherwise return vpTree
  else
    vpTree(vpFromPushedvp(node),
           vpListFromNode(node))
}

# Obtain the current viewport tree
# Either from the current location in the tree down
# or ALL of the tree
current.vpTree <- function(all=TRUE) {
  cpvp <- grid.Call(C_currentViewport)
  moving <- all && vpDepth() > 0
  if (moving) {
    savedpath <- current.vpPath()
    upViewport(0, recording=FALSE)
    cpvp <- grid.Call(C_currentViewport)
  }
  tree <- vpTreeFromNode(cpvp)
  if (moving) {
    downViewport(savedpath, recording=FALSE)
  }
  tree
}

current.transform <- function() {
    grid.Call(C_currentViewport)$trans
}

current.rotation <- function() {
    grid.Call(C_currentViewport)$rotation
}

# Call this function if you want the graphics device erased or moved
# on to a new page.  High-level plotting functions should call this.
# NOTE however, that if you write a function which calls grid.newpage,
# you should provide an argument to allow people to turn it off
# so that they can use your function within a parent viewport
# (rather than the whole device) if they want to.
grid.newpage <- function(recording=TRUE,
                         clearGroups=TRUE) {
    if (length(as.logical(clearGroups)) != 1)
        stop("Invalid 'clearGroups' argument")
    for (fun in getHook("before.grid.newpage"))  {
        if(is.character(fun)) fun <- get(fun)
        try(fun())
    }
    # NOTE that we do NOT do grid.Call here because we have to do
    # things slightly differently if grid.newpage is the first grid operation
    # on a new device
    .Call(C_newpagerecording)
    .Call(C_initGPar)
    .Call(C_newpage)
    .Call(C_clearDefinitions, as.logical(clearGroups))
    .Call(C_initViewportStack)
    if (recording) {
        .Call(C_initDisplayList)
        grDevices:::recordPalette()
        for (fun in getHook("grid.newpage"))  {
            if(is.character(fun)) fun <- get(fun)
            try(fun())
        }
    }
    invisible()
}

###########
# DISPLAY LIST FUNCTIONS
###########

# Keep a list of all drawing operations (since last grid.newpage()) so
# that we can redraw upon edit.

inc.display.list <- function() {
  display.list <- grid.Call(C_getDisplayList)
  dl.index <- grid.Call(C_getDLindex)
  dl.index <- dl.index + 1
  n <- length(display.list)
  # The " - 1" below is because dl.index is now stored internally
  # so is a C-style zero-based index rather than an R-style
  # 1-based index
  if (dl.index > (n - 1)) {
    temp <- display.list
    display.list <- vector("list", n + 100L)
    display.list[1L:n] <- temp
  }
  grid.Call(C_setDisplayList, display.list)
  grid.Call(C_setDLindex, as.integer(dl.index))
}

# This will either ...
#   (i) turn on AND INITIALISE the display list or ...
#   (ii) turn off AND ERASE the display list
grid.display.list <- function(on=TRUE) {
    grid.Call(C_setDLon, as.logical(on))
    if (on) {
        .Call(C_initDisplayList)
    } else {
        .Call(C_setDisplayList, NULL)
        .Call(C_setDLindex, 0L)
    }
}

record <- function(x) {
  if (grid.Call(C_getDLon))
    UseMethod("record")
}

# When there is a pop.viewport, the number of viewports popped
# gets put on the display list
record.default <- function(x) {
  if (!is.numeric(x))
    stop("invalid object inserted on the display list")
  grid.Call(C_setDLelt, x)
  inc.display.list()
}

record.grob <- function(x) {
  grid.Call(C_setDLelt, x)
  inc.display.list()
}

record.viewport <- function(x) {
  grid.Call(C_setDLelt, x)
  inc.display.list()
}

record.vpPath <- function(x) {
  grid.Call(C_setDLelt, x)
  inc.display.list()
}

# This controls whether grid is using the graphics engine's display list
engine.display.list <- function(on=TRUE) {
  grid.Call(C_setEngineDLon, as.logical(on))
}

# Rerun the grid DL
grid.refresh <- function() {
  draw.all()
}

# Call a function on each element of the grid display list
# AND replace the element with the result
# This is blood-curdlingly dangerous for the state of the
# display list
# Two token efforts at safety are made:
#   - generate all of the new elements first THEN assign them all
#     (so if there is an error in generating any one element
#      you don't end up with a trashed display list)
#   - check that the new element is either NULL or the same
#     class as the element it is replacing
grid.DLapply <- function(FUN, ...) {
    FUN <- match.fun(FUN)
    # Traverse DL and do something to each entry
#    gridDL <- grid.Call(C_getDisplayList)
    gridDLindex <- grid.Call(C_getDLindex)
    newDL <- vector("list", gridDLindex)
    for (i in 1:(gridDLindex - 1)) {
        elt <- grid.Call(C_getDLelt, i)
        newElt <- FUN(elt, ...)
        if (!(is.null(newElt) || inherits(newElt, class(elt))))
            stop("invalid modification of the display list")
        newDL[[i]] <- newElt
    }
    for (i in 1:(gridDLindex - 1)) {
        grid.Call(C_setDLindex, i)
        grid.Call(C_setDLelt, newDL[[i]])
    }
    grid.Call(C_setDLindex, gridDLindex)
}

# Wrapper for .Call and .Call.graphics
# Used to make sure that grid-specific initialisation occurs just before
# the first grid graphics output OR the first querying of grid state
# (on the current device)
# The general rule is you should use these rather than .Call or
# .Call.graphics unless you have a good reason and you know what
# you are doing -- this will be a bit of overkill, but is for safety
grid.Call <- function(fnname, ...) {
  .Call(C_gridDirty)
  .Call(dontCheck(fnname), ...)  # skip code analysis checks, keep runtime checks
}

grid.Call.graphics <- function(fnname, ...) {
  # Only record graphics operations on the graphics engine's display
  # list if the engineDLon flag is set
  engineDLon <- grid.Call(C_getEngineDLon)
  if (engineDLon) {
    # NOTE that we need a .Call.graphics(C_gridDirty) so that
    # the first thing on the engine display list is a dirty
    # operation;  this is necessary in case the display list is
    # played on another device (e.g., via replayPlot() or dev.copy())
    .Call.graphics(C_gridDirty)
    result <- .Call.graphics(dontCheck(fnname), ...)
  } else {
    .Call(C_gridDirty)
    result <- .Call(dontCheck(fnname), ...)
  }
  result
}

# A call to recordGraphics() outside of [pre|post]drawDetails methods
# will not record the expr on the grid DL.
# If a user REALLY wants to call recordGraphics(), they should use
# grid.record() instead
drawDetails.recordedGrob <- function(x, recording) {
  eval(x$expr, x$list, getNamespace("grid"))
}

grid.record <- function(expr, list,
                        name=NULL, gp=NULL, vp=NULL) {
  grid.draw(grob(expr=substitute(expr), list=list,
                 name=name, gp=gp, vp=vp, cl="recordedGrob"))
}

recordGrob <- function(expr, list,
                       name=NULL, gp=NULL, vp=NULL) {
  grob(expr=substitute(expr), list=list,
       name=name, gp=gp, vp=vp, cl="recordedGrob")
}

# Must only generate a grob, not modify drawing context
makeContent.delayedgrob <- function(x) {
    grob <- eval(x$expr, x$list, getNamespace("grid"))
    if (is.grob(grob)) {
        children <- gList(grob)
    } else if (is.gList(grob)) {
        children <- grob
    } else {
        stop("'expr' must return a grob or gList")
    }
    x <- setChildren(x, children)
    x
}

grid.delay <- function(expr, list,
                       name=NULL, gp=NULL, vp=NULL) {
    grid.draw(gTree(expr=substitute(expr), list=list,
                    name=name, gp=gp, vp=vp, cl="delayedgrob"))
}

delayGrob <- function(expr, list,
                      name=NULL, gp=NULL, vp=NULL) {
    gTree(expr=substitute(expr), list=list,
          name=name, gp=gp, vp=vp, cl="delayedgrob")
}

#  File src/library/grid/R/grob.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2016 The R Core Team
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

######################################
# Grid graphical objects
#######################################

################
# CLASS DEFN
################
# A "virtual" class "gDesc" underlies both "grob" and "gPath"

initGrobAutoName <- function() {
  index <- 0
  function(prefix="GRID", suffix="GROB") {
    index <<- index + 1
    paste(prefix, suffix, index, sep=".")
  }
}

grobAutoName <- initGrobAutoName()

# Function for user to call to get "autogenerated" grob name
grobName <- function(grob=NULL, prefix="GRID") {
    if (is.null(grob))
        grobAutoName(prefix)
    else {
        if (!is.grob(grob))
            stop("invalid 'grob' argument")
        else
            grobAutoName(prefix, class(grob)[1L])
    }
}

################
# CLASS DEFN
################
# A grob has a name, a gp, and a vp
# grob inherits from gDesc
checkvpSlot <- function(vp) {
  # vp can be a viewport, a viewport name, or a viewport path
  if (!is.null(vp))
    if (!inherits(vp, "viewport") &&
        !inherits(vp, "vpPath") &&
        !is.character(vp))
      stop("invalid 'vp' slot")
  # For interactive use, allow user to specify
  # vpPath directly (i.e., w/o calling vpPath)
  if (is.character(vp))
    vp <- vpPath(vp)
  vp
}

checkNameSlot <- function(x) {
  # Supply a default name if one is not given
  if (is.null(x$name))
    grobAutoName(suffix=class(x)[1L])
  else
    as.character(x$name)
}

checkgpSlot <- function(gp) {
  # gp must be a gpar
  if (!is.null(gp))
    if (!inherits(gp, "gpar"))
      stop("invalid 'gp' slot")
}

validDetails <- function(x) {
  UseMethod("validDetails")
}

validDetails.grob <- function(x) {
  x
}

validGrob <- function(x, ...) {
  UseMethod("validGrob")
}

validGrob.grob <- function(x, ...) {
  # Validate class-specific slots
  x <- validDetails(x)
  # Validate standard grob slots
  x$name <- checkNameSlot(x)
  checkgpSlot(x$gp)
  if (!is.null(x$vp))
    x$vp <- checkvpSlot(x$vp)
  return(x)
}

# This actually creates a new class derived from grob
# and returns an instance of that new class, all in one step
grob <- function(..., name=NULL, gp=NULL, vp=NULL, cl=NULL) {
  g <- list(..., name=name, gp=gp, vp=vp)
  if (!is.null(cl) &&
      !is.character(cl))
    stop("invalid 'grob' class")
  class(g) <- c(cl, "grob", "gDesc")
  validGrob(g)
}

grid.grob <- function(list.struct, cl=NULL, draw=TRUE) .Defunct("grob")

is.grob <- function(x) {
  inherits(x, "grob")
}

as.character.grob <- function(x, ...) {
  paste0(class(x)[1L], "[", x$name, "]")
}

print.grob <- function(x, ...) {
  cat(as.character(x), "\n")
  invisible(x)
}

################
# gPath CLASS DEFN
################
# gPath is a concatenated list of names specifying a path to a grob
# Functions for creating "paths" of viewport names

gPathFromVector <- function(names) {
  if (any(bad <- !is.character(names)))
      stop(ngettext(sum(bad), "invalid grob name", "invalid grob names"),
           domain = NA)
  # Break out any embedded .grid.pathSep's
  names <- unlist(strsplit(names, .grid.pathSep))
  n <- length(names)
  if (n < 1L)
    stop("a 'grob' path must contain at least one 'grob' name")
  path <- list(path = if (n==1) NULL else
               paste(names[1L:(n-1)], collapse = .grid.pathSep),
               name = names[n], n = n)
  class(path) <- c("gPath", "path")
  path
}

gPath <- function(...) {
  names <- c(...)
  gPathFromVector(names)
}

################
# gList CLASS DEFN
################
# Just a list of grobs
okGListelt <- function(x) {
  is.grob(x) || is.null(x) || is.gList(x)
}

is.gList <- function(x) {
    inherits(x, "gList")
}

as.gList <- function(x) {
    if (is.null(x)) {
        result <- list()
        class(result) <- "gList"
    } else if (is.grob(x)) {
        result <- list(x)
        class(result) <- "gList"
    } else if (is.gList(x)) {
        result <- x
    } else {
        stop("unable to coerce to \"gList\"")
    }
    result
}

gList <- function(...) {
    gl <- list(...)
    if (length(gl) == 0L ||
        all(sapply(gl, okGListelt, simplify=TRUE))) {
        # Ensure gList is "flat"
        # Don't want gList containing gList ...
        if (!all(sapply(gl, is.grob)))
            gl <- do.call("c", lapply(gl, as.gList))
        class(gl) <- c("gList")
        return(gl)
    } else {
        stop("only 'grobs' allowed in \"gList\"")
    }
}

addToGList <- function(x, gList) {
  UseMethod("addToGList")
}

addToGList.default <- function(x, gList) {
  if (is.null(x))
    gList
  else
    stop("invalid element to add to \"gList\"")
}

addToGList.grob <- function(x, gList) {
  if (is.null(gList))
    gList(x)
  else {
    gList[[length(gList) + 1L]] <- x
    return(gList)
  }
}

addToGList.gList <- function(x, gList) {
  gl <- c(gList, x)
  class(gl) <- "gList"
  return(gl)
}

as.character.gList <- function(x, ...) {
  paste0("(", paste(lapply(x, as.character), collapse=", "), ")")
}

print.gList <- function(x, ...) {
  cat(as.character(x), "\n")
  invisible(x)
}

`[.gList` <- function(x, index, ...) {
    cl <- class(x)
    result <- "["(unclass(x), index, ...)
    class(result) <- cl
    result
}

################
# gTree CLASS DEFN
################
# gTree extends grob
# A gTree has additional children slot
childName <- function(x) {
  x$name
}

setChildren <- function(x, children) {
  if (!inherits(x, "gTree"))
    stop("can only set 'children' for a \"gTree\"")
  if (!is.null(children) &&
      !inherits(children, "gList"))
    stop("'children' must be a \"gList\"")
  # Thin out NULL children
  if (!is.null(children)) {
    cl <- class(children)
    children <- children[!sapply(children, is.null)]
    class(children) <- cl
  }
  if (length(children)) {
    x$children <- children
    childNames <- sapply(children, childName)
    names(x$children) <- childNames
    x$childrenOrder <- childNames
  } else {
    x$children <- gList()
    x$childrenOrder <- character()
  }
  x
}

childNames <- function(gTree) {
  if (!inherits(gTree, "gTree"))
    stop("it is only valid to get 'children' from a \"gTree\"")
  gTree$childrenOrder
}

validGrob.gTree <- function(x, childrenvp, ...) {
  # Validate class-specific slots
  x <- validDetails(x)
  # Validate standard grob slots
  x$name <- checkNameSlot(x)
  checkgpSlot(x$gp)
  if (!is.null(x$vp))
    x$vp <- checkvpSlot(x$vp)
  # Only add childrenvp here so that gTree slots can
  # be validated before childrenvp get made
  # (making of childrenvp and children likely to depend
  #  on gTree slots)
  if (!is.null(childrenvp))
    x$childrenvp <- checkvpSlot(childrenvp)
  return(x)
}

gTree <- function(..., name=NULL, gp=NULL, vp=NULL,
                  children=NULL, childrenvp=NULL,
                  cl=NULL) {
  gt <- list(..., name=name, gp=gp, vp=vp)
  if (!is.null(cl) &&
      !is.character(cl))
    stop("invalid \"gTree\" class")
  class(gt) <- c(cl, "gTree", "grob", "gDesc")
  gt <- validGrob(gt, childrenvp)
  gt <- setChildren(gt, children)
  return(gt)
}

# A basic gTree that is JUST a collection of grobs
# (simply interface to gTree)
grobTree <- function(..., name=NULL, gp=NULL, vp=NULL,
                     childrenvp=NULL, cl=NULL) {
    gTree(children=gList(...),
          name=name, gp=gp, vp=vp,
          childrenvp=childrenvp, cl=cl)
}

################
# Getting just the names of the top-level grobs on the DL
################
getName <- function(elt) {
  if (inherits(elt, "grob"))
    elt$name
  else
    ""
}

getNames <- function() {
  dl <- grid.Call(C_getDisplayList)[1L:grid.Call(C_getDLindex)]
  names <- sapply(dl, getName)
  names[nzchar(names)]
}

################
# Getting/adding/removing/editing (children of [children of ...]) a gTree
################

# NOTE:  In order to cut down on repeated code, some of these
# (i.e., all but get and set) are inefficient and call get/set
# to do their work.  If speed becomes an issue, may have to
# revert to individual support for each function with highly
# repetitive code

# Get a grob from the display list
grid.get <- function(gPath, strict=FALSE, grep=FALSE, global=FALSE,
                     allDevices=FALSE) {
  if (allDevices)
    stop("'allDevices' not yet implemented")
  if (is.character(gPath))
    gPath <- gPath(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  getDLfromGPath(gPath, strict, grep, global)
}

# Just different defaults to grid.get for convenience
# Justified by usage patterns of Hadley Wickham
grid.gget <- function(..., grep=TRUE, global=TRUE) {
    grid.get(..., grep=grep, global=global)
}

# Get a child (of a child, of a child, ...) of a grob
getGrob <- function(gTree, gPath, strict=FALSE,
                    grep=FALSE, global=FALSE) {
  if (!inherits(gTree, "gTree"))
    stop("it is only valid to get a child from a \"gTree\"")
  if (is.character(gPath))
    gPath <- gPath(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (depth(gPath) == 1 && strict) {
    gTree$children[[gPath$name]]
  } else {
    if (!is.logical(grep))
      stop("invalid 'grep' value")
    grep <- rep(grep, length.out=depth(gPath))
    getGTree(gTree, NULL, gPath, strict, grep, global)
  }
}

# Set a grob on the display list
# nor is it valid to specify a global destination (i.e., no global arg)
grid.set <- function(gPath, newGrob, strict=FALSE, grep=FALSE,
                     redraw=TRUE) {
  if (is.character(gPath))
    gPath <- gPath(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  result <- setDLfromGPath(gPath, newGrob, strict, grep)
  # result$index will be non-zero if matched the gPath
  if (result$index) {
    # Get the current DL index
    dl.index <- grid.Call(C_getDLindex)
    # Destructively modify the DL elt
    grid.Call(C_setDLindex, as.integer(result$index))
    grid.Call(C_setDLelt, result$grob)
    # Reset the DL index
    grid.Call(C_setDLindex, as.integer(dl.index))
    if (redraw)
      draw.all()
  } else {
    stop("'gPath' does not specify a valid child")
  }
}

# Set a grob
# nor is it valid to specify a global destination (i.e., no global arg)
setGrob <- function(gTree, gPath, newGrob, strict=FALSE, grep=FALSE) {
  if (!inherits(gTree, "gTree"))
    stop("it is only valid to set a child of a \"gTree\"")
  if (!inherits(newGrob, "grob"))
    stop("it is only valid to set a 'grob' as child of a \"gTree\"")
  if (is.character(gPath))
    gPath <- gPath(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  if (depth(gPath) == 1 && strict) {
    # gPath must specify an existing child
    if (old.pos <- nameMatch(gPath$name, gTree$childrenOrder, grep)) {
      # newGrob name must match existing name
      if (match(gTree$childrenOrder[old.pos], newGrob$name, nomatch=0L)) {
        gTree$children[[newGrob$name]] <- newGrob
      } else {
          stop(gettextf("New 'grob' name (%s) does not match 'gPath' (%s)",
                        newGrob$name, gPath), domain = NA)
      }
    } else {
        stop("'gPath' does not specify a valid child")
    }
  } else {
    gTree <- setGTree(gTree, NULL, gPath, newGrob, strict, grep)
    if (is.null(gTree))
      stop("'gPath' does not specify a valid child")
  }
  gTree
}

# Add a grob to a grob on the display list
grid.add <- function(gPath, child, strict=FALSE,
                     grep=FALSE, global=FALSE, allDevices=FALSE,
                     redraw=TRUE) {
  if (allDevices)
    stop("'allDevices' not yet implemented")
  if (is.character(gPath))
    gPath <- gPath(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  addDLfromGPath(gPath, child, strict, grep, global, redraw)
}

# Add a grob to a gTree (or a child of a (child of a ...) gTree)
addGrob <- function(gTree, child, gPath=NULL, strict=FALSE,
                    grep=FALSE, global=FALSE, warn=TRUE) {
    if (!inherits(child, "grob"))
        stop("it is only valid to add a 'grob' to a \"gTree\"")
    if (is.null(gPath)) {
        addToGTree(gTree, child)
    } else {
        if (is.character(gPath))
            gPath <- gPath(gPath)
        # Only makes sense to specify a gPath for a gTree
        if (!inherits(gTree, "gTree"))
            stop("it is only valid to add a child to a \"gTree\"")
        if (!is.logical(grep))
            stop("invalid 'grep' value")
        grep <- rep(grep, length.out=depth(gPath))
        # result will be NULL if no match
        result <- addGTree(gTree, child, NULL, gPath, strict, grep, global)
        if (is.null(result)) {
            if (warn)
                warning(gettextf("'gPath' (%s) not found",
                                 as.character(gPath)),
                        domain = NA)
            gTree
        } else {
            result
        }
    }
}

# Remove a grob (or child of ...) from the display list
grid.remove <- function(gPath, warn=TRUE, strict=FALSE,
                        grep=FALSE, global=FALSE, allDevices=FALSE,
                        redraw=TRUE) {
  if (allDevices)
    stop("'allDevices' not yet implemented")
  if (is.character(gPath))
    gPath <- gPath(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  if (depth(gPath) == 1) {
    removeNameFromDL(gPath$name, strict, grep, global, warn, redraw)
  } else {
    name <- gPath$name
    gPath <- gPath(gPath$path)
    greppath <- grep[-length(grep)]
    grepname <- grep[length(grep)]
    removeDLFromGPath(gPath, name, strict, greppath, grepname,
                      global, warn, redraw)
  }
}

# Just different defaults to grid.remove for convenience
# Justified by usage patterns of Hadley Wickham
grid.gremove <- function(..., grep=TRUE, global=TRUE) {
    grid.remove(..., grep=grep, global=global)
}

# Remove a child from a (child of ...) gTree
removeGrob <- function(gTree, gPath, strict=FALSE,
                       grep=FALSE, global=FALSE, warn=TRUE) {
    if (!inherits(gTree, "gTree"))
        stop("it is only valid to remove a child from a \"gTree\"")
    if (is.character(gPath))
        gPath <- gPath(gPath)
    if (!inherits(gPath, "gPath"))
        stop("invalid 'gPath'")
    if (!is.logical(grep))
        stop("invalid 'grep' value")
    grep <- rep(grep, length.out=depth(gPath))
    if (depth(gPath) == 1) {
        # result will be NULL if no match
        result <- removeName(gTree, gPath$name, strict, grep, global, warn)
    } else {
        name <- gPath$name
        gPath <- gPath(gPath$path)
        greppath <- grep[-length(grep)]
        grepname <- grep[length(grep)]
        # result will be NULL if no match
        result <- removeGTree(gTree, name, NULL, gPath, strict,
                              greppath, grepname, global, warn)
    }
    if (is.null(result)) {
        if (warn)
            warning(gettextf("'gPath' (%s) not found", as.character(gPath)),
                    domain = NA)
        gTree
    } else {
        result
    }
}

# Edit a grob on the display list
grid.edit <- function(gPath, ..., strict=FALSE,
                      grep=FALSE, global=FALSE, allDevices=FALSE,
                      redraw=TRUE) {
  if (allDevices)
    stop("'allDevices' not yet implemented")
  if (is.character(gPath))
    gPath <- gPath(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  specs <- list(...)
  editDLfromGPath(gPath, specs, strict, grep, global, redraw)
}

# Just different defaults to grid.edit for convenience
# Justified by usage patterns of Hadley Wickham
grid.gedit <- function(..., grep=TRUE, global=TRUE) {
    grid.edit(..., grep=grep, global=global)
}

# Edit a (child of a ...) grob
editGrob <- function(grob, gPath=NULL, ..., strict=FALSE,
                     grep=FALSE, global=FALSE, warn=TRUE) {
    specs <- list(...)
    if (is.null(gPath)) {
        editThisGrob(grob, specs)
    } else {
        if (is.character(gPath))
            gPath <- gPath(gPath)
        # Only makes sense to specify a gPath for a gTree
        if (!inherits(grob, "gTree"))
            stop("it is only valid to edit a child of a \"gTree\"")
        if (!is.logical(grep))
            stop("invalid 'grep' value")
        grep <- rep(grep, length.out=depth(gPath))
        # result will be NULL if no match
        result <- editGTree(grob, specs, NULL, gPath, strict, grep, global)
        if (is.null(result)) {
            if (warn)
                warning(gettextf("'gPath' (%s) not found",
                                 as.character(gPath)),
                        domain = NA)
            grob
        } else {
            result
        }
    }
}

#########
# Generic "hook" to allow customised action on edit
#########
editDetails <- function(x, specs) {
  UseMethod("editDetails")
}

editDetails.default <- function(x, specs) {
  # Do nothing BUT return object being edited
  x
}

editDetails.gTree <- function(x, specs) {
  # Disallow editing children or childrenOrder slots directly
  if (any(specs %in% c("children", "childrenOrder")))
    stop("it is invalid to directly edit the 'children' or 'childrenOrder' slot")
  x
}

#########
# Helper functions for getting/adding/removing/editing grobs
#
# ASSUME down here that the grep argument has been replicated
# up to the length of the gPath argument
#########

# Find a "match" between a path$name and a grob$name
nameMatch <- function(pathName, grobName, grep) {
  if (grep) {
    pos <- grep(pathName, grobName)
    (length(pos) && pos == 1)
  } else {
    match(pathName, grobName, nomatch=0L)
  }
}

# Return the position of path$name in vector of names
# Return FALSE if not found
# If grep=TRUE, the answer may be a vector!
namePos <- function(pathName, names, grep) {
  if (grep) {
    pos <- grep(pathName, names)
    if (length(pos) == 0L)
      pos <- FALSE
  } else {
    pos <- match(pathName, names, nomatch=0L)
  }
  pos
}

partialPathMatch <- function(pathsofar, path, strict=FALSE, grep) {
  if (strict) {
    if (!any(grep))
      length(grep(paste0("^", pathsofar), path)) > 0L
    else {
      pathSoFarElts <- explode(pathsofar)
      pathElts <- explode(path)
      ok <- TRUE
      npsfe <- length(pathSoFarElts)
      index <- 1
      while (ok && index <= npsfe) {
        if (grep[index])
          ok <- (grep(pathSoFarElts[index], pathElts[index]) == 1)
        else
          ok <- match(pathSoFarElts[index], pathElts[index], nomatch=0L)
        index <- index + 1
      }
      ok
    }
  } else {
    # If we're not doing strict matching then anything from a full
    # path match to absolutely no match means a partial match
    # (i.e., keep looking)
    TRUE
  }
}

fullPathMatch <- function(pathsofar, gPath, strict, grep) {
  if (is.null(pathsofar))
    match <- (depth(gPath) == 1)
  else {
    path <- gPath$path
    if (!any(grep))
      if (strict)
        match <- match(pathsofar, path, nomatch=0L)
      else
        match <- (length(grep(paste0(path, "$"), pathsofar)) > 0L)
    else {
      pathSoFarElts <- explode(pathsofar)
      pathElts <- explode(path)
      npsfe <- length(pathSoFarElts)
      npe <- length(pathElts)
      if (npe > npsfe) {
        match <- FALSE
      } else {
        match <- TRUE
        index <- 1
        if (strict) {# pathSoFar same length as gPath
        } else {# pathSoFar could be longer than gPath
          pathSoFarElts <- pathSoFarElts[(npsfe - npe + 1):npsfe]
        }
        while (match && index <= npe) {
          if (grep[index])
            match <- (length(grep(pathElts[index], pathSoFarElts[index])) > 0L)
          else
            match <- match(pathSoFarElts[index], pathElts[index], nomatch = 0L)
          index <- index + 1
        }
      }
    }
  }
  match
}

#####
##### Get support
#####

# Add a grob to a result
growResult <- function(result, x) {
  UseMethod("growResult")
}

# Should only be when result is NULL
growResult.default <- function(result, x) {
  if (!is.null(result))
    stop("invalid 'result'")
  x
}

growResult.grob <- function(result, x) {
  if (is.grob(x))
    gList(result, x)
  else
    # x should be a gList
    addToGList(result, x)
}

growResult.gList <- function(result, x) {
  addToGList(x, result)
}

# A gPath may specify the child of a gTree
# (or the child of a child of a gTree, or ...)
getGrobFromGPath <- function(grob, pathsofar, gPath, strict,
                             grep, global) {
  UseMethod("getGrobFromGPath")
}

# If it's not a grob then fail
# Handles case when traversing DL
getGrobFromGPath.default <- function(grob, pathsofar, gPath, strict,
                                     grep, global) {
  NULL
}

getGrobFromGPath.grob <- function(grob, pathsofar, gPath, strict,
                                  grep, global) {
  if (depth(gPath) > 1)
    NULL
  else {
    if (nameMatch(gPath$name, grob$name, grep))
      grob
    else
      NULL
  }
}

getGTree <- function(gTree, pathsofar, gPath, strict, grep, global) {
  # Try to find pathsofar at start of gPath
  # NOTE: may be called directly with pathsofar=NULL
  if (is.null(pathsofar) ||
      (!strict && depth(gPath) == 1) ||
      partialPathMatch(pathsofar, gPath$path, strict, grep)) {
    found <- FALSE
    index <- 1
    grob <- NULL
    # Search children for match
    while (index <= length(gTree$childrenOrder) &&
           (!found || global)) {
      childName <- gTree$childrenOrder[index]
      child <- gTree$children[[childName]]
      # Special case when strict is FALSE and depth(gPath) is 1
      # Just check for gPath$name amongst children and recurse if no match
      if (!strict && depth(gPath) == 1) {
        if (nameMatch(gPath$name, childName, grep)) {
          grob <- growResult(grob, child)
          found <- TRUE
        } else {
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- getGrobFromGPath(child, newpathsofar,
                                                    gPath, strict,
                                                    grep, global))) {
            grob <- growResult(grob, newChild)
            found <- TRUE
          }
        }
      } else {
        # Only check for match with child if have full match with pathsofar
        # If it's a complete match, look for gPath$name amongst child
        # NOTE: may be called directly with pathsofar=NULL
        if (fullPathMatch(pathsofar, gPath, strict, grep)) {
          if (nameMatch(gPath$name, childName, grep[depth(gPath)])) {
            grob <- growResult(grob, child)
            found <- TRUE
          }
        # Otherwise recurse down child
        } else {
          # NOTE: may be called directly with pathsofar=NULL
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- getGrobFromGPath(child, newpathsofar,
                                                    gPath, strict,
                                                    grep, global))) {
            grob <- growResult(grob, newChild)
            found <- TRUE
          }
        }
      }
      index <- index + 1
    }
    if (found)
      grob
    else
      NULL
  } else {
    NULL
  }
}

getGrobFromGPath.gTree <- function(grob, pathsofar, gPath, strict,
                                   grep, global) {
  if (depth(gPath) == 1) {
    if (nameMatch(gPath$name, grob$name, grep))
      grob
    else
      if (strict)
        NULL
      else
        getGTree(grob,
                 if (is.null(pathsofar)) grob$name else pathsofar,
                 gPath, strict, grep, global)
  } else {
    getGTree(grob,
             if (is.null(pathsofar)) grob$name else pathsofar,
             gPath, strict, grep, global)
  }
}

getDLfromGPath <- function(gPath, strict, grep, global) {
  dl.index <- grid.Call(C_getDLindex)
  result <- NULL
  index <- 1
  while (index < dl.index &&
         (is.null(result) || global)) {
    grob <- getGrobFromGPath(grid.Call(C_getDLelt,
                                       as.integer(index)),
                             NULL, gPath, strict,
                             grep, global)
    if (!is.null(grob))
      result <- growResult(result, grob)
    index <- index + 1
  }
  result
}

#####
##### Set support
#####
# A gPath may specify the child of a gTree
# (or the child of a child of a gTree, or ...)
setGrobFromGPath <- function(grob, pathsofar, gPath, newGrob, strict, grep) {
  UseMethod("setGrobFromGPath")
}

# Ignore DL elements which are not grobs
setGrobFromGPath.default <- function(grob, pathsofar, gPath, newGrob,
                                     strict, grep) {
  NULL
}

setGrobFromGPath.grob <- function(grob, pathsofar, gPath, newGrob,
                                  strict, grep) {
  if (depth(gPath) > 1)
    NULL
  else {
    if (nameMatch(gPath$name, grob$name, grep))
      if (match(grob$name, newGrob$name, nomatch=0L))
        newGrob
      else
        NULL
    else
      NULL
  }
}

# Try to match gPath in gTree children
# Return NULL if cant' find match
# Return modified gTree if can find match
setGTree <- function(gTree, pathsofar, gPath, newGrob, strict, grep) {
  # Try to find pathsofar at start of gPath
  # NOTE: may be called directly with pathsofar=NULL
  if (is.null(pathsofar) ||
      (!strict && depth(gPath) == 1) ||
      partialPathMatch(pathsofar, gPath$path, strict, grep)) {
    found <- FALSE
    index <- 1
    # Search children for match
    while (index <= length(gTree$childrenOrder) && !found) {
      childName <- gTree$childrenOrder[index]
      child <- gTree$children[[childName]]
      # Special case when strict is FALSE and depth(gPath) is 1
      # Just check for gPath$name amongst children and recurse if no match
      if (!strict && depth(gPath) == 1) {
        if (nameMatch(gPath$name, childName, grep)) {
          if (match(childName, newGrob$name, nomatch=0L)) {
            gTree$children[[newGrob$name]] <- newGrob
            found <- TRUE
          } else {
            stop("the new 'grob' must have the same name as the old 'grob'")
          }
        } else {
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- setGrobFromGPath(child, newpathsofar,
                                                    gPath, newGrob,
                                                    strict, grep))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      } else {
        # Only check for match with child if have full match with pathsofar
        # If it's a complete match, look for gPath$name amongst child
        # NOTE: may be called directly with pathsofar=NULL
        if (fullPathMatch(pathsofar, gPath, strict, grep)) {
          if (nameMatch(gPath$name, childName, grep[depth(gPath)])) {
            if (match(childName, newGrob$name, nomatch=0L)) {
                gTree$children[[newGrob$name]] <- newGrob
                found <- TRUE
            } else {
                stop("the new 'grob' must have the same name as the old 'grob'")
            }
          }
        # Otherwise recurse down child
        } else {
          # NOTE: may be called directly with pathsofar=NULL
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- setGrobFromGPath(child, newpathsofar,
                                                    gPath, newGrob,
                                                    strict, grep))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      }
      index <- index + 1
    }
    if (found)
      gTree
    else
      NULL
  } else {
    NULL
  }
}

setGrobFromGPath.gTree <- function(grob, pathsofar, gPath, newGrob,
                                   strict, grep) {
  if (depth(gPath) == 1) {
    if (nameMatch(gPath$name, grob$name, grep))
      if (match(grob$name, newGrob$name, nomatch=0L))
        newGrob
      else
        stop("the new 'grob' must have the same name as the old 'grob'")
    else
      if (strict)
        NULL
      else
        setGTree(grob,
                 if (is.null(pathsofar)) grob$name else pathsofar,
                 gPath, newGrob, strict, grep)
  } else {
    setGTree(grob,
             # Initialise pathsofar if first time through
             if (is.null(pathsofar)) grob$name else pathsofar,
             gPath, newGrob, strict, grep)
  }
}

setDLfromGPath <- function(gPath, newGrob, strict, grep) {
  dl.index <- grid.Call(C_getDLindex)
  index <- 1
  result <- list(index=0, grob=NULL)
  while (index < dl.index &&
         result$index == 0) {
    result$grob <- setGrobFromGPath(grid.Call(C_getDLelt,
                                              as.integer(index)),
                                    NULL, gPath, newGrob, strict, grep)
    if (!is.null(result$grob))
      result$index <- index
    index <- index + 1
  }
  result
}

#####
##### Edit support
#####
editThisGrob <- function(grob, specs) {
  for (i in names(specs))
    if (nzchar(i))
      # Handle gp as special case
      if (match(i, "gp", nomatch=0))
        # Handle NULL as special case
        if (is.null(specs[[i]]))
          grob[i] <- list(gp=NULL)
        else
          grob$gp <- mod.gpar(grob$gp, specs$gp)
      # If there is no slot with the argument name, just ignore that argument
      else if (match(i, names(grob), nomatch=0))
        # Handle NULL as special case
        if (is.null(specs[[i]]))
          grob[i] <- eval(substitute(list(i=NULL)))
        else
          grob[[i]] <- specs[[i]]
      else
        warning(gettextf("slot '%s' not found", i), domain = NA)
  # Check grob slots are ok before trying to do anything with them
  # in editDetails
  # grob$childrenvp may be non-NULL for a gTree
  grob <- validGrob(grob, grob$childrenvp)
  editDetails(grob, specs)
}

# A gPath may specify the child of a gTree
# (or the child of a child of a gTree, or ...)
editGrobFromGPath <- function(grob, specs, pathsofar, gPath, strict,
                              grep, global) {
  UseMethod("editGrobFromGPath")
}

# If it's not a grob then fail
# Handles case when traversing DL
editGrobFromGPath.default <- function(grob, specs,
                                      pathsofar, gPath, strict,
                                      grep, global) {
  NULL
}

editGrobFromGPath.grob <- function(grob, specs,
                                   pathsofar, gPath, strict,
                                   grep, global) {
  if (depth(gPath) > 1)
    NULL
  else {
    if (nameMatch(gPath$name, grob$name, grep))
      editThisGrob(grob, specs)
    else
      NULL
  }
}

editGTree <- function(gTree, specs, pathsofar, gPath, strict,
                      grep, global) {
  # Try to find pathsofar at start of gPath
  # NOTE: may be called directly with pathsofar=NULL
  if (is.null(pathsofar) ||
      (!strict && depth(gPath) == 1) ||
      partialPathMatch(pathsofar, gPath$path, strict, grep)) {
    found <- FALSE
    index <- 1
    # Search children for match
    while (index <= length(gTree$childrenOrder) &&
           (!found || global)) {
      childName <- gTree$childrenOrder[index]
      child <- gTree$children[[childName]]
      # Special case when strict is FALSE and depth(gPath) is 1
      # Just check for gPath$name amongst children and recurse if no match
      if (!strict && depth(gPath) == 1) {
        if (nameMatch(gPath$name, childName, grep)) {
          gTree$children[[childName]] <- editThisGrob(child, specs)
          found <- TRUE
        } else {
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- editGrobFromGPath(child, specs,
                                                     newpathsofar,
                                                     gPath, strict,
                                                     grep, global))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      } else {
        # Only check for match with child if have full match with pathsofar
        # If it's a complete match, look for gPath$name amongst child
        # NOTE: may be called directly with pathsofar=NULL
        if (fullPathMatch(pathsofar, gPath, strict, grep)) {
          if (nameMatch(gPath$name, childName, grep[depth(gPath)])) {
            gTree$children[[childName]] <- editThisGrob(child, specs)
            found <- TRUE
          }
        # Otherwise recurse down child
        } else {
          # NOTE: may be called directly with pathsofar=NULL
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- editGrobFromGPath(child, specs,
                                                     newpathsofar,
                                                     gPath, strict,
                                                     grep, global))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      }
      index <- index + 1
    }
    if (found)
      gTree
    else
      NULL
  } else {
    NULL
  }
}

editGrobFromGPath.gTree <- function(grob, specs,
                                    pathsofar, gPath, strict,
                                    grep, global) {
  if (depth(gPath) == 1) {
    if (nameMatch(gPath$name, grob$name, grep))
      editThisGrob(grob, specs)
    else
      if (strict)
        NULL
      else
        editGTree(grob, specs,
                  if (is.null(pathsofar)) grob$name else pathsofar,
                  gPath, strict, grep, global)
  } else {
    editGTree(grob, specs,
              if (is.null(pathsofar)) grob$name else pathsofar,
              gPath, strict, grep, global)
  }
}

editDLfromGPath <- function(gPath, specs, strict, grep, global, redraw) {
  dl.index <- grid.Call(C_getDLindex)
  index <- 1
  grob <- NULL
  found <- FALSE
  while (index < dl.index &&
         (is.null(grob) || global)) {
    grob <- editGrobFromGPath(grid.Call(C_getDLelt,
                                        as.integer(index)),
                              specs,
                              NULL, gPath, strict, grep, global)
    if (!is.null(grob)) {
      # Destructively modify the DL elt
      grid.Call(C_setDLindex, as.integer(index))
      grid.Call(C_setDLelt, grob)
      # Reset the DL index
      grid.Call(C_setDLindex, as.integer(dl.index))
      found <- TRUE
    }
    index <- index + 1
  }
  if (!found)
    stop(gettextf("'gPath' (%s) not found", as.character(gPath)), domain = NA)
  else if (redraw)
    draw.all()
}

#####
##### Add support
#####

# Assume that child is a grob
addToGTree <- function(gTree, child) {
  if (!inherits(gTree, "gTree"))
    stop("it is only valid to add a child to a \"gTree\"")
  gTree$children[[child$name]] <- child
  # Handle case where child name already exists (so will be overwritten)
  if (old.pos <- match(child$name, gTree$childrenOrder, nomatch=0))
    gTree$childrenOrder <- gTree$childrenOrder[-old.pos]
  gTree$childrenOrder <- c(gTree$childrenOrder, child$name)
  gTree
}

# A gPath may specify the child of a gTree
# (or the child of a child of a gTree, or ...)
addGrobFromGPath <- function(grob, child, pathsofar, gPath, strict,
                             grep, global) {
  UseMethod("addGrobFromGPath")
}

# If it's not a grob then fail
# Handles case when traversing DL
addGrobFromGPath.default <- function(grob, child,
                                     pathsofar, gPath, strict,
                                     grep, global) {
  NULL
}

# If no match then fail
# If match then error!
addGrobFromGPath.grob <- function(grob, child,
                                  pathsofar, gPath, strict,
                                  grep, global) {
  if (depth(gPath) > 1)
    NULL
  else {
    if (nameMatch(gPath$name, grob$name, grep))
      stop("it is only valid to add a child to a \"gTree\"")
    else
      NULL
  }
}

# In this function, the grob being added is called "grob"
# (in all others it is called "child"
addGTree <- function(gTree, grob, pathsofar, gPath, strict,
                     grep, global) {
  # Try to find pathsofar at start of gPath
  # NOTE: may be called directly with pathsofar=NULL
  if (is.null(pathsofar) ||
      (!strict && depth(gPath) == 1) ||
      partialPathMatch(pathsofar, gPath$path, strict, grep)) {
    found <- FALSE
    index <- 1
    # Search children for match
    while (index <= length(gTree$childrenOrder) &&
           (!found || global)) {
      childName <- gTree$childrenOrder[index]
      child <- gTree$children[[childName]]
      # Special case when strict is FALSE and depth(gPath) is 1
      # Just check for gPath$name amongst children and recurse if no match
      if (!strict && depth(gPath) == 1) {
        if (nameMatch(gPath$name, childName, grep)) {
          gTree$children[[childName]] <- addToGTree(child, grob)
          found <- TRUE
        } else {
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- addGrobFromGPath(child, grob,
                                                    newpathsofar,
                                                    gPath, strict,
                                                    grep, global))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      } else {
        # Only check for match with child if have full match with pathsofar
        # If it's a complete match, look for gPath$name amongst child
        # NOTE: may be called directly with pathsofar=NULL
        if (fullPathMatch(pathsofar, gPath, strict, grep)) {
          if (nameMatch(gPath$name, childName, grep[depth(gPath)])) {
            gTree$children[[childName]] <- addToGTree(child, grob)
            found <- TRUE
          }
        # Otherwise recurse down child
        } else {
          # NOTE: may be called directly with pathsofar=NULL
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- addGrobFromGPath(child, grob,
                                                    newpathsofar,
                                                    gPath, strict,
                                                    grep, global))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      }
      index <- index + 1
    }
    if (found)
      gTree
    else
      NULL
  } else {
    NULL
  }
}

addGrobFromGPath.gTree <- function(grob, child,
                                   pathsofar, gPath, strict,
                                   grep, global) {
  if (depth(gPath) == 1) {
    if (nameMatch(gPath$name, grob$name, grep))
      addToGTree(grob, child)
    else
      if (strict)
        NULL
      else
        addGTree(grob, child,
                 if (is.null(pathsofar)) grob$name else pathsofar,
                 gPath, strict, grep, global)
  } else {
    addGTree(grob, child,
             if (is.null(pathsofar)) grob$name else pathsofar,
             gPath, strict, grep, global)
  }
}

addDLfromGPath <- function(gPath, child, strict, grep, global, redraw) {
  dl.index <- grid.Call(C_getDLindex)
  index <- 1
  grob <- NULL
  found <- FALSE
  while (index < dl.index &&
         (is.null(grob) || global)) {
    grob <- addGrobFromGPath(grid.Call(C_getDLelt,
                                       as.integer(index)),
                             child,
                             NULL, gPath, strict, grep, global)
    if (!is.null(grob)) {
      # Destructively modify the DL elt
      grid.Call(C_setDLindex, as.integer(index))
      grid.Call(C_setDLelt, grob)
      # Reset the DL index
      grid.Call(C_setDLindex, as.integer(dl.index))
      found <- TRUE
    }
    index <- index + 1
  }
  if (!found)
    stop(gettextf("'gPath' (%s) not found", gPath), domain = NA)
  else if (redraw)
    draw.all()
}

#####
##### Remove support
#####

removeFromGTree <- function(gTree, name, grep) {
  if (!inherits(gTree, "gTree"))
    stop("it is only valid to remove a child from a \"gTree\"")
  if (grep) {
    old.pos <- grep(name, gTree$childrenOrder)
    if (length(old.pos) == 0L)
      old.pos <- 0
  } else {
    old.pos <- match(name, gTree$childrenOrder, nomatch=0)
  }
  if (old.pos > 0) {
    # name might be a regexp so use real name
    gTree$children[[gTree$childrenOrder[old.pos]]] <- NULL
    gTree$childrenOrder <- gTree$childrenOrder[-old.pos]
    gTree
  } else {
    NULL
  }
}

# A gPath may specify the child of a gTree
# (or the child of a child of a gTree, or ...)
removeGrobFromGPath <- function(grob, name, pathsofar, gPath, strict,
                                grep, grepname, global, warn) {
  UseMethod("removeGrobFromGPath")
}

# If it's not a grob then fail
# Handles case when traversing DL
removeGrobFromGPath.default <- function(grob, name,
                                        pathsofar, gPath, strict,
                                        grep, grepname, global, warn) {
  NULL
}

# ALWAYS fail
# (either no match or match but grob has no children!)
removeGrobFromGPath.grob <- function(grob, name,
                                     pathsofar, gPath, strict,
                                     grep, grepname, global, warn) {
  NULL
}

removeGTree <- function(gTree, name, pathsofar, gPath, strict,
                        grep, grepname, global, warn) {
  # Try to find pathsofar at start of gPath
  # NOTE: may be called directly with pathsofar=NULL
  if (is.null(pathsofar) ||
      (!strict && depth(gPath) == 1) ||
      partialPathMatch(pathsofar, gPath$path, strict, grep)) {
    found <- FALSE
    index <- 1
    # Search children for match
    while (index <= length(gTree$childrenOrder) &&
           (!found || global)) {
      childName <- gTree$childrenOrder[index]
      child <- gTree$children[[childName]]
      # Special case when strict is FALSE and depth(gPath) is 1
      # Just check for gPath$name amongst children and recurse if no match
      if (!strict && depth(gPath) == 1) {
        # NOTE: child has to be a gTree if we hope to find a child in it!
        if (inherits(child, "gTree") &&
            nameMatch(gPath$name, childName, grep)) {
          newchild <- removeFromGTree(child, name, grepname)
          if (!is.null(newchild)) {
            gTree$children[[childName]] <- newchild
            found <- TRUE
          }
        } else {
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- removeGrobFromGPath(child, name,
                                                       newpathsofar,
                                                       gPath, strict,
                                                       grep, grepname,
                                                       global, warn))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      } else {
        # Only check for match with child if have full match with pathsofar
        # If it's a complete match, look for gPath$name amongst child
        # NOTE: may be called directly with pathsofar=NULL
        if (fullPathMatch(pathsofar, gPath, strict, grep)) {
          # NOTE: child has to be a gTree if we hope to find a child in it!
          if (inherits(child, "gTree") &&
              nameMatch(gPath$name, childName, grep[depth(gPath)])) {
            newchild <- removeFromGTree(child, name, grepname)
            if (!is.null(newchild)) {
              gTree$children[[childName]] <- newchild
              found <- TRUE
            }
          }
        # Otherwise recurse down child
        } else {
          # NOTE: may be called directly with pathsofar=NULL
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- removeGrobFromGPath(child, name,
                                                       newpathsofar,
                                                       gPath, strict,
                                                       grep, grepname,
                                                       global, warn))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      }
      index <- index + 1
    }
    if (found)
      gTree
    else
      NULL
  } else {
    NULL
  }
}

removeGrobFromGPath.gTree <- function(grob, name,
                                      pathsofar, gPath, strict,
                                      grep, grepname, global, warn) {
  if (depth(gPath) == 1) {
    if (nameMatch(gPath$name, grob$name, grep))
      removeFromGTree(grob, name, grepname)
    else
      if (strict)
        NULL
      else
        removeGTree(grob, name,
                    if (is.null(pathsofar)) grob$name else pathsofar,
                    gPath, strict, grep, grepname, global, warn)
  } else {
    removeGTree(grob, name,
                if (is.null(pathsofar)) grob$name else pathsofar,
                gPath, strict, grep, grepname, global, warn)
  }
}

removeDLFromGPath <- function(gPath, name, strict, grep, grepname, global,
                              warn, redraw) {
  dl.index <- grid.Call(C_getDLindex)
  index <- 1
  grob <- NULL
  found <- FALSE
  while (index < dl.index &&
         (is.null(grob) || global)) {
    grob <- removeGrobFromGPath(grid.Call(C_getDLelt, as.integer(index)),
                                name,
                                NULL, gPath, strict, grep, grepname,
                                global, warn)
    if (!is.null(grob)) {
      # Destructively modify the DL elt
      grid.Call(C_setDLindex, as.integer(index))
      grid.Call(C_setDLelt, grob)
      # Reset the DL index
      grid.Call(C_setDLindex, as.integer(dl.index))
      found <- TRUE
    }
    index <- index + 1
  }
  if (!found)
    stop(gettextf("gPath (%s) not found",
                  paste(gPath, name, sep=.grid.pathSep)),
                  domain = NA)
  else if (redraw)
    draw.all()
}

#####
##### Remove NAME support
#####

# NEVER called when strict=TRUE
removeGrobFromName <- function(grob, name, grep, global, warn) {
  UseMethod("removeGrobFromName")
}

removeGrobFromName.grob <- function(grob, name, grep, global, warn) {
  NULL
}

# For a gTree, just recurse straight back to removeName
removeGrobFromName.gTree <- function(grob, name, grep, global, warn) {
    removeName(grob, name, FALSE, grep, global, warn)
}

removeName <- function(gTree, name, strict, grep, global, warn) {
  found <- FALSE
  index <- 1
  # Search children for match
  while (index <= length(gTree$childrenOrder) &&
         (!found || global)) {
    childName <- gTree$childrenOrder[index]
    child <- gTree$children[[childName]]
    # Just check child name and recurse if no match
    if (nameMatch(name, childName, grep)) {
      # name might be a regexp, so get real name
      gTree$children[[gTree$childrenOrder[index]]] <- NULL
      gTree$childrenOrder <- gTree$childrenOrder[-index]
      found <- TRUE
      # If deleted the child, do NOT increase index!
    } else if (strict) {
      NULL
      index <- index + 1
    } else {
      if (!is.null(newChild <- removeGrobFromName(child, name,
                                                  grep, global, warn))) {
        gTree$children[[childName]] <- newChild
        found <- TRUE
      }
      index <- index + 1
    }
  }
  if (found)
    gTree
  else
    NULL
}

removeNameFromDL <- function(name, strict, grep, global, warn, redraw) {
  dl.index <- grid.Call(C_getDLindex)
  index <- 1
  grob <- NULL
  found <- FALSE
  while (index < dl.index &&
         (is.null(grob) || global)) {
    grob <- grid.Call(C_getDLelt, as.integer(index))
    if (inherits(grob, "grob")) {
      # If match top-level grob, remove it from DL
      if (nameMatch(name, grob$name, grep)) {
        # Destructively modify the DL elt
        grid.Call(C_setDLindex, as.integer(index))
        grid.Call(C_setDLelt, NULL)
        # Reset the DL index
        grid.Call(C_setDLindex, as.integer(dl.index))
        found <- TRUE
      # Otherwise search down it for match
      } else {
        if (!strict) {
          grob <- removeGrobFromName(grob, name, grep, global, warn)
          if (!is.null(grob)) {
            # Destructively modify the DL elt
            grid.Call(C_setDLindex, as.integer(index))
            grid.Call(C_setDLelt, grob)
            # Reset the DL index
            grid.Call(C_setDLindex, as.integer(dl.index))
            found <- TRUE
          }
        }
      }
    } else {
      grob <- NULL
    }
    index <- index + 1
  }
  if (!found) {
    if (warn)
        stop(gettextf("gPath (%s) not found", name), domain = NA)
  } else if (redraw)
    draw.all()
}

################
# Finding a grob from a grob name
################
findgrob <- function(x, name) {
  UseMethod("findgrob")
}

findgrob.default <- function(x, name) {
  NULL
}

findgrob.grob <- function(x, name) {
  if (match(name, x$name, nomatch=0L))
    x
  else
    NULL
}

findGrobinDL <- function(name) {
  dl.index <- grid.Call(C_getDLindex)
  result <- NULL
  index <- 1
  while (index < dl.index && is.null(result)) {
    result <- findgrob(grid.Call(C_getDLelt, as.integer(index)), name)
    index <- index + 1
  }
  if (is.null(result))
    stop(gettextf("grob '%s' not found", name), domain = NA)
  result
}

findGrobinChildren <- function(name, children) {
  nc <- length(children)
  result <- NULL
  index <- 1
  while (index <= nc && is.null(result)) {
    result <- findgrob(children[[index]], name)
    index <- index + 1
  }
  if (is.null(result))
    stop(gettextf("grob '%s' not found", name), domain = NA)
  result
}

################
# grid.draw
################
# Use generic function "draw" rather than generic function "print"
# because want graphics functions to produce graphics output
# without having to be evaluated at the command-line AND without having
# to necessarily produce a single graphical object as the return value
# (i.e., so that simple procedural code can be written just for its
# side-effects).
# For example, so that the following code will draw
# a rectangle AND a line:
#   temp <- function() { grid.lines(); grid.rect() }
#   temp()
grid.draw <- function(x, recording=TRUE) {
    # If 'x' is NULL, draw nothing
    if (!is.null(x))
        UseMethod("grid.draw")
}

grid.draw.viewport <- function(x, recording) {
  pushViewport(x, recording=FALSE)
}

grid.draw.vpPath <- function(x, recording) {
  # Assumes strict=FALSE, BUT in order to get onto
  # display list it must have worked => strict same as non-strict
  downViewport(x, recording=FALSE)
}

grid.draw.pop <- function(x, recording) {
  popViewport(x, recording=FALSE)
}

grid.draw.up <- function(x, recording) {
  upViewport(x, recording=FALSE)
}

pushgrobvp <- function(vp) {
  UseMethod("pushgrobvp")
}

pushgrobvp.viewport <- function(vp) {
  pushViewport(vp, recording=FALSE)
}

pushgrobvp.vpPath <- function(vp) {
  downViewport(vp, strict=TRUE, recording=FALSE)
}

popgrobvp <- function(vp) {
  UseMethod("popgrobvp")
}

popgrobvp.viewport <- function(vp) {
  # NOTE that the grob's vp may be a vpStack/List/Tree
  upViewport(depth(vp), recording=FALSE)
}

popgrobvp.vpPath <- function(vp) {
  upViewport(depth(vp), recording=FALSE)
}

preDraw <- function(x) {
  UseMethod("preDraw")
}

pushvpgp <- function(x) {
  if (!is.null(x$vp))
    pushgrobvp(x$vp)
  if (!is.null(x$gp)) {
    set.gpar(x$gp, grob=x)
  }
}

makeContext <- function(x) {
    UseMethod("makeContext")
}

makeContext.default <- function(x) {
    x
}

makeContent <- function(x) {
    UseMethod("makeContent")
}

makeContent.default <- function(x) {
    x
}

preDraw.grob <- function(x) {
    # Allow customisation of x$vp
    x <- makeContext(x)
    # automatically push/pop the viewport and set/unset the gpar
    pushvpgp(x)
    preDrawDetails(x)
    x
}

preDraw.gTree <- function(x) {
    # Allow customisation of x$vp (and x$childrenvp)
    x <- makeContext(x)
    # Make this gTree the "current grob" for evaluation of
    # grobwidth/height units via gPath
    # Do this as a .Call.graphics to get it onto the base display list
    grid.Call.graphics(C_setCurrentGrob, x)
    # automatically push/pop the viewport
    pushvpgp(x)
    # Push then "up" childrenvp
    if (!is.null(x$childrenvp)) {
        # Save any x$gp gpar settings
        tempgp <- grid.Call(C_getGPar)
        pushViewport(x$childrenvp, recording=FALSE)
        upViewport(depth(x$childrenvp), recording=FALSE)
        # reset the x$gp gpar settings
        # The upViewport above may have overwritten them with
        # the previous vp$gp settings
        grid.Call.graphics(C_setGPar, tempgp)
    }
    preDrawDetails(x)
    x
}

postDraw <- function(x) {
    UseMethod("postDraw")
}

postDraw.grob <- function(x) {
    postDrawDetails(x)
    if (!is.null(x$vp))
        popgrobvp(x$vp)
}

drawGrob <- function(x) {
    # Temporarily turn off the grid DL so that
    # nested calls to drawing code do not get recorded
    dlon <- grid.Call(C_setDLon, FALSE)
    # If get error or user-interrupt, need to reset state
    # Need to turn grid DL back on (if it was on)
    on.exit(grid.Call(C_setDLon, dlon))
    # Save current gpar
    tempgpar <- grid.Call(C_getGPar)
    # If get error or user-interrupt, need to reset state
    # Need to restore current grob (gtree predraw sets current grob)
    # Need to restore gpar settings (set by gtree itself and/or its vp)
    # This does not need to be a grid.Call.graphics() because
    # we are nested within a recordGraphics()
    # Do not call set.gpar because set.gpar accumulates cex
    on.exit(grid.Call(C_setGPar, tempgpar), add=TRUE)
    # Setting up the drawing context may involve modifying the grob
    # (typically only x$vp) but the modified grob is needed for postDraw()
    ## R_GE_DEBUG print(x)
    ## R_GE_DEBUG print(paste0("pre  preDraw: ", current.viewport()))
    x <- preDraw(x)
    ## R_GE_DEBUG print(paste0("post preDraw: ", current.viewport()))
    # Allow customisation of x
    # (should only return a basic grob that has a drawDetails()
    #  method, otherwise nothing will be drawn)
    x <- makeContent(x)
    ## For pattern fill resolution, attach the built grob to gp$fill
    recordGrobForPatternResolution(x)
    # Do any class-specific drawing
    drawDetails(x, recording=FALSE)
    ## R_GE_DEBUG print(paste0("pre  postDraw: ", current.viewport()))
    postDraw(x)
    ## R_GE_DEBUG print(paste0("post postDraw: ", current.viewport()))
}

grid.draw.grob <- function(x, recording=TRUE) {
    engineDLon <- grid.Call(C_getEngineDLon)
    if (engineDLon)
        recordGraphics(drawGrob(x),
                       list(x=x),
                       getNamespace("grid"))
    else
        drawGrob(x)
    if (recording)
        record(x)
    invisible()
}

drawGList <- function(x) {
    # DO NOT turn off grid DL.
    # A top-level gList does not itself go on the DL,
    # but its children do.
    # A gList which is part of some other grob (e.g., children
    # of a gTree) will be "protected" by the gTree
    # turning off the DL.
    lapply(x, grid.draw)
}

grid.draw.gList <- function(x, recording=TRUE) {
    engineDLon <- grid.Call(C_getEngineDLon)
    if (engineDLon)
        recordGraphics(drawGList(x),
                       list(x=x),
                       getNamespace("grid"))
    else
        drawGList(x)
    invisible()
}

drawGTree <- function(x) {
    # Temporarily turn off the grid DL so that
    # nested calls to drawing code do not get recorded
    dlon <- grid.Call(C_setDLon, FALSE)
    # If get error or user-interrupt, need to reset state
    # Need to turn grid DL back on (if it was on)
    on.exit(grid.Call(C_setDLon, dlon))
    # Save current grob and current gpar
    tempgrob <- grid.Call(C_getCurrentGrob)
    tempgpar <- grid.Call(C_getGPar)
    # If get error or user-interrupt, need to reset state
    # Need to restore current grob (gtree predraw sets current grob)
    # Need to restore gpar settings (set by gtree itself and/or its vp)
    # This does not need to be a grid.Call.graphics() because
    # we are nested within a recordGraphics()
    # Do not call set.gpar because set.gpar accumulates cex
    on.exit({ grid.Call(C_setGPar, tempgpar)
              grid.Call(C_setCurrentGrob, tempgrob)
            }, add=TRUE)
    # Setting up the drawing context may involve modifying the grob
    # (typically only x$vp) but the modified grob is needed for postDraw()
    ## R_GE_DEBUG print(x)
    ## R_GE_DEBUG print(paste0("pre  preDraw: ", current.viewport()))
    x <- preDraw(x)
    ## R_GE_DEBUG print(paste0("post preDraw: ", current.viewport()))
    # Allow customisation of x (should be confined to x$children)
    x <- makeContent(x)
    ## For pattern fill resolution, attach the built grob to gp$fill
    recordGTreeForPatternResolution(x)
    # Do any class-specific drawing
    drawDetails(x, recording=FALSE)
    # Draw all children IN THE RIGHT ORDER
    for (i in x$childrenOrder)
      grid.draw(x$children[[i]], recording=FALSE)
    ## R_GE_DEBUG print(paste0("pre  postDraw: ", current.viewport()))
    postDraw(x)
    ## R_GE_DEBUG print(paste0("post postDraw: ", current.viewport()))
}

grid.draw.gTree <- function(x, recording=TRUE) {
    engineDLon <- grid.Call(C_getEngineDLon)
    if (engineDLon)
        recordGraphics(drawGTree(x),
                       list(x=x),
                       getNamespace("grid"))
    else
        drawGTree(x)
    if (recording)
        record(x)
    invisible()
}

draw.all <- function() {
    grid.newpage(recording=FALSE)
    dl.index <- grid.Call(C_getDLindex)
    if (dl.index > 1)
        # Start at 2 because first element is viewport[ROOT]
        for (i in 2:dl.index) {
            grid.draw(grid.Call(C_getDLelt, as.integer(i - 1)),
                      recording=FALSE)
        }
}

draw.details <- function(x, recording) {
    .Defunct("drawDetails")
}

preDrawDetails <- function(x) {
    UseMethod("preDrawDetails")
}

preDrawDetails.grob <- function(x) {
}

postDrawDetails <- function(x) {
    UseMethod("postDrawDetails")
}

postDrawDetails.grob <- function(x) {
}

drawDetails <- function(x, recording) {
    UseMethod("drawDetails")
}

drawDetails.grob <- function(x, recording) {
}

grid.copy <- function(grob) {
    warning("this function is redundant and will disappear in future versions",
            domain = NA)
    grob
}

################################
# Flattening a grob

forceGrob <- function(x) {
    UseMethod("forceGrob")
}

# The default action is to leave 'x' untouched
# BUT it is also necessary to enforce the drawing context
# for viewports and vpPaths
forceGrob.default <- function(x) {
    grid.draw(x, recording=FALSE)
    x
}

# This allows 'x' to be modified, but may not
# change 'x' at all
forceGrob.grob <- function(x) {
    # Copy of the original object to allow a "revert"
    originalX <- x
    # Same set up as drawGrob()
    dlon <- grid.Call(C_setDLon, FALSE)
    on.exit(grid.Call(C_setDLon, dlon))
    tempgpar <- grid.Call(C_getGPar)
    on.exit(grid.Call(C_setGPar, tempgpar), add=TRUE)
    # Same drawing context set up as drawGrob()
    # including enforcing the drawing context
    x <- preDraw(x)
    # Same drawing content set up as drawGrob() ...
    x <- makeContent(x)
    # BUT NO DRAWING
    # Same context clean up as drawGrob()
    postDraw(x)
    # If 'x' has not changed, just return original 'x'
    # Also, do not bother with saving original
    # If 'x' has changed ...
    if (!identical(x, originalX)) {
        # Store the original object to allow a "revert"
        x$.ORIGINAL <- originalX
        # Return the 'x' that would have been drawn
        # This will typically be a standard R primitive
        # (which do not have makeContext() or makeContent()
        #  methods, only drawDetails())
        # BUT ot be safe add "forcedgrob" class so that subsequent
        # draws will NOT run makeContext() or makeContent()
        # methods
        class(x) <- c("forcedgrob", class(x))
    }
    x
}

# This allows 'x' to be modified, but may not
# change 'x' at all
forceGrob.gTree <- function(x) {
    # Copy of the original object to allow a "revert"
    originalX <- x
    # Same set up as drawGTree()
    dlon <- grid.Call(C_setDLon, FALSE)
    on.exit(grid.Call(C_setDLon, dlon))
    tempgrob <- grid.Call(C_getCurrentGrob)
    tempgpar <- grid.Call(C_getGPar)
    on.exit({ grid.Call(C_setGPar, tempgpar)
              grid.Call(C_setCurrentGrob, tempgrob)
            }, add=TRUE)
    # Same drawing context set up as drawGTree(),
    # including enforcing the drawing context
    x <- preDraw(x)
    # Same drawing content set up as drawGTree() ...
    x <- makeContent(x)
    # Ensure that children are also forced
    x$children <- do.call("gList", lapply(x$children, forceGrob))
    # BUT NO DRAWING
    # Same context clean up as drawGTree()
    postDraw(x)
    # If 'x' has changed ...
    if (!identical(x, originalX)) {
        # Store the original object to allow a "revert"
        x$.ORIGINAL <- originalX
        # Return the 'x' that would have been drawn
        # This will typically be a vanilla gTree with children to draw
        # (which will not have makeContext() or makeContent() methods)
        # BUT to be safe add "forcedgrob" class so that subsequent
        # draws will NOT run makeContext() or makeContent()
        # methods
        class(x) <- c("forcedgrob", class(x))
    }
    x
}

# A "forcedgrob" does NOT modify context or content at
# drawing time
makeContext.forcedgrob <- function(x) x

makeContent.forcedgrob <- function(x) x

grid.force <- function(x, ...) {
    UseMethod("grid.force")
}

grid.force.default <- function(x, redraw = FALSE, ...) {
    if (!missing(x))
        stop("Invalid force target")
    # Must upViewport(0) otherwise you risk running the display
    # list from something other than the ROOT viewport
    oldcontext <- upViewport(0, recording=FALSE)
    dl.index <- grid.Call(C_getDLindex)
    if (dl.index > 1) {
        # Start at 2 because first element is viewport[ROOT]
        for (i in 2:dl.index) {
            grid.Call(C_setDLindex, as.integer(i - 1))
            grid.Call(C_setDLelt,
                      forceGrob(grid.Call(C_getDLelt, as.integer(i - 1))))
        }
        grid.Call(C_setDLindex, dl.index)
    }
    if (redraw) {
        draw.all()
    }
    # Try to go back to original context
    if (length(oldcontext)) {
        seekViewport(oldcontext, recording=FALSE)
    }
}

grid.force.grob <- function(x, draw = FALSE, ...) {
    fx <- forceGrob(x)
    if (draw)
        grid.draw(fx)
    fx
}

grid.force.character <- function(x, ...) {
    grid.force(gPath(x), ...)
}

grid.force.gPath <- function(x,
                             strict=FALSE, grep=FALSE, global=FALSE,
                             redraw = FALSE, ...) {
    # Use viewports=TRUE so that get vpPaths in result
    paths <- grid.grep(x, viewports = TRUE,
                       strict = strict, grep = grep, global = global)
    f <- function(path, ...) {
        # Only force grobs or gTrees
        # (might have vpPaths because we said grid.grep(viewports=TRUE))
        if (!inherits(path, "gPath")) return()
        target <- grid.get(path, strict=TRUE)
        vpPath <- attr(path, "vpPath")
        depth <- 0
        if (nchar(vpPath))
            depth <- downViewport(vpPath, recording=FALSE)
        forcedgrob <- forceGrob(target, ...)
        if (depth > 0)
            upViewport(depth, recording=FALSE)
        grid.set(path, strict=TRUE, forcedgrob)
    }
    if (length(paths)) {
        # To get the force happening in the correct context ...
        oldcontext <- upViewport(0, recording=FALSE)
        if (global) {
            lapply(paths, f, ...)
        } else {
            f(paths, ...)
        }
        if (redraw) {
            draw.all()
        }
        # Try to go back to original context
        if (length(oldcontext))
            seekViewport(oldcontext, recording=FALSE)
    }
    invisible()
}

revert <- function(x) {
    UseMethod("revert")
}

revert.default <- function(x) {
    x
}

# Only need to revert "forcedgrob"s
revert.forcedgrob <- function(x) {
    x$.ORIGINAL
}

# No need for recursion for gTree because if top-level grob
# changed its children then top-level grob will have retained
# revert version of its entire self (including children)

# NOTE that things will get much trickier if allow
# grid.revert(gPath = ...)
grid.revert <- function(x, ...) {
    UseMethod("grid.revert")
}

grid.revert.default <- function(x, redraw=FALSE, ...) {
    if (!missing(x))
        stop("Invalid revert target")
    dl.index <- grid.Call(C_getDLindex)
    if (dl.index > 1) {
        # Start at 2 because first element is viewport[ROOT]
        for (i in 2:dl.index) {
            grid.Call(C_setDLindex, as.integer(i - 1))
            grid.Call(C_setDLelt,
                      revert(grid.Call(C_getDLelt, as.integer(i - 1))))
        }
        grid.Call(C_setDLindex, dl.index)
    }
    if (redraw) {
        draw.all()
    }
}

grid.revert.grob <- function(x, draw=FALSE, ...) {
    rx <- revert(x)
    if (draw) {
        grid.draw(x)
    }
    rx
}

grid.revert.character <- function(x, ...) {
    grid.revert(gPath(x), ...)
}

grid.revert.gPath <- function(x,
                              strict=FALSE, grep=FALSE, global=FALSE,
                              redraw = FALSE, ...) {
    paths <- grid.grep(x, strict = strict, grep = grep, global = global)
    f <- function(path, ...) {
        grid.set(path, strict=TRUE,
                 revert(grid.get(path, strict=TRUE), ...))
    }
    if (length(paths)) {
        if (global) {
            lapply(paths, f, ...)
        } else {
            f(paths, ...)
        }
        if (redraw) {
            draw.all()
        }
    }
    invisible()
}

###############################
# Reordering grobs

# Reorder the children of a gTree
# Order may be specified as a character vector
#   Character vector MUST name existing children
# Order may be specified as a numeric vector
#   (which makes it easy to say something like
#    "make last child the first child")
#   Numeric vector MUST be within range 1:numChildren
# Only unique order values used
# Any children NOT specified by order are appended to
#   front or back of order (depending on 'front' argument)
# Order is ALWAYS back-to-front
reorderGrob <- function(x, order, back=TRUE) {
    if (!inherits(x, "gTree"))
        stop("can only reorder 'children' for a \"gTree\"")
    order <- unique(order)
    oldOrder <- x$childrenOrder
    N <- length(oldOrder)
    if (is.character(order)) {
        # Convert to numeric
        order <- match(order, x$childrenOrder)
    }
    if (is.numeric(order)) {
        if (any(!is.finite(order)) ||
            !(all(order %in% 1:N))) {
            stop("Invalid 'order'")
        }
        if (back) {
            newOrder <- c(x$childrenOrder[order],
                          x$childrenOrder[-order])
        } else {
            newOrder <- c(x$childrenOrder[-order],
                          x$childrenOrder[order])
        }
    }
    x$childrenOrder <- newOrder
    x
}

# Reorder the children of a gTree on the display list
# (identified by a gPath)
# NOTE that it is possible for this operation to produce a grob
# that no longer draws (because it relies on another grob that
# used to be drawn before it, e.g., when the width of grob "b"
# is calculated from the width of grob "a")
# Do NOT allow reordering of grobs on the display list
# (it is not even clear what should happen in terms of reordering
#  grobs mixed with viewports PLUS the potential for ending up with
#  something that will not draw is pretty high)
# IF you want to reorder the grobs on the DL, do a grid.grab()
# first and then reorder the children of the resulting gTree
grid.reorder <- function(gPath, order, back=TRUE, grep=FALSE, redraw=TRUE) {
    grob <- grid.get(gPath, grep=grep)
    grid.set(gPath, reorderGrob(grob, order, back=back),
             grep=grep, redraw=redraw)
}


## Composite one or more grobs, using a compositing operator,
## THEN draw the result
## (normal drawing just composites each grob [in fact each shape!] with OVER)

##########################
## Support functions

finaliseGroup <- function(x) {
    source <- function() {
        cvp <- current.viewport()
        ## Push a viewport (with the current viewport layout and scales)
        ## WITH NO MASK to ensure that a group begins with no soft mask.
        ##
        ## "Before execution of the transparency group XObject's
        ## content stream, the current blend mode in the graphics
        ## state is initialized to Normal, the current stroking and
        ## nonstroking alpha constants to 1.0, and the current soft
        ## mask to None."
        ##
        ## https://www.adobe.com/content/dam/acom/en/devnet/pdf/pdfs/pdf_reference_archives/PDFReference.pdf">PDF</a>
        ##
        ## Current viewport layout and scales are preserved so that
        ## locations and dimensions of viewports and grobs in x$src
        ## retain their original meaning.
        ##
        ## Justification of the current viewport must also be preserved
        ## so that transformation on group use is calculated correctly.
        hjust <- resolveHJust(cvp$just, cvp$hjust)
        vjust <- resolveVJust(cvp$just, cvp$vjust)
        pushViewport(viewport(hjust, vjust,
                              just=c(hjust, vjust),
                              mask="none",
                              layout=cvp$layout,
                              xscale=cvp$xscale, yscale=cvp$yscale),
                     recording=FALSE)
        grid.draw(x$src, recording=FALSE)
        popViewport(recording=FALSE)
    }
    if (is.grob(x$dst)) {
        destination <- function() {
            cvp <- current.viewport()
            hjust <- resolveHJust(cvp$just, cvp$hjust)
            vjust <- resolveVJust(cvp$just, cvp$vjust)
            pushViewport(viewport(hjust, vjust,
                                  just=c(hjust, vjust),
                                  mask="none",
                                  layout=cvp$layout,
                                  xscale=cvp$xscale, yscale=cvp$yscale),
                         recording=FALSE)
            grid.draw(x$dst, recording=FALSE)
            popViewport(recording=FALSE)
        }
    } else { ## NULL ("transparent") destination
        destination <- NULL
    }
    list(src=source, op=x$op, dst=destination)
}

## group mappings are stored in 'grid' state in component 18 (17 zero-based)
## (see grid.h)
groupIndex <- 18

## Record group definition (in 'grid' state)
recordGroup <- function(x, ref) {
    devState <- get(".GRID.STATE", envir=.GridEvalEnv)[[dev.cur() - 1]]
    devStateGroups <- devState[[groupIndex]]
    cvp <- current.viewport()
    ct <- current.transform()
    if (x$coords) {
        ## Points of group definition, in inches, on the device
        pts <- groupPoints(x, TRUE)
        closedPoints <- transformCoords(pts, ct)
        pts <- groupPoints(x, FALSE)
        openPoints <- transformCoords(pts, ct)
    } else {
        closedPoints <- emptyGrobCoords(x$name)
        openPoints <- emptyGrobCoords(x$name)
    }
    group <- list(ref=ref,
                  ## Record location, size, angle for re-use in
                  ## different viewport
                  xy=deviceLoc(unit(resolveHJust(cvp$just, cvp$hjust), "npc"),
                               unit(resolveVJust(cvp$just, cvp$vjust), "npc"),
                               valueOnly=TRUE, device=TRUE),
                  xyin=deviceLoc(unit(resolveHJust(cvp$just, cvp$hjust), "npc"),
                                 unit(resolveVJust(cvp$just, cvp$vjust), "npc"),
                                 valueOnly=TRUE, device=FALSE),
                  wh=c(convertX(unit(1, "npc"), "in", valueOnly=TRUE),
                       convertY(unit(1, "npc"), "in", valueOnly=TRUE)),
                  r=current.rotation(),
                  closedPoints=closedPoints,
                  openPoints=openPoints)
    if (is.null(devStateGroups)) {
        grps <- list(group)
        names(grps) <- x$name
        devStateGroups <- grps
    } else {
        devStateGroups[[x$name]] <- group
    }
    grid.Call.graphics(C_setGridState,
                       as.integer(groupIndex - 1),
                       devStateGroups)
}

lookupGroup <- function(name) {
    ## Resolve group definition to ref (using 'grid' state)
    devState <- get(".GRID.STATE", envir=.GridEvalEnv)[[dev.cur() - 1]]
    devState[[groupIndex]][[name]]
}

##########################
## Transforms (for grid.use())

groupTranslate <- function(dx=0, dy=0) {
    translate <- diag(3)
    translate[3, 1] <- dx
    translate[3, 2] <- dy
    translate
}

defnTranslate <- function(group, inverse=FALSE, device=TRUE) {
    if (device) {
        xy <- group$xy
    } else {
        xy <- group$xyin
    }
    if (inverse) {
        groupTranslate(-xy$x, -xy$y)
    } else {
        groupTranslate(xy$x, xy$y)
    }
}

useTranslate <- function(inverse=FALSE, device=TRUE) {
    cvp <- current.viewport()
    xy <- deviceLoc(unit(resolveHJust(cvp$just, cvp$hjust), "npc"),
                    unit(resolveVJust(cvp$just, cvp$vjust), "npc"),
                    valueOnly=TRUE, device=device)
    if (inverse) {
        groupTranslate(-xy$x, -xy$y)
    } else {
        groupTranslate(xy$x, xy$y)
    }
}

viewportTranslate <- function(group, device=TRUE) {
    defnTranslate(group, inverse=TRUE, device=device) %*%
        useTranslate(device=device)
}

groupRotate <- function(r=0, device=TRUE) {
    ## Account for devices that have origin at top-left
    if (device && !.devUp()) r <- -r
    rotate <- diag(3)
    theta <- r/180*pi
    costheta <- cos(theta)
    sintheta <- sin(theta)
    rotate[1, 1] <- costheta
    rotate[1, 2] <- sintheta
    rotate[2, 1] <- -sintheta
    rotate[2, 2] <- costheta
    rotate    
}

defnRotate <- function(group, inverse=FALSE, device=TRUE) {
    if (inverse) {
        groupRotate(-group$r, device)
    } else {
        groupRotate(group$r, device)
    }
}

useRotate <- function(inverse=FALSE, device=TRUE) {
    r <- current.rotation()
    if (inverse) {
        groupRotate(-r, device)
    } else {
        groupRotate(r, device)
    }    
}

viewportRotate <- function(group, device=TRUE) {
    defnTranslate(group, inverse=TRUE, device=device) %*%
        defnRotate(group, inverse=TRUE, device=device) %*%
        useRotate(device=device) %*%
        ## NOTE: NOT useTranslate() because we are ONLY rotating
        defnTranslate(group, device=device)
}

groupScale <- function(sx=1, sy=1) {
    scale <- diag(3)
    scale[1, 1] <- sx
    scale[2, 2] <- sy
    scale    
}

defnScale <- function(group, inverse=FALSE) {
    if (inverse) {
        groupScale(1/group$wh[1], 1/group$wh[2])
    } else {
        groupScale(group$wh[1], group$wh[2])
    }
}

useScale <- function(inverse=FALSE) {
    wh <- c(convertX(unit(1, "npc"), "in", valueOnly=TRUE),
            convertY(unit(1, "npc"), "in", valueOnly=TRUE))
    if (inverse) {
        groupScale(1/wh[1], 1/wh[2])
    } else {
        groupScale(wh[1], wh[2])
    }
}

viewportScale <- function(group, device=TRUE) {
    defnTranslate(group, inverse=TRUE, device=device) %*%
        defnScale(group, inverse=TRUE) %*%
        useScale() %*%
        ## NOTE: NOT useTranslate() because we are ONLY scaling
        defnTranslate(group, device=device)
}

groupShear <- function(sx=0, sy=0) {
    shear <- diag(3)
    shear[1, 2] <- sy
    shear[2, 1] <- sx
    shear        
}

groupFlip <- function(flipX=FALSE, flipY=FALSE) {
    flip <- diag(3)
    if (flipX)
        flip[1, 1] <- -1
    if (flipY)
        flip[2, 2] <- -1
    flip
}

viewportTransform <- function(group,
                              shear=groupShear(),
                              flip=groupFlip(),
                              device=TRUE) {
    ## Account for devices that have origin at top-left
    if (device && !.devUp()) {
        shear[1, 2] <- -shear[1, 2]
        shear[2, 1] <- -shear[2, 1]
    }
    defnTranslate(group, inverse=TRUE, device=device) %*%
        defnRotate(group, inverse=TRUE, device=device) %*%
        defnScale(group, inverse=TRUE) %*%
        flip %*%
        useScale() %*%
        shear %*%
        useRotate(device=device) %*%
        useTranslate(device=device)
}

##########################
## Simple interface:  define and use group in one
drawDetails.GridGroup <- function(x, recording) {
    grp <- finaliseGroup(x)
    ref <- .defineGroup(grp$src, grp$op, grp$dst)
    ## Record group to allow later reuse
    recordGroup(x, ref)
    if (is.null(ref))
        warning("Group definition failed")
    else 
        .useGroup(ref, NULL)
}

groupGrob <- function(src,
                      op = "over",
                      dst = NULL,
                      coords = TRUE,
                      name = NULL, gp=gpar(), vp=NULL) {
    
    if (!is.grob(src))
        stop("Invalid source")
    if (!(is.grob(dst) || is.null(dst)))
        stop("Invalid destination")
    ## Check valid 'op'
    .opIndex(op)
    group <- gTree(src=src, op=op, dst=dst, coords=coords,
                   name=name, gp=gp, vp=vp, cl="GridGroup")
    group
}

grid.group <- function(src,
                       op = "over",
                       dst = NULL,
                       coords = TRUE, 
                       name = NULL, gp=gpar(), vp=NULL) {
    grid.draw(groupGrob(src, op, dst, coords, name, gp, vp))
}

##########################
## More complex interface:  separate define group and use group

drawDetails.GridDefine <- function(x, recording) {
    group <- finaliseGroup(x)
    ref <- .defineGroup(group$src, group$op, group$dst)
    recordGroup(x, ref)
}
    
defineGrob <- function(src,
                       op = "over",
                       dst = NULL,
                       coords = TRUE,
                       name = NULL, gp=gpar(), vp=NULL) {
    if (!is.grob(src))
        stop("Invalid source")
    if (!(is.grob(dst) || is.null(dst)))
        stop("Invalid destination")
    ## Check valid 'op'
    .opIndex(op)
    group <- gTree(src=src, op=op, dst=dst, coords=coords,
                   name=name, gp=gp, vp=vp, cl="GridDefine")
    group
}

grid.define <- function(src,
                        op = "over",
                        dst = NULL,
                        coords = TRUE,
                        name = NULL, gp=gpar(), vp=NULL) {
    grid.draw(defineGrob(src, op, dst, coords, name, gp, vp))
}

drawDetails.GridUse <- function(x, recording) {
    group <- lookupGroup(x$group)
    if (is.null(group))
        warning(paste0("Unknown group: ", x$group))
    else {
        transform <- x$transform(group, device=TRUE)
        if (!is.matrix(transform) ||
            !is.numeric(transform) ||
            !all(dim(transform) == 3) ||
            transform[1, 3] != 0 ||
            transform[2, 3] != 0 ||
            transform[3, 3] != 1) {
            warning("Invalid transform (nothing drawn)")
            return()
        }
        .useGroup(group$ref, transform)
    }
}

useGrob <- function(group, transform=viewportTransform,
                    name=NULL, gp=gpar(), vp=NULL) {
    if (!is.function(transform))
        stop("Invalid transform")
    grp <- gTree(group=as.character(group), transform=transform,
                 name=name, gp=gp, vp=vp, cl="GridUse")
    grp
}

grid.use <- function(group, transform=viewportTransform,
                     name=NULL, gp=gpar(), vp=NULL) {
    grid.draw(useGrob(group, transform, name, gp, vp))
}
                        
################################
## Other grob methods

groupPoints <- function(x, closed, ...) {
    if (is.null(x$dst))
        children <- gList(x$src)
    else
        children <- gList(x$src, x$dst)
    grobPoints(gTree(children=children), closed, ...)
}
    
grobCoords.GridGroup <- function(x, closed, ...) {
    if (is.null(x$dst))
        children <- gList(x$src)
    else
        children <- gList(x$src, x$dst)
    grobCoords(gTree(children=children, gp=x$gp, vp=x$vp),
               closed, ...)
}

## NOTE that we still create a gTree so that grobPoints.gTree(),
## via grobPoints.gList(), will still call grobCoords() on the
## "child" src and dst
grobPoints.GridGroup <- function(x, closed, ...) {
    groupPoints(x, closed, ...)
}

## A group definition does not draw anything BUT coords
## may be needed to resolve pattern fill
## AND group definition may be a sibling of a group use
## in which case the group use needs definition coordinates
## BY DEFAULT, coordinates are recorded (but that can be turned
## of by setting x$coords=FALSE) but NOT returned (but that
## can be turned on by setting return=TRUE)
definePoints <- function(x, closed, return=TRUE, ...) {
    ## Record definition so GridUse can find it,
    ## including calculating coordinates
    ## (no actual drawing required)
    recordGroup(x, "ignored")
    if (return) {
        ## Access the stored coordinates
        ## (for the case when the definition itself needs them,
        ##  e.g., for bbox to resolve pattern fill)
        group <- lookupGroup(x$name)
        if (closed)
            pts <- group$closedPoints
        else
            pts <- group$openPoints
        ## Coordinates recorded relative to device, so transform back to
        ## relative to current viewport
        ct <- current.transform()
        transformCoords(pts, solve(ct))
    } else {
        emptyGrobCoords(x$name)
    }
}

grobCoords.GridDefineChild <- function(x, closed, return=FALSE, ...) {
    definePoints(x, closed, return, ...)
}

grobCoords.GridDefine <- function(x, closed, return=FALSE, ...) {
    ## Create gTree to get automatic enforcement of 'gp' and 'vp'
    ## BUT avoid infinite loop by setting temporary class
    ## on GridUse grob
    class(x) <- c("GridDefineChild", class(x))
    grobCoords(gTree(children=gList(x), gp=x$gp, vp=x$vp),
               closed, return, ...)
}

grobPoints.GridDefine <- function(x, closed, return=FALSE, ...) {
    definePoints(x, closed, return, ...)
}

## A group use retrieves points from group definition
## and applies transformation to them
transformCoords <- function(coords, transform) {
    if (isEmptyCoords(coords)) return(coords)
    UseMethod("transformCoords")
}

transformCoords.GridCoords <- function(coords, transform) {
    new <- cbind(coords$x, coords$y, 1) %*% transform
    gridCoords(new[, 1], new[, 2])
}

transformCoords.GridGrobCoords <- function(coords, transform) {
    new <- lapply(coords, transformCoords, transform)
    gridGrobCoords(new, attr(coords, "name"), attr(coords, "rule"))
}

transformCoords.GridGTreeCoords <- function(coords, transform) {
    new <- lapply(coords, transformCoords, transform)
    gridGTreeCoords(new, attr(coords, "name"))
}

usePoints <- function(x, closed, ...) {
    group <- lookupGroup(x$group)
    if (is.null(group)) {
        warning(paste0("Unknown group: ", x$group))
        emptyGrobCoords(x$name)
    } else {
        transform <- x$transform(group, device=FALSE)
        if (!is.matrix(transform) ||
            !is.numeric(transform) ||
            !all(dim(transform) == 3) ||
            transform[1, 3] != 0 ||
            transform[2, 3] != 0 ||
            transform[3, 3] != 1) {
            warning("Invalid transform")
            emptyGrobCoords(x$name)
        }
        ## Apply inverse of viewport transform
        ## (because grobCoords() are supposed to be relative to
        ##  current viewport NOT device)
        transform <- transform %*% solve(current.transform())
        if (closed) {
            transformCoords(group$closedPoints, transform)
        } else {
            transformCoords(group$openPoints, transform)
        }
    }
}

grobCoords.GridUseChild <- function(x, closed, ...) {
    usePoints(x, closed, ...)
}

grobCoords.GridUse <- function(x, closed, ...) {
    ## Create gTree to get automatic enforcement of 'gp' and 'vp'
    ## BUT avoid infinite loop by setting temporary class
    ## on GridUse grob
    class(x) <- c("GridUseChild", class(x))
    grobCoords(gTree(children=gList(x), gp=x$gp, vp=x$vp),
               closed, ...)
}

grobPoints.GridUse <- function(x, closed, ...) {
    usePoints(x, closed, ...)
}
#  File src/library/grid/R/highlevel.R
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

######################################
## Example applications of grid    #
######################################

grid.strip <- function(label="whatever", range.full=c(0, 1),
                   range.thumb=c(.3, .6),
                   fill="#FFBF00", thumb="#FF8000",
                   vp=NULL) {
  diff.full <- diff(range.full)
  diff.thumb <- diff(range.thumb)
  if (!is.null(vp))
    pushViewport(vp)
  grid.rect(gp=gpar(col=NULL, fill=fill))
  grid.rect((range.thumb[1L] - range.full[1L])/diff.full, 0,
            diff.thumb/diff.full, 1,
            just=c("left", "bottom"),
            gp=gpar(col=NULL, fill=thumb))
  grid.text(as.character(label))
  if (!is.null(vp))
    popViewport()
}

grid.panel <- function(x = stats::runif(10), y = stats::runif(10),
                   zrange = c(0, 1), zbin = stats::runif(2),
                   xscale = extendrange(x),
                   yscale = extendrange(y),
                   axis.left = TRUE, axis.left.label = TRUE,
                   axis.right = FALSE, axis.right.label = TRUE,
                   axis.bottom = TRUE, axis.bottom.label = TRUE,
                   axis.top = FALSE, axis.top.label = TRUE,
                   vp=NULL) {
  if (!is.null(vp))
    pushViewport(vp)
  temp.vp <- viewport(layout=grid.layout(2, 1,
                         heights=unit(c(1, 1), c("lines", "null"))))
  pushViewport(temp.vp)
  strip.vp <- viewport(layout.pos.row=1, layout.pos.col=1,
                        xscale=xscale)
  pushViewport(strip.vp)
  grid.strip(range.full=zrange, range.thumb=zbin)
  grid.rect()
  if (axis.top)
    grid.xaxis(main=FALSE, label=axis.top.label)
  popViewport()
  plot.vp <- viewport(layout.pos.row=2, layout.pos.col=1,
                       xscale=xscale, yscale=yscale)
  pushViewport(plot.vp)
  grid.grill()
  grid.points(x, y, gp=gpar(col="blue"))
  grid.rect()
  if (axis.left)
    grid.yaxis(label=axis.left.label)
  if (axis.right)
    grid.yaxis(main=FALSE, label=axis.right.label)
  if (axis.bottom)
    grid.xaxis(label=axis.bottom.label)
  popViewport(2)
  if (!is.null(vp))
    popViewport()
  invisible(list(strip.vp = strip.vp, plot.vp = plot.vp))
}

grid.multipanel <- function(x = stats::runif(90), y = stats::runif(90),
                            z = stats::runif(90),
                            nplots = 9, nrow = 5, ncol = 2,
                            newpage = TRUE, vp = NULL)
{
    if (newpage)
        grid.newpage()
    if (!is.null(vp))
        pushViewport(vp)
    stopifnot(nplots >= 1)
    if((missing(nrow) || missing(ncol)) && !missing(nplots)) {
        ## determine 'smart' default ones
        rowcol <- grDevices::n2mfrow(nplots)
        nrow <- rowcol[1L]
        ncol <- rowcol[2L]
    }
    temp.vp <- viewport(layout = grid.layout(nrow, ncol))
    pushViewport(temp.vp)
    xscale <- extendrange(x)
    yscale <- extendrange(y)
    breaks <- seq.int(min(z), max(z), length.out = nplots + 1)
    for (i in 1L:nplots) {
        col <- (i - 1) %% ncol + 1
        row <- (i - 1) %/% ncol + 1
        panel.vp <- viewport(layout.pos.row = row,
                             layout.pos.col = col)
        panelx <- x[z >= breaks[i] & z <= breaks[i+1]]
        panely <- y[z >= breaks[i] & z <= breaks[i+1]]
        grid.panel(panelx, panely, range(z), c(breaks[i], breaks[i+1]),
                   xscale, yscale,
                   axis.left = (col == 1),
                   axis.right = (col == ncol || i == nplots),
                   axis.bottom = (row == nrow),
                   axis.top = (row == 1),
                   axis.left.label = is.even(row),
                   axis.right.label = is.odd(row),
                   axis.bottom.label = is.even(col),
                   axis.top.label = is.odd(col),
                   vp = panel.vp)
    }
    grid.text("Compression Ratio", unit(.5, "npc"), unit(-4, "lines"),
              gp = gpar(fontsize = 20),
              just = "center", rot = 0)
    grid.text("NOx (micrograms/J)", unit(-4, "lines"), unit(.5, "npc"),
              gp = gpar(fontsize = 20),
              just = "centre", rot = 90)
    popViewport()
    if (!is.null(vp))
        popViewport()
}

grid.show.layout <- function(l, newpage=TRUE, vp.ex=0.8,
                             bg="light grey",
                             cell.border="blue", cell.fill="light blue",
                             cell.label=TRUE, label.col="blue",
                             unit.col="red", vp=NULL, ...) {
  if (!is.layout(l))
    stop("'l' must be a layout")
  if (newpage)
    grid.newpage()
  if (!is.null(vp))
    pushViewport(vp)
  grid.rect(gp=gpar(col=NULL, fill=bg))
  vp.mid <- viewport(0.5, 0.5, vp.ex, vp.ex, layout=l)
  pushViewport(vp.mid)
  grid.rect(gp=gpar(fill="white"))
  gp.red <- gpar(col=unit.col)
  for (i in 1L:l$nrow)
    for (j in 1L:l$ncol) {
      vp.inner <- viewport(layout.pos.row=i, layout.pos.col=j)
      pushViewport(vp.inner)
      grid.rect(gp=gpar(col=cell.border, fill=cell.fill))
      if (cell.label)
        grid.text(paste0("(", i, ", ", j, ")"), gp=gpar(col=label.col))
      if (j==1)
        # recycle heights if necessary
        grid.text(format("["(l$heights, i, top=FALSE), ...), gp=gp.red,
              just=c("right", "centre"),
              x=unit(-.05, "inches"), y=unit(.5, "npc"), rot=0)
      if (i==l$nrow)
        # recycle widths if necessary
        grid.text(format("["(l$widths, j, top=FALSE), ...), gp=gp.red,
              just=c("centre", "top"),
              x=unit(.5, "npc"), y=unit(-.05, "inches"), rot=0)
      if (j==l$ncol)
        # recycle heights if necessary
        grid.text(format("["(l$heights, i, top=FALSE), ...), gp=gp.red,
              just=c("left", "centre"),
              x=unit(1, "npc") + unit(.05, "inches"), y=unit(.5, "npc"),
              rot=0)
      if (i==1)
        # recycle widths if necessary
        grid.text(format("["(l$widths, j, top=FALSE), ...), gp=gp.red,
              just=c("centre", "bottom"),
              x=unit(.5, "npc"), y=unit(1, "npc") + unit(.05, "inches"),
              rot=0)
      popViewport()
    }
  popViewport()
  if (!is.null(vp))
    popViewport()
  ## return the viewport used to represent the parent viewport
  invisible(vp.mid)
}

grid.show.viewport <- function(v, parent.layout=NULL, newpage=TRUE, vp.ex=0.8,
                               border.fill="light grey",
                               vp.col="blue", vp.fill="light blue",
                               scale.col="red",
                               vp=NULL)
{
    ## if the viewport has a non-NULL layout.pos.row or layout.pos.col
    ## AND the viewport has a parent AND the parent has a layout
    ## represent the location of the viewport in the parent's layout ...
    if ((!is.null(v$layout.pos.row) || !is.null(v$layout.pos.col)) &&
        !is.null(parent.layout)) {
        if (!is.null(vp))
            pushViewport(vp)
        vp.mid <- grid.show.layout(parent.layout, vp.ex=vp.ex,
                                   cell.border="grey", cell.fill="white",
                                   cell.label=FALSE, newpage=newpage)
        pushViewport(vp.mid)
        pushViewport(v)
        gp.red <- gpar(col=scale.col)
        grid.rect(gp=gpar(col="blue", fill="light blue"))
        at <- grid.pretty(v$xscale)
        grid.xaxis(at=c(min(at), max(at)), gp=gp.red)
        at <- grid.pretty(v$yscale)
        grid.yaxis(at=c(min(at), max(at)), gp=gp.red)
        popViewport(2)
        if (!is.null(vp))
            popViewport()
    } else {
        if (newpage)
            grid.newpage()
        if (!is.null(vp))
            pushViewport(vp)
        grid.rect(gp=gpar(col=NULL, fill=border.fill))
        ## generate a viewport within the "top" viewport (vp) to represent the
        ## parent viewport of the viewport we are "show"ing (v).
        ## This is so that annotations at the edges of the
        ## parent viewport will be at least partially visible
        vp.mid <- viewport(0.5, 0.5, vp.ex, vp.ex)
        pushViewport(vp.mid)
        grid.rect(gp=gpar(fill="white"))
        x <- v$x
        y <- v$y
        w <- v$width
        h <- v$height
        pushViewport(v)
        grid.rect(gp=gpar(col=vp.col, fill=vp.fill))
        ## represent the "native" scale
        gp.red <- gpar(col=scale.col)
        at <- grid.pretty(v$xscale)
        grid.xaxis(at=c(min(at), max(at)), gp=gp.red)
        at <- grid.pretty(v$yscale)
        grid.yaxis(at=c(min(at), max(at)), gp=gp.red)
        grid.text(as.character(w), gp=gp.red,
                  just=c("centre", "bottom"),
                  x=unit(.5, "npc"), y=unit(1, "npc") + unit(.05, "inches"))
        grid.text(as.character(h), gp=gp.red,
                  just=c("left", "centre"),
                  x=unit(1, "npc") + unit(.05, "inches"), y=unit(.5, "npc"))
        popViewport()
        ## annotate the location and dimensions of the viewport
        grid.lines(unit.c(x, x), unit.c(unit(0, "npc"), y),
                   gp=gpar(col=scale.col, lty="dashed"))
        grid.lines(unit.c(unit(0, "npc"), x), unit.c(y, y),
                   gp=gpar(col=scale.col, lty="dashed"))
        grid.text(as.character(x), gp=gp.red,
                  just=c("centre", "top"),
                  x=x, y=unit(-.05, "inches"))
        grid.text(as.character(y), gp=gp.red,
                  just=c("bottom"),
                  x=unit(-.05, "inches"), y=y, rot=90)
        popViewport()
        if (!is.null(vp))
            popViewport()
    }
}

## old grid.legend <-
function(pch, labels, frame=TRUE,
                        hgap=unit(0.5, "lines"), vgap=unit(0.5, "lines"),
                        default.units="lines",
                        gp=gpar(), draw=TRUE,
                        vp=NULL) {
  ## Type checking on arguments
  labels <- as.character(labels)
  nkeys <- length(labels)
  if (length(pch) != nkeys)
    stop("'pch' and 'labels' not the same length")
  if (!is.unit(hgap))
    hgap <- unit(hgap, default.units)
  if (length(hgap) != 1)
    stop("'hgap' must be single unit")
  if (!is.unit(vgap))
    vgap <- unit(vgap, default.units)
  if (length(vgap) != 1)
    stop("'vgap' must be single unit")
  gf <- grid.frame(layout=grid.layout(nkeys, 2), vp=vp, gp=gp, draw=FALSE)
  for (i in 1L:nkeys) {
    if (i==1) {
      symbol.border <- unit.c(vgap, hgap, vgap, hgap)
      text.border <- unit.c(vgap, unit(0, "npc"), vgap, hgap)
    }
    else {
      symbol.border <- unit.c(vgap, hgap, unit(0, "npc"), hgap)
      text.border <- unit.c(vgap, unit(0, "npc"), unit(0, "npc"), hgap)
    }
    grid.pack(gf, grid.points(.5, .5, pch=pch[i], draw=FALSE),
              col=1, row=i, border=symbol.border,
              width=unit(1, "lines"), height=unit(1, "lines"),
              force.width=TRUE, draw=FALSE)
    grid.pack(gf, grid.text(labels[i], x=0, y=.5, just=c("left", "centre"),
                            draw=FALSE),
              col=2, row=i, border=text.border, draw=FALSE)
  }
  if (draw)
    grid.draw(gf)
  gf
}


legendGrob <-
    function(labels, nrow, ncol, byrow=FALSE,
	     do.lines = has.lty || has.lwd, lines.first=TRUE,
	     hgap=unit(1, "lines"), vgap=unit(1, "lines"),
	     default.units="lines",
	     pch, gp=gpar(), vp=NULL)
{
    ## Type checking on arguments; labels: character, symbol or expression:
    labels <- as.graphicsAnnot(labels)
    labels <- if(is.character(labels)) as.list(labels) else as.expression(labels)
    nkeys <- if(is.call(labels)) 1 else length(labels)
    if(nkeys == 0) return(nullGrob(vp=vp))
    if (!is.unit(hgap))
	hgap <- unit(hgap, default.units)
    if (length(hgap) != 1) stop("'hgap' must be single unit")
    if (!is.unit(vgap))
	vgap <- unit(vgap, default.units)
    if (length(vgap) != 1) stop("'vgap' must be single unit")
    ## nrow, ncol
    miss.nrow <- missing(nrow)
    miss.ncol <- missing(ncol)
    if(miss.nrow && miss.ncol) {ncol <- 1; nrow <- nkeys} # defaults to 1-column legend
    else if( miss.nrow && !miss.ncol) nrow <- ceiling(nkeys / ncol)
    else if(!miss.nrow &&  miss.ncol) ncol <- ceiling(nkeys / nrow)
    if(nrow < 1) stop("'nrow' must be >= 1")
    if(ncol < 1) stop("'ncol' must be >= 1")
    if(nrow * ncol < nkeys)
        stop("nrow * ncol < #{legend labels}")
    ## pch, gp
    if(has.pch <- !missing(pch) && length(pch) > 0) pch <- rep_len(pch, nkeys)
    if(doGP <- length(nmgp <- names(gp)) > 0) {
	if(has.lty  <-  "lty" %in% nmgp) gp$lty  <- rep_len(gp$lty, nkeys)
	if(has.lwd  <-  "lwd" %in% nmgp) gp$lwd  <- rep_len(gp$lwd, nkeys)
	if(has.col  <-  "col" %in% nmgp) gp$col  <- rep_len(gp$col,  nkeys)
	if(has.fill <- "fill" %in% nmgp) gp$fill <- rep_len(gp$fill, nkeys)
    } else {
	gpi <- gp
	if(missing(do.lines)) do.lines <- FALSE
    }

    ## main
    u0 <- unit(0, "npc")
    u1 <- unit(1, "char")
    ord <- if(lines.first) 1:2 else 2:1
    fg <- frameGrob(vp = vp)	  # set up basic frame grob (for packing)
    for (i in seq_len(nkeys)) {
	if(doGP) {
	    gpi <- gp
	    if(has.lty)	 gpi$lty <- gp$lty[i]
	    if(has.lwd)	 gpi$lwd <- gp$lwd[i]
	    if(has.col)	 gpi$col <- gp$col[i]
	    if(has.fill) gpi$fill<- gp$fill[i]
	}
	if(byrow) {
	    ci <- 1+ (i-1) %%  ncol
	    ri <- 1+ (i-1) %/% ncol
	} else {
	    ci <- 1+ (i-1) %/% nrow
	    ri <- 1+ (i-1) %%  nrow
	}
	## borders; unit.c creates a 4-vector of borders (bottom, left, top, right)
	vg <- if(ri != nrow) vgap else u0
	symbol.border <- unit.c(vg, u0, u0, 0.5 * hgap)
	text.border   <- unit.c(vg, u0, u0, if(ci != ncol) hgap else u0)

	## points/lines grob:
	plGrob <- if(has.pch && do.lines)
	    gTree(children = gList(linesGrob (0:1, 0.5, gp=gpi),
		  pointsGrob(0.5, 0.5, default.units="npc", pch=pch[i], gp=gpi))[ord])
	else if(has.pch) pointsGrob(0.5, 0.5, default.units="npc", pch=pch[i], gp=gpi)
	else if(do.lines) linesGrob(0:1, 0.5, gp=gpi)
	else nullGrob() # should not happen...
	fg <- packGrob(fg, plGrob,
		       col = 2*ci-1, row = ri, border = symbol.border,
		       width = u1, height = u1, force.width = TRUE)
	## text grob: add the labels
	gpi. <- gpi
	gpi.$col <- "black" # maybe needs its own 'gp' in the long run (?)
	fg <- packGrob(fg, textGrob(labels[[i]], x = 0, y = 0.5,
				    just = c("left", "centre"), gp=gpi.),
		       col = 2*ci, row = ri, border = text.border)
    }
    fg
}

grid.legend <- function(..., draw=TRUE)
{
    g <- legendGrob(...)# will error out if '...' has nonsense
    if (draw)
	grid.draw(g)
    invisible(g)
}

## Just a wrapper for a sample series of grid commands
grid.plot.and.legend <- function() {
  grid.newpage()
  top.vp <- viewport(width=0.8, height=0.8)
  pushViewport(top.vp)
  x <- stats::runif(10)
  y1 <- stats::runif(10)
  y2 <- stats::runif(10)
  pch <- 1L:3
  labels <- c("Girls", "Boys", "Other")
  lf <- frameGrob()
  plot <- gTree(children=gList(rectGrob(),
                  pointsGrob(x, y1, pch=1),
                  pointsGrob(x, y2, pch=2),
                  xaxisGrob(),
                  yaxisGrob()))
  lf <- packGrob(lf, plot)
  lf <- packGrob(lf, grid.legend(labels, pch=pch, draw=FALSE),
                 height=unit(1,"null"), side="right")
  grid.draw(lf)
}

#  File src/library/grid/R/interactive.R
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


grid.locator <- function(unit="native") {
  location <- c(grid.Call(C_locator), 1)
  if (is.na(location[1L]))
    invisible(NULL)
  else {
    transform <- solve(current.transform())
    location <- (location %*% transform)
    # The inverse viewport transform is from device coordinates into
    # inches relative to the current viewport
    location <- unit(location/location[3L], "inches")
    list(x=convertX(location[1L], unit),
         y=convertY(location[2L], unit))
  }
}

#  File src/library/grid/R/just.R
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

# NOTE: the order of the strings in these conversion functions must
# match the order of the enums in ../src/lattice.h
# NOTE: the result of match() is an integer, but subtracting 1 converts
# to real => have to convert back to integer for passing to C code

# If the user specifies two values, the first is horizontal
# justification and the second is vertical

# If the user specifies only one value, use the following
# conversion table to give a second default value
#
# bottom  -->  centre, bottom
# left    -->  left,   centre
# right   -->  right,  centre
# top     -->  centre, top
# centre  -->  centre, centre

valid.charjust <- function(just) {
  if (length(just) == 1) {
    # single value may be any valid just
    just <- as.integer(match(just[1L], c("left", "right", "bottom", "top",
                                        "centre", "center")) - 1)
    if (anyNA(just))
      stop("invalid justification")
  } else if (length(just) > 1) {
    # first value must be one of "left", "right", "centre", or "center"
    just[1L] <- as.integer(match(just[1L], c("left", "right", "bottom", "top",
                                           "centre", "center")) - 1)
    if (!(just[1L] %in% c(0, 1, 4, 5)))
      stop("invalid horizontal justification")
    # second value must be one of "bottom", "top", "centre", or "center"
    just[2L] <- as.integer(match(just[2L], c("left", "right", "bottom", "top",
                                           "centre", "center")) - 1)
    if (!(just[2L] %in% c(2, 3, 4, 5)))
      stop("invalid vertical justification")
    just <- as.integer(just)
  }
  # Extend to length 2 if necessary
  if (length(just) < 2) {
    if (length(just) == 0)
      just <- c(4, 4)
    else
      just <- switch (just[1L] + 1,
                      c(0, 4), # left
                      c(1, 4), # right
                      c(4, 2), # bottom
                      c(4, 3), # top
                      c(4, 4), # centre
                      c(4, 4)) # center
  }
  # Convert to numeric
  just <- c(switch(just[1L] + 1, 0, 1, NA, NA, 0.5, 0.5),
            switch(just[2L] + 1, NA, NA, 0, 1, 0.5, 0.5))
  # Final paranoid check
  if (anyNA(just))
    stop("invalid justification")
  just
}

valid.numjust <- function(just) {
  if (length(just) == 0) {
    c(0.5, 0.5)
  } else {
    if (length(just) < 2) {
      c(just, 0.5)
    } else {
      just
    }
  }
}

valid.just <- function(just) {
  if (is.character(just))
    valid.charjust(just)
  else {
    valid.numjust(as.numeric(just))
  }
}

resolveHJust <- function(just, hjust) {
  if (is.null(hjust) || length(hjust) == 0)
    valid.just(just)[1L]
  else
    hjust
}

resolveVJust <- function(just, vjust) {
  if (is.null(vjust) || length(vjust) == 0)
    valid.just(just)[2L]
  else
    vjust
}
#  File src/library/grid/R/layout.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2016 The R Core Team
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


is.layout <- function(l) {
  inherits(l, "layout")
}

# FIXME:  The internal C code now does a lot of recycling of
# unit values, units, and data.  Can some/most/all of the
# recycling stuff below be removed ?
valid.layout <- function(nrow, ncol, widths, heights, respect, just) {
  nrow <- as.integer(nrow)
  ncol <- as.integer(ncol)
  # make sure we're dealing with a unit object
  if (!is.logical(respect)) {
    respect <- as.matrix(respect)
    if (!is.matrix(respect) || any(dim(respect) != c(nrow, ncol)))
      stop("'respect' must be logical or an 'nrow' by 'ncol' matrix")
    }
  if (is.matrix(respect)) {
    respect.mat <- matrix(as.integer(respect),
                          dim(respect)[1L],
                          dim(respect)[2L])
    respect <- 2
  }
  else respect.mat <- matrix(0L, nrow, ncol)

  valid.just <- valid.just(just)
  l <- list(nrow = nrow, ncol = ncol,
            widths = widths, heights = heights,
            respect = respect, valid.respect=as.integer(respect),
            respect.mat = respect.mat,
            just=just, valid.just=valid.just)
  class(l) <- "layout"
  l
}

layout.torture <- function() {
  top.vp <- viewport(y=0, height=unit(1, "npc") - unit(1.5, "lines"),
                     just=c("centre", "bottom"))
  do.label <- function(label) {
    grid.rect(y=1, height=unit(1.5, "lines"),
              just=c("center", "top"))
    grid.text(label,
              y=unit(1, "npc") - unit(1, "lines"),
              gp=gpar(font=2))
  }
  # 1 = all relative widths and heights
  grid.show.layout(grid.layout(3,2), vp=top.vp)
  do.label("All dimensions relative -- no respect")
  # (1) with full respect
  grid.show.layout(grid.layout(3,2, respect=TRUE), vp=top.vp)
  do.label("All dimensions relative -- full respect")
  # (1) with partial respect
  grid.show.layout(grid.layout(3,2,respect=matrix(c(1,0,0,0,0,0), 3L, 2L, TRUE)),
                   vp=top.vp)
  do.label("All dimensions relative -- only top-left cell respected")
  # (1) with slightly weirder partial respect
  grid.show.layout(grid.layout(3,2,respect=matrix(c(1,0,0,0,0,1), 3L, 2L, TRUE)),
                   vp=top.vp)
  do.label("All relative -- top-left, bottom-right respected")
  # 2 = combination of absolute and relative widths and heights
  grid.show.layout(grid.layout(2, 3,
                       widths=unit(c(2,4,1), c("null", "cm", "null")),
                       heights=unit(c(6,4), c("cm", "null"))),
                   vp=top.vp)
  do.label("Absolute and relative -- no respect")
  # (2) with full respect
  grid.show.layout(grid.layout(2, 3,
                       widths=unit(c(2,4,1), c("null", "cm", "null")),
                       heights=unit(c(6,4), c("cm", "null")), respect=TRUE),
                   vp=top.vp)
  do.label("Absolute and relative -- full respect")
  # (2) with partial respect
  grid.show.layout(grid.layout(2, 3,
                       widths=unit(c(2,4,1), c("null", "cm", "null")),
                       heights=unit(c(6,4), c("cm", "null")),
                       respect=matrix(c(0,0,0,0,0,1), 2L, 3L, TRUE)),
                   vp=top.vp)
  do.label("Absolute and relative -- bottom-right respected")
}

# Return the region allocated by the layout of the current viewport
layoutRegion <- function(layout.pos.row=1, layout.pos.col=1) {
  region <- grid.Call(C_layoutRegion,
                      # This conversion matches the vailidity check in
                      # valid.viewport()
                      if (is.null(layout.pos.row)) layout.pos.row
                      else as.integer(rep(layout.pos.row, length.out=2)),
                      if (is.null(layout.pos.col)) layout.pos.col
                      else as.integer(rep(layout.pos.col, length.out=2)))
  list(left=unit(region[1L], "npc"),
       bottom=unit(region[2L], "npc"),
       width=unit(region[3L], "npc"),
       height=unit(region[4L], "npc"))
}

####################
# Accessors
####################

layout.nrow <- function(lay) {
  lay$nrow
}

layout.ncol <- function(lay) {
  lay$ncol
}

layout.widths <- function(lay) {
  lay$widths
}

layout.heights <- function(lay) {
  lay$heights
}

layout.respect <- function(lay) {
  switch(lay$respect + 1,
         FALSE,
         TRUE,
         lay$respect.mat)
}

####################
# Public constructor function
####################
grid.layout <- function (nrow = 1, ncol = 1,
                         widths = unit(rep_len(1, ncol), "null"),
                         heights = unit(rep_len(1, nrow), "null"),
                         default.units = "null",
                         respect = FALSE,
                         just="centre")
{
  if (!is.unit(widths))
    widths <- unit(widths, default.units)
  if (!is.unit(heights))
    heights <- unit(heights, default.units)
  valid.layout(nrow, ncol, widths, heights, respect, just)
}

####################
# Utility Functions
####################

dim.layout <- function(x) {
    c(x$nrow, x$ncol)
}
#  File src/library/grid/R/ls.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2016 The R Core Team
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


# Code for listing objects in various grid "namespaces"
# (gTrees, vpTrees, and the grid display list)

# Return a "gridListing" object,
# ... either ...
# "gridVectorListing", which is just character vector,
#     "grobListing", or "vpListing", or "vpNameListing", or
#     "vpPopListing", or "vpUpListing",
# ... or ...
# "gridListListing", which is list of "gridListing" objects,
#      "grobListListing", or "vpListListing", ...
# ... or ...
# "gridTreeListing", which is list of parent "gridVectorListing" object
#                    plus children "gridListing" object,
#      "gTreeListing", or "vpTreeListing", or "vpNameTreeListing"
#      (vpStack or vpTree produces a "vpTreeListing").
#      (vpPath [depth > 1] produces a "vpNameTreeListing").
#
# "vpListListing", and all "gridTreeListing" objects have a "depth" attribute

# The print method will print these in some format, but by having
# a separate object, others can capture the result and format the
# printing themselves.

grid.ls <- function(x=NULL, grobs=TRUE, viewports=FALSE, fullNames=FALSE,
                    recursive=TRUE, print=TRUE, flatten=TRUE, ...) {
    # If 'x' is NULL, list the grobs on the DL
    if (is.null(x)) {
        listing <- gridListDL(grobs=grobs, viewports=viewports,
                              fullNames=fullNames, recursive=recursive)
    } else {
        listing <- gridList(x, grobs=grobs, viewports=viewports,
                            fullNames=fullNames, recursive=recursive)
    }
    if (flatten) {
        listing <- flattenListing(listing)
    }
    if (is.logical(print)) {
        if (print) {
            print(listing)
        }
    } else if (is.function(print)) {
        print(listing, ...)
    } else {
        stop("invalid 'print' argument")
    }
    invisible(listing)
}

gridListDL <- function(x, grobs=TRUE, viewports=FALSE,
                       fullNames=FALSE, recursive=TRUE) {
    if (is.null(dev.list())) {
        result <- list(gridList(NULL))
    } else {
        display.list <- grid.Call(C_getDisplayList)
        dl.index <- grid.Call(C_getDLindex)
        result <- lapply(display.list[1L:dl.index], gridList,
                         grobs=grobs, viewports=viewports,
                         fullNames=fullNames, recursive=recursive)
        names(result) <- NULL
    }
    class(result) <- c("gridListListing", "gridListing")
    result
}

gridList <- function(x, grobs=TRUE, viewports=FALSE,
                     fullNames=FALSE, recursive=TRUE) {
    UseMethod("gridList")
}

gridList.default <- function(x, grobs=TRUE, viewports=FALSE,
                             fullNames=FALSE, recursive=TRUE) {
    if (is.null(x)) {
        # This handles empty slots in the display list
        result <- character()
        class(result) <- "gridListing"
    } else {
        stop("invalid object in 'listing'")
    }
    result
}

# Grob methods
gridList.grob <- function(x, grobs=TRUE, viewports=FALSE,
                          fullNames=FALSE, recursive=TRUE) {
    if (grobs) {
        if (fullNames) {
            result <- as.character(x)
        } else {
            result <- x$name
        }
        class(result) <- c("grobListing", "gridVectorListing", "gridListing")
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    if (viewports) {
        # Call makeContext() to get x$vp at drawing time
        x <- makeContext(x)
    }
    if (viewports && !is.null(x$vp)) {
        # Bit dodgy this bit
        # Emulates an "upViewport" on the DL
        n <- depth(x$vp)
        class(n) <- "up"
        result <- list(gridList(x$vp,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive),
                       result,
                       gridList(n,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive))
        class(result) <- c("gridListListing", "gridListing")
    }
    result
}

gridList.gList <- function(x, grobs=TRUE, viewports=FALSE,
                           fullNames=FALSE, recursive=TRUE) {
    # Allow for grobs=FALSE but viewports=TRUE
    if (grobs || viewports) {
        if (length(x) == 0L) {
            result <- character()
            class(result) <- "gridListing"
        } else {
            result <- lapply(x, gridList,
                             grobs=grobs, viewports=viewports,
                             fullNames=fullNames, recursive=recursive)
            class(result) <- c("gListListing", "gridListListing",
                               "gridListing")
        }
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

gridList.gTree <- function(x, grobs=TRUE, viewports=FALSE,
                           fullNames=FALSE, recursive=TRUE) {
    if (fullNames) {
        name <- as.character(x)
    } else {
        name <- x$name
    }
    class(name) <- c("grobListing", "gridVectorListing", "gridListing")
    if (viewports) {
        # Call makeContext() to get x$vp and x$childrenvp at drawing time
        x <- makeContext(x)
    }
    if (recursive) {
        # Allow for grobs=FALSE but viewports=TRUE
        result <- gridList(x$children[x$childrenOrder],
                           grobs=grobs, viewports=viewports,
                           fullNames=fullNames, recursive=recursive)
        if (viewports && !is.null(x$childrenvp)) {
            # Bit dodgy this bit
            # Emulates an "upViewport" on the DL
            n <- depth(x$childrenvp)
            class(n) <- "up"
            result <- list(gridList(x$childrenvp,
                                    grobs=grobs, viewports=viewports,
                                    fullNames=fullNames,
                                    recursive=recursive),
                           gridList(n,
                                    grobs=grobs, viewports=viewports,
                                    fullNames=fullNames,
                                    recursive=recursive),
                           result)
            class(result) <- c("gridListListing", "gridListing")
        }
        if (grobs) {
            result <- list(parent=name,
                           children=result)
            class(result) <- c("gTreeListing", "gridTreeListing",
                               "gridListing")
        } else if (!viewports) {
            result <- character()
            class(result) <- "gridListing"
        }
    } else {
        if (grobs) {
            result <- name
        } else {
            result <- character()
            class(result) <- "gridListing"
        }
    }
    if (viewports && !is.null(x$vp)) {
        # Bit dodgy this bit
        # Emulates an "upViewport" on the DL
        n <- depth(x$vp)
        class(n) <- "up"
        result <- list(gridList(x$vp,
                                grobs=grobs, viewports=viewports,
                                fullNames=fullNames,
                                recursive=recursive),
                       result,
                       gridList(n,
                                grobs=grobs, viewports=viewports,
                                fullNames=fullNames,
                                recursive=recursive))
        class(result) <- c("gridListListing", "gridListing")
    }
    result
}

# Viewport methods
gridList.viewport <- function(x, grobs=TRUE, viewports=FALSE,
                              fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        if (fullNames) {
            result <- as.character(x)
        } else {
            result <- x$name
        }
        class(result) <- c("vpListing", "gridVectorListing", "gridListing")
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

# ... are arugments to gridList
listvpListElement <- function(x, ...) {
    n <- depth(x)
    class(n) <- "up"
    result <- list(gridList(x, ...),
                   gridList(n, ...))
    class(result) <- c("gridListListing", "gridListing")
    result
}

gridList.vpList <- function(x, grobs=TRUE, viewports=FALSE,
                            fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        if (length(x) == 0L) {
            result <- character()
            class(result) <- "gridListing"
        } else if (length(x) == 1L) {
            result <- gridList(x[[1L]],
                              grobs=grobs, viewports=viewports,
                              fullNames=fullNames,
                              recursive=recursive)
        } else {
            result <- c(lapply(x[-length(x)], listvpListElement,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive),
                        list(gridList(x[[length(x)]],
                                     grobs=grobs, viewports=viewports,
                                     fullNames=fullNames,
                                     recursive=recursive)))
            attr(result, "depth") <- depth(x[[length(x)]])
            class(result) <- c("vpListListing", "gridListListing",
                               "gridListing")
        }
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

gridList.vpStack <- function(x, grobs=TRUE, viewports=FALSE,
                             fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        if (length(x) == 0L) {
            result <- character()
            class(result) <- "gridListing"
        } else if (length(x) == 1L || !recursive) {
            result <- gridList(x[[1L]],
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames, recursive=recursive)
        } else {
            theRest <- x[-1L]
            class(theRest) <- "vpStack"
            result <- gridList(theRest,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive)
            result <- list(parent=gridList(x[[1L]],
                             grobs=grobs, viewports=viewports,
                             fullNames=fullNames,
                             recursive=recursive),
                           children=result)
            attr(result, "depth") <- depth(x)
            class(result) <- c("vpTreeListing", "gridTreeListing",
                               "gridListing")
        }
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

gridList.vpTree <- function(x, grobs=TRUE, viewports=FALSE,
                            fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        if (recursive) {
            result <- gridList(x$children,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames, recursive=recursive)
            # Parent can only be a plain viewport
            result <- list(parent=gridList(x$parent,
                             grobs=grobs, viewports=viewports,
                             fullNames=fullNames,
                             recursive=recursive),
                           children=result)
            attr(result, "depth") <- depth(x$children) + 1
            class(result) <- c("vpTreeListing", "gridTreeListing",
                               "gridListing")
        } else {
            result <- gridList(x$parent,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames, recursive=recursive)
        }
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

# This handles downViewports in the display list
gridList.vpPath <- function(x, grobs=TRUE, viewports=FALSE,
                            fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        # Have to account for top-level downViewports that are
        # non-strict (i.e., they could navigate down quite a long way)
        # In particular, when the vpPath navigates down more
        # levels than there are names in the vpPath
        recordedDepth <- attr(x, "depth")
        if (!is.null(recordedDepth) && recordedDepth != depth(x)) {
            # In this case, need to prepend a fake path on the front
            # so that subsequent upViewport()s will work
            x <- vpPathFromVector(c(rep("...", recordedDepth - depth(x)),
                                    explode(as.character(x))))
        }
        # This would be simpler if paths were kept as vectors
        # but that redesign is a bit of an undertaking
        if (depth(x) == 1) {
            if (fullNames) {
                result <- paste0("downViewport[", x$name, "]")
            } else {
                result <- x$name
            }
            class(result) <- c("vpNameListing", "gridVectorListing",
                               "gridListing")
        } else if (depth(x) == 2) {
            result <- gridList(vpPath(x$name),
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive)
            result <- list(parent=gridList(vpPath(x$path),
                             grobs=grobs, viewports=viewports,
                             fullNames=fullNames,
                             recursive=recursive),
                           children=result)
            attr(result, "depth") <- depth(x)
            # Inherit updateVPDepth and updateVPPath methods
            # from vpTreeListing
            class(result) <- c("vpNameTreeListing", "vpTreeListing",
                               "gridTreeListing", "gridListing")
        } else {
            path <- explode(x$path)
            result <- gridList(vpPathFromVector(c(path[-1L], x$name)),
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive)
            result <- list(parent=gridList(vpPath(path[1L]),
                             grobs=grobs, viewports=viewports,
                             fullNames=fullNames,
                             recursive=recursive),
                           children=result)
            attr(result, "depth") <- depth(x)
            # Inherit updateVPDepth and updateVPPath methods
            # from vpTreeListing
            class(result) <- c("vpNameTreeListing", "vpTreeListing",
                               "gridTreeListing", "gridListing")
        }
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

# This handles popViewports in the display list
gridList.pop <- function(x, grobs=TRUE, viewports=FALSE,
                         fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        result <- as.character(x)
        if (fullNames) {
            result <- paste0("popViewport[", result, "]")
        }
        class(result) <- c("vpPopListing", "gridVectorListing", "gridListing")
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

# This handles upViewports in the display list
gridList.up <- function(x, grobs=TRUE, viewports=FALSE,
                        fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        result <- as.character(x)
        if (fullNames) {
            result <- paste0("upViewport[", result, "]")
        }
        class(result) <- c("vpUpListing", "gridVectorListing", "gridListing")
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

######################
# flatten methods for gridListing objects
######################

incDepth <- function(depth, n=1) {
    depth + n
}

decrDepth <- function(depth, x) {
    n <- as.numeric(gsub("^.+\\[", "",
                         gsub("\\]$", "",
                              as.character(x))))
    depth - n
}

# updateDepth modifies depth from sibling to sibling
# (flatListing methods take care of parent to child updates of depth)
updateGDepth <- function(x, gdepth) {
    UseMethod("updateGDepth")
}

updateGDepth.default <- function(x, gdepth) {
    gdepth
}

updateVPDepth <- function(x, vpdepth) {
    UseMethod("updateVPDepth")
}

updateVPDepth.default <- function(x, vpdepth) {
    vpdepth
}

updateVPDepth.vpListing <- function(x, vpdepth) {
    incDepth(vpdepth)
}

updateVPDepth.vpNameListing <- function(x, vpdepth) {
    incDepth(vpdepth)
}

updateVPDepth.vpListListing <- function(x, vpdepth) {
    incDepth(vpdepth, attr(x, "depth"))
}

updateVPDepth.vpUpListing <- function(x, vpdepth) {
    decrDepth(vpdepth, x)
}

updateVPDepth.vpPopListing <- function(x, vpdepth) {
    decrDepth(vpdepth, x)
}

updateVPDepth.vpTreeListing <- function(x, vpdepth) {
    incDepth(vpdepth, attr(x, "depth"))
}

incPath <- function(oldpath, addition) {
    if (nchar(oldpath) > 0) {
        paste0(oldpath, .grid.pathSep, as.character(addition))
    } else {
        as.character(addition)
    }
}

decrPath <- function(oldpath, x) {
    bits <- strsplit(oldpath, .grid.pathSep)[[1L]]
    n <- as.numeric(gsub("^.+\\[", "",
                         gsub("\\]$", "",
                              as.character(x))))
    if ((m <- (length(bits) - n)) == 0L) {
        ""
    } else {
	paste(bits[seq_len(m)], collapse=.grid.pathSep)
    }
}

updateGPath <- function(x, gpath) {
    UseMethod("updateGPath")
}

updateGPath.default <- function(x, gpath) {
    gpath
}

updateVPPath <- function(x, vppath) {
    UseMethod("updateVPPath")
}

updateVPPath.default <- function(x, vppath) {
    vppath
}

updateVPPath.vpListing <- function(x, vppath) {
    incPath(vppath, x)
}

updateVPPath.vpNameListing <- function(x, vppath) {
    incPath(vppath, x)
}

updateVPPath.vpListListing <- function(x, vppath) {
    incPath(vppath, x[[length(x)]])
}

updateVPPath.vpUpListing <- function(x, vppath) {
    decrPath(vppath, x)
}

updateVPPath.vpPopListing <- function(x, vppath) {
    decrPath(vppath, x)
}

updateVPPath.vpTreeListing <- function(x, vppath) {
    incPath(vppath,
            paste0(updateVPPath(x$parent, ""), .grid.pathSep,
                   updateVPPath(x$children, "")))
}

flatListing <- function(x, gDepth=0, vpDepth=0, gPath="", vpPath="") {
    UseMethod("flatListing")
}

flatListing.gridListing <- function(x, gDepth=0, vpDepth=0,
                                    gPath="", vpPath="") {
    if (length(x)) {
        list(name=as.character(x),
             gDepth=gDepth,
             vpDepth=vpDepth,
             gPath=gPath,
             vpPath=vpPath,
             type=class(x)[1L])
    } else {
        list(name=character(),
             gDepth=numeric(),
             vpDepth=numeric(),
             gPath=character(),
             vpPath=character(),
             type=character())
    }
}

flatListing.gTreeListing <- function(x, gDepth=0, vpDepth=0,
                                     gPath="", vpPath="") {
    # Increase gDepth and gPath
    flatChildren <- flatListing(x$children, incDepth(gDepth, 1), vpDepth,
                                incPath(gPath, x$parent), vpPath)
    list(name=c(as.character(x$parent), flatChildren$name),
         gDepth=c(gDepth, flatChildren$gDepth),
         vpDepth=c(vpDepth, flatChildren$vpDepth),
         gPath=c(gPath, flatChildren$gPath),
         vpPath=c(vpPath, flatChildren$vpPath),
         type=c(class(x)[1L], flatChildren$type))
}

OLDflatListing.vpTreeListing <- function(x, gDepth=0, vpDepth=0,
                                      gPath="", vpPath="") {
    # Increase vpDepth and vpPath
    flatChildren <- flatListing(x$children, gDepth, incDepth(vpDepth, 1),
                                gPath, incPath(vpPath, x$parent))
    list(name=c(as.character(x$parent), flatChildren$name),
         gDepth=c(gDepth, flatChildren$gDepth),
         vpDepth=c(vpDepth, flatChildren$vpDepth),
         gPath=c(gPath, flatChildren$gPath),
         vpPath=c(vpPath, flatChildren$vpPath),
         type=c(class(x)[1L], flatChildren$type))
}

flatListing.vpTreeListing <- function(x, gDepth=0, vpDepth=0,
                                      gPath="", vpPath="") {
    flatParent <- flatListing(x$parent, gDepth, vpDepth,
                              gPath, vpPath)
    depth <- attr(x$parent, "depth")
    if (is.null(depth)) {
        depth <- 1
    }
    # Increase vpDepth and vpPath
    flatChildren <- flatListing(x$children, gDepth, incDepth(vpDepth, depth),
                                gPath, updateVPPath(x$parent, vpPath))
    list(name=c(flatParent$name, flatChildren$name),
         gDepth=c(flatParent$gDepth, flatChildren$gDepth),
         vpDepth=c(flatParent$vpDepth, flatChildren$vpDepth),
         gPath=c(flatParent$gPath, flatChildren$gPath),
         vpPath=c(flatParent$vpPath, flatChildren$vpPath),
         type=c(flatParent$type, flatChildren$type))
}

flatListing.vpNameTreeListing <- function(x, gDepth=0, vpDepth=0,
                                      gPath="", vpPath="") {
    # Increase vpDepth and vpPath
    flatChildren <- flatListing(x$children, gDepth, incDepth(vpDepth, 1),
                                gPath, incPath(vpPath, x$parent))
    list(name=c(as.character(x$parent), flatChildren$name),
         gDepth=c(gDepth, flatChildren$gDepth),
         vpDepth=c(vpDepth, flatChildren$vpDepth),
         gPath=c(gPath, flatChildren$gPath),
         vpPath=c(vpPath, flatChildren$vpPath),
         type=c(class(x)[1L], flatChildren$type))
}

flatListing.gridListListing <- function(x, gDepth=0, vpDepth=0,
                                        gPath="", vpPath="") {
    n <- length(x)
    listListing <- list(name=character(),
                        gDepth=numeric(),
                        vpDepth=numeric(),
                        gPath=character(),
                        vpPath=character(),
                        type=character())
    for (i in 1L:n) {
        componentListing <- flatListing(x[[i]], gDepth, vpDepth,
                                        gPath, vpPath)
        listListing$name <- c(listListing$name,
                              componentListing$name)
        listListing$gDepth <- c(listListing$gDepth,
                                componentListing$gDepth)
        listListing$vpDepth <- c(listListing$vpDepth,
                                 componentListing$vpDepth)
        listListing$gPath <- c(listListing$gPath,
                               componentListing$gPath)
        listListing$vpPath <- c(listListing$vpPath,
                                componentListing$vpPath)
        listListing$type <- c(listListing$type,
                              componentListing$type)
        gPath <- updateGPath(x[[i]], gPath)
        vpPath <- updateVPPath(x[[i]], vpPath)
        gDepth <- updateGDepth(x[[i]], gDepth)
        vpDepth <- updateVPDepth(x[[i]], vpDepth)
    }
    listListing
}

flattenListing <- function(x) {
    listing <- flatListing(x)
    class(listing) <- "flatGridListing"
    listing
}

print.flatGridListing <- function(x, ...) {
    nestedListing(x, ...)
    invisible(x)
}

######################
# Print functions for flatGridListings
######################

nestedListing <- function(x, gindent="  ", vpindent=gindent) {

    makePrefix <- function(indent, depth) {
        indents <- rep(indent, length(depth))
        indents <- mapply(rep, indents, depth)
        sapply(indents, paste, collapse="")
    }

    if (!inherits(x, "flatGridListing"))
        stop("invalid listing")
    cat(paste0(makePrefix(gindent, x$gDepth),
               makePrefix(vpindent, x$vpDepth),
               x$name),
        sep = "\n")
}

pathListing <- function(x, gvpSep=" | ", gAlign=TRUE) {

    appendToPrefix <- function(path, name) {
        emptyPath <- nchar(path) == 0
        ifelse(emptyPath,
               name,
               paste(path, name, sep = .grid.pathSep))
    }

    padPrefix <- function(path, maxLen) {
        paste0(path, strrep(" ", maxLen - nchar(path)))
    }

    if (!inherits(x, "flatGridListing"))
        stop("invalid 'listing'")
    vpListings <- seq_along(x$name) %in% grep("^vp", x$type)
    paths <- x$vpPath
    # Only if viewport listings
    if (sum(vpListings) > 0) {
        paths[vpListings] <- appendToPrefix(paths[vpListings],
                                            x$name[vpListings])
        # If viewports are shown, then allow extra space before grobs
        maxLen <- max(nchar(paths[vpListings]))
    }
    else
	maxLen <- max(nchar(paths))

    # Only if grob listings
    if (sum(!vpListings) > 0) {
        if (gAlign) {
            paths[!vpListings] <- padPrefix(paths[!vpListings], maxLen)
        }
        paths[!vpListings] <- paste0(paths[!vpListings],
				     gvpSep,
				     appendToPrefix(x$gPath[!vpListings],
						    x$name[!vpListings]))
    }
    cat(paths, sep = "\n")
}

grobPathListing <- function(x, ...) {
    subset <- grep("^g", x$type)
    if (length(subset)) {
        cl <- class(x)
        subListing <- lapply(x, `[`, subset)
        class(subListing) <- cl
        pathListing(subListing, ...)
    }
}

# Tidy up the vpPath from grid.ls() to remove ROOT if it is there
clean <- function(paths) {
    sapply(lapply(paths,
                  function(x) {
                      pieces <- explode(x)
                      if (length(pieces) && pieces[1] == "ROOT")
                          pieces <- pieces[-1]
                      pieces
                  }),
           function(x) {
               if (length(x))
                   as.character(vpPath(x))
               else ""
           })
}

# Given a gPath, return complete grob paths that match from the display list
grid.grep <- function(path, x = NULL, grobs = TRUE, viewports = FALSE,
                      strict = FALSE, grep = FALSE, global = FALSE,
                      no.match = character(), vpPath = viewports) {
    if (!inherits(path, "gPath"))
        path <- gPath(path)
    depth <- depth(path)
    grep <- rep(grep, length.out = depth)

    # Get each piece of the path as a sequential char vector
    pathPieces <- explode(path)

    if (is.null(x)) {
        dl <- grid.ls(grobs=grobs,
                      viewports=viewports || vpPath,
                      print = FALSE)
    } else {
        dl <- grid.ls(x,
                      grobs=grobs,
                      viewports=viewports || vpPath,
                      print = FALSE)
    }
    if (!length(dl$name))
        return(no.match)
    # Only keep vpListing and grobListing
    names <- names(dl)
    dl <- lapply(dl,
                 function(x) {
                     if (viewports) {
                         keep <- dl$type == "vpListing" |
                             dl$type == "grobListing" |
                             dl$type == "gTreeListing"
                     } else {
                         keep <- dl$type == "grobListing" |
                             dl$type == "gTreeListing"                         
                     }
                     x[keep]
                 })
    names(dl) <- names
    # "depth" is vpDepth for vpListing and gDepth for grobListing
    # "path" is gPath for vpListing and vpPath for grobListing
    if (is.null(x)) {
        # (remove "ROOT" from path and depth)
        dl$depth <- ifelse(dl$type == "vpListing", dl$vpDepth - 1, dl$gDepth)
        dl$path <- ifelse(dl$type == "vpListing", clean(dl$vpPath), dl$gPath)
    } else {
        dl$depth <- ifelse(dl$type == "vpListing", dl$vpDepth, dl$gDepth)
        dl$path <- ifelse(dl$type == "vpListing", dl$vpPath, dl$gPath)
    }
    # Limit our search only to grobs whose depth matches ours
    # For not strict, we're only looking at the grob names, so all
    # depths apply.
    matchingDepths <- if (! strict) which((dl$depth + 1) >= depth)
                      else which((dl$depth + 1) == depth)
    if (!length(matchingDepths))
        return(no.match)

    nMatches <- 0
    searchMatches <- vector("list", length(matchingDepths))
    # For each name of the correct path length
    for (i in matchingDepths) {
        dlPathPieces <-
            if (dl$depth[i] > 0)
                c(explode(dl$path[i]), dl$name[i])
            else
                dl$name[i]
        matches <- logical(depth)
        if (!strict) {
            # NOTE that we already know that the dlPath is AT LEAST as long
            # as the path
            depthOffset <- 0
            while (depthOffset + depth <= dl$depth[i] + 1 &&
                   !all(matches)) {
                for (j in 1:depth) {
                    matches[j] <-
                        if (grep[j])
                            grepl(pathPieces[j], dlPathPieces[depthOffset + j])
                        else
                            pathPieces[j] == dlPathPieces[depthOffset + j]
                }
                depthOffset <- depthOffset + 1
            }
        } else {
            # Check whether we need to grep this level or not, attempt match
            # NOTE that we already know that path and dlPath are same length
            for (j in 1:depth) {
                matches[j] <-
                    if (grep[j])
                        grepl(pathPieces[j], dlPathPieces[j])
                    else
                        pathPieces[j] == dlPathPieces[j]
            }
        }
        # We have found a grob
        if (all(matches)) {
            if (!global) {
                # Returning early to avoid further searching
                if (dl$type[i] == "vpListing") {
                    result <- do.call("vpPath", list(dlPathPieces))
                } else {
                    result <- do.call("gPath", list(dlPathPieces))
                    if (vpPath) {
                        attr(result, "vpPath") <- clean(dl$vpPath[i])
                    } else {
                        attr(result, "vpPath") <- ""
                    }
                }
                return(result)
            } else {
                nMatches <- nMatches + 1
                if (dl$type[i] == "vpListing") {
                    result <- do.call("vpPath",
                                      list(dlPathPieces))
                } else {
                    result <- do.call("gPath",
                                      list(dlPathPieces))
                    if (vpPath) {
                        attr(result, "vpPath") <- clean(dl$vpPath[i])
                    } else {
                        attr(result, "vpPath") <- ""
                    }
                }
                searchMatches[[nMatches]] <- result
            }
        }
    }

    if (!nMatches)
        return(no.match)

    # We may have allocated a list too large earlier,
    # subset to only matching results
    searchMatches <- searchMatches[1:nMatches]

    return(searchMatches)
}

createMask <- function(mask, type="alpha") {
    force(mask)
    maskFun <- function() {
        grid.draw(mask, recording=FALSE)
    }
    result <- list(f=.mask(maskFun, type), ref=NULL)
    class(result) <- "GridMask"
    result
}   

isMask <- function(x) {
    inherits(x, "GridMask")
}

## "resolve" masks
resolveMask <- function(mask) {
    UseMethod("resolveMask")
}

resolveMask.NULL <- function(mask) {
    .setMask(NULL, NULL)
    NULL
}

resolveMask.GridMask <- function(mask) {
    ref <- .setMask(mask$f, mask$ref)
    resolvedMask(mask, ref)
}

resolvedMask <- function(mask, ref) {
    UseMethod("resolvedMask")
}

resolvedMask.GridMask <- function(mask, ref) {
    mask$ref <- ref
    class(mask) <- c("GridResolvedMask", class(mask))
    mask
}    

resolvedMask.GridResolvedMask <- function(mask, ref) {
    mask$ref <- ref
    mask
}    

unresolveMask <- function(mask) {
    UseMethod("unresolveMask")
}
    
## Unresolved masks just pass through
unresolveMask.GridMask <- function(mask) {
    mask
}

unresolveMask.GridResolvedMask <- function(mask) {
    result <- list(f=mask$f, ref=NULL)
    class(result) <- "GridMask"
    result
}

## User interface
as.mask <- function(x, type=c("alpha", "luminance")) {
    if (!is.grob(x))
        stop("Only a grob can be converted to a mask")
    createMask(x, match.arg(type))
}

# API to access detailed text metric info
#
#  Copyright (C) 1995-2012 The R Core Team

# This first function does NOT return a "unit" object
# It is just access to font metric info in the calling context
# (similar to the convert*() functions, with corresponding caveats on use)

grid.textMetric <- function(string) {

}

# It should be possible to define units like "strascent" and "strdescent"
#  File src/library/grid/R/origin.R
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

valid.origin <- function(origin) {
  origin <- as.integer(match(origin,
                             c("bottom.left", "top.left",
                               "bottom.right", "top.right")) - 1)
  if (anyNA(origin))
    stop("invalid 'origin'")
  origin
}

origin.left <- function(origin) {
  switch (origin,
          bottom.left = TRUE,
          bottom.right = FALSE,
          top.left = TRUE,
          top.right = FALSE)
}

origin.right <- function(origin) {
  switch (origin,
          bottom.left = FALSE,
          bottom.right = TRUE,
          top.left = FALSE,
          top.right = TRUE)
}

origin.bottom <- function(origin) {
  switch (origin,
          bottom.left = TRUE,
          bottom.right = TRUE,
          top.left = FALSE,
          top.right = FALSE)
}

origin.top <- function(origin) {
  switch (origin,
          bottom.left = FALSE,
          bottom.right = FALSE,
          top.left = TRUE,
          top.right = TRUE)
}

swap.origin.horizontal <- function(origin) {
  switch (origin,
          bottom.left = "bottom.right",
          bottom.right = "bottom.left",
          top.left = "top.right",
          top.right = "top.left")
}

swap.origin.vertical <- function(origin) {
  switch (origin,
          bottom.left = "top.left",
          bottom.right = "top.right",
          top.left = "bottom.left",
          top.right = "bottom.right")
}

## Mark a grob as something to draw as a single "path"
as.path <- function(x, gp=gpar(), rule=c("winding", "evenodd")) {
    if (!is.grob(x))
        stop("Only a grob can be converted to a path")
    path <- list(grob=x, gp=gp, rule=match.arg(rule))
    class(path) <- "GridPath"
    path
}


## Stroke the outline of a path defined by a grob
drawDetails.GridStroke <- function(x, recording) {
    path <- function() {
        grid.draw(x$path, recording=FALSE)
    }
    grid.Call.graphics(C_stroke, path)
}

strokeGrob <- function(x, ...) {
    UseMethod("strokeGrob")
}

strokeGrob.grob <- function(x, name=NULL, gp=gpar(), vp=NULL, ...) {
    ## grob() rather than gTree() because a "path" only produces
    ## a single shape
    ## (e.g., pattern fills should ALWAYS be resolved)
    stroke <- grob(path=x, name=name, gp=gp, vp=vp, cl="GridStroke")
    stroke
}

strokeGrob.GridPath <- function(x, name=NULL, vp=NULL, ...) {
    stroke <- grob(path=x$grob, name=name, gp=x$gp, vp=vp, cl="GridStroke")
    stroke
}

grid.stroke <- function(...) {
    grid.draw(strokeGrob(...))
}

## Fill the outline of a path defined by a grob
drawDetails.GridFill <- function(x, recording) {
    path <- function() {
        grid.draw(x$path, recording=FALSE)
    }
    grid.Call.graphics(C_fill, path, .ruleIndex(x$rule))
}

fillGrob <- function(x, ...) {
    UseMethod("fillGrob")
}

fillGrob.grob <- function(x, rule=c("winding", "evenodd"),
                          name=NULL, gp=gpar(), vp=NULL, ...) {
    fill <- grob(path=x, rule=match.arg(rule),
                 name=name, gp=gp, vp=vp, cl="GridFill")
    fill
}

fillGrob.GridPath <- function(x, name=NULL, vp=NULL, ...) {
    fill <- grob(path=x$grob, rule=x$rule,
                 name=name, gp=x$gp, vp=vp, cl="GridFill")
    fill
}

grid.fill <- function(...) {
    grid.draw(fillGrob(...))
}



## Stroke and fill the outline of a path defined by a grob
drawDetails.GridFillStroke <- function(x, recording) {
    path <- function() {
        grid.draw(x$path, recording=FALSE)
    }
    grid.Call.graphics(C_fillStroke, path, .ruleIndex(x$rule))
}

fillStrokeGrob <- function(x, ...) {
    UseMethod("fillStrokeGrob")
}

fillStrokeGrob.grob <- function(x, rule=c("winding", "evenodd"),
                                name=NULL, gp=gpar(), vp=NULL, ...) {
    fillStroke <- grob(path=x, rule=match.arg(rule),
                       name=name, gp=gp, vp=vp, cl="GridFillStroke")
    fillStroke
}

fillStrokeGrob.GridPath <- function(x, name=NULL, vp=NULL, ...) {
    fillStroke <- grob(path=x$grob, rule=x$rule,
                       name=name, gp=x$gp, vp=vp, cl="GridFillStroke")
    fillStroke
}

grid.fillStroke <- function(...) {
    grid.draw(fillStrokeGrob(...))
}

################################
## Other grob methods

flattenCoords <- function(coords, name, rule) {
    UseMethod("flattenCoords")
}

flattenCoords.GridGrobCoords <- function(coords, name, rule) {
    names(coords) <- rep("1", length(coords))
    attr(coords, "rule") <- rule
    coords
}

flattenCoords.GridGTreeCoords <- function(coords, name, rule) {
    childCoords <- lapply(coords, flattenCoords, name, rule)
    coords <- do.call(c, childCoords)
    names(coords) <- rep("1", length(coords))
    gridGrobCoords(coords, name, rule)
}   

## NOTE, we need to create a gTree to include any 'gp' and 'vp' settings
## (and cannot just call grobCoords(x) because, although 'x' is a gTree
##  is has no children - the "child" is x$grob)
grobCoords.GridStroke <- function(x, closed, ...) {
    if (closed) {
        emptyGrobCoords(x$name)
    } else {
        coords <- grobCoords(gTree(children=gList(x$path), gp=x$gp, vp=x$vp),
                             closed, ...)
        flattenCoords(coords, x$name, NULL)
    }
}

## NOTE, like grobPoints.gList(), we need to call grobCoords()
## on the "child" grob so that it can perform any relevant set up
grobPoints.GridStroke <- function(x, closed, ...) {
    if (closed) {
        emptyGrobCoords(x$name)
    } else {
        coords <- grobCoords(x$path, closed, ...)
        flattenCoords(coords, x$name, NULL)
    }
}

grobCoords.GridFill <- function(x, closed, ...) {
    if (closed) {
        coords <- grobCoords(gTree(children=gList(x$path), gp=x$gp, vp=x$vp),
                             closed, ...)
        flattenCoords(coords, x$name, x$rule)
    } else {
        emptyGrobCoords(x$name)
    }
}

grobPoints.GridFill <- function(x, closed, ...) {
    if (closed) {
        coords <- grobCoords(x$path, closed, ...)
        flattenCoords(coords, x$name, x$rule)
    } else {
        emptyGrobCoords(x$name)
    }
}

grobCoords.GridFillStroke <- function(x, closed, ...) {
    coords <- grobCoords(gTree(children=gList(x$path), gp=x$gp, vp=x$vp),
                         closed, ...)
    flattenCoords(coords, x$name, x$rule)    
}

grobPoints.GridFillStroke <- function(x, closed, ...) {
    coords <- grobCoords(x$path, closed, ...)
    flattenCoords(coords, x$name, x$rule)
}


## Create R objects defining patterns

is.pattern <- function(x) {
    inherits(x, "GridPattern")
}

is.patternList <- function(x) {
    inherits(x, "GridPatternList")
}

linearGradient <- function(colours = c("black", "white"),
                           stops = seq(0, 1, length.out = length(colours)),
                           x1 = unit(0, "npc"), y1 = unit(0, "npc"),  
                           x2 = unit(1, "npc"), y2 = unit(1, "npc"),
                           default.units = "npc",
                           extend = c("pad", "repeat", "reflect", "none"),
                           group = TRUE) {

    nstops <- max(length(colours), length(stops))
    if (nstops < 1)
        stop("colours and stops must be at least length 1")
    colours <- rep(colours, length.out = nstops)
    stops <- rep(stops, length.out = nstops)

    if (! is.unit(x1))
        x1 <- unit(x1, default.units)
    if (! is.unit(x2))
        x2 <- unit(x2, default.units)
    if (! is.unit(y1))
        y1 <- unit(y1, default.units)
    if (! is.unit(y2))
        y2 <- unit(y2, default.units)

    if (length(x1) != 1 || length(x2) != 1 ||
        length(y1) != 1 || length(y2) != 1)
        stop("x1, y1, x2, and y2 must all be length 1")

    grad <- list(x1 = x1, y1 = y1,
                 x2 = x2, y2 = y2,
                 stops = as.numeric(stops), colours = colours,
                 extend = match.arg(extend), group = as.logical(group))
    class(grad) <- c("GridLinearGradient", "GridPattern")
    grad
}

radialGradient <- function(colours = c("black", "white"),
                           stops = seq(0, 1, length.out = length(colours)),
                           cx1 = unit(.5, "npc"), cy1 = unit(.5, "npc"),
                           r1 = unit(0, "npc"),
                           cx2 = unit(.5, "npc"), cy2 = unit(.5, "npc"),
                           r2 = unit(.5, "npc"),
                           default.units = "npc",
                           extend = c("pad", "repeat", "reflect", "none"),
                           group = TRUE) {

    nstops <- max(length(colours), length(stops))
    if (nstops < 1)
        stop("colours and stops must be at least length 1")
    colours <- rep(colours, length.out = nstops)
    stops <- rep(stops, length.out = nstops)

    if (!is.unit(cx1))
        cx1 <- unit(cx1, default.units)
    if (!is.unit(cy1))
        cy1 <- unit(cy1, default.units)
    if (!is.unit(r1))
        r1 <- unit(r1, default.units)
    if (!is.unit(cx2))
        cx2 <- unit(cx2, default.units)
    if (!is.unit(cy2))
        cy2 <- unit(cy2, default.units)
    if (!is.unit(r2))
        r2 <- unit(r2, default.units)

    if (length(cx1) != 1 || length(cx2) != 1 ||
        length(cy1) != 1 || length(cy2) != 1 ||
        length(r1) != 1  || length(r2) != 1)
        stop("cx1, cy1, cx2, cy2, r1, and r2 must all be length 1")

    grad <- list(cx1 = cx1, cy1 = cy1, r1=r1,
                 cx2 = cx2, cy2 = cy2, r2=r2,
                 stops = as.numeric(stops), colours = colours,
                 extend = match.arg(extend), group = as.logical(group))
    class(grad) <- c("GridRadialGradient", "GridPattern")
    grad
}

## Wrap the pattern grob in a gTree with "initial" 'gp' settings
## for the grob to inherit
## (we are particularly concerned about the grob inheriting the
##  fill from its parent, which would mean infinite recursion)
## AND wrap that gTree in a function that draws it.
pattern <- function(grob,
                    x = 0.5, y = 0.5, width = 1, height = 1,
                    default.units = "npc",
                    just="centre", hjust=NULL, vjust=NULL,
                    extend = c("pad", "repeat", "reflect", "none"),
                    gp = gpar(fill="transparent"),
                    group = TRUE) {

    if (! is.unit(x))
        x <- unit(x, default.units)
    if (! is.unit(y))
        y <- unit(y, default.units)
    if (! is.unit(width))
        width <- unit(width, default.units)
    if (! is.unit(height))
        height <- unit(height, default.units)
    hjust = resolveHJust(just, hjust)
    vjust = resolveVJust(just, vjust)

    if (length(x) != 1 || length(y) != 1 ||
        length(width) != 1 || length(height) != 1)
        stop("x, y, width, and height must all be length 1")

    force(grob)
    if (!is.grob(grob))
        stop("Pattern must be based on grob")
    
    ## Do NOT want x$gp$fill to be NULL because that would mean
    ## that 'x' inherits its fill from the grob that it is
    ## filling, which means infinite recursion
    if (is.null(gp$fill)) {
        gp$fill <- "transparent"
        warning("Missing pattern fill has been set to transparent")
    }
    patternFun <- function() {
        grid.draw(gTree(children=gList(grob), gp=gp), recording=FALSE)
    }
    pat <- list(f=patternFun,
                x=x, y=y, width=width, height=height,
                hjust=hjust, vjust=vjust,
                extend = match.arg(extend), group = as.logical(group))
    class(pat) <- c("GridTilingPattern", "GridPattern")
    pat
}

################################################################################
## Pattern resolution

## If CURRENT gp$fill is a pattern then need to attach
## "built" grob (post makeContent() call) to gp$fill for (subsequent)
## resolution of the pattern
## NOTE that this is NOT just for grob$gp$fill because inherited
## gp$fill may be an unresolved pattern
## NOTE that can indiscriminately attach grob to both resolved
## and unresolved patterns because resolved patterns will just ignore
## attached grob
recordGrobForPatternResolution <- function(x) {
    gpar <- grid.Call(C_getGPar)
    if (is.pattern(gpar$fill)) {
        attr(gpar$fill, "grob") <- x
        class(gpar$fill) <- c("GridGrobPattern", class(gpar$fill))
        grid.Call(C_setGPar, gpar)
    } else if (is.patternList(gpar$fill)) {
        attr(gpar$fill, "grob") <- x
        class(gpar$fill) <- c("GridGrobPatternList", class(gpar$fill))
        grid.Call(C_setGPar, gpar)
    }
}

## If gTree has a pattern in gp$fill AND gp$fill$group then attach
## "built" gTree (post makeContent() call) to gp$fill for resolution
## of the pattern AND resolve the pattern
## (if gTree has a LIST of patterns in gp$fill, this may resolve
##  *some* patterns in the list)
## NOTE that this IS with gp$fill because inherited fills pass through
## (see below)
## NOTE that patterns within pattern list that are unresolved
## will get grob attached to list and get resolved in drawing of
## gTree children
## NOTE that inherited patterns SHOULD pass through untouched
## (if group is TRUE inherited pattern should already be resolved
##  and if group is FALSE inherited pattern should be passed through)
recordGTreeForPatternResolution <- function(x) {
    if ((is.pattern(x$gp$fill) && x$gp$fill$group) ||
        is.patternList(x$gp$fill)) { 
        gpar <- grid.Call(C_getGPar)
        attr(gpar$fill, "grob") <- x
        resolvedFill <- resolveFill(gpar$fill, 1)
        ## Resolution may generate NULL (e.g., if gTree has nothing to fill)
        if (is.null(resolvedFill)) {
            gpar$fill <- "transparent"
        } else {
            gpar$fill <- resolvedFill
        }
        grid.Call(C_setGPar, gpar)
    }
}

resolvedPattern <- function(pattern, ref) {
    pattern$ref <- ref
    class(pattern) <- c("GridResolvedPattern", class(pattern))
    pattern
}

## Called when drawing a grob
resolveFill <- function(fill, ...) {
    UseMethod("resolveFill")
}

## Simple fills include an R colour (integer or string) or NA
## These just pass through
resolveFill.default <- function(fill, ...) {
    fill
}

## A pattern fill that has already been resolved
resolveFill.GridResolvedPattern <- function(fill, ...) {
    fill
}

## A pattern fill that needs resolving
## (a grid::GridPattern)
## This will handle viewports (with a single pattern)
resolveFill.GridPattern <- function(fill, ...) {
    resolvePattern(fill)
}

## This will handle viewports (with a list of patterns)
resolveFill.GridPatternList <- function(fill, ...) {
    ## NOTE that some patterns may be resolved, but others may not
    ## (so we cannot mark the entire list as resolved)
    resolvedPatterns <- lapply(fill,
                               function(x) {
                                   if (x$group)
                                       resolvePattern(x)
                                   else
                                       x
                               })
    class(resolvedPatterns) <- class(fill)
    resolvedPatterns
}

## This will handle grobs (with a single pattern)
resolveFill.GridGrobPattern <- function(fill, index=1, ...) {
    ## All predrawing has been done
    ## ('return' is just for GridDefine grobs)
    pts <- grobPoints(attr(fill, "grob"), closed=TRUE, return=TRUE)
    if (!isEmptyCoords(pts)) {
        if (fill$group || length(pts) == 1) {
            ## Pattern is relative to bounding box of all shapes
            bbox <- coordsBBox(pts)
        } else {
            ## Pattern is relative to bounding box of individual shapes
            if (index > length(pts)) {
                warning("grob drawing produces more shapes than grob coords
(recycling coords)")
                index <- (index - 1) %% length(pts) + 1
            }
            ## Individual shape may consist of more than one set of
            ## coordinates (e.g., single path consists of distinct shapes)
            shapeIndex <- names(pts) %in% index
            ## Fallback if 'pts' does not have names to identify shapes
            if (!any(shapeIndex)) {
                shapeIndex <- index
            }
            bbox <- coordsBBox(pts, shapeIndex)
        }
        ## Temporary viewport for calculations, so do NOT record on grid DL
        ## Also, ensure NO mask and NO clip
        ## (at least initially) for resolution of pattern
        ## Also, set fill to "transparent"
        ## (to avoid this viewport picking up the fill being resolved)
        pushViewport(viewport(bbox$left, bbox$bottom, bbox$width, bbox$height,
                              default.units="in",
                              just=c("left", "bottom"),
                              clip="off", mask="none",
                              gp=gpar(fill="transparent")),
                     recording=FALSE)
        pattern <- resolvePattern(fill)
        popViewport(recording=FALSE)
        pattern
    } else {
        warning("Pattern fill applied to object with no inside")
        ## Set fill to transparent
        "transparent"
    }
}

## This will handle grobs (with a list of patterns)
resolveFill.GridGrobPatternList <- function(fill, ...) {
    ## All predrawing has been done
    pts <- grobPoints(attr(fill, "grob"), closed=TRUE, return=TRUE)
    if (!isEmptyCoords(pts)) {
        resolvedFills <- vector("list", length(pts))
        for (i in seq_along(pts)) {
            ## Recycle patterns if necessary
            which <- (i - 1) %% length(fill) + 1
            ## Only resolve if not already resolved
            if (inherits(fill[[which]], "GridResolvedPattern")) {
                resolvedFills[[i]] <- fill[[which]]
                next
            }
            if (fill[[which]]$group || length(pts) == 1) {
                ## Pattern is relative to bounding box of all shapes
                bbox <- coordsBBox(pts)
            } else {
                ## Pattern is relative to bounding box of individual shapes
                ## Individual shape may consist of more than one set of
                ## coordinates (e.g., single path consists of distinct shapes)
                shapeIndex <- names(pts) %in% i
                ## Fallback if 'pts' does not have names to identify shapes
                if (!any(shapeIndex)) {
                    shapeIndex <- i
                }
                bbox <- coordsBBox(pts, shapeIndex)
            }
            ## Temporary viewport for calculations, so do NOT record on grid DL
            ## Also, ensure NO mask and NO clip
            ## (at least initially) for resolution of pattern
            ## Also, set fill to "transparent"
            ## (to avoid this viewport picking up the fill being resolved)
            pushViewport(viewport(bbox$left, bbox$bottom,
                                  bbox$width, bbox$height,
                                  default.units="in",
                                  just=c("left", "bottom"),
                                  clip="off", mask="none",
                                  gp=gpar(fill="transparent")),
                         recording=FALSE)
            pattern <- resolvePattern(fill[[which]])
            popViewport(recording=FALSE)
            resolvedFills[[i]] <- pattern
        }
        class(resolvedFills) <- c("GridResolvedPatternList",
                                  class(fill))
        resolvedFills
    } else {
        warning("Pattern fill applied to object with no inside")
        ## Set fill to transparent
        "transparent"
    }
}

resolveGTreeFill <- function(fill, pts) {
    if (!isEmptyCoords(pts)) {
        ## Pattern is relative to bounding box of all shapes
        bbox <- coordsBBox(pts)
        ## Temporary viewport for calculations, so do NOT record on grid DL
        ## Also, ensure NO mask and NO clip
        ## (at least initially) for resolution of pattern
        ## Also, set fill to "transparent"
        ## (to avoid this viewport picking up the fill being resolved)
        pushViewport(viewport(bbox$left, bbox$bottom,
                              bbox$width, bbox$height,
                              default.units="in",
                              just=c("left", "bottom"),
                              clip="off", mask="none",
                              gp=gpar(fill="transparent")),
                     recording=FALSE)
        pattern <- resolvePattern(fill)
        popViewport(recording=FALSE)
        pattern
    } 
}

## This will handle gTrees (with a single pattern)
## This should ONLY be called if fill$group is TRUE
## (see recordGTreeForPatternResolution())
resolveFill.GridGTreePattern <- function(fill, index=1, ...) {
    ## All predrawing has been done
    pts <- grobPoints(attr(fill, "grob"), closed=TRUE, return=TRUE)
    resolveGTreeFill(fill, pts)
}

## This will handle gTrees (with a list of patterns)
## For each fill ...
##   if fill$group resolve relative to gTree bbox
##   else pass through
## (similar to resolveFill.GridPatternList() for viewports)
resolveFill.GridGTreePatternList <- function(fill, ...) {
    ## All predrawing has been done
    pts <- grobPoints(attr(fill, "grob"), closed=TRUE, return=TRUE)
    resolveOneFill <- function(x) {
        if (x$group) {
            resolveGTreeFill(x, pts)
        } else {
            x
        }
    }
    resolvedFills <- lapply(fill, resolveOneFill)
    class(resolvedFills) <- class(fill)
    resolvedFills
}

resolvePattern <- function(pattern) {
    UseMethod("resolvePattern")
}

resolvePattern.GridLinearGradient <- function(pattern) {
    p1 <- deviceLoc(pattern$x1, pattern$y1, valueOnly=TRUE, device=TRUE)
    p2 <- deviceLoc(pattern$x2, pattern$y2, valueOnly=TRUE, device=TRUE)
    index <- .setPattern(.linearGradientPattern(pattern$colours,
                                                pattern$stops,
                                                p1$x, p1$y, p2$x, p2$y,
                                                extend=pattern$extend))
    resolvedPattern(pattern, index)
}

resolvePattern.GridRadialGradient <- function(pattern) {
    c1 <- deviceLoc(pattern$cx1, pattern$cy1, valueOnly=TRUE, device=TRUE)
    r1 <- min(sqrt(sum(unlist(deviceDim(unit(0, "in"), pattern$r1,
                                        valueOnly=TRUE, device=TRUE))^2)),
              sqrt(sum(unlist(deviceDim(pattern$r1, unit(0, "in"), 
                                        valueOnly=TRUE, device=TRUE))^2)))
    c2 <- deviceLoc(pattern$cx2, pattern$cy2, valueOnly=TRUE, device=TRUE)
    r2 <- min(sqrt(sum(unlist(deviceDim(unit(0, "in"), pattern$r2,
                                        valueOnly=TRUE, device=TRUE))^2)),
              sqrt(sum(unlist(deviceDim(pattern$r2, unit(0, "in"), 
                                        valueOnly=TRUE, device=TRUE))^2)))
    index <- .setPattern(.radialGradientPattern(pattern$colours,
                                                pattern$stops,
                                                c1$x, c1$y, r1,
                                                c2$x, c2$y, r2,
                                                extend=pattern$extend))
    resolvedPattern(pattern, index)
}

resolvePattern.GridTilingPattern <- function(pattern) {
    xy <- deviceLoc(pattern$x, pattern$y, valueOnly=TRUE, device=TRUE)
    wh <- deviceDim(pattern$width, pattern$height, valueOnly=TRUE, device=TRUE)
    left <- xy$x - pattern$hjust*wh$w
    bottom <- xy$y - pattern$vjust*wh$h
    index <- .setPattern(.tilingPattern(pattern$f,
                                        left, bottom, wh$w, wh$h,
                                        extend=pattern$extend))
    resolvedPattern(pattern, index)
}

## Used when "grab"ing the display list to "demote"
## a resolved pattern
unresolveFill <- function(fill) {
    UseMethod("unresolveFill")
}

## Simple fills include an R colour (integer or string) or NA
## These just pass through
unresolveFill.default <- function(fill) {
    fill
}

unresolveFill.GridPattern <- function(fill) {
    unresolvePattern(fill)
}

unresolvePattern <- function(pattern) {
    UseMethod("unresolvePattern")
}
    
## Unresolved patterns just pass through
unresolvePattern.GridPattern <- function(pattern) {
    pattern
}

unresolvePattern.GridResolvedPattern <- function(pattern) {
    pattern$ref <- NULL
    class(pattern) <-
        class(pattern)[!(class(pattern) %in% "GridResolvedPattern")]
    pattern
}

#  File src/library/grid/R/primitives.R
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


# Function that creates a description of an arrow head
# to add to a line
arrow <- function(angle=30, length=unit(0.25, "inches"),
                  ends="last", type="open") {
    angle <- as.numeric(angle)
    if (!is.unit(length))
        stop("'length' must be a 'unit' object")
    ends <- as.integer(match(ends, c("first", "last", "both")))
    type <- as.integer(match(type, c("open", "closed")))
    if (anyNA(ends) || anyNA(type) ||
        length(ends) == 0 || length(type) == 0)
        stop("invalid 'ends' or 'type' argument")
    a <- list(angle=angle, length=length,
              ends=ends, type=type)
    class(a) <- "arrow"
    a
}

length.arrow <- function(x) {
    max(do.call("max", lapply(x, length)),
                length(x$length))
}

rep.arrow <- function(x, ...) {
    maxn <- length(x)
    newa <- list(angle=rep(x$angle, length.out=maxn),
                 length=rep(x$length, length.out=maxn),
                 ends=rep(x$ends, length.out=maxn),
                 type=rep(x$type, length.out=maxn))
    newa <- lapply(newa, rep, ...)
    class(newa) <- "arrow"
    newa
}

# Method for subsetting "arrow" objects
`[.arrow` <- function(x, index, ...) {
    if (length(index) == 0 ||
        (is.logical(index) && sum(index) == 0))
        return(NULL)
    maxn <- length(x)
    newa <- list(angle=rep(x$angle, length.out=maxn),
                 length=rep(x$length, length.out=maxn),
                 ends=rep(x$ends, length.out=maxn),
                 type=rep(x$type, length.out=maxn))
    newa <- lapply(X = newa, FUN = `[`, index, ...)
    class(newa) <- "arrow"
    newa
}

str.arrow <- function(object, ...) {
  o <- oldClass(object)
  oldClass(object) <- setdiff(o, "arrow")
  str(object)
}

######################################
# move-to and line-to primitives
######################################
validDetails.move.to <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  # Make sure that x and y are of length 1
  if (length(x$x) > 1 || length(x$y) > 1)
    stop("'x' and 'y' must have length 1")
  x
}

drawDetails.move.to <- function(x, recording=TRUE) {
  grid.Call.graphics(C_moveTo, x$x, x$y)
}

moveToGrob <- function(x=0, y=0,
                       default.units="npc",
                       name=NULL, vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y,
       name=name, vp=vp, cl="move.to")
}

grid.move.to <- function(x=0, y=0,
                         default.units="npc",
                         name=NULL, draw=TRUE, vp=NULL) {
  mtg <- moveToGrob(x=x, y=y, default.units=default.units,
                    name=name, vp=vp)
  if (draw)
    grid.draw(mtg)
  invisible(mtg)
}

validDetails.line.to <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  # Make sure that x and y are of length 1
  if (length(x$x) > 1 || length(x$y) > 1)
    stop("'x' and 'y' must have length 1")
  if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
      stop("invalid 'arrow' argument")
  x
}

drawDetails.line.to <- function(x, recording=TRUE) {
  grid.Call.graphics(C_lineTo, x$x, x$y, x$arrow)
}

lineToGrob <- function(x=1, y=1,
                       default.units="npc",
                       arrow=NULL,
                       name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y, arrow=arrow,
       name=name, gp=gp, vp=vp, cl="line.to")
}

grid.line.to <- function(x=1, y=1,
                         default.units="npc",
                         arrow=NULL,
                         name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  ltg <- lineToGrob(x=x, y=y, default.units=default.units, arrow=arrow,
                    name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(ltg)
  invisible(ltg)
}

######################################
# LINES primitive
######################################
validDetails.lines <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
      stop("invalid 'arrow' argument")
  x
}

drawDetails.lines <- function(x, recording=TRUE) {
    grid.Call.graphics(C_lines, x$x, x$y,
                       list(as.integer(1L:max(length(x$x), length(x$y)))),
                       x$arrow)
}

xDetails.lines <- function(x, theta) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.lines <- function(x, theta) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.lines <- function(x) {
  bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.lines <- function(x) {
  bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

linesGrob <- function(x=unit(c(0, 1), "npc"),
                      y=unit(c(0, 1), "npc"),
                      default.units="npc",
                      arrow=NULL,
                      name=NULL, gp=gpar(), vp=NULL) {
  # Allow user to specify unitless vector;  add default units
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y,
       arrow=arrow, name=name, gp=gp, vp=vp, cl="lines")
}

grid.lines <- function(x=unit(c(0, 1), "npc"),
                       y=unit(c(0, 1), "npc"),
                       default.units="npc",
                       arrow=NULL,
                       name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  lg <- linesGrob(x=x, y=y,
                  default.units=default.units, arrow=arrow,
                  name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(lg)
  invisible(lg)
}

######################################
# POLYLINES primitive
######################################
# Very similar to LINES primitive, but allows
# multiple polylines via 'id' and 'id.lengths' args
# as per POLYGON primitive
validDetails.polyline <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
      stop("'x' and 'y' must be units")
  if (!is.null(x$id) && !is.null(x$id.lengths))
      stop("it is invalid to specify both 'id' and 'id.lengths'")
  if (length(x$x) != length(x$y))
      stop("'x' and 'y' must be same length")
  if (!is.null(x$id) && (length(x$id) != length(x$x)))
      stop("'x' and 'y' and 'id' must all be same length")
  if (!is.null(x$id))
      x$id <- as.integer(x$id)
  if (!is.null(x$id.lengths) && (sum(x$id.lengths) != length(x$x)))
      stop("'x' and 'y' and 'id.lengths' must specify same overall length")
  if (!is.null(x$id.lengths))
      x$id.lengths <- as.integer(x$id.lengths)
  if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
      stop("invalid 'arrow' argument")
  x
}

drawDetails.polyline <- function(x, recording=TRUE) {
    if (is.null(x$id) && is.null(x$id.lengths))
        grid.Call.graphics(C_lines, x$x, x$y,
                           list(as.integer(seq_along(x$x))),
                           x$arrow)
    else {
        if (is.null(x$id)) {
            n <- length(x$id.lengths)
            id <- rep(1L:n, x$id.lengths)
        } else {
            n <- length(unique(x$id))
            id <- x$id
        }
        index <- split(as.integer(seq_along(x$x)), id)
        grid.Call.graphics(C_lines, x$x, x$y, index, x$arrow)
    }
}

xDetails.polyline <- function(x, theta) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.polyline <- function(x, theta) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.polyline <- function(x) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[3L], "inches")
}

heightDetails.polyline <- function(x) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[4L], "inches")
}

polylineGrob <- function(x=unit(c(0, 1), "npc"),
                         y=unit(c(0, 1), "npc"),
                         id=NULL, id.lengths=NULL,
                         default.units="npc",
                         arrow=NULL,
                         name=NULL, gp=gpar(), vp=NULL) {
    # Allow user to specify unitless vector;  add default units
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    grob(x=x, y=y, id=id, id.lengths=id.lengths,
         arrow=arrow, name=name, gp=gp, vp=vp, cl="polyline")
}

grid.polyline <- function(...) {
    grid.draw(polylineGrob(...))
}

######################################
# SEGMENTS primitive
######################################
validDetails.segments <- function(x) {
  if (!is.unit(x$x0) || !is.unit(x$x1) ||
      !is.unit(x$y0) || !is.unit(x$y1))
    stop("'x0', 'y0', 'x1', and 'y1' must be units")
  if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
      stop("invalid 'arrow' argument")
  x
}

drawDetails.segments <- function(x, recording=TRUE) {
  grid.Call.graphics(C_segments, x$x0, x$y0, x$x1, x$y1, x$arrow)
}

segmentBounds <- function(x, theta) {
    n <- max(length(x$x0), length(x$x1),
             length(x$y0), length(x$y1))
    x0 <- rep(x$x0, length.out=n)
    x1 <- rep(x$x1, length.out=n)
    y0 <- rep(x$y0, length.out=n)
    y1 <- rep(x$y1, length.out=n)
    grid.Call(C_locnBounds, unit.c(x0, x1), unit.c(y0, y1), theta)
}

xDetails.segments <- function(x, theta) {
    bounds <- segmentBounds(x, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.segments <- function(x, theta) {
    bounds <- segmentBounds(x, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.segments <- function(x) {
    bounds <- segmentBounds(x, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[3L], "inches")
}

heightDetails.segments <- function(x) {
    bounds <- segmentBounds(x, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[4L], "inches")
}

segmentsGrob <- function(x0=unit(0, "npc"), y0=unit(0, "npc"),
                         x1=unit(1, "npc"), y1=unit(1, "npc"),
                         default.units="npc",
                         arrow=NULL,
                         name=NULL, gp=gpar(), vp=NULL) {
  # Allow user to specify unitless vector;  add default units
  if (!is.unit(x0))
    x0 <- unit(x0, default.units)
  if (!is.unit(x1))
    x1 <- unit(x1, default.units)
  if (!is.unit(y0))
    y0 <- unit(y0, default.units)
  if (!is.unit(y1))
    y1 <- unit(y1, default.units)
  grob(x0=x0, y0=y0, x1=x1, y1=y1, arrow=arrow, name=name, gp=gp, vp=vp,
       cl="segments")
}

grid.segments <- function(x0=unit(0, "npc"), y0=unit(0, "npc"),
                          x1=unit(1, "npc"), y1=unit(1, "npc"),
                          default.units="npc",
                          arrow=NULL,
                          name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  sg <- segmentsGrob(x0=x0, y0=y0, x1=x1, y1=y1,
                     default.units=default.units,
                     arrow=arrow,
                     name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(sg)
  invisible(sg)
}

######################################
# ARROWS primitive
######################################

# Superceded by 'arrow' arg to line-drawing primitives
# which contains an "arrow" object
validDetails.arrows <- function(x) {
  if ((!is.null(x$x) && !is.unit(x$x)) ||
      (!is.null(x$y) && !is.unit(x$y)))
    stop("'x' and 'y' must be units or NULL")
  if (!is.unit(x$length))
    stop("'length' must be a 'unit' object")
  x$ends <- as.integer(match(x$ends, c("first", "last", "both")))
  x$type <- as.integer(match(x$type, c("open", "closed")))
  if (anyNA(x$ends) || anyNA(x$type))
    stop("invalid 'ends' or 'type' argument")
  x
}

drawDetails.arrows <- function(x, recording=TRUE) {
  if (is.null(x$x)) { # y should be null too
    if (!is.null(x$y))
      stop("corrupt 'arrows' object")
    lineThing <- getGrob(x, childNames(x))
    # This could be done via method dispatch, but that really
    # seemed like overkill
    # OTOH, this is NOT user-extensible
    # AND the code for, e.g., "lines" is not located with
    # the other grid.lines code so changes there are unlikely
    # to propagate to here (e.g., add an id arg to grid.lines?
    if (inherits(lineThing, "line.to")) {
      x1 <- NULL
      x2 <- lineThing$x
      y1 <- NULL
      y2 <- lineThing$y
      xnm1 <- NULL
      xn <- lineThing$x
      ynm1 <- NULL
      yn <- lineThing$y
    } else if (inherits(lineThing, "lines")) {
      # x or y may be recycled
      n <- max(length(lineThing$x),
               length(lineThing$y))
      xx <- rep(lineThing$x, length.out=2)
      x1 <- xx[1L]
      x2 <- xx[2L]
      xx <- rep(lineThing$x, length.out=n)
      xnm1 <- xx[n - 1]
      xn <- xx[n]
      yy <- rep(lineThing$y, length.out=2)
      y1 <- yy[1L]
      y2 <- yy[2L]
      yy <- rep(lineThing$y, length.out=n)
      ynm1 <- yy[n - 1]
      yn <- yy[n]
    } else { # inherits(lineThing, "segments")
      x1 <- lineThing$x0
      x2 <- lineThing$x1
      xnm1 <- lineThing$x0
      xn <- lineThing$x1
      y1 <- lineThing$y0
      y2 <- lineThing$y1
      ynm1 <- lineThing$y0
      yn <- lineThing$y1
    }
  } else {
    # x or y may be recycled
    n <- max(length(x$x), length(x$y))
    xx <- rep(x$x, length.out=2)
    x1 <- xx[1L]
    x2 <- xx[2L]
    xx <- rep(x$x, length.out=n)
    xnm1 <- xx[n - 1]
    xn <- xx[n]
    yy <- rep(x$y, length.out=2)
    y1 <- yy[1L]
    y2 <- yy[2L]
    yy <- rep(x$y, length.out=n)
    ynm1 <- yy[n - 1]
    yn <- yy[n]
    grid.Call.graphics(C_lines, x$x, x$y,
                       list(as.integer(1L:n)),
                       NULL)
  }
  grid.Call.graphics(C_arrows, x1, x2, xnm1, xn, y1, y2, ynm1, yn,
                     x$angle, x$length, x$ends, x$type)
}

widthDetails.arrows <- function(x) {
  if (is.null(x$x)) { # y should be null too
    if (!is.null(x$y))
      stop("corrupt 'arrows' object")
    lineThing <- getGrob(x, childNames(x))
    widthDetails(lineThing)
  } else {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
      unit(0, "inches")
    else
      unit(bounds[3L], "inches")
  }
}

heightDetails.arrows <- function(x) {
  if (is.null(x$x)) { # y should be null too
    if (!is.null(x$y))
      stop("corrupt 'arrows' object")
    lineThing <- getGrob(x, childNames(x))
    heightDetails(lineThing)
  } else {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
      unit(0, "inches")
    else
      unit(bounds[4L], "inches")
  }
}

arrowsGrob <- function(x=c(0.25, 0.75), y=0.5,
                       default.units="npc",
                       grob=NULL,
                       angle=30, length=unit(0.25, "inches"),
                       ends="last", type="open",
                       name=NULL, gp=gpar(), vp=NULL) {
    .Defunct(msg="'arrowsGrob' is defunct; use 'arrow' arguments to line drawing functions")
}

grid.arrows <- function(x=c(0.25, 0.75), y=0.5,
                        default.units="npc",
                        grob=NULL,
                        angle=30, length=unit(0.25, "inches"),
                        ends="last", type="open",
                        name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
    .Defunct(msg="'grid.arrows' is defunct; use 'arrow' arguments to line drawing functions")
}

######################################
# POLYGON primitive
######################################

validDetails.polygon <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  if (!is.null(x$id) && !is.null(x$id.lengths))
    stop("it is invalid to specify both 'id' and 'id.lengths'")
  if (length(x$x) != length(x$y))
    stop("'x' and 'y' must be same length")
  if (!is.null(x$id) && (length(x$id) != length(x$x)))
    stop("'x' and 'y' and 'id' must all be same length")
  if (!is.null(x$id))
    x$id <- as.integer(x$id)
  if (!is.null(x$id.lengths) && (sum(x$id.lengths) != length(x$x)))
    stop("'x' and 'y' and 'id.lengths' must specify same overall length")
  if (!is.null(x$id.lengths))
    x$id.lengths <- as.integer(x$id.lengths)
  x
}

drawDetails.polygon <- function(x, recording=TRUE) {
  if (is.null(x$id) && is.null(x$id.lengths))
    grid.Call.graphics(C_polygon, x$x, x$y,
                       list(as.integer(seq_along(x$x))))
  else {
    if (is.null(x$id)) {
      n <- length(x$id.lengths)
      id <- rep(1L:n, x$id.lengths)
    } else {
      n <- length(unique(x$id))
      id <- x$id
    }
    index <- split(as.integer(seq_along(x$x)), id)
    grid.Call.graphics(C_polygon, x$x, x$y, index)
  }
}

xDetails.polygon <- function(x, theta) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.polygon <- function(x, theta) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.polygon <- function(x) {
  bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.polygon <- function(x) {
  bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

polygonGrob <- function(x=c(0, 0.5, 1, 0.5), y=c(0.5, 1, 0.5, 0),
                        id=NULL, id.lengths=NULL,
                        default.units="npc",
                        name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y, id=id,
       id.lengths=id.lengths,
       name=name, gp=gp, vp=vp, cl="polygon")
}

grid.polygon <- function(x=c(0, 0.5, 1, 0.5), y=c(0.5, 1, 0.5, 0),
                         id=NULL, id.lengths=NULL,
                         default.units="npc",
                         name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  pg <- polygonGrob(x=x, y=y, id=id, id.lengths=id.lengths,
                    default.units=default.units,
                    name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(pg)
  invisible(pg)
}

######################################
# PATH primitive
######################################

validDetails.pathgrob <- function(x) {
    if (!is.unit(x$x) || !is.unit(x$y))
        stop("'x' and 'y' must be units")
    if (!is.null(x$id) && !is.null(x$id.lengths))
        stop("it is invalid to specify both 'id' and 'id.lengths'")
    if (length(x$x) != length(x$y))
        stop("'x' and 'y' must be same length")
    if (!is.null(x$id) && (length(x$id) != length(x$x)))
        stop("'x' and 'y' and 'id' must all be same length")
    if (!is.null(x$id))
        x$id <- as.integer(x$id)
    if (!is.null(x$pathId))
    	x$pathId <- as.integer(x$pathId)
    if (!is.null(x$id.lengths) && (sum(x$id.lengths) != length(x$x)))
        stop("'x' and 'y' and 'id.lengths' must specify same overall length")
    if (!is.null(x$pathId.lengths) && (sum(x$pathId.lengths) != length(x$x)))
    	stop("'x' and 'y' and 'pathId.lengths' must specify same overall length")
    if (!is.null(x$id.lengths))
        x$id.lengths <- as.integer(x$id.lengths)
    if (!is.null(x$pathId.lengths))
    	x$pathId.lengths <- as.integer(x$pathId.lengths)
    x
}

xDetails.pathgrob <- function(x, theta) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.pathgrob <- function(x, theta) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.pathgrob <- function(x) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[3L], "inches")
}

heightDetails.pathgrob <- function(x) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[4L], "inches")
}


drawDetails.pathgrob <- function(x, recording=TRUE) {
    hasMultiple <- !(is.null(x$pathId) && is.null(x$pathId.lengths))
    if (hasMultiple) {
        if (is.null(x$pathId)) {
            n <- length(x$pathId.lengths)
            pathId <- rep(1L:n, x$pathId.lengths)
        } else {
            pathId <- x$pathId
        }
    }
    if (is.null(x$id) && is.null(x$id.lengths)) {
        if (hasMultiple) {
            grid.Call.graphics(C_polygon, x$x, x$y,
                               split(as.integer(seq_along(x$x)), pathId))
        } else {
            grid.Call.graphics(C_polygon, x$x, x$y,
                               list(as.integer(seq_along(x$x))))
        }
    } else {
        if (is.null(x$id)) {
            n <- length(x$id.lengths)
            id <- rep(1L:n, x$id.lengths)
        } else {
            n <- length(unique(x$id))
            id <- x$id
        }
        if (hasMultiple) {
            index <- mapply(split,
                            x=split(as.integer(seq_along(x$x)), pathId), 
                            f=split(id, pathId),
                            SIMPLIFY = FALSE, USE.NAMES = FALSE)
        } else {
            index <- list(split(as.integer(seq_along(x$x)), id))
        }
        grid.Call.graphics(C_path, x$x, x$y, index,
                           switch(x$rule, winding=1L, evenodd=0L))
    }
}

pathGrob <- function(x, y,
                     id=NULL, id.lengths=NULL,
                     pathId=NULL, pathId.lengths=NULL,
                     rule="winding",
                     default.units="npc",
                     name=NULL, gp=gpar(), vp=NULL) {
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    grob(x=x, y=y, id=id, id.lengths=id.lengths,
         pathId=pathId, pathId.lengths=pathId.lengths,
         rule=rule,
         name=name, gp=gp, vp=vp, cl="pathgrob")
}

grid.path <- function(...) {
  grid.draw(pathGrob(...))
}

######################################
# XSPLINE primitive
######################################

validDetails.xspline <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("x and y must be units")
  if (!is.null(x$id) && !is.null(x$id.lengths))
    stop("it is invalid to specify both 'id' and 'id.lengths'")
  nx <- length(x$x)
  ny <- length(x$y)
  if (nx != ny)
    stop("'x' and 'y' must be same length")
  if (!is.null(x$id) && (length(x$id) != nx))
    stop("'x' and 'y' and 'id' must all be same length")
  if (!is.null(x$id))
    x$id <- as.integer(x$id)
  if (!is.null(x$id.lengths) && (sum(x$id.lengths) != nx))
    stop("'x' and 'y' and 'id.lengths' must specify same overall length")
  if (!is.null(x$id.lengths))
    x$id.lengths <- as.integer(x$id.lengths)
  if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
      stop("invalid 'arrow' argument")
  if (any(x$shape < -1 | x$shape > 1))
    stop("'shape' must be between -1 and 1")
  x$open <- as.logical(x$open)
  # Force all first and last shapes to be 0 for open xsplines
  if (x$open) {
      x$shape <- rep(x$shape, length.out=nx)
      # Watch out for id or id.length!
      index <- xsplineIndex(x)
      first <- sapply(index, min)
      last <- sapply(index, max)
      x$shape[c(first, last)] <- 0
  }
  x
}

xsplineIndex <- function(x) {
  if (is.null(x$id) && is.null(x$id.lengths))
      list(as.integer(seq_along(x$x)))
  else {
    if (is.null(x$id)) {
      n <- length(x$id.lengths)
      id <- rep(1L:n, x$id.lengths)
    } else {
      n <- length(unique(x$id))
      id <- x$id
    }
    split(as.integer(seq_along(x$x)), id)
  }
}

drawDetails.xspline <- function(x, recording=TRUE) {
    grid.Call.graphics(C_xspline, x$x, x$y, x$shape, x$open, x$arrow,
                       x$repEnds, xsplineIndex(x))
}

xDetails.xspline <- function(x, theta) {
  bounds <- grid.Call(C_xsplineBounds, x$x, x$y, x$shape, x$open, x$arrow,
                      x$repEnds, xsplineIndex(x), theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[1L], "inches")
}

yDetails.xspline <- function(x, theta) {
  bounds <- grid.Call(C_xsplineBounds, x$x, x$y, x$shape, x$open, x$arrow,
                      x$repEnds, xsplineIndex(x), theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[2L], "inches")
}

widthDetails.xspline <- function(x) {
  bounds <- grid.Call(C_xsplineBounds, x$x, x$y, x$shape, x$open, x$arrow,
                      x$repEnds, list(as.integer(seq_along(x$x))), 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.xspline <- function(x) {
  bounds <- grid.Call(C_xsplineBounds, x$x, x$y, x$shape, x$open, x$arrow,
                      x$repEnds, list(as.integer(seq_along(x$x))), 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

xsplineGrob <- function(x=c(0, 0.5, 1, 0.5), y=c(0.5, 1, 0.5, 0),
                        id=NULL, id.lengths=NULL,
                        default.units="npc",
                        shape=0, open=TRUE, arrow=NULL, repEnds=TRUE,
                        name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y, shape=shape, open=open,
       id=id, id.lengths=id.lengths, arrow=arrow, repEnds=repEnds,
       name=name, gp=gp, vp=vp, cl="xspline")
}

grid.xspline <- function(...) {
  grid.draw(xsplineGrob(...))
}

xsplinePoints <- function(x) {
    # Mimic drawGrob() to ensure x$vp and x$gp enforced
    dlon <- grid.Call(C_setDLon, FALSE)
    on.exit(grid.Call(C_setDLon, dlon))
    tempgpar <- grid.Call(C_getGPar)
    on.exit(grid.Call(C_setGPar, tempgpar), add=TRUE)
    preDraw(x)
    # Raw pts in dev coords
    devPoints <- grid.Call(C_xsplinePoints,
                           x$x, x$y, x$shape, x$open, x$arrow,
                           x$repEnds, xsplineIndex(x), 0)
    postDraw(x)
    # Convert to units in inches
    unitPoints <- lapply(devPoints,
                         function(x) {
                             names(x) <- c("x", "y")
                             x$x <- unit(x$x, "inches")
                             x$y <- unit(x$y, "inches")
                             x
                         })
    if (length(unitPoints) == 1)
        unitPoints <- unitPoints[[1]]
    unitPoints
}

######################################
# BEZIER primitive
######################################

# A bezier grob that works of a (not-100% accurate) approximation
# using X-splines

# X-Spline approx to Bezier
Ms <- 1/6*rbind(c(1, 4, 1, 0),
                c(-3, 0, 3, 0),
                c(3, -6, 3, 0),
                c(-1, 3, -3, 1))
Msinv <- solve(Ms)
# Bezier control matrix
Mb <- rbind(c(1, 0, 0, 0),
            c(-3, 3, 0, 0),
            c(3, -6, 3, 0),
            c(-1, 3, -3, 1))

splinePoints <- function(xb, yb, idIndex) {
    xs <- unlist(lapply(idIndex,
                        function(i) {
                            Msinv %*% Mb %*% xb[i]
                        }))
    ys <- unlist(lapply(idIndex,
                        function(i) {
                            Msinv %*% Mb %*% yb[i]
                        }))
    list(x=xs, y=ys)
}

splinegrob <- function(x) {
    xx <- convertX(x$x, "inches", valueOnly=TRUE)
    yy <- convertY(x$y, "inches", valueOnly=TRUE)
    sp <- splinePoints(xx, yy, xsplineIndex(x))
    xsplineGrob(sp$x, sp$y, default.units="inches",
                id=x$id, id.lengths=x$id.lengths,
                shape=1, repEnds=FALSE,
                arrow=x$arrow, name=x$name,
                gp=x$gp, vp=x$vp)
}

validDetails.beziergrob <- function(x) {
    if (!is.unit(x$x) ||
        !is.unit(x$y))
        stop("x and y must be units")
    if (!is.null(x$id) && !is.null(x$id.lengths))
        stop("it is invalid to specify both 'id' and 'id.lengths'")
    nx <- length(x$x)
    ny <- length(x$y)
    if (nx != ny)
        stop("'x' and 'y' must be same length")
    if (!is.null(x$id) && (length(x$id) != nx))
        stop("'x' and 'y' and 'id' must all be same length")
    if (!is.null(x$id))
        x$id <- as.integer(x$id)
    if (!is.null(x$id.lengths) && (sum(x$id.lengths) != nx))
        stop("'x' and 'y' and 'id.lengths' must specify same overall length")
    if (!is.null(x$id.lengths))
        x$id.lengths <- as.integer(x$id.lengths)
    if (is.null(x$id) && is.null(x$id.lengths)) {
        if (length(x$x) != 4L)
            stop("must have exactly 4 control points")
    } else {
        if (is.null(x$id)) {
            n <- length(x$id.lengths)
            id <- rep(1L:n, x$id.lengths)
        } else {
            id <- x$id
        }
        xper <- split(x$x, id)
        if (any(lengths(xper) != 4L))
            stop("must have exactly 4 control points per Bezier curve")
    }
    if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
        stop("invalid 'arrow' argument")
    x
}

makeContent.beziergrob <- function(x) {
    splinegrob(x)
}

xDetails.beziergrob <- function(x, theta) {
    xDetails(splinegrob(x), theta)
}

yDetails.beziergrob <- function(x, theta) {
    yDetails(splinegrob(x), theta)
}

widthDetails.beziergrob <- function(x) {
    widthDetails(splinegrob(x))
}

heightDetails.beziergrob <- function(x) {
    heightDetails(splinegrob(x))
}

bezierGrob <- function(x=c(0, 0.5, 1, 0.5), y=c(0.5, 1, 0.5, 0),
                       id=NULL, id.lengths=NULL,
                       default.units="npc", arrow=NULL,
                       name=NULL, gp=gpar(), vp=NULL) {
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    grob(x=x, y=y,
         id=id, id.lengths=id.lengths, arrow=arrow,
         name=name, gp=gp, vp=vp, cl="beziergrob")
}

grid.bezier <- function(...) {
    grid.draw(bezierGrob(...))
}

bezierPoints <- function(x) {
    sg <- splinegrob(x)
    # splinegrob() does not make use of x$vp
    sg$vp <- x$vp
    xsplinePoints(sg)
}


######################################
# CIRCLE primitive
######################################

validDetails.circle <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y) ||
      !is.unit(x$r))
    stop("'x', 'y', and 'r' must be units")
  x
}

drawDetails.circle <- function(x, recording=TRUE) {
  grid.Call.graphics(C_circle, x$x, x$y, x$r)
}

xDetails.circle <- function(x, theta) {
  bounds <- grid.Call(C_circleBounds, x$x, x$y, x$r, theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[1L], "inches")
}

yDetails.circle <- function(x, theta) {
  bounds <- grid.Call(C_circleBounds, x$x, x$y, x$r, theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[2L], "inches")
}

widthDetails.circle <- function(x) {
  bounds <- grid.Call(C_circleBounds, x$x, x$y, x$r, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.circle <- function(x) {
  bounds <- grid.Call(C_circleBounds, x$x, x$y, x$r, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

circleGrob <- function(x=0.5, y=0.5, r=0.5,
                       default.units="npc",
                       name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  if (!is.unit(r))
    r <- unit(r, default.units)
  grob(x=x, y=y, r=r, name=name, gp=gp, vp=vp, cl="circle")
}

grid.circle <- function(x=0.5, y=0.5, r=0.5,
                        default.units="npc",
                        name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  cg <- circleGrob(x=x, y=y, r=r,
                   default.units=default.units,
                   name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(cg)
  invisible(cg)
}

######################################
# RECT primitive
######################################
validDetails.rect <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y) ||
      !is.unit(x$width) ||
      !is.unit(x$height))
    stop("'x', 'y', 'width', and 'height' must be units")
  valid.just(x$just)
  if (!is.null(x$hjust))
    x$hjust <- as.numeric(x$hjust)
  if (!is.null(x$vjust))
    x$vjust <- as.numeric(x$vjust)
  x
}

drawDetails.rect <- function(x, recording=TRUE) {
    grid.Call.graphics(C_rect, x$x, x$y, x$width, x$height,
                       resolveHJust(x$just, x$hjust),
                       resolveVJust(x$just, x$vjust))
}

xDetails.rect <- function(x, theta) {
  bounds <- grid.Call(C_rectBounds, x$x, x$y, x$width, x$height,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[1L], "inches")
}

yDetails.rect <- function(x, theta) {
  bounds <- grid.Call(C_rectBounds, x$x, x$y, x$width, x$height,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[2L], "inches")
}

widthDetails.rect <- function(x) {
  bounds <- grid.Call(C_rectBounds, x$x, x$y, x$width, x$height,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.rect <- function(x) {
  bounds <- grid.Call(C_rectBounds, x$x, x$y, x$width, x$height,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

rectGrob <- function(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                     width=unit(1, "npc"), height=unit(1, "npc"),
                     just="centre", hjust=NULL, vjust=NULL,
                     default.units="npc",
                     name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  if (!is.unit(width))
    width <- unit(width, default.units)
  if (!is.unit(height))
    height <- unit(height, default.units)
  grob(x=x, y=y, width=width, height=height, just=just,
       hjust=hjust, vjust=vjust,
       name=name, gp=gp, vp=vp, cl="rect")
}

grid.rect <- function(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                      width=unit(1, "npc"), height=unit(1, "npc"),
                      just="centre", hjust=NULL, vjust=NULL,
                      default.units="npc",
                      name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  rg <- rectGrob(x=x, y=y, width=width, height=height, just=just,
                 hjust=hjust, vjust=vjust,
                 default.units=default.units,
                 name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(rg)
  invisible(rg)
}

######################################
# RASTER primitive
######################################

validDetails.rastergrob <- function(x) {
    if (!(is.raster(x$raster) || inherits(x$raster, "nativeRaster")))
        x$raster <- as.raster(x$raster)
    if (!is.unit(x$x) ||
        !is.unit(x$y) ||
        (!is.null(x$width) && !is.unit(x$width)) ||
        (!is.null(x$height) && !is.unit(x$height)))
        stop("'x', 'y', 'width', and 'height' must be units")
    valid.just(x$just)
    if (!is.null(x$hjust))
        x$hjust <- as.numeric(x$hjust)
    if (!is.null(x$vjust))
        x$vjust <- as.numeric(x$vjust)
    x
}

resolveRasterSize <- function(x) {
    if (is.null(x$width)) {
        if (is.null(x$height)) {
            rasterRatio <- dim(x$raster)[1]/dim(x$raster)[2]
            vpWidth <- convertWidth(unit(1, "npc"), "inches", valueOnly=TRUE)
            vpHeight <- convertHeight(unit(1, "npc"), "inches", valueOnly=TRUE)
            vpRatio <- vpHeight/vpWidth
            if (rasterRatio > vpRatio) {
                x$height <- unit(vpHeight, "inches")
                x$width <- unit(vpHeight*dim(x$raster)[2]/dim(x$raster)[1],
                                "inches")
            } else {
                x$width <- unit(vpWidth, "inches")
                x$height <- unit(vpWidth*dim(x$raster)[1]/dim(x$raster)[2],
                                 "inches")
            }
        } else {
            h <- convertHeight(x$height, "inches", valueOnly=TRUE)
            x$width <- unit(h*dim(x$raster)[2]/dim(x$raster)[1],
                            "inches")
        }
    } else {
        if (is.null(x$height)) {
            w <- convertWidth(x$width, "inches", valueOnly=TRUE)
            x$height <- unit(w*dim(x$raster)[1]/dim(x$raster)[2],
                             "inches")
        }
    }
    x
}

drawDetails.rastergrob <- function(x, recording=TRUE) {
    # At this point resolve NULL width/height based on
    # image dimensions
    x <- resolveRasterSize(x)
    if (is.null(x$width)) {
        if (is.null(x$height)) {
            rasterRatio <- dim(x$raster)[1]/dim(x$raster)[2]
            vpWidth <- convertWidth(unit(1, "npc"), "inches", valueOnly=TRUE)
            vpHeight <- convertHeight(unit(1, "npc"), "inches", valueOnly=TRUE)
            vpRatio <- vpHeight/vpWidth
            if (rasterRatio > vpRatio) {
                x$height <- unit(vpHeight, "inches")
                x$width <- unit(vpHeight*dim(x$raster)[2]/dim(x$raster)[1],
                                "inches")
            } else {
                x$width <- unit(vpWidth, "inches")
                x$height <- unit(vpWidth*dim(x$raster)[1]/dim(x$raster)[2],
                                 "inches")
            }
        } else {
            h <- convertHeight(x$height, "inches", valueOnly=TRUE)
            x$width <- unit(h*dim(x$raster)[2]/dim(x$raster)[1],
                            "inches")
        }
    } else {
        if (is.null(x$height)) {
            w <- convertWidth(x$width, "inches", valueOnly=TRUE)
            x$height <- unit(w*dim(x$raster)[1]/dim(x$raster)[2],
                             "inches")
        }
    }
    grid.Call.graphics(C_raster, x$raster,
                       x$x, x$y, x$width, x$height,
                       resolveHJust(x$just, x$hjust),
                       resolveVJust(x$just, x$vjust),
                       x$interpolate)
}

xDetails.rastergrob <- function(x, theta) {
    x <- resolveRasterSize(x)
    bounds <- grid.Call(C_rectBounds, x$x, x$y, x$width, x$height,
                        resolveHJust(x$just, x$hjust),
                        resolveVJust(x$just, x$vjust),
                        theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.rastergrob <- function(x, theta) {
    x <- resolveRasterSize(x)
    bounds <- grid.Call(C_rectBounds, x$x, x$y, x$width, x$height,
                        resolveHJust(x$just, x$hjust),
                        resolveVJust(x$just, x$vjust),
                        theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.rastergrob <- function(x) {
    x <- resolveRasterSize(x)
    bounds <- grid.Call(C_rectBounds, x$x, x$y, x$width, x$height,
                        resolveHJust(x$just, x$hjust),
                        resolveVJust(x$just, x$vjust),
                        0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[3L], "inches")
}

heightDetails.rastergrob <- function(x) {
    x <- resolveRasterSize(x)
    bounds <- grid.Call(C_rectBounds, x$x, x$y, x$width, x$height,
                        resolveHJust(x$just, x$hjust),
                        resolveVJust(x$just, x$vjust),
                        0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[4L], "inches")
}

rasterGrob <- function(image,
                       x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                       width=NULL, height=NULL,
                       just="centre", hjust=NULL, vjust=NULL,
                       interpolate=TRUE,
                       default.units="npc",
                       name=NULL, gp=gpar(), vp=NULL) {

    if (inherits(image, "nativeRaster"))
        raster <- image
    else
        raster <- as.raster(image)
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    if (!is.null(width) && !is.unit(width))
        width <- unit(width, default.units)
    if (!is.null(height) && !is.unit(height))
        height <- unit(height, default.units)
    grob(raster=raster, x=x, y=y, width=width, height=height, just=just,
         hjust=hjust, vjust=vjust, interpolate=interpolate,
         name=name, gp=gp, vp=vp, cl="rastergrob")
}

grid.raster <- function(image,
                        x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                        width=NULL, height=NULL,
                        just="centre", hjust=NULL, vjust=NULL,
                        interpolate=TRUE,
                        default.units="npc",
                        name=NULL, gp=gpar(), vp=NULL) {
    rg <- rasterGrob(image,
                     x=x, y=y, width=width, height=height, just=just,
                     hjust=hjust, vjust=vjust, interpolate=interpolate,
                     default.units=default.units,
                     name=name, gp=gp, vp=vp)
    grid.draw(rg)
}

######################################
# TEXT primitive
######################################
validDetails.text <- function(x) {
  if (!is.language(x$label))
    x$label <- as.character(x$label)
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  x$rot <- as.numeric(x$rot)
  if (!all(is.finite(x$rot)) || length(x$rot) == 0)
    stop("invalid 'rot' value")
  valid.just(x$just)
  if (!is.null(x$hjust))
    x$hjust <- as.numeric(x$hjust)
  if (!is.null(x$vjust))
    x$vjust <- as.numeric(x$vjust)
  x$check.overlap <- as.logical(x$check.overlap)
  x
}

drawDetails.text <- function(x, recording=TRUE) {
  grid.Call.graphics(C_text, as.graphicsAnnot(x$label),
                     x$x, x$y,
                     resolveHJust(x$just, x$hjust),
                     resolveVJust(x$just, x$vjust),
                     x$rot, x$check.overlap)
}

xDetails.text <- function(x, theta) {
  bounds <- grid.Call(C_textBounds, as.graphicsAnnot(x$label),
                      x$x, x$y,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      x$rot, theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[1L], "inches")
}

yDetails.text <- function(x, theta) {
  bounds <- grid.Call(C_textBounds, as.graphicsAnnot(x$label),
                      x$x, x$y,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      x$rot, theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[2L], "inches")
}

widthDetails.text <- function(x) {
  bounds <- grid.Call(C_textBounds, as.graphicsAnnot(x$label),
                      x$x, x$y,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      x$rot, as.numeric(NA))
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.text <- function(x) {
  bounds <- grid.Call(C_textBounds, as.graphicsAnnot(x$label),
                      x$x, x$y,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      x$rot, as.numeric(NA))
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

ascentDetails.text <- function(x) {
    if (length(x$label) == 1) {
        metrics <- grid.Call(C_stringMetric, as.graphicsAnnot(x$label))
        unit(metrics[[1]], "inches")
    } else {
        heightDetails(x)
    }
}

descentDetails.text <- function(x) {
    if (length(x$label) == 1) {
        metrics <- grid.Call(C_stringMetric, as.graphicsAnnot(x$label))
        unit(metrics[[2]], "inches")
    } else {
        unit(0, "inches")
    }
}

textGrob <- function(label, x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                     just="centre", hjust=NULL, vjust=NULL,
                     rot=0, check.overlap=FALSE,
                     default.units="npc",
                     name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(label=label, x=x, y=y, just=just, hjust=hjust, vjust=vjust,
       rot=rot, check.overlap=check.overlap,
       name=name, gp=gp, vp=vp, cl="text")
}

grid.text <- function(label, x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                      just="centre", hjust=NULL, vjust=NULL,
                      rot=0, check.overlap=FALSE,
                      default.units="npc",
                      name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  tg <- textGrob(label=label, x=x, y=y, just=just,
                 hjust=hjust, vjust=vjust, rot=rot,
                 check.overlap=check.overlap,
                 default.units=default.units,
                 name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(tg)
  invisible(tg)
}

######################################
# POINTS primitive
######################################
valid.pch <- function(pch) {
  if (length(pch) == 0L)
    stop("zero-length 'pch'")
  if (is.null(pch))
    pch <- 1L
  else if (!is.character(pch))
    pch <- as.integer(pch)
  pch
}

validDetails.points <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y) ||
      !is.unit(x$size))
    stop("'x', 'y' and 'size' must be units")
  if (length(x$x) != length(x$y))
    stop("'x' and 'y' must be 'unit' objects and have the same length")
  x$pch <- valid.pch(x$pch)
  x
}

drawDetails.points <- function(x, recording=TRUE) {
  grid.Call.graphics(C_points, x$x, x$y, x$pch, x$size)
}

# FIXME:  does not take into account the size of the symbols
xDetails.points <- function(x, theta) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.points <- function(x, theta) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.points <- function(x) {
  bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.points <- function(x) {
  bounds <- grid.Call(C_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

pointsGrob <- function(x=stats::runif(10),
                       y=stats::runif(10),
                       pch=1, size=unit(1, "char"),
                       default.units="native",
                       name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y, pch=pch, size=size,
       name=name, gp=gp, vp=vp, cl="points")
}

grid.points <- function(x=stats::runif(10),
                        y=stats::runif(10),
                        pch=1, size=unit(1, "char"),
                        default.units="native",
                        name=NULL, gp=gpar(),
                        draw=TRUE, vp=NULL) {
  pg <- pointsGrob(x=x, y=y, pch=pch, size=size,
                   default.units=default.units,
                   name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(pg)
  invisible(pg)
}

######################################
# CLIP primitive
######################################
validDetails.clip <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y) ||
      !is.unit(x$width) ||
      !is.unit(x$height))
    stop("'x', 'y', 'width', and 'height' must be units")
  if (length(x$x) > 1 || length(x$y) > 1 ||
      length(x$width) > 1 || length(x$height) > 1)
    stop("'x', 'y', 'width', and 'height' must all be units of length 1")
  valid.just(x$just)
  if (!is.null(x$hjust))
    x$hjust <- as.numeric(x$hjust)
  if (!is.null(x$vjust))
    x$vjust <- as.numeric(x$vjust)
  x
}

drawDetails.clip <- function(x, recording=TRUE) {
  grid.Call.graphics(C_clip, x$x, x$y, x$width, x$height,
                     resolveHJust(x$just, x$hjust),
                     resolveVJust(x$just, x$vjust))
}

clipGrob <- function(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                     width=unit(1, "npc"), height=unit(1, "npc"),
                     just="centre", hjust=NULL, vjust=NULL,
                     default.units="npc",
                     name=NULL, vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  if (!is.unit(width))
    width <- unit(width, default.units)
  if (!is.unit(height))
    height <- unit(height, default.units)
  grob(x=x, y=y, width=width, height=height, just=just,
       hjust=hjust, vjust=vjust,
       name=name, vp=vp, cl="clip")
}

grid.clip <- function(...) {
  grid.draw(clipGrob(...))
}


######################################
# NULL primitive
######################################

validDetails.null <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  if (length(x$x) > 1 || length(x$y) > 1)
    stop("'x' and 'y' must all be units of length 1")
  x
}

drawDetails.null <- function(x, recording=TRUE) {
    # Deliberate null op.
    # NOTE: nothing will go on the graphics engine DL
    # This is ok I think because these grobs are only
    # useful on the grid DL (for other grid code to query
    # their size or location).
}

xDetails.null <- function(x, theta) {
    bounds <- grid.Call(C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.null <- function(x, theta) {
    bounds <- grid.Call( C_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

# Deliberately ZERO
widthDetails.null <- function(x) {
    unit(0, "inches")
}

heightDetails.null <- function(x) {
    unit(0, "inches")
}

# A grob with GUARANTEED zero-width
# also GUARANTEED NOT to draw anything
nullGrob <- function(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                     default.units="npc",
                     name=NULL, vp=NULL) {
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    grob(x=x, y=y, name=name, vp=vp, cl="null")
}

# Convenient way to get nullGrob on the grid display list
grid.null <- function(...) {
    grid.draw(nullGrob(...))
}
#  File src/library/grid/R/roundrect.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2016 The R Core Team
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


# Good idea to choose r as absolute unit or "snpc"
roundrectGrob <- function(x=0.5, y=0.5, width=1, height=1,
                          default.units="npc",
                          r=unit(0.1, "snpc"),
                          just="centre",
                          name=NULL, gp=NULL, vp=NULL) {
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    if (!is.unit(width))
        width <- unit(width, default.units)
    if (!is.unit(height))
        height <- unit(height, default.units)
    grob(x=x, y=y, width=width, height=height, r=r, just=just,
         name=name, gp=gp, vp=vp, cl="roundrect")
}

grid.roundrect <- function(...) {
  grid.draw(roundrectGrob(...))
}

validDetails.roundrect <- function(x) {
    if (!is.unit(x$x) ||
        !is.unit(x$y) ||
        !is.unit(x$width) ||
        !is.unit(x$height))
        stop("'x', 'y', 'width', and 'height' must be units")
    if (!is.unit(x$r))
        stop("'r' must be a 'unit' object")
    valid.just(x$just)
    # Make sure that x and y are of length 1
    if (length(x$x) != 1 || length(x$y) != 1 ||
        length(x$width) != 1 || length(x$height) != 1)
        stop("'x', 'y', 'width', and 'height' must have length 1")
    x
}

makeContext.roundrect <- function(x) {
    rrvp <- viewport(x$x, x$y, x$width, x$height, just=x$just,
                     name="rrvp")
    if (!is.null(x$vp)) {
        x$vp <- vpStack(x$vp, rrvp)
    } else {
        x$vp <- rrvp
    }
    x
}

# x, y, is the real corner
roundCorner <- function(num, x, y, r) {
  n <- 10*4
  t <- seq(0, 2*pi, length.out=n)
  cost <- cos(t)
  sint <- sin(t)
  if (num == 1) {
    xc <- x + r
    yc <- y + r
    subset <- (n/2):(3*n/4)
  } else if (num == 2) {
    xc <- x + r
    yc <- y - r
    subset <- (n/4):(n/2)
  } else if (num == 3) {
    xc <- x - r
    yc <- y - r
    subset <- 1L:(n/4)
  } else if (num == 4) {
    xc <- x - r
    yc <- y + r
    subset <- (3*n/4):n
  }
  list(x=xc + (cost*r)[subset], y=yc + (sint*r)[subset])
}

rrpoints <- function(x) {
  left <- 0
  bottom <- 0
  right <- convertX(unit(1, "npc"), "inches", valueOnly=TRUE)
  top <- convertY(unit(1, "npc"), "inches", valueOnly=TRUE)
  r <- min(convertWidth(x$r, "inches", valueOnly=TRUE),
           convertHeight(x$r, "inches", valueOnly=TRUE))
  corner1 <- roundCorner(1, left, bottom, r)
  corner2 <- roundCorner(2, left, top, r)
  corner3 <- roundCorner(3, right, top, r)
  corner4 <- roundCorner(4, right, bottom, r)
  xx <- unit(c(left + r, right - r, corner4$x,
               right, right, corner3$x,
               right - r, left + r, corner2$x,
               left, left, corner1$x),
             "inches")
  yy <- unit(c(bottom, bottom, corner4$y,
               bottom + r, top - r, corner3$y,
               top, top, corner2$y,
               top - r, bottom + r, corner1$y),
             "inches")
  list(x=xx, y=yy)
}

makeContent.roundrect <- function(x) {
    boundary <- rrpoints(x)
    polygonGrob(boundary$x, boundary$y,
                name=x$name, gp=x$gp, vp=x$vp)
}

xDetails.roundrect <- function(x, theta) {
    boundary <- rrpoints(x)
    bounds <- grid.Call(C_locnBounds, boundary$x, boundary$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.roundrect <- function(x, theta) {
    boundary <- rrpoints(x)
    bounds <- grid.Call(C_locnBounds, boundary$x, boundary$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.roundrect <- function(x) {
    boundary <- rrpoints(x)
    bounds <- grid.Call(C_locnBounds, boundary$x, boundary$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[3L], "inches")
}

heightDetails.roundrect <- function(x) {
    boundary <- rrpoints(x)
    bounds <- grid.Call(C_locnBounds, boundary$x, boundary$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[4L], "inches")
}



#  File src/library/grid/R/size.R
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

# These functions are used to evaluate "grobwidth" and
# "grobheight" units.
# They are usually called from within the C code
# (specifically, from within unit.c)
# It should be noted that they only give the width/height
# of the grob in the current drawing context
# (i.e., evaluating the width/height in another context
#  will not necessarily give the same result)

# The C code to evaluate "grobwidth" and "grobheight" calls
# the preDrawDetails and postDrawDetails generics before and
# after the call to width/height() to allow for complex grobs which
# construct their own viewports.

#########
# X locations on edge
#########

xDetails <- function(x, theta) {
  UseMethod("xDetails")
}

xDetails.default <- function(x, theta) {
  unit(0.5, "npc")
}

#########
# Y locations on edge
#########

yDetails <- function(x, theta) {
  UseMethod("yDetails")
}

yDetails.default <- function(x, theta) {
  unit(0.5, "npc")
}

#########
# WIDTHS
#########

# We are doing this in R code to provide generics like widthDetails
# so that users can customise the behaviour for complex grobs by
# writing their own (R code!) methods
width <- function(x) {
  widthDetails(x)
}

widthDetails <- function(x) {
  UseMethod("widthDetails", x)
}

widthDetails.default <- function(x) {
  unit(1, "null")
}

#########
# HEIGHTS
#########
height <- function(x) {
  heightDetails(x)
}

heightDetails <- function(x) {
  UseMethod("heightDetails", x)
}

heightDetails.default <- function(x) {
  unit(1, "null")
}

ascentDetails <- function(x) {
    UseMethod("ascentDetails", x)
}

ascentDetails.default <- heightDetails.default

ascentDetails.grob <- function(x) {
    heightDetails(x)
}

descentDetails <- function(x) {
    UseMethod("descentDetails", x)
}

descentDetails.default <- function(x) {
    unit(0, "inches")
}

#########
# Some functions that might be useful for determining the sizes
# of your grobs
#########

# Dimensions which depend on the parent context EITHER don't make
# sense (e.g., no good to have the parent width depend on the child's
# width unit(1, "grobwidth", <child>), which depends on the parent's
# width unit(.1, "npc"), ...) OR are slightly ambiguous
# (e.g., gf <- grid.frame(); grid.pack(gf, grid.rect(width=unit(.1, "npc")))
# makes the area allocated to the rectangle .1 of the frame area, but
# then the rectangle only occupies .1 of _that_ allocated area;  my head
# hurts !).  The first sort will actually lead to infinite loops so
# watch out for that;  the second sort I just don't want to have to deal with.
#
# On the other hand, dimensions which do not depend on the parent context
# are much easier to deal with (e.g., "inches", "cm", "lines", ...)
#
# So this function takes a unit and returns absolute values
# untouched and replaces other values with unit(1, "null")

absolute.size <- function(unit) {
  absolute.units(unit)
}


################################################################################
## Render typeset glyphs

validDetails.glyphgrob <- function(x) {
    if (!inherits(x$glyphInfo, "RGlyphInfo"))
        stop("Invalid glyph info")
    if (!is.unit(x$x) || !is.unit(x$y))
        stop("'x' and 'y' must be units")
    ## Make sure that x and y are of length > 0
    if (length(x$x) < 1 || length(x$y) < 1)
        stop("'x' and 'y' must have length > 0")
    if (!inherits(x$hjust, "GlyphJust") || !inherits(x$vjust, "GlyphJust"))
        stop("'hjust' and 'vjust' must be glyphJust() values")
    if (length(x$hjust) != 1 || length(x$vjust) != 1)
        stop("'hjust' and 'vjust' must have length 1")
    x
}

glyphHOffset <- function(glyphInfo, hjust) {
    gx <- convertWidth(unit(glyphInfo$glyphs$x, "bigpts"), "in", valueOnly=TRUE)
    glyphWidth <- glyphInfo$width
    width <- convertWidth(unit(glyphWidth, "bigpts"), "in",
                          valueOnly=TRUE)
    names(width) <- names(glyphWidth)
    hAnchor <- glyphInfo$hAnchor
    if (is.numeric(hjust)) {
        justName <- names(hjust)
        if (is.null(justName)) {
            gx - convertWidth(unit(hAnchor["left"], "bigpts"), "in",
                              valueOnly=TRUE) -
                hjust*width[1]
        } else {
            if (!justName %in% names(width)) {
                warning("Unknown width; using first width")
                gx - convertWidth(unit(hAnchor["left"], "bigpts"), "in",
                                  valueOnly=TRUE) -
                    hjust*width[1]
            } else {
                anchor <- glyphWidthLeft(glyphWidth, justName)
                gx - convertWidth(unit(hAnchor[anchor], "bigpts"), "in",
                              valueOnly=TRUE) -
                    hjust*width[justName]
            }
        }
    } else {
        if (!hjust %in% names(hAnchor)) {
            warning("Unknown anchor; using left justification")
            gx - convertWidth(unit(hAnchor["left"], "bigpts"), "in",
                              valueOnly=TRUE)
        } else {
            gx - convertWidth(unit(hAnchor[hjust], "bigpts"), "in",
                              valueOnly=TRUE)
        }
    }
}

glyphHJust <- function(x, glyphInfo, hjust) {
    x <- convertX(x, "in", valueOnly=TRUE)
    x + glyphHOffset(glyphInfo, hjust)
}

glyphVOffset <- function(glyphInfo, vjust) {
    gy <- convertHeight(unit(glyphInfo$glyphs$y, "bigpts"), "in",
                        valueOnly=TRUE)
    glyphHeight <- glyphInfo$height
    height <- convertHeight(unit(glyphHeight, "bigpts"), "in",
                            valueOnly=TRUE)
    names(height) <- names(glyphHeight)
    vAnchor <- glyphInfo$vAnchor
    if (is.numeric(vjust)) {
        justName <- names(vjust)
        if (is.null(justName)) {
            gy - convertHeight(unit(vAnchor["bottom"], "bigpts"), "in",
                               valueOnly=TRUE) -
                vjust*height[1]
        } else {
            if (!justName %in% names(height)) {
                warning("Unknown height; using first height")
                gy - convertHeight(unit(vAnchor["bottom"], "bigpts"), "in",
                                   valueOnly=TRUE) -
                    vjust*height[1]
            } else {
                anchor <- glyphHeightBottom(glyphHeight, justName)
                gy - convertHeight(unit(vAnchor[anchor], "bigpts"), "in",
                                   valueOnly=TRUE) -
                    vjust*height[justName]
            }
        }
    } else {
        if (!vjust %in% names(vAnchor)) {
            warning("Unknown anchor; using bottom justification")
            gy - convertHeight(unit(vAnchor["bottom"], "bigpts"), "in",
                               valueOnly=TRUE)
        } else {
            gy - convertHeight(unit(vAnchor[vjust], "bigpts"), "in",
                               valueOnly=TRUE)
        }
    }
}

glyphVJust <- function(y, glyphInfo, vjust) {
    y <- convertY(y, "in", valueOnly=TRUE)
    y + glyphVOffset(glyphInfo, vjust)
}

drawDetails.glyphgrob <- function(x, recording=TRUE) {
    ## Calculate runs of glyphs
    fontstring <- unlist(do.call(paste,
                                 c(x$glyphInfo$glyphs[c("font", "size",
                                                        "colour")],
                                   list(sep=":"))))
    runs <- rle(fontstring)
    ## Calculate final glyph positions
    gx <- unit(glyphHJust(x$x, x$glyphInfo, x$hjust), "in")
    gy <- unit(glyphVJust(x$y, x$glyphInfo, x$vjust), "in")
    ## Replace NA colours with current gp$col
    naCol <- is.na(x$glyphInfo$glyphs$colour)
    if (any(naCol))
        x$glyphInfo$glyphs$colour[naCol] <- get.gpar("col")$col[1]
    ## Call dev->glyph() for each run of glyphs
    grid.Call.graphics(C_glyph,
                       as.integer(runs$lengths),
                       x$glyphInfo, gx, gy)
}

glyphRect <- function(x) {
    gx <- glyphHJust(x$x, x$glyphInfo, x$hjust)
    gy <- glyphVJust(x$y, x$glyphInfo, x$vjust)
    w <- convertWidth(unit(x$glyphInfo$width[1], "bigpts"),
                      "in", valueOnly=TRUE)
    h <- convertHeight(unit(x$glyphInfo$height[1], "bigpts"),
                       "in", valueOnly=TRUE)
    left <- min(gx)
    bottom <- min(gy)
    rectGrob(left, bottom, w, h, default.units="in", just=c("left", "bottom"))
}

xDetails.glyphgrob <- function(x, theta) {
    xDetails(glyphRect(x), theta)
}

yDetails.glyphgrob <- function(x, theta) {
    yDetails(glyphRect(x), theta)
}

widthDetails.glyphgrob <- function(x) {
    widthDetails(glyphRect(x))
}

heightDetails.glyphgrob <- function(x) {
    heightDetails(glyphRect(x))
}

xDetails.glyphgrob <- function(x, theta) {
    gx <- glyphHJust(x$x, x$glyphInfo, x$hjust)
    gy <- glyphVJust(x$y, x$glyphInfo, x$vjust)
    w <- convertWidth(unit(x$glyphInfo$width[1], "bigpts"),
                      "in", valueOnly=TRUE)
    h <- convertHeight(unit(x$glyphInfo$height[1], "bigpts"),
                       "in", valueOnly=TRUE)
    left <- min(gx)
    bottom <- min(gy)
    xDetails(rectGrob(left, bottom, w, h, default.units="in",
                      just=c("left", "bottom")),
             theta)
}

grobPoints.glyphgrob <- function(x, closed=TRUE, ...) {
    if (closed) {
        gx <- glyphHJust(x$x, x$glyphInfo, x$hjust)
        gy <- glyphVJust(x$y, x$glyphInfo, x$vjust)
        w <- convertWidth(unit(x$glyphInfo$width[1], "bigpts"),
                          "in", valueOnly=TRUE)
        h <- convertHeight(unit(x$glyphInfo$height[1], "bigpts"),
                           "in", valueOnly=TRUE)
        left <- min(gx)
        bottom <- min(gy)
        right <- left + w
        top <- bottom + h
        gridGrobCoords(list("1"=gridCoords(x=c(left, left, right, right),
                                           y=c(bottom, top, top, bottom))),
                       x$name)
    } else {
        emptyGrobCoords(x$name)
    }    
}

glyphGrob <- function(glyphInfo,
                      x=.5, y=.5, default.units="npc",
                      hjust="centre", vjust="centre",
                      gp=gpar(), vp=NULL, name=NULL) {
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    grob(glyphInfo=glyphInfo, x=x, y=y,
         hjust=glyphJust(hjust), vjust=glyphJust(vjust),
         gp=gp, vp=vp, name=name,
         cl="glyphgrob")    
}

grid.glyph <- function(...) {
    grid.draw(glyphGrob(...))
}
#  File src/library/grid/R/unit.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2016 The R Core Team
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


# Create an object of class "unit"
# Simple units are of the form 'unit(1, "cm")' or 'unit(1L:3, "cm")' or
# 'unit(c(1,3,6), c("cm", "inch", "npc"))'
# More complicated units are of the form 'unit(1, "string", "a string")'
# or 'unit(1, "grob", a.grob)'
unit <- function(x, units, data = NULL) {
    # Old units passed to a grob will get passed on to unit() and need to be 
    # upgraded instead of converted to the default unit
    if (is.unit(x)) return(upgradeUnit(x))

    x <- as.numeric(x)
    units <- as.character(units)
    if (length(x) == 0 || length(units) == 0)
        stop("'x' and 'units' must have length > 0")
    if (is.null(data)) {
        data <- list(NULL)
    } else if (is.language(data)) {
        data <- list(as.expression(data))
    } else if (is.character(data) || is.grob(data) || inherits(data, "gPath")) {
    	data <- list(data)
    }
    .Call(C_constructUnits, x, data, units)
}

single_unit <- function(x, data, valid_units) {
    `class<-`(list(
        list(
            x,
            data,
            valid_units
        )
    ), c("unit", "unit_v2"))
}
grid.convert <- function(x, unitTo, axisFrom="x", typeFrom="location",
                         axisTo=axisFrom, typeTo=typeFrom,
                         valueOnly=FALSE) {
    .Defunct("convertUnit")
}

convertUnit <- function(x, unitTo, axisFrom="x", typeFrom="location",
                        axisTo=axisFrom, typeTo=typeFrom,
                        valueOnly=FALSE) {
  whatfrom <- match(axisFrom, c("x", "y")) - 1L +
    2L*(match(typeFrom, c("location", "dimension")) - 1L)
  whatto <- match(axisTo, c("x", "y")) - 1L +
    2L*(match(typeTo, c("location", "dimension")) - 1L)
  x <- upgradeUnit(x)
  if (is.na(whatfrom) || is.na(whatto))
    stop("invalid 'axis' or 'type'")
  value <- grid.Call(C_convert, x, as.integer(whatfrom),
                 as.integer(whatto), valid.units(unitTo))
  if (!valueOnly)
    unit(value, unitTo)
  else
    value
}

grid.convertX <- function(x, unitTo, valueOnly=FALSE) {
  .Defunct("convertX")
}

convertX <- function(x, unitTo, valueOnly=FALSE) {
  convertUnit(x, unitTo, "x", "location", "x", "location",
              valueOnly=valueOnly)
}

grid.convertY <- function(x, unitTo, valueOnly=FALSE) {
  .Defunct("convertY")
}

convertY <- function(x, unitTo, valueOnly=FALSE) {
  convertUnit(x, unitTo, "y", "location", "y", "location",
              valueOnly=valueOnly)
}

grid.convertWidth <- function(x, unitTo, valueOnly=FALSE) {
  .Defunct("convertWidth")
}

convertWidth <- function(x, unitTo, valueOnly=FALSE) {
  convertUnit(x, unitTo, "x", "dimension", "x", "dimension",
              valueOnly=valueOnly)
}

grid.convertHeight <- function(x, unitTo, valueOnly=FALSE) {
  .Defunct("convertHeight")
}

convertHeight <- function(x, unitTo, valueOnly=FALSE) {
  convertUnit(x, unitTo, "y", "dimension", "y", "dimension",
              valueOnly=valueOnly)
}

convertNative <- function(unit, dimension="x", type="location") {
  .Defunct("convertUnit")
}

deviceLoc <- function(x, y, valueOnly=FALSE, device=FALSE) {
    if (!(is.unit(x) && length(x) > 0 &&
          is.unit(y) && length(y) > 0 &&
          is.logical(device) && length(device) == 1 && !is.na(device)))
        stop("'x' and 'y' must be units and 'device' must be logical")
    result <- grid.Call(C_devLoc, x, y, as.logical(device))
    names(result) <- c("x", "y")
    if (valueOnly) {
        result
    } else {
        if (device) {
            list(x=unit(result$x, "native"), y=unit(result$y, "native"))
        } else {
            list(x=unit(result$x, "in"), y=unit(result$y, "in"))
        }
    }
}

deviceDim <- function(w, h, valueOnly=FALSE, device=FALSE) {
    if (!(is.unit(w) && length(w) > 0 &&
          is.unit(h) && length(h) > 0 &&
          is.logical(device) && length(device) == 1 && !is.na(device)))
        stop("'w' and 'h' must be units and 'device' must be logical")
    result <- grid.Call(C_devDim, w, h, as.logical(device))
    names(result) <- c("w", "h")
    if (valueOnly) {
        result
    } else {
        if (device) {
            list(w=unit(result$w, "native"), h=unit(result$h, "native"))
        } else {
            list(w=unit(result$w, "in"), h=unit(result$h, "in"))
        }
    }
}

# This is like the "convert" functions:  it evaluates units (immediately)
# in the current context
calcStringMetric <- function(text) {
    # .Call rather than .Call.graphics because it is a one-off calculation
    metric <- grid.Call(C_stringMetric, text)
    names(metric) <- c("ascent", "descent", "width")
    metric
}

# NOTE: the order of the strings in these conversion functions must
# match the order of the enums in ../src/grid.h
# AND in ../src/unit.c (see UnitTable)
# NOTE: ../src/unit.c also allows some pseudonyms (e.g., "in" for "inches")
.grid.unit.list <- c("npc", "cm", "inches", "lines",
                     "native", "null", "snpc", "mm",
                     "points", "picas", "bigpts",
                     "dida", "cicero", "scaledpts",
                     "strwidth", "strheight",
                     "strascent", "strdescent",
                     "vplayoutwidth", "vplayoutheight", "char",
                     "grobx", "groby", "grobwidth", "grobheight",
                     "grobascent", "grobdescent",
                     "mylines", "mychar", "mystrwidth", "mystrheight",
                     "sum", "min", "max")

valid.units <- function(units) {
  .Call(C_validUnits, units)
}

# Printing, formating, and coercion
unitDesc <- function(x, format = FALSE, ...) {
    amount <- if (format) format(x[[1]], ...) else x[[1]]
    unit <- units[as.character(x[[3]])]
    if (unit %in% c('sum', 'min', 'max')) {
        paste0(if (amount == 1) '' else paste0(amount, '*'),
               unit, '(',
               paste(lapply(unclass(x[[2]]), unitDesc, format = format, ...),
                     collapse = ', '),
               ')')
    } else {
        paste0(amount, unit)
    }
}

unitType <- function(x, recurse=FALSE) {
    x <- upgradeUnit(x)
    if (is.simpleUnit(x)) {
        names <- rep_len(units[[as.character(attr(x, "unit"))]], length(x))
        if (recurse) {
            unit <- as.list(names)
            names(unit) <- names
            unit
        } else {
            names
        }
    } else {
        unit <- lapply(unclass(x), `[[`, 3)
        names <- unlist(units[as.character(unit)], use.names=FALSE)
        if (recurse) {
            sub <- names %in% c("sum", "min", "max")
            if (any(sub)) {
                unit[sub] <- lapply(unclass(x)[sub],
                                    function(u) unitType(u[[2]], recurse))
            }
            if (any(!sub)) {
                unit[!sub] <- names[!sub]
            }
            names(unit) <- names
            unit
        } else {
            names
        }
    }
}

as.character.unit <- function(x, ...) {
    x <- upgradeUnit(x) # guard against old unit
    vapply(unclass(as.unit(x)), unitDesc, character(1))
}
as.double.unit <- function(x, ...) {
    x <- upgradeUnit(x) # guard against old unit
    vapply(unclass(x), `[[`, numeric(1), 1L)
}
as.vector.unit <- as.double.unit
format.unit <- function(x, ...) {
    x <- upgradeUnit(x) # guard against old unit
    vapply(unclass(as.unit(x)), unitDesc, character(1), format = TRUE, ...)
}
print.unit <- function(x, ...) {
    print(as.character(x), quote = FALSE, ...)
    invisible(x)
}
as.double.simpleUnit <- function(x, ...) as.double(unclass(x), ...)
as.vector.simpleUnit <- function(x, ...) as.double(unclass(x), ...)

upgradeUnit <- function(x) {
    if (is.newUnit(x)) return(x)
    UseMethod("upgradeUnit")
}
upgradeUnit.unit <- function(x) {
    unit(unclass(x), attr(x, "unit"), attr(x, 'data'))
}
upgradeUnit.unit.list <- function(x) {
    do.call(unit.c, lapply(unclass(x), upgradeUnit))
}
upgradeUnit.unit.arithmetic <- function(x) {
    fun <- .subset2(x, "fname")
    arg1 <- .subset2(x, "arg1")
    if (inherits(arg1, "unit")) arg1 <- upgradeUnit(arg1)
    arg2 <- .subset2(x, "arg2")
    if (inherits(arg2, "unit")) arg2 <- upgradeUnit(arg2)
    
    do.call(fun, list(arg1, arg2))
}
upgradeUnit.default <- function(x) {
    stop("Not a unit object")
}
is.unit <- function(x) {
    inherits(x, 'unit')
}
is.newUnit <- function(x) {
    inherits(x, 'unit_v2')
}
is.simpleUnit <- function(x) inherits(x, 'simpleUnit')
identicalUnits <- function(x) .Call(C_conformingUnits, x)

as.unit <- function(x) {
	.Call(C_asUnit, x)
}

str.unit <- function(object, ...) {
    object <- upgradeUnit(object)
    object <- unclass(as.unit(object))
    for (i in seq_along(object)) {
      unit <- object[[i]]
      cat('[[', i, ']] Amount: ', unit[[1]], '; Unit: ', 
          units[[as.character(unit[[3]])]], '; Data: ', 
          if (is.null(unit[[2]])) 
              'none' 
          else 
              paste(as.character(unit[[2]]), collapse = ", "), '\n', sep = '')
    }
}
#########################
# UNIT ARITHMETIC STUFF
#########################

Summary.unit <- function(..., na.rm=FALSE) {
    units <- list(...)
    units <- units[!vapply(units, is.null, logical(1))]
    ok <- switch(.Generic, "sum" = 201L, "min" = 202L, "max" = 203L, 0L)
    if (ok == 0)
        stop(gettextf("'Summary' function '%s' not meaningful for units",
                      .Generic), domain = NA)
    # Optimise for simple units
    identicalSimple <- identicalUnits(units)
    if (!is.null(identicalSimple)) {
        res <- switch(
            .Generic,
            "sum" = sum(unlist(units)),
            "min" = min(unlist(units)),
            "max" = max(unlist(units)),
        )
        return(`attributes<-`(res, list(
            class = c("simpleUnit", "unit", "unit_v2"), 
            unit = identicalSimple
        )))
    }
    # NOTE that this call to unit.c makes sure that arg1 is
    # a single unit object
    x <- unlist(lapply(units, as.unit), recursive = FALSE)
    class(x) <- c('unit', 'unit_v2')
    matchUnits <- .Call(C_matchUnit, x, ok)
    x <- unclass(x)
    nMatches <- length(matchUnits)
    
    if (nMatches != 0) {
        data <- lapply(x, .subset2, 2L)
        amount <- vapply(x, .subset2, numeric(1), 1L)[matchUnits]
        amount <- rep(amount, lengths(data[matchUnits]))
        matchData <- unclass(unlist(data[matchUnits], recursive = FALSE))
        for (i in seq_along(amount)) {
            if (amount[i] != 1) 
                matchData[[i]][[1]] <- matchData[[i]][[1]] * amount[i]
        }
        if (nMatches == length(x)) {
            data <- matchData
        } else {
            data <- c(x[-matchUnits], matchData)
        }
    } else {
        data <- x
    }
    single_unit(1, `class<-`(data, c('unit', 'unit_v2')), ok)
}
Ops.unit <- function(e1, e2) {
    ok <- switch(.Generic, "+"=TRUE, "-"=TRUE, "*"=TRUE, "/"=TRUE, FALSE)
    if (!ok)
        stop(gettextf("operator '%s' not meaningful for units", .Generic),
             domain = NA)
    # Unary
    if (missing(e2)) {
        if (.Generic %in% c('*', '/')) 
            stop("'*' or '/' cannot be used as a unary operator")
        if (.Generic == '-') {
            if (is.simpleUnit(e1)) {
                attr <- attributes(e1)
                e1 <- -as.vector(e1)
                attributes(e1) <- attr
            } else {
                e1 <- upgradeUnit(e1)  # guard against old unit
                e1 <- .Call(C_flipUnits, e1)
            }
        }
        return(e1)
    }
    # Multiply
    if (.Generic %in% c("*", "/")) {
        # can only multiply a unit by a scalar
        if (nzchar(.Method[1L])) {
            if (nzchar(.Method[2L])) stop("only one operand may be a unit")
            if (!is.numeric(e2)) stop("non-unit operand must be numeric")
            unit <- e1
            value <- e2
        } else {
            if (!is.numeric(e1)) stop("non-unit operand must be numeric")
            if (.Generic == "/") stop("can't divide with a unit")
            unit <- e2
            value <- e1
        }
        if (.Generic == "/") value <- 1 / value
        if (is.simpleUnit(unit)) {
            attr <- attributes(unit)
            unit <- value * as.vector(unit)
            attributes(unit) <- attr
        } else {
            unit <- upgradeUnit(unit)  # guard against old unit
            unit <- .Call(C_multUnits, unit, value)
        }
        return(unit)
    }
    # Add and sub remains
    if (!nzchar(.Method[1L]) && !nzchar(.Method[2L])) {
        stop("both operands must be units")
    }
    if ((is.simpleUnit(e1) && is.simpleUnit(e2)) && (attr(e1, 'unit') == attr(e2, 'unit'))) {
        attr <- attributes(e1)
        unit <- switch(
            .Generic, 
            "-" = as.vector(e1) - as.vector(e2), 
            "+" = as.vector(e1) + as.vector(e2)
        )
        return(`attributes<-`(unit, attr))
    }
    # Convert subtraction to addition
    if (.Generic == '-') {
        e2 <- -e2
    }
    .Call(C_addUnits, as.unit(e1), as.unit(e2))
}

unit.pmin <- function(...) {
    pSummary(..., op = 'min')
}

unit.pmax <- function(...) {
    pSummary(..., op = 'max')
}

unit.psum <- function(...) {
    pSummary(..., op = 'sum')
}

pSummary <- function(..., op) {
    units <- list(...)
    units <- units[lengths(units) != 0]
    # optimisation for simple units
    identicalSimple <- identicalUnits(units)
    if (!is.null(identicalSimple)) {
        res <- switch(
            op,
            "sum" = Reduce(`+`, lapply(units, unclass)),
            "min" = do.call(pmin, lapply(units, unclass)),
            "max" = do.call(pmax, lapply(units, unclass))
        )
        return(`attributes<-`(res, list(
            class = c("simpleUnit", "unit", "unit_v2"), 
            unit = identicalSimple
        )))
    }
    op <- switch(op, sum = 201L, min = 202L, max = 203L)
    .Call(C_summaryUnits, units, op)
}

#########################
# Unit subsetting
#########################

## 'top' argument retained to avoid breaking any uses from it
## that are hang-overs from old unit implementation

`[.unit` <- function(x, index, ..., top = TRUE) {
    x <- upgradeUnit(x) # guard against old unit
    attr <- attributes(x)
    x <- unclass(x)
    n <- length(x)
    if (is.numeric(index) && any(index > n)) {
        if (top) stop('index out of bounds ("unit" subsetting)', call. = FALSE)
        index <- (seq_len(n)[index] - 1L) %% n + 1L
    }
    x <- x[index]
    if (length(x) == 0) {
        stop('Cannot create zero-length unit vector ("unit" subsetting)', call. = FALSE)
    }
    if (!is.null(attr$names)) attr$names <- attr$names[index]
    `attributes<-`(x, attr)
}
`[[.unit` <- function(x, index, ...) {
    if (length(index) != 1) {
        stop("index must be of length 1", call. = FALSE)
    }
    x[index]
}
`[<-.unit` <- function(x, i, value) {
    x <- upgradeUnit(x) # guard against old unit
    attr <- attributes(x)
    simpleResult <- FALSE
    if (is.simpleUnit(x)) {
        if (!(is.simpleUnit(value) && attr(x, 'unit') == attr(value, 'unit'))) {
            x <- as.unit(x)
            value <- as.unit(value)
        } else {
            simpleResult <- TRUE
        }
    } else {
        value <- as.unit(value)
    }
    x <- unclass(x)
    x[i] <- value
    if (simpleResult) {
        attributes(x) <- attr
    } else {
        class(x) <- c("unit", "unit_v2")
    }
    x
}
`[[<-.unit` <- function(x, i, value) {
    if (length(i) != 1) {
        stop("index must be of length 1", call. = FALSE)
    }
    if (length(value) != 1) {
        stop("replacement must be of length 1", call. = FALSE)
    }
    x[i] <- value
    x
}
#########################
# "c"ombining unit objects
#########################

# NOTE that I have not written methods for c()
# because method dispatch occurs on the first argument to
# "c" so c(unit(...), ...) would come here, but c(whatever, unit(...), ...)
# would go who-knows-where.
# A particularly nasty example is:  c(1, unit(1, "npc")) which will
# produce the same result as c(1, 1)
# Same problem for trying to control c(<unit>, <unit.arithmetic>)
# versus c(<unit.arithmetic>, <unit>), etc ...

# If any arguments are unit.arithmetic or unit.list, then the result will be
# unit.list

unit.c <- function(..., check = TRUE) {
    x <- list(...)
    identicalSimple <- identicalUnits(x)
    if (!is.null(identicalSimple)) {
        `attributes<-`(unlist(x), list(class = c('simpleUnit', 'unit', 'unit_v2'), unit = identicalSimple))
    } else {
        `class<-`(unlist(lapply(x, as.unit), recursive = FALSE), c('unit', 'unit_v2'))
    }
}

#########################
# rep'ing unit objects
#########################

rep.unit <- function(x, times = 1, length.out = NA, each = 1, ...) {
    index <- rep(seq_along(x), times = times, length.out = length.out, each = each)
    x[index]
}

# Vestige from when rep() was not generic
unit.rep <- function (x, ...)
{
  warning("'unit.rep' has been deprecated in favour of a unit method for the generic rep function", domain = NA)
  rep(x, ...)
}

#########################
# Length of unit objects
#########################

# Vestige of when length was not generic and a custom length method was needed
unit.length <- function(unit) {
   warning("'unit.length' has been deprecated in favour of a unit method for the generic length function", domain = NA)
   length(unit)
}

#########################
# Convenience functions
#########################

stringWidth <- function(string) {
    n <- length(string)
    if (is.language(string)) {
        string <- as.expression(string)
        data <- vector("list", n)
        for (i in 1L:n)
            data[[i]] <- string[i]
    } else {
        data <- as.list(as.character(string))
    }
    unit(rep_len(1, n), "strwidth", data=data)
}

stringHeight <- function(string) {
    n <- length(string)
    if (is.language(string)) {
        string <- as.expression(string)
        data <- vector("list", n)
        for (i in 1L:n)
            data[[i]] <- string[i]
    } else {
        data <- as.list(as.character(string))
    }
    unit(rep_len(1, n), "strheight", data=data)
}

stringAscent <- function(string) {
    n <- length(string)
    if (is.language(string)) {
        string <- as.expression(string)
        data <- vector("list", n)
        for (i in 1L:n)
            data[[i]] <- string[i]
    } else {
        data <- as.list(as.character(string))
    }
    unit(rep_len(1, n), "strascent", data=data)
}

stringDescent <- function(string) {
    n <- length(string)
    if (is.language(string)) {
        string <- as.expression(string)
        data <- vector("list", n)
        for (i in 1L:n)
            data[[i]] <- string[i]
    } else {
        data <- as.list(as.character(string))
    }
    unit(rep_len(1, n), "strdescent", data=data)
}

convertTheta <- function(theta) {
    if (is.character(theta))
        # Allow some aliases for common angles
        switch(theta,
               east=0,
               north=90,
               west=180,
               south=270,
               stop("invalid 'theta'"))
    else
        # Ensure theta in [0, 360)
        theta <- as.numeric(theta) %% 360
}

# grobX
grobX <- function(x, theta) {
    UseMethod("grobX", x)
}

grobX.grob <- function(x, theta) {
  unit(convertTheta(theta), "grobx", data=x)
}

grobX.gList <- function(x, theta) {
  unit(rep(convertTheta(theta), length(x)), "grobx", data=x)
}

grobX.gPath <- function(x, theta) {
  unit(convertTheta(theta), "grobx", data=x)
}

grobX.default <- function(x, theta) {
  unit(convertTheta(theta), "grobx", data=gPath(as.character(x)))
}

# grobY
grobY <- function(x, theta) {
    UseMethod("grobY", x)
}

grobY.grob <- function(x, theta) {
  unit(convertTheta(theta), "groby", data=x)
}

grobY.gList <- function(x, theta) {
  unit(rep(convertTheta(theta), length(x)), "groby", data=x)
}

grobY.gPath <- function(x, theta) {
  unit(convertTheta(theta), "groby", data=x)
}

grobY.default <- function(x, theta) {
  unit(convertTheta(theta), "groby", data=gPath(as.character(x)))
}

# grobWidth
grobWidth <- function(x) {
  UseMethod("grobWidth")
}

grobWidth.grob <- function(x) {
  unit(1, "grobwidth", data=x)
}

grobWidth.gList <- function(x) {
  unit(rep_len(1, length(x)), "grobwidth", data=x)
}

grobWidth.gPath <- function(x) {
  unit(1, "grobwidth", data=x)
}

grobWidth.default <- function(x) {
  unit(1, "grobwidth", data=gPath(as.character(x)))
}

# grobHeight
grobHeight <- function(x) {
  UseMethod("grobHeight")
}

grobHeight.grob <- function(x) {
  unit(1, "grobheight", data=x)
}

grobHeight.gList <- function(x) {
  unit(rep_len(1, length(x)), "grobheight", data=x)
}

grobHeight.gPath <- function(x) {
  unit(1, "grobheight", data=x)
}

grobHeight.default <- function(x) {
  unit(1, "grobheight", data=gPath(as.character(x)))
}

# grobAscent
grobAscent <- function(x) {
  UseMethod("grobAscent")
}

grobAscent.grob <- function(x) {
  unit(1, "grobascent", data=x)
}

grobAscent.gList <- function(x) {
  unit(rep_len(1, length(x)), "grobascent", data=x)
}

grobAscent.gPath <- function(x) {
  unit(1, "grobascent", data=x)
}

grobAscent.default <- function(x) {
  unit(1, "grobascent", data=gPath(as.character(x)))
}

# grobDescent
grobDescent <- function(x) {
  UseMethod("grobDescent")
}

grobDescent.grob <- function(x) {
  unit(1, "grobdescent", data=x)
}

grobDescent.gList <- function(x) {
  unit(rep_len(1, length(x)), "grobdescent", data=x)
}

grobDescent.gPath <- function(x) {
  unit(1, "grobdescent", data=x)
}

grobDescent.default <- function(x) {
  unit(1, "grobdescent", data=gPath(as.character(x)))
}

#########################
# Function to decide which values in a unit are "absolute" (do not depend
# on parent's drawing context or size)
#########################

absolute.units <- function(unit) {
    .Call(C_absoluteUnits, unit)
}

# Lookup table for unit ids
# This table MUST correspond to the enumeration in grid.h
units <- list(
    '0' = "npc",
    '1' = "cm",
    '2' = "inches",
    '3' = "lines",
    '4' = "native",
    '5' = "null",
    '6' = "snpc",
    '7' = "mm",
    '8' = "points",
    '9' = "picas",
    '10' = "bigpts",
    '11' = "dida",
    '12' = "cicero",
    '13' = "scaledpts",
    '14' = "strwidth",
    '15' = "strheight",
    '16' = "strascent",
    '17' = "strdescent",
    '18' = "char",
    '19' = "grobx",
    '20' = "groby",
    '21' = "grobwidth",
    '22' = "grobheight",
    '23' = "grobascent",
    '24' = "grobdescent",

    '103' = "mylines",
    '104' = "mychar",
    '105' = "mystrwidth",
    '106' = "mystrheight",

    '201' = "sum",
    '202' = "min",
    '203' = "max",

    '1001' = "centimetre",
    '1001' = "centimetres",
    '1001' = "centimeter",
    '1001' = "centimeters",
    '1002' = "in",
    '1002' = "inch",
    '1003' = "line",
    '1007' = "millimetre",
    '1007' = "millimetres",
    '1007' = "millimeter",
    '1007' = "millimeters",
    '1008' = "point",
    '1008' = "pt"
)

#  File src/library/grid/R/util.R
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


# Define a convenience function that is easy to call from C code
grid.top.level.vp <- function() {
    vp <- viewport(clip=TRUE, mask=FALSE, name="ROOT")
    # Must mimic order of events in push.vp.viewport()
    # In particular, add 'gpar' and 'parentgpar' IN THE RIGHT ORDER
    # before calling pushedvp()
    vp$parentgpar <- gpar()
    vp$gpar <- gpar()
    pushedvp(vp)
}

# An accessor for getting at the grid global state structure
# to make debugging easier for me;  all I have to type is grid:::STATE()
STATE <- function() {
  get(".GRID.STATE", envir=.GridEvalEnv)
}

is.even <- function(x) x %% 2 == 0

is.odd <- function(x) !is.even(x)


grid.pretty <- function(range, n = 5L) {
  if (!is.numeric(range))
      stop("'range' must be numeric")
  .Call(C_pretty2, range, n)
}

#  File src/library/grid/R/viewport.R
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


initvpAutoName <- function() {
  index <- 0
  function() {
    index <<- index + 1
    paste0("GRID.VP.", index)
  }
}

vpAutoName <- initvpAutoName()

vpObject <- function(x, y, width, height, just,
                     gp, clip, mask,
                     xscale, yscale, angle,
                     layout, layout.pos.row, layout.pos.col,
                     name) {
    
    vp <- list(x = x, y = y, width = width, height = height,
               justification = just,
               gp = gp,
               clip = clip,
               xscale = xscale,
               yscale = yscale,
               angle = angle,
               layout = layout,
               layout.pos.row = layout.pos.row,
               layout.pos.col = layout.pos.col,
               valid.just = valid.just(just),
               valid.pos.row = layout.pos.row,
               valid.pos.col = layout.pos.col,
               name = name,
               ## A whole lot of blank slots that pushedvp() fills in
               parentgpar = NULL,
               gpar = NULL,
               trans = NULL,
               widths = NULL,
               heights = NULL,
               width.cm = NULL,
               height.cm = NULL,
               rotation = NULL,
               cliprect = NULL,
               parent = NULL,
               children = NULL,
               devwidth = NULL,
               devheight = NULL,
               clippath = NULL,
               ## Some viewport slots that were added later on
               ## (pairs of 'vp' and 'pushedvp' slots)
               mask = mask,
               resolvedmask = NULL)
    class(vp) <- "viewport"
    vp
}

# NOTE: The order of the elements in viewports and pushedvps are
# VERY IMPORTANT because the C code accesses them using constant
# indices (i.e., if you change the order here the world will end!
valid.viewport <- function(x, y, width, height, just,
                           gp, clip, mask,
                           xscale, yscale, angle,
                           layout, layout.pos.row, layout.pos.col,
                           name) {
  if (length(x) > 1 || length(y) > 1 ||
      length(width) > 1 || length(height) > 1)
    stop("'x', 'y', 'width', and 'height' must all be units of length 1")
  if (!is.gpar(gp))
    stop("invalid 'gp' value")
  if (!is.logical(clip)) {
      if (is.grob(clip)) {
          clip <- createClipPath(as.path(clip))
      } else if (inherits(clip, "GridPath")) {
          clip <- createClipPath(clip)
      } else {
          clip <- switch(as.character(clip),
                         on=TRUE,
                         off=NA,
                         inherit=FALSE,
                         stop("invalid 'clip' value"))
      }
  }
  if (!is.logical(mask)) {
      if (isMask(mask)) {
          ## Do nothing
      } else if (is.grob(mask)) {
          mask <- createMask(mask)
      } else {
          mask <- switch(as.character(mask),
                         inherit=TRUE,
                         none=FALSE,
                         stop("invalid 'mask' value"))
      }
  }
  # Ensure both 'xscale' and 'yscale' are numeric (brute force defense)
  xscale <- as.numeric(xscale)
  yscale <- as.numeric(yscale)
  if (!is.numeric(xscale) || length(xscale) != 2 ||
      any(!is.finite(xscale)) || diff(xscale) == 0)
    stop("invalid 'xscale' in viewport")
  if (!is.numeric(yscale) || length(yscale) != 2 ||
      any(!is.finite(yscale)) || diff(yscale) == 0)
    stop("invalid 'yscale' in viewport")
  if (!is.numeric(angle) || length(angle) != 1 ||
      !is.finite(angle))
    stop("invalid 'angle' in viewport")
  if (!(is.null(layout) || is.layout(layout)))
    stop("invalid 'layout' in viewport")
  if (!is.null(layout.pos.row)) {
    layout.pos.row <- as.integer(range(layout.pos.row))
    if (any(!is.finite(layout.pos.row)))
      stop("invalid 'layout.pos.row' in viewport")
  }
  if (!is.null(layout.pos.col)) {
    layout.pos.col <- as.integer(range(layout.pos.col))
    if (any(!is.finite(layout.pos.col)))
      stop("invalid 'layout.pos.col' in viewport")
  }
  # If name is NULL then we give it a default
  # Otherwise it should be a valid R name
  if (is.null(name))
    name <- vpAutoName()
  # Put all the valid things first so that are found quicker
  ## Order is VERY important
  vpObject(x, y, width, height, valid.just(just),
           gp, clip, mask, xscale, yscale,
           angle, layout, layout.pos.row, layout.pos.col,
           name)
}

# When a viewport is pushed, an internal copy is stored along
# with plenty of additional information relevant to the state
# at the time of being pushed (this is all used to return to this
# viewport without having to repush it)
pushedvp <- function(vp) {
    # NOTE that this function is only called from C code:
    # either directly from L_setviewport() or indirectly from initVP()
    # via grid.top.level.vp()
    # vp$gpar and vp$parentgpar are both set previously in push.vp.viewport()
    pvp <- vp
    ## Children of this pushedvp will be stored
    ## in an environment
    pvp$children = new.env(hash=TRUE, parent=baseenv())
    ## Initial value of 0 means that the viewport will
    ## be pushed "properly" the first time, calculating
    ## transformations, etc ...
    pvp$devwidthcm <- 0
    pvp$devheightcm <- 0
    class(pvp) <- c("pushedvp", class(vp))
    pvp
}

vpFromPushedvp <- function(pvp) {
    ## Unresolve any resolved fills
    if (!is.null(pvp$gp$fill)) {
        pvp$gp$fill <- unresolveFill(pvp$gp$fill)
    }
    ## Unresolve any clip paths or masks
    if (isClipPath(pvp$clip)) {
        pvp$clip <- unresolveClipPath(pvp$clip)
    }
    if (isMask(pvp$mask)) {
        pvp$mask <- unresolveMask(pvp$mask)
    }
    ## Only keep non-pushedvp content
    with(unclass(pvp),
         vpObject(x, y, width, height, justification,
                  gp, clip, mask, xscale, yscale,
                  angle, layout, layout.pos.row, layout.pos.col,
                  name))
}

as.character.viewport <- function(x, ...) {
  paste0("viewport[", x$name, "]")
}

as.character.vpList <- function(x, ...) {
  paste0("(", paste(vapply(x, as.character, ""), collapse=", "), ")")
}

as.character.vpStack <- function(x, ...) {
  paste(vapply(x, as.character, ""), collapse="->")
}

as.character.vpTree <- function(x, ...) {
  paste(x$parent, x$children, sep="->")
}

print.viewport <- function(x, ...) {
  cat(as.character(x), "\n")
  invisible(x)
}

width.details.viewport <- function(x) {
  absolute.size(x$width)
}

height.details.viewport <- function(x) {
  absolute.size(x$height)
}

# How many "levels" in viewport object
depth <- function(x, ...) {
  UseMethod("depth")
}

depth.viewport <- function(x, ...) {
  1
}

depth.vpList <- function(x, ...) {
  # When pushed, the last element of the vpList is pushed last
  # so we are left whereever that leaves us
  depth(x[[length(x)]], ...)
}

depth.vpStack <- function(x, ...) {
  # Elements in the stack may be vpStacks or vpLists or vpTrees
  # so need to sum all the depths
  sum(sapply(x, depth, ..., simplify=TRUE))
}

depth.vpTree <- function(x, ...) {
  # When pushed, the last element of the vpTree$children is
  # pushed last so we are left wherever that leaves us
  depth(x$parent, ...) + depth(x$children[[length(x$children)]], ...)
}

depth.path <- function(x, ...) {
  x$n
}

####################
# Accessors
####################

viewport.layout <- function(vp) {
  vp$layout
}

viewport.transform <- function(vp) {
    .Defunct("current.transform")
}

####################
# Public Constructor
####################
viewport <- function(x = unit(0.5, "npc"),
                     y = unit(0.5, "npc"),
                     width = unit(1, "npc"),
                     height = unit(1, "npc"),
                     default.units = "npc",
                     just = "centre",
                     gp = gpar(),
                     clip = "inherit",
                     mask = "inherit", # or "none" or grob
                     # FIXME: scales are only linear at the moment
                     xscale = c(0, 1),
                     yscale = c(0, 1),
                     angle = 0,
                     # Layout for arranging children of this viewport
                     layout = NULL,
                     # Position of this viewport in parent's layout
                     layout.pos.row = NULL,
                     layout.pos.col = NULL,
                     # This is down here to avoid breaking
                     # existing code
                     name=NULL) {
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    if (!is.unit(width))
        width <- unit(width, default.units)
    if (!is.unit(height))
        height <- unit(height, default.units)
    valid.viewport(x, y, width, height, just,
                   gp, clip, mask, xscale, yscale, angle,
                   layout, layout.pos.row, layout.pos.col, name)
}

is.viewport <- function(vp) {
  inherits(vp, "viewport")
}

#############
# Some classes derived from viewport
#############

viewportorpath <- function(x) {
    is.viewport(x) || inherits(x, "vpPath")
}

vpListFromList <- function(vps) {
  if (all(sapply(vps, viewportorpath, simplify=TRUE))) {
    class(vps) <- c("vpList", "viewport")
    vps
  } else {
    stop("only viewports allowed in 'vpList'")
  }
}

# Viewports will be pushed in parallel
vpList <- function(...) {
  vps <- list(...)
  vpListFromList(vps)
}

# Viewports will be pushed in series
vpStack <- function(...) {
  vps <- list(...)
  if (all(sapply(vps, viewportorpath, simplify=TRUE))) {
    class(vps) <- c("vpStack", "viewport")
    vps
  } else {
    stop("only viewports allowed in 'vpStack'")
  }
}

# Viewports will be pushed as a tree
vpTree <- function(parent, children) {
  if (viewportorpath(parent) && inherits(children, "vpList")) {
    tree <- list(parent=parent, children=children)
    class(tree) <- c("vpTree", "viewport")
    tree
  } else {
    stop("'parent' must be a viewport and 'children' must be a 'vpList' in 'vpTree'")
  }
}

# A function for setting all gpars for vpStack/List/Tree
# Used in size.R
setvpgpar <- function(vp) {
  UseMethod("setvpgpar")
}

setvpgpar.viewport <- function(vp) {
  if (!is.null(vp$gp))
    set.gpar(vp$gp)
}

setvpgpar.vpStack <- function(vp) {
  lapply(vp, setvpgpar)
}

setvpgpar.vpList <- function(vp) {
  setvpgpar(vp[[length(vp)]])
}

setvpgpar.vpTree <- function(vp) {
  setvpgpar(vp$parent)
  setvpgpar(vp$children)
}

#############
# Functions for creating "paths" of viewport names
#############
.grid.pathSep <- "::"

vpPathFromVector <- function(names) {
  if (any(bad <- !is.character(names)))
      stop(ngettext(sum(bad),
                    "invalid viewport name",
                    "invalid viewport names"),
           domain = NA)
  names <- unlist(strsplit(names, .grid.pathSep))
  n <- length(names)
  if (n < 1)
    stop("a viewport path must contain at least one viewport name")
  path <- list(path=if (n==1) NULL else
               paste(names[seq_len(n-1L)], collapse=.grid.pathSep),
               name=names[n],
               n=n)
  class(path) <- c("vpPath", "path")
  path
}

vpPath <- function(...) {
  names <- c(...)
  vpPathFromVector(names)
}

as.character.path <- function(x, ...) {
  if (x$n == 1)
    x$name
  else
    paste(x$path, x$name, sep=.grid.pathSep)
}

print.path <- function(x, ...) {
  cat(as.character(x), "\n")
  invisible(x)
}

`[.vpPath` <- function(x, index, ...) {
  names <- unlist(strsplit(as.character(x), .grid.pathSep))[index]
  vpPathFromVector(names)
}

# Explode path$path
explode <- function(x) {
    UseMethod("explode")
}

explode.character <- function(x) {
    unlist(strsplit(x, .grid.pathSep))
}

explode.path <- function(x) {
  if (x$n == 1)
    x$name
  else
    c(explode(x$path), x$name)
}


#############
# Some handy viewport functions
#############

# Create a viewport with margins given in number of lines
plotViewport <- function(margins=c(5.1, 4.1, 4.1, 2.1), ...) {
  margins <- rep(as.numeric(margins), length.out=4)
  viewport(x=unit(margins[2L], "lines"),
           width=unit(1, "npc") - unit(sum(margins[c(2,4)]), "lines"),
           y=unit(margins[1L], "lines"),
           height=unit(1, "npc") - unit(sum(margins[c(1,3)]), "lines"),
           just=c("left", "bottom"),
           ...)
}

# Create a viewport from data
# If xscale not specified then determine from x
# If yscale not specified then determine from y
dataViewport <- function(xData = NULL, yData = NULL,
                         xscale = NULL, yscale = NULL, extension = 0.05, ...)
{
    extension <- rep(extension, length.out = 2)
    if (is.null(xscale)) {
        if (is.null(xData))
            stop("must specify at least one of 'x' or 'xscale'")
        xscale <- extendrange(xData, f = extension[1L])
    }
    if (is.null(yscale)) {
        if (is.null(yData))
            stop("must specify at least one of 'y' or 'yscale'")
        yscale <- extendrange(yData, f = extension[2L])
    }
    viewport(xscale = xscale, yscale = yscale, ...)
}

editViewport <- function(vp=current.viewport(), ...) {
    edits <- list(...)
    vp <- vpFromPushedvp(vp)
    vp[names(edits)] <- edits
    valid.viewport(vp$x, vp$y, vp$width, vp$height, vp$just,
                   vp$gp, vp$clip, vp$mask, vp$xscale, vp$yscale, vp$angle,
                   vp$layout, vp$layout.pos.row, vp$layout.pos.col, vp$name)
}
#  File src/library/grid/R/zzz.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2016 The R Core Team
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

## environment used for evaluation in the C code
## assigned here to protect from GC, but otherwise unused at R level
.GridEvalEnv <- new.env()

# This should be the only grid global variable(?)
# It contains the list of state structures corresponding to the
# state for each device.
# The state structures are stored in here so that they do not
# get garbage collected.
assign(".GRID.STATE", vector("list", 64L), envir = .GridEvalEnv)
## 64 comes from the maximum number of R devices allowed to be open at
## one time, see R_MaxDevices in Graphics.h.

.noGenerics <- TRUE

utils::globalVariables(c("n", "vp", "path"))

.onLoad <- function(libname, pkgname)
{
    ## want eval in C code to see unexported objects
    environment(.GridEvalEnv) <- asNamespace("grid")
    .Call(C_initGrid, .GridEvalEnv)
    .grid.loaded <<- TRUE
}

.onUnload <- function(libpath)
{
    if (.grid.loaded) {
        ## Kill all existing devices to avoid replay
        ## of display list which tries to run grid code
        ## Not very friendly to other registered graphics systems
        ## but its safety first for now
        if(length(.Devices) > 1L)
            warning("shutting down all devices when unloading 'grid' namespace",
                    call. = FALSE)
        graphics.off()
        .Call(C_killGrid)
    }
    library.dynam.unload("grid", libpath)
}

## .gridplot.hook <- function()
## {
##     pushViewport(viewport(width=unit(1, "npc") - unit(1, "lines"),
## 			  x=0, just="left"))
##     grid.text(paste("help(", ..nameEx, ")"),
## 	      x=unit(1, "npc") + unit(0.5, "lines"),
## 	      y=unit(0.8, "npc"), rot=90,
## 	      gp=gpar(col="orchid"))
## }
