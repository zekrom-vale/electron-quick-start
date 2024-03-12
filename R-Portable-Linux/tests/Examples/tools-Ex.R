pkgname <- "tools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('tools')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CRANtools")
### * CRANtools

flush(stderr()); flush(stdout())

### Name: CRANtools
### Title: CRAN Package Repository Tools
### Aliases: CRAN_package_db CRAN_check_results CRAN_check_details
###   CRAN_check_issues summarize_CRAN_check_status R_CRAN_WEB R_CRAN_SRC

### ** Examples


cleanEx()
nameEx("HTMLheader")
### * HTMLheader

flush(stderr()); flush(stdout())

### Name: HTMLheader
### Title: Generate a Standard HTML Header for R Help
### Aliases: HTMLheader
### Keywords: utilities documentation

### ** Examples

cat(HTMLheader("This is a sample header"), sep="\n")



cleanEx()
nameEx("Rd2HTML")
### * Rd2HTML

flush(stderr()); flush(stdout())

### Name: Rd2HTML
### Title: Rd Converters
### Aliases: Rd2txt Rd2HTML Rd2ex Rd2latex
### Keywords: documentation

### ** Examples
cleanEx()
nameEx("Rd2txt_options")
### * Rd2txt_options

flush(stderr()); flush(stdout())

### Name: Rd2txt_options
### Title: Set Formatting Options for Text Help
### Aliases: Rd2txt_options
### Keywords: documentation

### ** Examples




cleanEx()
nameEx("Rdutils")
### * Rdutils

flush(stderr()); flush(stdout())

### Name: Rdutils
### Title: Rd Utilities
### Aliases: Rd_db
### Keywords: utilities documentation

### ** Examples


cleanEx()
nameEx("assertCondition")
### * assertCondition

flush(stderr()); flush(stdout())

### Name: assertCondition
### Title: Asserting Error Conditions
### Aliases: assertCondition assertWarning assertError
### Keywords: programming error

### ** Examples

  assertError(sqrt("abc"))
  assertWarning(matrix(1:8, 4,3))

  assertCondition( ""-1 ) # ok, any condition would satisfy this

try( assertCondition(sqrt(2), "warning") )
## .. Failed to get warning in evaluating sqrt(2)
     assertCondition(sqrt("abc"), "error")   # ok
try( assertCondition(sqrt("abc"), "warning") )# -> error: had no warning
     assertCondition(sqrt("abc"), "error")
  ## identical to assertError() call above

assertCondition(matrix(1:5, 2,3), "warning")
try( assertCondition(matrix(1:8, 4,3), "error") )
## .. Failed to get expected error ....

## either warning or worse:
assertCondition(matrix(1:8, 4,3), "error","warning") # OK
assertCondition(matrix(1:8, 4, 3), "warning") # OK

## when both are signalled:
ff <- function() { warning("my warning"); stop("my error") }
    assertCondition(ff(), "warning")
## but assertWarning does not allow an error to follow
try(assertWarning(ff()))
    assertCondition(ff(), "error")          # ok
assertCondition(ff(), "error", "warning") # ok (quietly, catching warning)

## assert that assertC..() does not assert [and use *one* argument only]
assertCondition( assertCondition(sqrt( 2   ), "warning") )
assertCondition( assertCondition(sqrt("abc"), "warning"), "error")
assertCondition( assertCondition(matrix(1:8, 4,3), "error"),
                "error")



cleanEx()
nameEx("bibstyle")
### * bibstyle

flush(stderr()); flush(stdout())

### Name: bibstyle
### Title: Select or Define a Bibliography Style
### Aliases: bibstyle getBibstyle
### Keywords: utilties documentation

### ** Examples

## Don't show: 
options(useFancyQuotes = FALSE)
## End(Don't show)
refs <-
c(bibentry(bibtype = "manual",
    title = "R: A Language and Environment for Statistical Computing",
    author = person("R Core Team"),
    organization = "R Foundation for Statistical Computing",
    address = "Vienna, Austria",
    year = 2013,
    url = "https://www.R-project.org"),
  bibentry(bibtype = "article",
    author = c(person(c("George", "E.", "P."), "Box"),
               person(c("David",  "R."),      "Cox")),
    year = 1964,
    title = "An Analysis of Transformations",
    journal = "Journal of the Royal Statistical Society, Series B",
    volume = 26, number = 2, pages = "211--243",
    doi = "10.1111/j.2517-6161.1964.tb00553.x"))

bibstyle("unsorted", sortKeys = function(refs) seq_along(refs),
    fmtPrefix = function(paper) paste0("[", paper$.index, "]"),
       .init = TRUE)
print(refs, .bibstyle = "unsorted")



cleanEx()
nameEx("buildVignettes")
### * buildVignettes

flush(stderr()); flush(stdout())

### Name: buildVignettes
### Title: List and Build Package Vignettes
### Aliases: buildVignettes pkgVignettes
### Keywords: utilities documentation

### ** Examples

gVigns <- pkgVignettes("grid")



cleanEx()
nameEx("charsets")
### * charsets

flush(stderr()); flush(stdout())

### Name: charsets
### Title: Conversion Tables between Character Sets
### Aliases: Adobe_glyphs charset_to_Unicode
### Keywords: datasets

### ** Examples

## find Adobe names for ISOLatin2 chars.
latin2 <- charset_to_Unicode[, "ISOLatin2"]
aUnicode <- as.hexmode(paste0("0x", Adobe_glyphs$unicode))
keep <- aUnicode %in% latin2
aUnicode <- aUnicode[keep]
aAdobe <- Adobe_glyphs[keep, 1]
## first match
aLatin2 <- aAdobe[match(latin2, aUnicode)]
## all matches
bLatin2 <- lapply(1:256, function(x) aAdobe[aUnicode == latin2[x]])
format(bLatin2, justify = "none")



cleanEx()
nameEx("checkFF")
### * checkFF

flush(stderr()); flush(stdout())

### Name: checkFF
### Title: Check Foreign Function Calls
### Aliases: checkFF print.checkFF
### Keywords: programming utilities

### ** Examples


cleanEx()
nameEx("checkPoFiles")
### * checkPoFiles

flush(stderr()); flush(stdout())

### Name: checkPoFiles
### Title: Check Translation Files for Inconsistent Format Strings
### Aliases: checkPoFile checkPoFiles
### Keywords: utilities

### ** Examples

## Not run: 
##D checkPoFiles("de", "/path/to/R/src/directory")
## End(Not run)



cleanEx()
nameEx("checkRdaFiles")
### * checkRdaFiles

flush(stderr()); flush(stdout())

### Name: checkRdaFiles
### Title: Report on Details of Saved Images or Re-saves them
### Aliases: checkRdaFiles resaveRdaFiles
### Keywords: utilities

### ** Examples
## Not run: 
##D ## from a package top-level source directory
##D paths <- sort(Sys.glob(c("data/*.rda", "data/*.RData")))
##D (res <- checkRdaFiles(paths))
##D ## pick out some that may need attention
##D bad <- is.na(res$ASCII) | res$ASCII | (res$size > 1e4 & res$compress == "none")
##D res[bad, ]
## End(Not run)


cleanEx()
nameEx("check_packages_in_dir")
### * check_packages_in_dir

flush(stderr()); flush(stdout())

### Name: check_packages_in_dir
### Title: Check Source Packages and Their Reverse Dependencies
### Aliases: check_packages_in_dir summarize_check_packages_in_dir_depends
###   summarize_check_packages_in_dir_results
###   summarize_check_packages_in_dir_timings check_packages_in_dir_changes
###   check_packages_in_dir_details
### Keywords: utilities

### ** Examples

## Not run: 
##D ## Check packages in dir without reverse dependencies:
##D check_packages_in_dir(dir)
##D ## Check packages in dir and their reverse dependencies using the
##D ## defaults (all repositories in getOption("repos"), all "strong"
##D ## reverse dependencies, no recursive reverse dependencies):
##D check_packages_in_dir(dir, reverse = list())
##D ## Check packages in dir with their reverse dependencies from CRAN,
##D ## using all strong reverse dependencies and reverse suggests:
##D check_packages_in_dir(dir,
##D                       reverse = list(repos = getOption("repos")["CRAN"],
##D                                      which = "most"))                   
##D ## Check packages in dir with their reverse dependencies from CRAN,
##D ## using '--as-cran' for the former but not the latter:
##D check_packages_in_dir(dir,
##D                       check_args = c("--as-cran", ""),
##D                       reverse = list(repos = getOption("repos")["CRAN"]))
## End(Not run)



cleanEx()
nameEx("delimMatch")
### * delimMatch

flush(stderr()); flush(stdout())

### Name: delimMatch
### Title: Delimited Pattern Matching
### Aliases: delimMatch
### Keywords: character

### ** Examples

x <- c("\\value{foo}", "function(bar)")
delimMatch(x)
delimMatch(x, c("(", ")"))



cleanEx()
nameEx("dependsOnPkgs")
### * dependsOnPkgs

flush(stderr()); flush(stdout())

### Name: dependsOnPkgs
### Title: Find Reverse Dependencies
### Aliases: dependsOnPkgs
### Keywords: utilities

### ** Examples


cleanEx()
nameEx("encoded")
### * encoded

flush(stderr()); flush(stdout())

### Name: encoded_text_to_latex
### Title: Translate non-ASCII Text to LaTeX Escapes
### Aliases: encoded_text_to_latex
### Keywords: utilities

### ** Examples

x <- "fa\xE7ile"
encoded_text_to_latex(x, "latin1")
## Not run: 
##D ## create a tex file to show the upper half of 8-bit charsets
##D x <- rawToChar(as.raw(160:255), multiple = TRUE)
##D (x <- matrix(x, ncol = 16, byrow = TRUE))
##D xx <- x
##D xx[] <- encoded_text_to_latex(x, "latin1") # or latin2 or latin9
##D xx <- apply(xx, 1, paste, collapse = "&")
##D con <- file("test-encoding.tex", "w")
##D header <- c(
##D "\\documentclass{article}",
##D "\\usepackage[T1]{fontenc}",
##D "\\usepackage{Rd}",
##D "\\begin{document}",
##D "\\HeaderA{test}{}{test}",
##D "\\begin{Details}\\relax",
##D "\\Tabular{cccccccccccccccc}{")
##D trailer <- c("}", "\\end{Details}", "\\end{document}")
##D writeLines(header, con)
##D writeLines(paste0(xx, "\\\\"), con)
##D writeLines(trailer, con)
##D close(con)
##D ## and some UTF_8 chars
##D x <- intToUtf8(as.integer(
##D     c(160:383,0x0192,0x02C6,0x02C7,0x02CA,0x02D8,
##D       0x02D9, 0x02DD, 0x200C, 0x2018, 0x2019, 0x201C,
##D       0x201D, 0x2020, 0x2022, 0x2026, 0x20AC)),
##D                multiple = TRUE)
##D x <- matrix(x, ncol = 16, byrow = TRUE)
##D xx <- x
##D xx[] <- encoded_text_to_latex(x, "UTF-8")
##D xx <- apply(xx, 1, paste, collapse = "&")
##D con <- file("test-utf8.tex", "w")
##D writeLines(header, con)
##D writeLines(paste(xx, "\\\\", sep = ""), con)
##D writeLines(trailer, con)
##D close(con)
## End(Not run)


cleanEx()
nameEx("fileutils")
### * fileutils

flush(stderr()); flush(stdout())

### Name: fileutils
### Title: File Utilities
### Aliases: file_ext file_path_as_absolute file_path_sans_ext
###   list_files_with_exts list_files_with_type
### Keywords: file

### ** Examples


cleanEx()
nameEx("find_gs_cmd")
### * find_gs_cmd

flush(stderr()); flush(stdout())

### Name: find_gs_cmd
### Title: Find a GhostScript Executable
### Aliases: find_gs_cmd R_GSCMD GSC

### ** Examples
## Not run: 
##D ## Suppose a Solaris system has GhostScript 9.00 on the path and
##D ## 9.07 in /opt/csw/bin.  Then one might set
##D Sys.setenv(R_GSCMD = "/opt/csw/bin/gs")
## End(Not run)


cleanEx()
nameEx("getVignetteInfo")
### * getVignetteInfo

flush(stderr()); flush(stdout())

### Name: getVignetteInfo
### Title: Get Information on Installed Vignettes
### Aliases: getVignetteInfo
### Keywords: utilities documentation

### ** Examples




cleanEx()
nameEx("loadRdMacros")
### * loadRdMacros

flush(stderr()); flush(stdout())

### Name: loadRdMacros
### Title: Load User-defined Rd Help System Macros
### Aliases: loadRdMacros loadPkgRdMacros
### Keywords: utilities documentation

### ** Examples

f <- tempfile()
writeLines(r"(
\newcommand{\Rlogo}{
  \if{html}{\figure{Rlogo.svg}{options: width=100 alt="R logo"}}
  \if{latex}{\figure{Rlogo.pdf}{options: width=0.5in}}
}
)", f)
m <- loadRdMacros(f)
ls(m)
ls(parent.env(parent.env(m)))
parse_Rd(textConnection(r"(\Rlogo)"), fragment = TRUE, macros = m)



cleanEx()
nameEx("makevars")
### * makevars

flush(stderr()); flush(stdout())

### Name: makevars
### Title: User and Site Compilation Variables
### Aliases: makevars_user makevars_site
### Keywords: utilities

### ** Examples
## Don't show: 
checkMV <- function(r)
  stopifnot(is.character(r),
            length(r) == 0 || (length(r) == 1 && file.exists(r)))
checkMV(makevars_user())
checkMV(makevars_site())
## End(Don't show)



cleanEx()
nameEx("md5sum")
### * md5sum

flush(stderr()); flush(stdout())

### Name: md5sum
### Title: Compute MD5 Checksums
### Aliases: md5sum
### Keywords: utilities

### ** Examples

as.vector(md5sum(dir(R.home(), pattern = "^COPY", full.names = TRUE)))



cleanEx()
nameEx("package_dependencies")
### * package_dependencies

flush(stderr()); flush(stdout())

### Name: package_dependencies
### Title: Computations on the Dependency Hierarchy of Packages
### Aliases: package_dependencies
### Keywords: utilities

### ** Examples




cleanEx()
nameEx("package_native_routine_registration_skeleton")
### * package_native_routine_registration_skeleton

flush(stderr()); flush(stdout())

### Name: package_native_routine_registration_skeleton
### Title: Write Skeleton for Adding Native Routine Registration to a
###   Package
### Aliases: package_native_routine_registration_skeleton

### ** Examples
## Not run: 
##D ## with a completed splines/DESCRIPTION file,
##D tools::package_native_routine_registration_skeleton('splines',,,FALSE)
##D ## produces
##D #include <R.h>
##D #include <Rinternals.h>
##D #include <stdlib.h> // for NULL
##D #include <R_ext/Rdynload.h>
##D 
##D /* FIXME: 
##D    Check these declarations against the C/Fortran source code.
##D */
##D 
##D /* .Call calls */
##D extern SEXP spline_basis(SEXP, SEXP, SEXP, SEXP);
##D extern SEXP spline_value(SEXP, SEXP, SEXP, SEXP, SEXP);
##D 
##D static const R_CallMethodDef CallEntries[] = {
##D     {"spline_basis", (DL_FUNC) &spline_basis, 4},
##D     {"spline_value", (DL_FUNC) &spline_value, 5},
##D     {NULL, NULL, 0}
##D };
##D 
##D void R_init_splines(DllInfo *dll)
##D {
##D     R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
##D     R_useDynamicSymbols(dll, FALSE);
##D }
## End(Not run)


cleanEx()
nameEx("parseLatex")
### * parseLatex

flush(stderr()); flush(stdout())

### Name: parseLatex
### Title: Experimental Functions to Work with LaTeX Code
### Aliases: parseLatex deparseLatex latexToUtf8
### Keywords: utilities documentation

### ** Examples


cleanEx()
nameEx("print.via.format")
### * print.via.format

flush(stderr()); flush(stdout())

### Name: .print.via.format
### Title: Printing Utilities
### Aliases: .print.via.format
### Keywords: utilities

### ** Examples

## The function is simply defined as
 function (x, ...) {
    writeLines(format(x, ...))
    invisible(x)
 }

## is used for simple print methods in R, and as prototype for new methods.



cleanEx()
nameEx("pskill")
### * pskill

flush(stderr()); flush(stdout())

### Name: pskill
### Title: Kill a Process
### Aliases: pskill SIGHUP SIGINT SIGQUIT SIGKILL SIGTERM SIGSTOP SIGTSTP
###   SIGCONT SIGCHLD SIGUSR1 SIGUSR2
### Keywords: utilities

### ** Examples
## Not run: 
##D pskill(c(237, 245), SIGKILL)
## End(Not run)


cleanEx()
nameEx("showNonASCII")
### * showNonASCII

flush(stderr()); flush(stdout())

### Name: showNonASCII
### Title: Pick Out Non-ASCII Characters
### Aliases: showNonASCII showNonASCIIfile
### Keywords: utilities

### ** Examples

out <- c(
"fa\xE7ile test of showNonASCII():",
"\\details{",
"   This is a good line",
"   This has an \xfcmlaut in it.",
"   OK again.",
"}")
f <- tempfile()
cat(out, file = f, sep = "\n")

## IGNORE_RDIFF_BEGIN
showNonASCIIfile(f)
## IGNORE_RDIFF_END
unlink(f)



cleanEx()
nameEx("toHTML")
### * toHTML

flush(stderr()); flush(stdout())

### Name: toHTML
### Title: Display an Object in HTML
### Aliases: toHTML toHTML.packageIQR toHTML.news_db
### Keywords: utilities documentation

### ** Examples

cat(toHTML(demo(package = "base")), sep = "\n")



cleanEx()
nameEx("undoc")
### * undoc

flush(stderr()); flush(stdout())

### Name: undoc
### Title: Find Undocumented Objects
### Aliases: undoc print.undoc
### Keywords: documentation

### ** Examples

undoc("tools")                  # Undocumented objects in 'tools'



cleanEx()
nameEx("updatePACKAGES")
### * updatePACKAGES

flush(stderr()); flush(stdout())

### Name: update_PACKAGES
### Title: Update Existing PACKAGES Files
### Aliases: update_PACKAGES
### Keywords: file utilities

### ** Examples

## Not run: 
##D write_PACKAGES("c:/myFolder/myRepository") # on Windows
##D update_PACKAGES("c:/myFolder/myRepository") # on Windows
##D write_PACKAGES("/pub/RWin/bin/windows/contrib/2.9",
##D type = "win.binary") # on Linux
##D update_PACKAGES("/pub/RWin/bin/windows/contrib/2.9",
##D type = "win.binary") # on Linux
## End(Not run)


cleanEx()
nameEx("userdir")
### * userdir

flush(stderr()); flush(stdout())

### Name: userdir
### Title: R User Directories
### Aliases: R_user_dir

### ** Examples

## IGNORE_RDIFF_BEGIN
R_user_dir("FOO", "cache")
## Create one, platform agnostically, must work if <normal> :
(Rdb <- R_user_dir("base"))
if(noD <- !dir.exists(Rdb)) # should work user specifically:
    dir.create(Rdb, recursive=TRUE)
stopifnot(dir.exists(Rdb)) # "everywhere"
dir(Rdb) # typically empty
if(noD) unlink(Rdb) # cleaning up
list.files(R_user_dir("grid"), full.names = TRUE)
## IGNORE_RDIFF_END



cleanEx()
nameEx("vignetteEngine")
### * vignetteEngine

flush(stderr()); flush(stdout())

### Name: vignetteEngine
### Title: Set or Get a Vignette Processing Engine
### Aliases: vignetteEngine
### Keywords: utilities documentation

### ** Examples

str(vignetteEngine("Sweave"))



cleanEx()
nameEx("vignetteInfo")
### * vignetteInfo

flush(stderr()); flush(stdout())

### Name: vignetteInfo
### Title: Basic Information about a Vignette
### Aliases: vignetteInfo
### Keywords: utilities

### ** Examples

  gridEx <- system.file("doc", "grid.Rnw", package = "grid")
  vi <- vignetteInfo(gridEx)
  str(vi)



cleanEx()
nameEx("writePACKAGES")
### * writePACKAGES

flush(stderr()); flush(stdout())

### Name: write_PACKAGES
### Title: Generate PACKAGES Files
### Aliases: write_PACKAGES
### Keywords: file utilities

### ** Examples

## Not run: 
##D write_PACKAGES("c:/myFolder/myRepository")  # on Windows
##D write_PACKAGES("/pub/RWin/bin/windows/contrib/2.9",
##D                type = "win.binary")  # on Linux
## End(Not run)


cleanEx()
nameEx("xgettext")
### * xgettext

flush(stderr()); flush(stdout())

### Name: xgettext
### Title: Extract Translatable Messages from R Files in a Package
### Aliases: xgettext xngettext xgettext2pot
### Keywords: utilities

### ** Examples
## Not run: 
##D ## in a source-directory build (not typical!) of R;
##D ## otherwise, download and unpack the R sources, and replace
##D ## R.home()  by  "<my_path_to_source_R>" :
##D xgettext(file.path(R.home(), "src", "library", "splines"))
## End(Not run)

## Create source package-like  <tmp>/R/foo.R  and get text from it:
tmpPkg <- tempdir()
tmpRDir <- file.path(tmpPkg, "R")
dir.create(tmpRDir, showWarnings = FALSE)
fnChar <- paste(sep = "\n",
  "foo <- function(x) {",
  "  if (x < -1)  stop('too small')",
  "  # messages unduplicated (not so for ngettext)",
  "  if (x < -.5) stop('too small')",
  "  if (x < 0) {",
  "    warning(",
  "      'sqrt(x) is', sqrt(as.complex(x)),",
  "      ', which may be too small'",
  "    )",
  "  }",
  "  # calls with domain=NA are skipped",
  "  if (x == 0) cat(gettext('x is 0!\n', domain=NA))",
  "  # gettext strings may be ignored due to 'outer' domain=NA",
  "  if (x > 10) stop('x is ', gettextf('%.2f', x), domain=NA)",
  "  x",
  "}")

writeLines(fnChar, con = file.path(tmpRDir, "foo.R"))

## [[1]] : suppressing (tmpfile) name to make example Rdiff-able
xgettext(tmpPkg, asCall=TRUE )[[1]] # default; shows  ' sqrt(as.complex(x)) '
xgettext(tmpPkg, asCall=FALSE)[[1]] # doesn't ; but then ' %.2f '

unlink(tmpRDir, recursive=TRUE)



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
