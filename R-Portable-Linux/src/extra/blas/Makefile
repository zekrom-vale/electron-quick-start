#
# ${R_HOME}/src/extra/blas/Makefile


srcdir = .
top_srcdir = ../../..

top_builddir = ../../..
subdir = src/extra/blas
R_HOME = $(top_builddir)

include $(top_builddir)/Makeconf
.SUFFIXES:
.SUFFIXES: .c .f .f90 .d .o
.f90.o:
	$(FC) $(ALL_FCFLAGS) -c  $< -o $@

ALL_CFLAGS = $(ALL_CFLAGS_LO)
ALL_FFLAGS = $(ALL_FFLAGS_LO)
ALL_FCFLAGS = $(ALL_FFLAGS_LO)

SOURCES = blas00.c blas.f cmplxblas.f blas2.f90 cmplxblas2.f90

Rblas_la = libRblas$(R_DYLIB_EXT)
#Rblas_la = libRblas.a
##  is used on macOS
## first for internal BLAS
Rblas_la_LIBADD =  $(FLIBS_IN_SO)
## then external one
Rblas_la_LIBADD0 = 

all: Makefile  R

Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status $(SOURCES_C)
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

R:
	@$(MAKE) $(Rblas_la)
	@$(MAKE) rhome="$(abs_top_builddir)" Rblas_install

blas2.o: blas2.f90
	$(FC) $(ALL_FCFLAGS) -c  $< -o $@
cmplxblas2.o: cmplxblas2.f90
	$(FC) $(ALL_FCFLAGS) -c  $< -o $@

blas_OBJS=blas.o blas2.o  cmplxblas.o cmplxblas2.o
$(Rblas_la): $(blas_OBJS)
	$(DYLIB_LINK) -o $(Rblas_la) $(blas_OBJS) $(Rblas_la_LIBADD)
	#$(AR) -cr $(Rblas_la) $(blas_OBJS)

#$(Rblas_la): blas00.o
#	$(DYLIB_LINK) -o $(Rblas_la) $^  $(Rblas_la_LIBADD0)

Rblas_install: $(Rblas_la)
	@$(MKINSTALLDIRS) "$(DESTDIR)$(Rexeclibdir)"
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(Rblas_la) \
	  "$(DESTDIR)$(Rexeclibdir)/$(Rblas_la)"

Rblas_install-strip: $(Rblas_la)
	@$(MKINSTALLDIRS) "$(DESTDIR)$(Rexeclibdir)"
	@$(INSTALL) $(Rblas_la) "$(DESTDIR)$(Rexeclibdir)/$(Rblas_la)"
	@if test -n "$(STRIP_SHARED_LIB)"; then \
	  $(STRIP_SHARED_LIB) "$(DESTDIR)$(Rexeclibdir)/$(Rblas_la)" ;\
	fi

mostlyclean: clean
clean:
	@-rm -f *.o $(Rblas_la)
distclean: clean
	@-rm -f Makefile
maintainer-clean: distclean

install:  Rblas_install
install-strip:  Rblas_install-strip

uninstall TAGS info dvi check:

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)
DISTFILES = Makefile.in Makefile.win README Rblas.def \
  $(SOURCES) $(HEADERS)

distdir: $(DISTFILES)
	@for f in $(DISTFILES); do \
	  test -f $(distdir)/$${f} \
	    || ln $(srcdir)/$${f} $(distdir)/$${f} 2>/dev/null \
	    || cp -p $(srcdir)/$${f} $(distdir)/$${f}; \
	done

