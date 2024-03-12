#
# ${R_HOME}/src/modules/lapack/Makefile


srcdir = .
top_srcdir = ../../..

top_builddir = ../../..
subdir = src/modules/lapack
R_HOME = $(top_builddir)

include $(top_builddir)/Makeconf
.SUFFIXES:
.SUFFIXES: .c .f .f90 .d .o
.f90.o:
	$(FC) $(ALL_FCFLAGS) -c  $< -o $@

SOURCES_C = Lapack.c # vecLibg95c.c
SOURCES_F = # vecLibg95f.f

## BUILD_LTO_CHECK_TRUE relies on using the internal LAPACK
DEPENDS = $(SOURCES_C:.c=.d)
SOURCES = $(SOURCES_C) $(SOURCES_F) $(SOURCES_F90)
OBJECTS = $(SOURCES_C:.c=.o) $(SOURCES_F:.f=.o) flexiblas.o
HEADERS = Lapack.h

## Sources for the -lRlapack.  On OSX this may be linked directly to
## vecLib/Accelerate, so needs the zdot* fix.
LIBSOURCES = dlamch.f dlapack.f \
   cmplx.f
LIBSOURCES90 = dlartg.f90 dlassq.f90 la_constants.f90 la_xisnan.f90 \
   zlartg.f90 zlassq.f90
LIBOBJECTS = $(LIBSOURCES:.f=.o) $(LIBSOURCES90:.f90=.o) \
  # vecLibg95f.o vecLibg95c.o

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)
DISTFILES = \
	LICENSE README Makefile.in Makefile.win \
	Rlapackrc.rc dllversion.rc \
	$(HEADERS) Lapack.c init_win.c \
	vecLibg95c.c vecLibg95f.f \
	cmplx.f dlapack.f dlamch.f \
	dlartg.f90 dlassq.f90 la_constants.f90 la_xisnan.f90 \
	zlartg.f90 zlassq.f90 \
	dble.txt cmplx.txt R.patch

lapack_la = lapack$(SHLIB_EXT)
lapack_la_OBJECTS = $(OBJECTS)
## Formerly used for passing '+s' to SHLIB_LINK on HP-UX.
lapack_la_LDFLAGS = 
lapack_la_LIBADD = $(LIBR1) # $(LIBINTL)

Rlapack_la = libRlapack$(R_DYLIB_EXT)
#Rlapack_la = libRlapack.a
Rlapack_la_OBJECTS = $(LIBOBJECTS)
## This used to skip FLIBS, but that relied on all the functions
## needed being loaded into R.bin/libR.so, and that need not be the
## case for a static libfortran or if a package is compiled under a
## different compiler.
Rlapack_la_LIBADD =  $(FLIBS_IN_SO) # $(LIBR)

ALL_CFLAGS = $(ALL_CFLAGS_LO)
ALL_FFLAGS = $(ALL_FFLAGS_LO)
ALL_FCFLAGS = $(ALL_FFLAGS_LO)

all: R

flexiblas.o: ../../main/flexiblas.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c  $< -o $@

Makefile: $(srcdir)/Makefile.in \
  $(top_builddir)/config.status \
  $(SOURCES)
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

Makedeps: Makefile $(DEPENDS)
	@cat $(DEPENDS) >> Makefile
	@touch $@

## treated separately: e.g. we need to force stores with 
## gcc >= 3.4.0 at least on i686
## (SHLIB_FFLAGS is apparently currently unused)
LTO = 
DLAMC_FFLAGS=$(FPICFLAGS) $(SHLIB_FFLAGS) -g -O2 -msse2 -mfpmath=sse $(LTO_FC)
dlamch.o: dlamch.f
	$(FC) $(DLAMC_FFLAGS) -c $< -o $@

la_constants.mod: la_constants.o
la_xisnan.mod: la_xisnan.o
dlartg.o dlassq.o zlartg.o zlassq.o la_xisnan.o: la_constants.o
dlassq.o zlassq.o: la_xisnan.o

R: Makedeps
	@$(MAKE) rhome="$(abs_top_builddir)" $(Rlapack_la) Rlapack_install
	@$(MAKE) $(lapack_la)

$(lapack_la): $(lapack_la_OBJECTS)
	$(SHLIB_LINK) -o $@ $(lapack_la_LDFLAGS) $(lapack_la_OBJECTS) $(lapack_la_LIBADD) -L"$(R_HOME)/lib$(R_ARCH)" -lRlapack -L"$(R_HOME)/lib$(R_ARCH)" -lRblas $(FLIBS_IN_SO)
## <FIXME>
## Not sure if we want to do this ...
	@$(MAKE) rhome="$(abs_top_builddir)" install
## </FIXME>

## Include BLAS here, as with (static) ATLAS that pulls all the
## BLAS routines into one place.
$(Rlapack_la): $(LIBOBJECTS)
	$(DYLIB_LINK) -o $@ $(LIBOBJECTS)  -L"$(R_HOME)/lib$(R_ARCH)" -lRblas $(Rlapack_la_LIBADD)
	#$(AR) -cr $@ $(LIBOBJECTS)

Rlapack_install: $(Rlapack_la) 
	@$(MKINSTALLDIRS) "$(DESTDIR)$(Rexeclibdir)"
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(Rlapack_la) \
	  "$(Rexeclibdir)/$(Rlapack_la)"

install: installdirs
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(lapack_la) "$(DESTDIR)$(Rexecmodulesdir)/$(lapack_la)"
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(Rlapack_la) "$(DESTDIR)$(Rexeclibdir)/$(Rlapack_la)"
installdirs:
	@$(MKINSTALLDIRS) "$(DESTDIR)$(Rexecmodulesdir)"
install-strip-Rlapack:
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(Rlapack_la) "$(DESTDIR)$(Rexeclibdir)/$(Rlapack_la)"
	@if test -n "$(STRIP_SHARED_LIB)"; then \
	  $(STRIP_SHARED_LIB) "$(DESTDIR)$(Rexeclibdir)/$(Rlapack_la)"; \
	fi
install-strip:  install-strip-Rlapack
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(lapack_la) "$(DESTDIR)$(Rexecmodulesdir)/$(lapack_la)"
	@if test -n "$(STRIP_SHARED_LIB)"; then \
	  $(STRIP_SHARED_LIB) "$(DESTDIR)$(Rexecmodulesdir)/$(lapack_la)"; \
	fi
uninstall:
	@rm -f "$(Rexecmodulesdir)/$(lapack_la)"

mostlyclean: clean
clean:
	@-rm -Rf .libs _libs
	@-rm -f Makedeps *.d *.o *.a *.lo *.la *.mod \
	  $(lapack_la) $(Rlapack_la)
distclean: clean
	@-rm -f Makefile
maintainer-clean: distclean

TAGS info dvi check:

distdir: $(DISTFILES)
	@for f in $(DISTFILES); do \
	  test -f $(distdir)/$${f} \
	    || ln $(srcdir)/$${f} $(distdir)/$${f} 2>/dev/null \
	    || cp -p $(srcdir)/$${f} $(distdir)/$${f}; \
	done

## Automagically generated dependencies:

Lapack.o: Lapack.c ../../../src/include/config.h \
 ../../../src/include/Defn.h ../../../src/include/R_ext/Complex.h \
 ../../../src/include/Rinternals.h ../../../src/include/R_ext/Arith.h \
 ../../../src/include/R_ext/libextern.h \
 ../../../src/include/R_ext/Boolean.h ../../../src/include/R_ext/Error.h \
 ../../../src/include/R_ext/Memory.h ../../../src/include/R_ext/Utils.h \
 ../../../src/include/R_ext/Print.h ../../../src/include/R_ext/Rdynload.h \
 ../../../src/include/Rinlinedfuns.h ../../../src/include/Errormsg.h \
 Lapack.h ../../../src/include/R_ext/RS.h ../../../src/include/Rconfig.h \
 ../../../src/include/R_ext/Lapack.h ../../../src/include/R_ext/BLAS.h \
 ../../../src/include/Rmodules/Rlapack.h
