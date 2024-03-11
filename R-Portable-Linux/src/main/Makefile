#
# ${R_HOME}/src/main/Makefile


srcdir = .
top_srcdir = ../..

top_builddir = ../..
subdir = src/main
# next is needed for shared BLAS
R_HOME=$(top_builddir)

include $(top_builddir)/Makeconf

SOURCES_C = \
	CommandLineArgs.c \
	Rdynload.c Renviron.c RNG.c \
	agrep.c altclasses.c altrep.c apply.c arithmetic.c array.c attrib.c \
	bind.c builtin.c \
	character.c clippath.c coerce.c \
	colors.c complex.c connections.c context.c cum.c \
	dcf.c datetime.c debug.c deparse.c devices.c \
	dotcode.c dounzip.c dstruct.c duplicate.c \
	edit.c engine.c envir.c errors.c eval.c \
	flexiblas.c format.c \
	gevents.c gram.c gram-ex.c graphics.c grep.c \
	identical.c inlined.c inspect.c internet.c iosupport.c \
	lapack.c list.c localecharset.c logic.c \
	machine.c main.c mapply.c mask.c match.c memory.c \
	names.c \
	objects.c options.c \
	paste.c patterns.c platform.c plot.c plot3d.c plotmath.c \
	print.c printarray.c printvector.c printutils.c qsort.c \
	radixsort.c random.c raw.c registration.c relop.c rlocale.c \
	saveload.c scan.c seq.c serialize.c sort.c source.c split.c \
	sprintf.c startup.c subassign.c subscript.c subset.c summary.c sysutils.c \
	times.c \
	unique.c util.c \
	version.c \
	g_alab_her.c g_cntrlify.c g_fontdb.c g_her_glyph.c

SOURCES_F =  xxxpr.f

## If the substitutes are needed, the corresponding objects are put by
## configure in  
EXTRA_SOURCES_C = \
	Rmain.c alloca.c mkdtemp.c strdup.c strncasecmp.c

DEPENDS = $(SOURCES_C:.c=.d) $(EXTRA_SOURCES_C:.c=.d)
SOURCES = $(SOURCES_C) $(SOURCES_F)
OBJECTS = $(SOURCES_C:.c=.o) $(SOURCES_F:.f=.o)  
HEADERS = \
	RBufferUtils.h Rcomplex.h Rstrptime.h \
	arithmetic.h \
	basedecl.h \
	contour-common.h \
	datetime.h \
	duplicate.h \
	gzio.h \
	machar.c \
	qsort-body.c \
	rlocale_data.h rlocale_tolower.h rlocale_toupper.h rlocale_widths.h \
	split-incl.c \
	unzip.h \
	valid_utf8.h \
	xspline.c \
	g_cntrlify.h g_control.h g_extern.h g_her_metr.h g_jis.h


distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)
DISTFILES = Makefile.in Makefile.win \
	$(HEADERS) \
	$(SOURCES_C) \
	$(EXTRA_SOURCES_C) \
	$(SOURCES_F) \
	gram.y

TRE_CPPFLAGS =  -I$(top_srcdir)/src/extra
XDR_CPPFLAGS = # -I$(top_srcdir)/src/extra/xdr
XDR_CPPFLAGS = -I/usr/include/tirpc

## platform.c needs $(CURL_CPPFLAGS).
ALL_CPPFLAGS = $(TRE_CPPFLAGS) $(XDR_CPPFLAGS) $(R_XTRA_CPPFLAGS) \
  $(CURL_CPPFLAGS) $(CPPFLAGS) -I$(top_srcdir)/src/nmath $(DEFS)
ALL_CFLAGS = $(ALL_CFLAGS_LO)
ALL_FFLAGS = $(ALL_FFLAGS_LO)
#ALL_FFLAGS = $(ALL_FFLAGS_PIE)


## use an explicit library: there might be an unsatisfactory -lintl around
R_TRE =  ../extra/tre/libtre.a
R_XDR = # ../extra/xdr/libxdr.a
R_LIBINTL = # ../extra/intl/libintl.a
R_TZONE = # ../extra/tzone/libtz.a

MAIN_LIBS = ../unix/libunix.a ../appl/libappl.a ../nmath/libnmath.a
MAIN_OBJS = `ls ../unix/*.o ../appl/*.o ../nmath/*.o`
EXTRA_STATIC_LIBS = $(R_TRE) $(R_XDR) $(R_LIBINTL) $(R_TZONE)
STATIC_LIBS = $(MAIN_LIBS) $(EXTRA_STATIC_LIBS)

EXTRA_LIBS = $(BLAS_LIBS) $(FLIBS)   $(READLINE_LIBS) $(LIBS)

R_binary = R.bin
R_bin_OBJECTS = Rmain.o #$(OBJECTS)
#R_bin_LDADD = $(MAIN_OBJS) $(EXTRA_STATIC_LIBS) $(EXTRA_LIBS)
## Linked against -lRblas because -lR is and otherwise ld complains.
R_bin_LDADD = -lR -lRblas
## This should depend on MAIN_OBJS not MAIN_LIBS, but we can't use that.
## There is also a dependence on libRblas when that is internal and static.
#R_bin_DEPENDENCIES = $(MAIN_LIBS) $(EXTRA_STATIC_LIBS)# $(top_builddir)/etc/R.exp

libR_la = libR$(R_DYLIB_EXT)
libR_la_OBJECTS = $(OBJECTS)
libR_la_LIBADD =  $(MAIN_OBJS) $(EXTRA_STATIC_LIBS) $(EXTRA_LIBS)  # -Wl,-bE:$(top_builddir)/etc/R.exp
libR_la_DEPENDENCIES = $(STATIC_LIBS) $(R_TZONE)  # $(top_builddir)/etc/R.exp

## The next is needed for macOS only at present
LIBR_LDFLAGS = 


all: R

Makefile: $(srcdir)/Makefile.in \
  $(top_builddir)/config.status \
  ../include/config.h \
  ../include/Rversion.h \
  $(SOURCES_C) $(EXTRA_SOURCES_C)
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

Makedeps: Makefile $(DEPENDS)
	@cat $(DEPENDS) >> Makefile
	@touch $@

## This target has been changed to ensure that R.bin and libR get
## installed in the build tree if necessary, even if the corresponding
## objects are not re-made.
## Note that dependencies in the test Makefiles used to depend on these objects,
## so copy-if-change is used.
R: Makedeps
	@$(MAKE) install-lib-local
#	@$(MAKE) install-static-local
	@$(MAKE) install-bin-local


## is this portable?  Documented as GNU extension.
../include/config.h ../include/Rversion.h:
	(cd $(@D); $(MAKE) $(@F))

$(R_binary): $(R_bin_OBJECTS) $(R_bin_DEPENDENCIES)
	$(MAIN_LINK) -o $@ $(R_bin_OBJECTS) $(R_bin_LDADD)

libR.a: $(OBJECTS) $(STATIC_LIBS)
	-@mkdir libs
	@(cd libs; for l in $(STATIC_LIBS); do $(AR) -x ../$$l; done)
	@rm -Rf $@
	$(AR) -cr $@ $(OBJECTS) libs/*o
	$(RANLIB) $@
	@rm -Rf libs

install-bin-local: $(R_binary)
	@$(MAKE) rhome="$(abs_top_builddir)" install-bin

install-static-local: libR.a
	@$(MAKE) rhome="$(abs_top_builddir)" install-static

$(top_builddir)/etc/R.exp: $(OBJECTS) $(MAIN_LIBS)
	@$(SHELL) $(top_srcdir)/tools/ldAIX4 -o $@ $(OBJECTS) $(MAIN_LIBS)

## is this portable?  Documented as GNU extension.
$(MAIN_LIBS):
	(cd $(@D); $(MAKE) $(@F))

../extra/intl/libintl.a:
	(cd $(@D); $(MAKE))


libR: $(libR_la)

$(libR_la): $(libR_la_OBJECTS) $(libR_la_DEPENDENCIES)
	$(DYLIB_LINK) $(LIBR_LDFLAGS) -o $@ $(libR_la_OBJECTS) $(libR_la_LIBADD)

install-lib-local: $(libR_la)
	@$(MAKE) rhome="$(abs_top_builddir)" install-lib


# suppress #line directives
YFLAGS=-l

$(srcdir)/gram.c: # $(srcdir)/gram.y
	@$(ECHO) "re-making gram.c"
	$(YACC) $(YFLAGS) $(srcdir)/gram.y
	$(SHELL) $(top_srcdir)/tools/move-if-change y.tab.c $(srcdir)/gram.c

install: installdirs
	@$(MAKE) install-bin
	@if test -f $(libR_la); then $(MAKE) install-lib; fi
#	@$(MAKE) install-static
## may not need $(rhome)/lib if not static/shared libR, external blas and lapack.
installdirs:
	@$(MKINSTALLDIRS) "$(DESTDIR)$(Rexecbindir2)"
install-bin: installdirs
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(R_binary) "$(DESTDIR)$(Rexecbindir2)/R"
install-lib: installdirs
	@$(MKINSTALLDIRS) "$(DESTDIR)$(Rexeclibdir)"
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(libR_la) "$(DESTDIR)$(Rexeclibdir)/$(libR_la)"
install-static: installdirs
	@$(MKINSTALLDIRS) "$(DESTDIR)$(Rexeclibdir)"
	@$(SHELL) $(top_srcdir)/tools/copy-if-change libR.a "$(DESTDIR)$(Rexeclibdir)/libR.a"
install-strip: installdirs
	@${INSTALL_PROGRAM} -s $(R_binary) "$(DESTDIR)$(Rexecbindir2)/R"
	@if test -f $(libR_la); then $(MAKE) install-lib; fi
	@if test -n "$(STRIP_SHARED_LIB)"; then \
	  if test -f $(libR_la); then $(STRIP_SHARED_LIB) "$(DESTDIR)$(Rexeclibdir)/$(libR_la)"; fi; \
	fi
#	@$(MAKE) install-strip-static
install-strip-static: installdirs
	@$(MAKE) install-static
	@if test -n "$(STRIP_STATIC_LIB)"; then \
	  if test -f libR.a; then $(STRIP_STATIC_LIB) "$(DESTDIR)$(Rexeclibdir)/libR.a"; fi; \
	fi
uninstall:
	@rm -f "$(DESTDIR)$(Rexecbindir)/exec/$(R_ARCH)/R" "$(DESTDIR)$(Rexecbindir)/R.bin"
	@rm -f "$(DESTDIR)$(Rexeclibdir)/libR$(R_DYLIB_EXT)"
	@rm -f "$(DESTDIR)$(Rexeclibdir)/libR.a"

mostlyclean: clean
clean:
	@-rm -f $(top_builddir)/etc/R.exp
	@-rm -Rf .libs _libs
	@-rm -f *core Makedeps *.d *.o *.lo *.la *$(R_DYLIB_EXT) \
	  $(R_binary) libR.a
distclean: clean
	@-rm -f Makefile
maintainer-clean: distclean
	@$(ECHO) "This command is intended for maintainers to use; it"
	@$(ECHO) "deletes files that may need special rules to rebuild"
	@-rm -f $(srcdir)/gram.c $(srcdir)/gramLatex.c $(srcdir)/gramRd.c

tags: TAGS
TAGS: $(SOURCES) $(EXTRA_SOURCES_C) $(HEADERS)
	etags $(SOURCES) $(EXTRA_SOURCES_C) $(HEADERS)

## Unused targets
info dvi check:

distdir: $(DISTFILES)
	@for f in $(DISTFILES); do \
	  test -f $(distdir)/$${f} \
	    || ln $(srcdir)/$${f} $(distdir)/$${f} 2>/dev/null \
	    || cp -p $(srcdir)/$${f} $(distdir)/$${f}; \
	done

## Automagically generated dependencies:

CommandLineArgs.o: CommandLineArgs.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/RStartup.h
Rdynload.o: Rdynload.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h \
 ../../src/include/Rmath.h ../../src/include/Rdynpriv.h
Renviron.o: Renviron.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Fileio.h \
 ../../src/include/Rinterface.h
RNG.o: RNG.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/Random.h \
 ../../src/nmath/nmath2.h
agrep.o: agrep.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h ../../src/extra/tre/tre.h \
 ../../src/extra/tre/tre-config.h
altclasses.o: altclasses.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/Altrep.h \
 ../../src/include/Print.h ../../src/include/Defn.h \
 ../../src/include/R_ext/PrtUtil.h ../../src/include/R_ext/Itermacros.h
altrep.o: altrep.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/R_ext/Altrep.h
apply.o: apply.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h
arithmetic.o: arithmetic.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h \
 ../../src/include/Rmath.h ../../src/include/R_ext/Itermacros.h \
 arithmetic.h
array.o: array.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rmath.h \
 ../../src/include/R_ext/RS.h ../../src/include/Rconfig.h \
 ../../src/include/R_ext/BLAS.h ../../src/include/R_ext/Itermacros.h \
 duplicate.h Rcomplex.h
attrib.o: attrib.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rmath.h
bind.o: bind.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/PrtUtil.h \
 ../../src/include/R_ext/Itermacros.h RBufferUtils.h duplicate.h
builtin.o: builtin.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Print.h \
 ../../src/include/Defn.h ../../src/include/R_ext/PrtUtil.h \
 ../../src/include/Fileio.h ../../src/include/Rconnections.h \
 ../../src/include/R_ext/Connections.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h
character.o: character.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h \
 ../../src/include/R_ext/RS.h ../../src/include/Rconfig.h \
 ../../src/include/R_ext/Itermacros.h ../../src/include/rlocale.h \
 RBufferUtils.h
clippath.o: clippath.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h
coerce.o: coerce.c ../../src/include/config.h ../../src/include/Parse.h \
 ../../src/include/IOStuff.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/R_ext/Parse.h ../../src/include/Internal.h \
 ../../src/include/Rmath.h ../../src/include/Print.h \
 ../../src/include/Defn.h ../../src/include/R_ext/PrtUtil.h \
 ../../src/include/R_ext/Itermacros.h
colors.o: colors.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h
complex.o: complex.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rmath.h arithmetic.h \
 Rcomplex.h ../../src/include/R_ext/Itermacros.h
connections.o: connections.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h \
 ../../src/include/Fileio.h ../../src/include/Rconnections.h \
 ../../src/include/R_ext/Connections.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h ../../src/include/R_ext/Riconv.h gzio.h
context.o: context.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h
cum.o: cum.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h
dcf.o: dcf.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rconnections.h \
 ../../src/include/R_ext/Connections.h ../../src/extra/tre/tre.h \
 ../../src/extra/tre/tre-config.h
datetime.o: datetime.c ../../src/include/config.h \
 ../../src/include/Rmath.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h Rstrptime.h ../../src/include/rlocale.h
debug.o: debug.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h
deparse.o: deparse.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Print.h \
 ../../src/include/Defn.h ../../src/include/R_ext/PrtUtil.h \
 ../../src/include/Fileio.h RBufferUtils.h \
 ../../src/include/Rconnections.h ../../src/include/R_ext/Connections.h
devices.o: devices.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Graphics.h \
 ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h ../../src/include/Rgraphics.h \
 ../../src/include/GraphicsBase.h
dotcode.o: dotcode.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rmath.h \
 ../../src/include/Rdynpriv.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h
dounzip.o: dounzip.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Fileio.h unzip.h ../../src/include/Rconnections.h \
 ../../src/include/R_ext/Connections.h
dstruct.o: dstruct.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h
duplicate.o: duplicate.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h duplicate.h
edit.o: edit.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Print.h \
 ../../src/include/Defn.h ../../src/include/R_ext/PrtUtil.h \
 ../../src/include/Fileio.h ../../src/include/Parse.h \
 ../../src/include/IOStuff.h ../../src/include/R_ext/Parse.h \
 ../../src/include/Rinterface.h ../../src/include/R_ext/RStartup.h
engine.o: engine.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h \
 ../../src/include/R_ext/Applic.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h ../../src/include/R_ext/BLAS.h \
 ../../src/include/Rmath.h ../../src/include/rlocale.h xspline.c
envir.o: envir.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/Callbacks.h
errors.o: errors.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Startup.h \
 ../../src/include/R_ext/RStartup.h ../../src/include/Rconnections.h \
 ../../src/include/R_ext/Connections.h ../../src/include/Rinterface.h \
 ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h ../../src/include/Rmath.h \
 ../../src/include/rlocale.h
eval.o: eval.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rinterface.h \
 ../../src/include/Fileio.h arithmetic.h ../../src/include/Rmath.h \
 ../../src/include/Rdynpriv.h ../../src/include/Parse.h \
 ../../src/include/IOStuff.h ../../src/include/R_ext/Parse.h
flexiblas.o: flexiblas.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h
format.o: format.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Rmath.h ../../src/include/Print.h \
 ../../src/include/Defn.h ../../src/include/R_ext/PrtUtil.h \
 ../../src/include/R_ext/Itermacros.h
gevents.o: gevents.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Rmath.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h
gram.o: gram.c ../../src/include/config.h ../../src/include/IOStuff.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Fileio.h \
 ../../src/include/Parse.h ../../src/include/R_ext/Parse.h \
 ../../src/include/rlocale.h ../../src/include/Rconnections.h \
 ../../src/include/R_ext/Connections.h
gram-ex.o: gram-ex.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h
graphics.o: graphics.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Graphics.h \
 ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h ../../src/include/Rgraphics.h \
 ../../src/include/Rmath.h
grep.o: grep.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h ../../src/extra/tre/tre.h \
 ../../src/extra/tre/tre-config.h
identical.o: identical.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h
inlined.o: inlined.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h
inspect.o: inspect.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h
internet.o: internet.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h \
 ../../src/include/Rconnections.h ../../src/include/R_ext/Connections.h \
 ../../src/include/Rdynpriv.h ../../src/include/Rmodules/Rinternet.h
iosupport.o: iosupport.c ../../src/include/config.h \
 ../../src/include/IOStuff.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h
lapack.o: lapack.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Rdynpriv.h ../../src/include/Rmodules/Rlapack.h
list.o: list.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h
localecharset.o: localecharset.c ../../src/include/config.h
logic.o: logic.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/Itermacros.h
machine.o: machine.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rinterface.h machar.c
main.o: main.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h \
 ../../src/include/Rinterface.h ../../src/include/IOStuff.h \
 ../../src/include/Fileio.h ../../src/include/Parse.h \
 ../../src/include/R_ext/Parse.h ../../src/include/Startup.h \
 ../../src/include/R_ext/RStartup.h ../../src/include/R_ext/Callbacks.h \
 ../../src/include/R_ext/RS.h ../../src/include/Rconfig.h
mapply.o: mapply.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h
mask.o: mask.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h
match.o: match.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h
memory.o: memory.c ../../src/include/config.h \
 ../../src/include/R_ext/RS.h ../../src/include/Rconfig.h \
 ../../src/include/R_ext/Print.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h \
 ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h \
 ../../src/include/R_ext/Rallocators.h ../../src/include/Rmath.h \
 ../../src/include/Print.h ../../src/include/Defn.h \
 ../../src/include/R_ext/PrtUtil.h RBufferUtils.h
names.o: names.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Print.h \
 ../../src/include/Defn.h ../../src/include/R_ext/PrtUtil.h arithmetic.h \
 ../../src/include/Rinterface.h
objects.o: objects.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h
options.o: options.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Print.h \
 ../../src/include/Defn.h ../../src/include/R_ext/PrtUtil.h
paste.o: paste.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Print.h \
 ../../src/include/Defn.h ../../src/include/R_ext/PrtUtil.h \
 RBufferUtils.h
patterns.o: patterns.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h
platform.o: platform.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h \
 ../../src/include/Rinterface.h RBufferUtils.h ../../src/include/Fileio.h \
 ../../src/extra/tre/tre.h ../../src/extra/tre/tre-config.h
plot.o: plot.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Graphics.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h ../../src/include/Rgraphics.h \
 ../../src/include/Print.h ../../src/include/Defn.h \
 ../../src/include/R_ext/PrtUtil.h ../../src/include/Rmath.h
plot3d.o: plot3d.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rmath.h \
 ../../src/include/Graphics.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h ../../src/include/Rgraphics.h \
 ../../src/include/Print.h ../../src/include/Defn.h \
 ../../src/include/R_ext/PrtUtil.h contour-common.h
plotmath.o: plotmath.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/rlocale.h \
 ../../src/include/Rmath.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h
print.o: print.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Print.h \
 ../../src/include/Defn.h ../../src/include/R_ext/PrtUtil.h \
 ../../src/include/Fileio.h ../../src/include/Rconnections.h \
 ../../src/include/R_ext/Connections.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h
printarray.o: printarray.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Print.h \
 ../../src/include/Defn.h ../../src/include/R_ext/PrtUtil.h
printvector.o: printvector.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Print.h \
 ../../src/include/Defn.h ../../src/include/R_ext/PrtUtil.h \
 ../../src/include/R_ext/Itermacros.h
printutils.o: printutils.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Rmath.h \
 ../../src/include/Print.h ../../src/include/Defn.h \
 ../../src/include/R_ext/PrtUtil.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h ../../src/include/Rconnections.h \
 ../../src/include/R_ext/Connections.h RBufferUtils.h \
 ../../src/include/rlocale.h
qsort.o: qsort.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rmath.h \
 ../../src/include/R_ext/RS.h ../../src/include/Rconfig.h qsort-body.c
radixsort.o: radixsort.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h
random.o: random.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/R_ext/Itermacros.h ../../src/include/R_ext/Random.h \
 ../../src/include/R_ext/RS.h ../../src/include/Rconfig.h \
 ../../src/include/Rmath.h ../../src/include/Internal.h
raw.o: raw.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h
registration.o: registration.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/Applic.h \
 ../../src/include/R_ext/RS.h ../../src/include/Rconfig.h \
 ../../src/include/R_ext/BLAS.h ../../src/include/R_ext/Linpack.h \
 ../../src/include/R_ext/Callbacks.h ../../src/include/Rdynpriv.h \
 basedecl.h
relop.o: relop.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rmath.h \
 ../../src/include/R_ext/Itermacros.h
rlocale.o: rlocale.c ../../src/include/config.h \
 ../../src/include/rlocale.h ../../src/include/R_ext/Riconv.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h rlocale_widths.h
saveload.o: saveload.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h \
 ../../src/include/Rinterface.h ../../src/include/Rmath.h \
 ../../src/include/Fileio.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h RBufferUtils.h \
 /usr/include/tirpc/rpc/types.h /usr/include/tirpc/netconfig.h \
 /usr/include/tirpc/rpc/xdr.h ../../src/include/Rconnections.h \
 ../../src/include/R_ext/Connections.h
scan.o: scan.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Fileio.h \
 ../../src/include/Rconnections.h ../../src/include/R_ext/Connections.h \
 ../../src/include/Print.h ../../src/include/Defn.h \
 ../../src/include/R_ext/PrtUtil.h ../../src/include/rlocale.h \
 RBufferUtils.h
seq.o: seq.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rmath.h \
 ../../src/include/R_ext/Itermacros.h RBufferUtils.h
serialize.o: serialize.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Rmath.h \
 ../../src/include/Fileio.h ../../src/include/Rversion.h \
 ../../src/include/R_ext/Riconv.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h /usr/include/tirpc/rpc/types.h \
 /usr/include/tirpc/netconfig.h /usr/include/tirpc/rpc/xdr.h \
 ../../src/include/Rconnections.h ../../src/include/R_ext/Connections.h
sort.o: sort.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rmath.h \
 ../../src/include/R_ext/RS.h ../../src/include/Rconfig.h \
 ../../src/include/R_ext/Itermacros.h
source.o: source.c ../../src/include/config.h ../../src/include/Parse.h \
 ../../src/include/IOStuff.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/R_ext/Parse.h ../../src/include/Internal.h \
 ../../src/include/Fileio.h ../../src/include/Rconnections.h \
 ../../src/include/R_ext/Connections.h
split.o: split.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/Itermacros.h \
 split-incl.c
sprintf.o: sprintf.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h RBufferUtils.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h
startup.o: startup.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Fileio.h ../../src/include/Startup.h \
 ../../src/include/R_ext/RStartup.h
subassign.o: subassign.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h \
 ../../src/include/R_ext/RS.h ../../src/include/Rconfig.h \
 ../../src/include/R_ext/Itermacros.h
subscript.o: subscript.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/Itermacros.h
subset.o: subset.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h
summary.o: summary.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/Itermacros.h \
 duplicate.h
sysutils.o: sysutils.c ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/Internal.h \
 ../../src/include/R_ext/Riconv.h ../../src/include/Rinterface.h \
 ../../src/include/rlocale.h RBufferUtils.h \
 ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h
times.o: times.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h
unique.o: unique.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/Altrep.h \
 ../../src/include/R_ext/Itermacros.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h ../../src/include/R_ext/Random.h
util.o: util.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/R_ext/Itermacros.h \
 ../../src/include/Print.h ../../src/include/Defn.h \
 ../../src/include/R_ext/PrtUtil.h ../../src/include/rlocale.h \
 ../../src/include/R_ext/Riconv.h valid_utf8.h \
 ../../src/include/R_ext/Applic.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h ../../src/include/R_ext/BLAS.h \
 ../../src/include/Rmath.h
version.o: version.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h \
 ../../src/include/Internal.h ../../src/include/Rversion.h
g_alab_her.o: g_alab_her.c g_extern.h ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h g_control.h g_her_metr.h
g_cntrlify.o: g_cntrlify.c g_extern.h ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h g_control.h g_cntrlify.h \
 g_jis.h
g_fontdb.o: g_fontdb.c g_extern.h ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h
g_her_glyph.o: g_her_glyph.c g_extern.h ../../src/include/config.h \
 ../../src/include/Defn.h ../../src/include/R_ext/Complex.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Print.h \
 ../../src/include/R_ext/Rdynload.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/GraphicsEngine.h \
 ../../src/include/R_ext/GraphicsDevice.h
Rmain.o: Rmain.c ../../src/include/Rinterface.h \
 ../../src/include/R_ext/Boolean.h
alloca.o: alloca.c ../../src/include/config.h
mkdtemp.o: mkdtemp.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/Rdynload.h \
 ../../src/include/Rinlinedfuns.h ../../src/include/Errormsg.h
strdup.o: strdup.c ../../src/include/config.h
strncasecmp.o: strncasecmp.c ../../src/include/config.h
