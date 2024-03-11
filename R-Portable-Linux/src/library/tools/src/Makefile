#
# ${R_HOME}/src/library/tools/src/Makefile

pkg = tools


srcdir = .
top_srcdir = ../../../..

top_builddir = ../../../..
subdir = src/library/$(pkg)/src
R_HOME = $(top_builddir)
R_SHARE_DIR = $(R_HOME)/share
R_INCLUDE_DIR = $(R_HOME)/include

SOURCES_C = text.c init.c Rmd5.c md5.c signals.c install.c getfmts.c http.c \
  gramLatex.c gramRd.c pdscan.c
DEPENDS = $(SOURCES_C:.c=.d)
OBJECTS = $(SOURCES_C:.c=.o)

# need config.h
PKG_CPPFLAGS = -I../../../include -I$(top_srcdir)/src/include \
  -DHAVE_CONFIG_H -I$(top_srcdir)/src/main
PKG_CFLAGS = $(C_VISIBILITY)

SHLIB = $(pkg).so

all: Makedeps
	@$(MAKE) shlib

Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

Makedeps: Makefile $(DEPENDS)
	@cat $(DEPENDS) >> Makefile
	@touch $@

shlib: $(SHLIB)
	@$(MKINSTALLDIRS) \
	  $(top_builddir)/library/$(pkg)/libs
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(SHLIB) \
	  $(top_builddir)/library/$(pkg)/libs/$(SHLIB)

include $(R_HOME)/etc${R_ARCH}/Makeconf
include $(top_srcdir)/share/make/shlib.mk
## override to use LTO for linking with --enable-lto=R
LTO = 

mostlyclean: clean
clean:
	@-rm -Rf .libs _libs
	@-rm -f Makedeps *.d *.o *$(SHLIB_EXT)
distclean: clean
	@-rm -f Makefile
maintainer-clean: distclean

# suppress #line directives
YFLAGS=-l

$(srcdir)/gramRd.c: # $(srcdir)/gramRd.y
	@$(ECHO) "re-making gramRd.c"
	$(YACC) $(YFLAGS) $(srcdir)/gramRd.y
	$(SHELL) $(top_srcdir)/tools/move-if-change y.tab.c $(srcdir)/gramRd.c	
# depend on gramRd.c to avoid race condition with parallel make
$(srcdir)/gramLatex.c: # $(srcdir)/gramLatex.y $(srcdir)/gramRd.c
	@$(ECHO) "re-making gramLatex.c"
	$(YACC) $(YFLAGS) $(srcdir)/gramLatex.y
	$(SHELL) $(top_srcdir)/tools/move-if-change y.tab.c $(srcdir)/gramLatex.c

## Automagically generated dependencies:
text.o: text.c ../../../../include/R.h ../../../../include/Rconfig.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Complex.h \
 ../../../../include/R_ext/Constants.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Print.h \
 ../../../../include/R_ext/Random.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/RS.h tools.h ../../../../include/Rinternals.h \
 ../../../../include/R_ext/Rdynload.h
init.o: init.c tools.h ../../../../include/R.h \
 ../../../../include/Rconfig.h ../../../../include/R_ext/Arith.h \
 ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Complex.h \
 ../../../../include/R_ext/Constants.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Print.h \
 ../../../../include/R_ext/Random.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/RS.h ../../../../include/Rinternals.h \
 ../../../../include/R_ext/Rdynload.h \
 ../../../../include/R_ext/Visibility.h
Rmd5.o: Rmd5.c ../../../include/config.h ../../../include/Defn.h \
 ../../../../include/R_ext/Complex.h ../../../../include/Rinternals.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/Print.h ../../../../include/R_ext/Rdynload.h \
 ../../../include/Rinlinedfuns.h ../../../include/Errormsg.h tools.h \
 ../../../../include/R.h ../../../../include/Rconfig.h \
 ../../../../include/R_ext/Constants.h ../../../../include/R_ext/Random.h \
 ../../../../include/R_ext/RS.h md5.h
md5.o: md5.c ../../../../include/Rconfig.h md5.h
signals.o: signals.c ../../../include/config.h tools.h \
 ../../../../include/R.h ../../../../include/Rconfig.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Complex.h \
 ../../../../include/R_ext/Constants.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Print.h \
 ../../../../include/R_ext/Random.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/RS.h ../../../../include/Rinternals.h \
 ../../../../include/R_ext/Rdynload.h
install.o: install.c ../../../include/config.h ../../../include/Defn.h \
 ../../../../include/R_ext/Complex.h ../../../../include/Rinternals.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/Print.h ../../../../include/R_ext/Rdynload.h \
 ../../../include/Rinlinedfuns.h ../../../include/Errormsg.h
getfmts.o: getfmts.c ../../../include/config.h ../../../include/Defn.h \
 ../../../../include/R_ext/Complex.h ../../../../include/Rinternals.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/Print.h ../../../../include/R_ext/Rdynload.h \
 ../../../include/Rinlinedfuns.h ../../../include/Errormsg.h \
 ../../../../src/main/RBufferUtils.h
http.o: http.c ../../../../include/Rinternals.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Complex.h \
 ../../../../include/R_ext/Error.h ../../../../include/R_ext/Memory.h \
 ../../../../include/R_ext/Utils.h ../../../../include/R_ext/Print.h \
 ../../../../include/R_ext/Rdynload.h ../../../../include/Rconfig.h \
 tools.h ../../../../include/R.h ../../../../include/R_ext/Constants.h \
 ../../../../include/R_ext/Random.h ../../../../include/R_ext/RS.h
gramLatex.o: gramLatex.c ../../../include/config.h \
 ../../../include/Defn.h ../../../../include/R_ext/Complex.h \
 ../../../../include/Rinternals.h ../../../../include/R_ext/Arith.h \
 ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/Print.h ../../../../include/R_ext/Rdynload.h \
 ../../../include/Rinlinedfuns.h ../../../include/Errormsg.h \
 ../../../include/Parse.h ../../../include/IOStuff.h \
 ../../../../include/R_ext/Parse.h ../../../../include/R_ext/RS.h \
 ../../../../include/Rconfig.h
gramRd.o: gramRd.c ../../../include/config.h ../../../include/Defn.h \
 ../../../../include/R_ext/Complex.h ../../../../include/Rinternals.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/Print.h ../../../../include/R_ext/Rdynload.h \
 ../../../include/Rinlinedfuns.h ../../../include/Errormsg.h \
 ../../../include/Parse.h ../../../include/IOStuff.h \
 ../../../../include/R_ext/Parse.h ../../../../include/R_ext/RS.h \
 ../../../../include/Rconfig.h ../../../../include/Rmath.h \
 ../../../include/Rconnections.h ../../../../include/R_ext/Connections.h
pdscan.o: pdscan.c ../../../../include/R.h ../../../../include/Rconfig.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Complex.h \
 ../../../../include/R_ext/Constants.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Print.h \
 ../../../../include/R_ext/Random.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/RS.h ../../../../include/Rdefines.h \
 ../../../../include/Rinternals.h ../../../../include/R_ext/Rdynload.h \
 tools.h
