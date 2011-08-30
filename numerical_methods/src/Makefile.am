## Process this file with automake to produce Makefile.in
## 

INCLUDES = 

AM_CFLAGS = \
	 -Wall \
	 -DG_LOG_DOMAIN=\"SLTOOLSC\" \
	 $(GLIB_CFLAGS) \
	 $(GSL_CFLAGS)  \
	 -I$(top_srcdir)

SUBDIRS = 

lib_LTLIBRARIES = libsltoolsc.la

libsltoolsc_la_SOURCES =          \
	hausdorff_distance.h

libsltoolsc_la_LIBADD = \
	$(GLIB_LIBS) \
	$(GSL_LIBS)

libsltoolsc_la_LDFLAGS = -version-info $(SLTOOLSC_LT_VERSION)

AM_CPPFLAGS = \
	$(GLIB_CFLAGS) \
	$(GSL_CFLAGS)

nobase_sltoolsc_include_HEADERS = \
	hausdorff_distance.h

sltoolsc_includedir = \
	$(pkgincludedir)

EXTRA_DIST = \
	$(nobase_sltoolsc_include_HEADERS) \
	build_cfg.h.in

-include $(INTROSPECTION_MAKEFILE)
INTROSPECTION_GIRS =
INTROSPECTION_SCANNER_ARGS = --add-include-path=$(srcdir) -I$(srcdir)/../ -I$(top_builddir) \
  --identifier-prefix=slc --warn-all --verbose

INTROSPECTION_COMPILER_ARGS = --includedir=$(srcdir)

if HAVE_INTROSPECTION

introspection_sources = perturbative_method.h $(libsltoolsc_la_SOURCES)

SLtoolsC-1.0.gir: libsltoolsc.la
SLtoolsC_1_0_gir_INCLUDES = GObject-2.0
SLtoolsC_1_0_gir_CFLAGS = $(INCLUDES)
SLtoolsC_1_0_gir_LIBS = libsltoolsc.la
SLtoolsC_1_0_gir_FILES = $(introspection_sources)
INTROSPECTION_GIRS += SLtoolsC-1.0.gir

girdir = $(datadir)/gir-1.0
gir_DATA = $(INTROSPECTION_GIRS)
  
typelibdir = $(libdir)/girepository-1.0
typelib_DATA = $(INTROSPECTION_GIRS:.gir=.typelib)
      
CLEANFILES = $(gir_DATA) $(typelib_DATA)

endif
