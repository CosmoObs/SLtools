## Comentario nos Makefiles.am sao com #

SUBDIRS = src/

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = sltoolsc.pc

ACLOCAL_AMFLAGS = -I m4
DISTCHECK_CONFIGURE_FLAGS=--enable-gtk-doc --enable-introspection

libsltoolscdocdir = ${docdir}-$(VERSION)
libsltoolscdoc_DATA = \
	README     \
	COPYING    \
	AUTHORS    \
	ChangeLog  \
	INSTALL    \
	NEWS

EXTRA_DIST = \
	$(libsltoolscdoc_DATA) \
	sltoolsc.pc.in

AM_CFLAGS = 

# Copy all the spec files. Of cource, only one is actually used.
dist-hook:
	for specfile in *.spec; do \
		if test -f $$specfile; then \
			cp -p $$specfile $(distdir); \
		fi \
	done
