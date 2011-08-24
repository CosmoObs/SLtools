dnl sltoolsc
dnl major: The most recent interface number that this library implements. 
dnl minor: The implementation number of the current interface. 
dnl micro: The difference between the newest and oldest interfaces that this library implements. In other words, the library implements all the interface numbers in the range from number current - age to current. 

m4_define([sltoolsc_major_version], [0])
m4_define([sltoolsc_minor_version], [0])
m4_define([sltoolsc_micro_version], [0])
m4_define([sltoolsc_version],
          [sltoolsc_major_version.sltoolsc_minor_version.sltoolsc_micro_version])
          
AC_INIT(sltoolsc, [sltoolsc_version])
AC_CONFIG_SRCDIR([perturbative_method/src/perturbative_method.h])
AC_CONFIG_MACRO_DIR([m4])

AC_PREREQ(2.62)
AM_INIT_AUTOMAKE([gnu -Wno-portability])

PACKAGE=sltoolsc
SLTOOLSC_MAJOR_VERSION=sltoolsc_major_version
SLTOOLSC_MINOR_VERSION=sltoolsc_minor_version
SLTOOLSC_MICRO_VERSION=sltoolsc_micro_version
SLTOOLSC_VERSION=$SLTOOLSC_MAJOR_VERSION.$SLTOOLSC_MINOR_VERSION.$SLTOOLSC_MICRO_VERSION
SLTOOLSC_LT_VERSION="$SLTOOLSC_MAJOR_VERSION:$SLTOOLSC_MINOR_VERSION:$SLTOOLSC_MICRO_VERSION"

AC_SUBST(SLTOOLSC_MAJOR_VERSION)
AC_SUBST(SLTOOLSC_MINOR_VERSION)
AC_SUBST(SLTOOLSC_MICRO_VERSION)
AC_SUBST(SLTOOLSC_VERSION)
AC_SUBST(SLTOOLSC_LT_VERSION)

AM_CONFIG_HEADER(config.h)

AM_MAINTAINER_MODE

dnl
dnl Gobject Introspection
dnl
dnl add the correct m4 directory location.
AC_CONFIG_MACRO_DIR([m4])
  
dnl any typos here, and you will end up with a message about HAVE_INTROSPECTION not being defined
GOBJECT_INTROSPECTION_CHECK([0.6.7])

dnl things required by automake
AC_PROG_MAKE_SET

dnl Check for which system.
AC_CANONICAL_HOST

dnl Checks for programs.
AC_LANG(C)
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL   
AC_PROG_LN_S
AC_PROG_CC_STDC
dnl LT_INIT([disable-shared])
dnl AC_DISABLE_SHARED
AC_PROG_LIBTOOL

dnl check for gtk-doc
GTK_DOC_CHECK([1.14],[--flavour no-tmpl])

dnl Set PACKAGE_DATA_DIR in config.h.
if test "x${prefix}" = "xNONE"; then
  AC_DEFINE_UNQUOTED(PACKAGE_DATA_DIR, "${ac_default_prefix}/share/${PACKAGE}-${VERSION}", [PACKAGE_DATA_DIR])
else
  AC_DEFINE_UNQUOTED(PACKAGE_DATA_DIR, "${prefix}/share/${PACKAGE}-${VERSION}", [PACKAGE_DATA_DIR])
fi

dnl Set PACKAGE_SOURCE_DIR in config.h.
packagesrcdir=`cd $srcdir && pwd`
AC_DEFINE_UNQUOTED(PACKAGE_SOURCE_DIR, "${packagesrcdir}",[PACKAGE_SOURCE_DIR])

dnl Dependencies
GLIB_REQUIRED_VERSION=2.20.0
GSL_REQUIRED_VERSION=1.13

AC_SUBST(GLIB_REQUIRED_VERSION)
AC_SUBST(GSL_REQUIRED_VERSION)

AC_DEFINE([GSL_RANGE_CHECK_OFF], [], [optimize gsl access functions])

AC_C_INLINE

if test "$ac_cv_c_inline" != no ; then
  AC_DEFINE([HAVE_INLINE], [1], [use inline functions in GSL])
  AC_SUBST(HAVE_INLINE)
fi

dnl ***************************************************************************
dnl Check for multi-thread.
dnl ***************************************************************************

AC_ARG_WITH([thread-pool-max], 
            AC_HELP_STRING([--with-thread-pool-max=N],[Number of threads in the pool]), 
            [thread_pool_max='yes'], [thread_pool_max='no'])

if test "x$thread_pool_max" = "xyes"; then
  AC_DEFINE_UNQUOTED([CL_THREAD_POOL_MAX],[$with_thread_pool_max], [Number of threads in the pool])
fi

dnl ***************************************************************************
dnl Check for math.h
dnl ***************************************************************************

AC_CHECK_HEADERS(math.h,, AC_MSG_ERROR([math.h header files missing]))
AC_CHECK_LIB([m],[cos])

dnl ***************************************************************************
dnl Check for glib.
dnl ***************************************************************************

PKG_CHECK_MODULES(GLIB, [glib-2.0 >= $GLIB_REQUIRED_VERSION gthread-2.0 >= $GLIB_REQUIRED_VERSION gobject-2.0 >= $GLIB_REQUIRED_VERSION])
AC_SUBST(GLIB_CFLAGS)
AC_SUBST(GLIB_LIBS)
TEST_SHARED_LIBS="$TEST_SHARED_LIBS $GLIB_LIBS"
TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS $GLIB_CFLAGS"

dnl ***************************************************************************
dnl Check for gsl
dnl ***************************************************************************

PKG_CHECK_MODULES(GSL, [gsl >= $GSL_REQUIRED_VERSION])
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
TEST_SHARED_LIBS="$TEST_SHARED_LIBS $GSL_LIBS"
TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS $GSL_CFLAGS"

PKG_CHECK_EXISTS([gsl >= $GSL_GLF_REQUIRED_VERSION], [AC_DEFINE_UNQUOTED([HAVE_GSL_GLF],[1], [gsl support fixed gauss legendre rules])], [])

dnl ***************************************************************************
dnl Appending libs and setting up TEST_SHARED
dnl ***************************************************************************

LIBS="$LIBS $APPEND_LIBS"
TEST_SHARED_LIBS="$TEST_SHARED_LIBS $LIBS"
TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS $CFLAGS"

dnl ***************************************************************************
dnl Setting up
dnl ***************************************************************************

AC_CONFIG_FILES([
sltoolsc.pc
Makefile
perturbative_method/src/Makefile
])

AC_OUTPUT
