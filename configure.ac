
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([eclib], [20150228], [john.cremona@gmail.com])
AM_INIT_AUTOMAKE([-Wall])
AC_MSG_NOTICE([Configuring eclib...])
AC_CONFIG_SRCDIR([libsrc])
#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_PROG_AR

# Initialise libtools:
LT_INIT

# Checks for programs:
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LIBTOOL

# Checks for system features and libraries, via macros defined in m4/ directory
AX_CXX_COMPILE_STDCXX_0X
# Uncomment this to check for C++-11 compliance (no C++-11 features are used)
#AX_CXX_COMPILE_STDCXX_11([],[optional])
AX_CXX_HEADER_TR1_UNORDERED_MAP
AX_CXX_HEADER_UNORDERED_MAP
# Boost 1.35.0+ required for Asio (1.25.0 for Thread)
AX_BOOST_BASE([1.35.0])
AX_BOOST_SYSTEM
AX_BOOST_ASIO
AX_BOOST_THREAD

if [test "$ax_cv_boost"        = "yes"] && 
   [test "$ax_cv_boost_system" = "yes"] &&
   [test "$ax_cv_boost_asio"   = "yes"] &&
   [test "$ax_cv_boost_thread" = "yes"]; then
   AC_DEFINE([ECLIB_MULTITHREAD],[1],[Define if multithreading available and requested])
fi

# Check for libraries
# The following AC_ARG_WITH section were written by Volker Braun

dnl Find PARI
AC_ARG_WITH(pari,
  [AS_HELP_STRING([--with-pari=<path>],
                  [prefix of PARI installation. e.g. /usr/local or /usr])],
  [],
  [AC_CHECK_LIB([pari], [pari_kernel_init],[ ],
    AC_MSG_ERROR([Pari library not found.
    Maybe you want to call configure with the --with-pari=<path> option?
    This tells configure where to find the PARI library and headers.
    e.g. --with-pari=/usr/local or --with-pari=/usr]),[ ])]
)
case "$with_pari" in
  ""|yes)
    PARI_LDFLAGS=""
    PARI_CFLAGS=""
    GP="gp"
    ;;
  no)
    AC_MSG_ERROR([Compiling without libpari impossible (but --without-pari given)])
    ;;
  *)
    PARI_LDFLAGS="-L${with_pari}/lib"
    PARI_CFLAGS="-I${with_pari}/include"
    GP="${with_pari}/bin/gp"
    ;;
esac
PARI_LIBS="-lpari"
AC_SUBST(PARI_LIBS)
AC_SUBST(PARI_CFLAGS)
AC_SUBST(PARI_LDFLAGS)
AC_SUBST(GP)

dnl Find the NTL Library
AC_ARG_WITH(ntl,
  [AS_HELP_STRING([--with-ntl=<path>],
                  [prefix of NTL installation. e.g. /usr/local or /usr])],
  [],
  []
)
case "$with_ntl" in
  ""|yes)
    NTL_LDFLAGS=""
    NTL_CFLAGS=""
    ;;
  no)
    AC_MSG_ERROR([Compiling without libntl impossible (but --without-ntl given)])
    ;;
  *)
    NTL_LDFLAGS="-L${with_ntl}/lib"
    NTL_CFLAGS="-I${with_ntl}/include"
    ;;
esac

AC_LANG_PUSH(C++)
CPPFLAGS_save="$CPPFLAGS"
LDFLAGS_save="$LDFLAGS"
CPPFLAGS="$CPPFLAGS $NTL_CFLAGS"
LDFLAGS="$LDFLAGS $NTL_LDFLAGS"
AC_CHECK_HEADER(NTL/ZZ.h, found_ntl=yes, found_ntl=no)
LDFLAGS="$LDFLAGS_save"
CPPFLAGS="$CPPFLAGS_save"
AC_LANG_POP(C++)

if test "x$found_ntl" = "xno"; then
  AC_MSG_ERROR([Could not find NTL, you might want to use --with-ntl=<path> to point to its location])
fi

NTL_LIBS="-lntl"
AC_SUBST(NTL_LIBS)
AC_SUBST(NTL_CFLAGS)
AC_SUBST(NTL_LDFLAGS)

dnl Find the FLINT Library and determine if it has the nmod_mat module
dnl (in FLINT from version 2.3)
dnl
dnl NB If your version of FLINT includes the optional hmod_mat module,
dnl then compile with FLINT_LEVEL=2 by either using configure as usual
dnl and then doing "make CXXFLAGS='-DFLINT_LEVEL=2'", or adding
dnl CXXFLAGS='-DFLINT_LEVEL=2' to the ./configure command.  In either
dnl case, ignore the warnings about FLINT_LEVEL being redfined.
dnl
dnl First set defaults for when no suitable FLINT is found
dnl
AC_ARG_WITH(flint,
  [AS_HELP_STRING([--with-flint=<path>],[prefix of FLINT installation.])],
  [  
    if test "$withval" = "no"; then
       want_flint=no
    elif test "$withval" = "yes"; then
        want_flint=yes
        AC_CHECK_LIB([flint], [fmpz_init, nmod_mat_rref],
        [FLINT_CFLAGS="-DFLINT_LEVEL=1";
         FLINT_LIBS="-lflint";
         AC_MSG_NOTICE([Building with FLINT support from system install])],
        [AC_MSG_NOTICE([No suitable FLINT Library (version >=2.3) found])],
        [])
    else
        want_flint=yes
        FLINT_CFLAGS="-DFLINT_LEVEL=1 -I$withval/include"
	FLINT_LDFLAGS="-L$withval/lib";
        FLINT_LIBS="-lflint"
    fi
    ],
    [want_flint=no]
)

if test "x$want_flint" = "xno"; then
       FLINT_CFLAGS="-DFLINT_LEVEL=0"
       FLINT_LDFLAGS=""
       FLINT_LIBS=""
fi

AC_SUBST(FLINT_LIBS)
AC_SUBST(FLINT_CFLAGS)
AC_SUBST(FLINT_LDFLAGS)

# Checks for header files.
AC_FUNC_ALLOCA
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT64_T
AC_C_RESTRICT
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor memmove memset pow sqrt strchr])

AC_ARG_ENABLE([allprogs],
   AS_HELP_STRING([--disable-allprogs],
                  [only build the library and mwrank, otherwise build all programs and tests]),
     [case "${enableval}" in
       yes) allprogs=true ;;
       no)  allprogs=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-allprogs]) ;;
     esac],[allprogs=true])
AM_CONDITIONAL([ALLPROGS], [test x$allprogs = xtrue])

AC_CONFIG_FILES([
  libsrc/Makefile
  libsrc/gpslave.cc
  tests/Makefile
  progs/Makefile
  doc/Makefile
  man/Makefile
  Makefile
])
AC_OUTPUT
