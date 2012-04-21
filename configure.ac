#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([eclib], [2012-04-21], [john.cremona@gmail.com])
AM_INIT_AUTOMAKE([-Wall])
AC_MSG_NOTICE([Configuring eclib...])
AC_CONFIG_SRCDIR([libsrc])
#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Initialise libtools:
LT_INIT

# Checks for programs:
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL

# Checks for libraries.

dnl Find GMP
AC_ARG_WITH(gmp,
  [AS_HELP_STRING([--with-gmp=<path>],
                  [prefix of gmp installation. e.g. /usr/local or /usr])],
  [GMP_PREFIX="$with_gmp"],
  [AC_CHECK_LIB([gmp], [__gmpz_init], ,
    AC_MSG_ERROR([gmp library not found.
    Maybe you want to call configure with the --with-gmp=<path> option?
    This tells configure where to find the gmp library and headers.
    e.g. --with-gmp=/usr/local or --with-gmp=/usr]))]
)
if test "x$GMP_PREFIX" != x ; then
  GMP_LIBS=" -L${GMP_PREFIX}/lib"
  GMP_CFLAGS=" -I${GMP_PREFIX}/include"
else
  GMP_LIBS=""
  GMP_CFLAGS=""
fi
GMP_LIBS="${GMP_LIBS} -lgmp"
AC_SUBST(GMP_LIBS)
AC_SUBST(GMP_CFLAGS)

# The following AC_ARG_WITH section were written by Volker Braun

dnl Find PARI
AC_ARG_WITH(pari,
  [AS_HELP_STRING([--with-pari=<path>],
                  [prefix of PARI installation. e.g. /usr/local or /usr])],
  [PARI_PREFIX="$with_pari"],
  [AC_CHECK_LIB([pari], [pari_kernel_init], ,
    AC_MSG_ERROR([Pari library not found.
    Maybe you want to call configure with the --with-pari=<path> option?
    This tells configure where to find the PARI library and headers.
    e.g. --with-pari=/usr/local or --with-pari=/usr]))]
)
if test "x$PARI_PREFIX" != x ; then
  PARI_LIBS=" -L${PARI_PREFIX}/lib"
  PARI_CFLAGS=" -I${PARI_PREFIX}/include"
else
  PARI_LIBS=""
  PARI_CFLAGS=""
fi
PARI_LIBS="${PARI_LIBS} -lpari"
AC_SUBST(PARI_LIBS)
AC_SUBST(PARI_CFLAGS)

dnl Find the NTL Library
AC_ARG_WITH(ntl,
  [AS_HELP_STRING([--with-ntl=<path>],
                  [prefix of NTL installation. e.g. /usr/local or /usr])],
  [NTL_PREFIX="$with_ntl"],
  [AC_CHECK_LIB([ntl], [_ntl_gmul], ,
    AC_MSG_ERROR([NTL Library not found.
    Maybe you want to call configure with the --with-ntl=<path> option?
    This tells configure where to find the NTL library and headers.
    e.g. --with-ntl=/usr/local or --with-ntl=/usr]))]
)
if test "x$NTL_PREFIX" != x ; then
  NTL_LIBS=" -L${NTL_PREFIX}/lib"
  NTL_CFLAGS="${NTL_ARITH} -I${NTL_PREFIX}/include"
else
  NTL_LIBS=""
  NTL_CFLAGS="${NTL_ARITH}"
fi
NTL_LIBS="${NTL_LIBS} -lgmp -lntl"
AC_SUBST(NTL_LIBS)
AC_SUBST(NTL_CFLAGS)

AM_CXXFLAGS="${NTL_CFLAGS} ${PARI_CFLAGS} ${GMP_CFLAGS} -DNTL_ALL -DUSE_PARI_FACTORING -DMETHOD=2 -DNEW_OP_ORDER"
AM_LDFLAGS="${NTL_LIBS} ${PARI_LIBS} ${GMP_LIBS}"

AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])

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
   AS_HELP_STRING([--enable-allprogs],
                  [build all programs and tests (default is to only build the library and mwrank)]),
     [case "${enableval}" in
       yes) allprogs=true ;;
       no)  allprogs=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-allprogs]) ;;
     esac],[allprogs=true])
AM_CONDITIONAL([ALLPROGS], [test x$allprogs = xtrue])

AC_CONFIG_FILES([
  libsrc/Makefile
  tests/Makefile
  Makefile
])
AC_OUTPUT
