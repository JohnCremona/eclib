#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([eclib-flint], [2012-10-5], [john.cremona@gmail.com])
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

# The following AC_ARG_WITH section were written by Volker Braun

dnl Find GMP
AC_ARG_WITH(gmp,
  [AS_HELP_STRING([--with-gmp=<path>],
                  [prefix of GMP installation. e.g. /usr/local or /usr])],
  [],
  [AC_CHECK_LIB([gmp], [__gmpz_add],[ ],
    AC_MSG_ERROR([GMP library not found.
    Maybe you want to call configure with the --with-gmp=<path> option?
    This tells configure where to find the GMP library and headers.
    e.g. --with-gmp=/usr/local or --with-gmp=/usr]))]
)
case "$with_gmp" in
  ""|yes)
    GMP_LIBS=""
    GMP_CFLAGS=""
    ;;
  no)
    AC_MSG_ERROR([Compiling without libgmp impossible (but --without-gmp given)])
    ;;
  *)
    GMP_LIBS="-L${with_gmp}/lib"
    GMP_CFLAGS="-I${with_gmp}/include"
    ;;
esac
GMP_LIBS="${GMP_LIBS} -lgmp"
AC_SUBST(GMP_LIBS)
AC_SUBST(GMP_CFLAGS)

dnl Find PARI
AC_ARG_WITH(pari,
  [AS_HELP_STRING([--with-pari=<path>],
                  [prefix of PARI installation. e.g. /usr/local or /usr])],
  [],
  [AC_CHECK_LIB([pari], [pari_kernel_init],[ ],
    AC_MSG_ERROR([Pari library not found.
    Maybe you want to call configure with the --with-pari=<path> option?
    This tells configure where to find the PARI library and headers.
    e.g. --with-pari=/usr/local or --with-pari=/usr]),[${GMP_LIBS}])]
)
case "$with_pari" in
  ""|yes)
    PARI_LIBS=""
    PARI_CFLAGS=""
    ;;
  no)
    AC_MSG_ERROR([Compiling without libpari impossible (but --without-pari given)])
    ;;
  *)
    PARI_LIBS="-L${with_pari}/lib"
    PARI_CFLAGS="-I${with_pari}/include"
    ;;
esac
PARI_LIBS="${PARI_LIBS} -lpari"
AC_SUBST(PARI_LIBS)
AC_SUBST(PARI_CFLAGS)

dnl Find the NTL Library
AC_ARG_WITH(ntl,
  [AS_HELP_STRING([--with-ntl=<path>],
                  [prefix of NTL installation. e.g. /usr/local or /usr])],
  [],
  [AC_CHECK_LIB([ntl], [_ntl_gmul],[ ],
    AC_MSG_ERROR([NTL Library not found.
    Maybe you want to call configure with the --with-ntl=<path> option?
    This tells configure where to find the NTL library and headers.
    e.g. --with-ntl=/usr/local or --with-ntl=/usr]),[${GMP_LIBS}])]
)
case "$with_ntl" in
  ""|yes)
    NTL_LIBS=""
    NTL_CFLAGS=""
    ;;
  no)
    AC_MSG_ERROR([Compiling without libntl impossible (but --without-ntl given)])
    ;;
  *)
    NTL_LIBS="-L${with_ntl}/lib"
    NTL_CFLAGS="-I${with_ntl}/include"
    ;;
esac
NTL_LIBS="${NTL_LIBS} -lntl"
AC_SUBST(NTL_LIBS)
AC_SUBST(NTL_CFLAGS)

dnl Find the FLINT Library
AC_ARG_WITH(flint,
  [AS_HELP_STRING([--with-flint=<path>],
                  [prefix of FLINT installation. e.g. /usr/local or /usr])],
  [],
  [AC_CHECK_LIB([flint], [fmpz_init, nmod_mat_rref_, hmod_mat_rref_],[ ],
    AC_MSG_ERROR([FLINT Library not found.
    Maybe you want to call configure with the --with-flint=<path> option?
    This tells configure where to find the FLINT library and headers.
    e.g. --with-flint=/usr/local or --with-flint=/usr]))]
)
case "$with_flint" in
  ""|yes)
    FLINT_LIBS=""
    FLINT_CFLAGS=""
    ;;
  no)
    AC_MSG_ERROR([Compiling without libflint impossible (but --without-flint given)])
    ;;
  *)
    FLINT_LIBS="${with_flint}/lib/libflint.so"
    FLINT_CFLAGS="-I${with_flint}/include"
    ;;
esac
#FLINT_LIBS="${FLINT_LIBS} -lflint"
AC_SUBST(FLINT_LIBS)
AC_SUBST(FLINT_CFLAGS)

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
  tests/Makefile
  progs/Makefile
  doc/Makefile
  man/Makefile
  Makefile
])
AC_OUTPUT
