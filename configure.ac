
AC_INIT([isolator],
        m4_esyscmd([./git-version-gen | tr -d '\n']),
        [dcjones@cs.washington.edu])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_MACRO_DIR([m4])

AM_PROG_AR
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX
LT_INIT

# Check for SIMD instruction sets.
AX_EXT

# Check for SIMD compiler intrinsics.
AC_CHECK_HEADERS([immintrin.h])

# Check for isnan
AC_MSG_CHECKING([for isnan])
AC_TRY_LINK([#include <math.h>
], [if (!isnan(3.14159)) isnan(2.7183);], ok=yes, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(HAVE_ISNAN,1,[Define if the isnan() function/macro is available.])
fi
AC_MSG_RESULT(${ok})

# Check for isinf
AC_MSG_CHECKING([for isinf])
AC_TRY_LINK([#include <math.h>
], [if (!isinf(3.14159)) isinf(2.7183);], ok=yes, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(HAVE_ISINF,1,[Define if the isinf() function/macro is available.])
fi
AC_MSG_RESULT(${ok})

# Zlib, for BAM's sake.
AX_CHECK_ZLIB

# Boost
BOOST_REQUIRE([1.48.0])
BOOST_THREADS
BOOST_SYSTEM
BOOST_TIMER

# GSL
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
AX_PATH_GSL

# Configure option to disable compliler optimization.
opt_CFLAGS="-Wall -Wextra -pedantic -g -O3 -DNDEBUG"
dbg_CFLAGS="-Wall -Wextra -pedantic -g -O0"

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging info (default is no)])],
              [], [enable_debug=no])

AS_IF([test "x$enable_debug" = xyes],
      [CFLAGS="$dbg_CFLAGS"],
      [CFLAGS="$opt_CFLAGS"])

CXXFLAGS=$CFLAGS
CFLAGS="$CFLAGS -std=gnu99"

# Don't crash on >4GB files.
AC_DEFINE([_FILE_OFFSET_BITS], [64],
          [Do not crash on >4GB files on 32bit machines.])

# Disable bounds checking in GSL.
AC_DEFINE([GSL_RANGE_CHECK_OFF], [],
          [Make some GSL functions a tiny bit faster.])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/gtf/Makefile
                 src/hat-trie/Makefile
                 src/samtools/Makefile
                 src/sqlite3/Makefile])

AM_CONFIG_HEADER([src/config.h])

AC_OUTPUT
