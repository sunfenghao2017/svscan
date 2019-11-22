# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([sver], [0.0.0], [wljlinksly@gmail.com])
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])
AC_PREFIX_DEFAULT(/usr/local/)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for typedefs, structures, and compiler characteristics.
AC_SUBST(AM_CXXFLAGS, "-std=c++11 -static -g -Wall -O3 -Wextra -Wno-unknown-pragmas -Wno-unused-function")
AC_SUBST(LDFLAGS, "$LDFLAGS -lhts -pthread -lz -lxlsxwriter -lbwa")
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
