## configure.ac
## tRNAscan-SE
##
## Settings for generating a configure file 

AC_PREREQ([2.69])
AC_INIT([tRNAscan-SE],[2.0.0],[pchan@soe.ucsc.edu])
AC_CONFIG_HEADERS([config.h])

AC_ARG_ENABLE(debugging,
  AS_HELP_STRING([--enable-debugging],[enable debugging, default: no]),
[
case "${enableval}" in
  yes) enable_debugging=yes ;;
  no)  enable_debugging=no ;;
  *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debugging]) ;;
esac],
[enable_debugging=no])

# Checks for programs.
AC_PROG_CC

if test "$enable_debugging" != "no"; then
  CFLAGS="-g -Wall"
else
  CFLAGS="-O3"
fi

# standard install program
AC_PROG_INSTALL

# Check for perl/get path
AC_PATH_PROG(PERL, perl)

# Check perl version
AC_DEFUN([AX_PROG_PERL_VERSION],
[AC_CACHE_CHECK([for Perl version $1 or later], [ax_cv_prog_perl_version],
    [AS_IF(["$PERL" -e 'require $1;' >/dev/null 2>&1], [$2], [$3])])])
AX_PROG_PERL_VERSION([5.0.0],[AC_MSG_RESULT(yes)],[AC_MSG_ERROR([not compatible with this version of perl], [1])])

# Checks for header files.
AC_FUNC_ALLOCA
AC_CHECK_HEADERS([limits.h memory.h stddef.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([bzero memset re_comp regcomp sqrt strcasecmp strchr strstr])

# use automake
AM_INIT_AUTOMAKE([foreign])

# use Makefile.in to generate Makefile
AC_CONFIG_FILES([Makefile])

AC_OUTPUT
