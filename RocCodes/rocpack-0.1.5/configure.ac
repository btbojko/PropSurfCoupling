AC_PREREQ([2.69])
AC_INIT([rocpack], [0.1.5], [guilherme@amadio.org])
AC_CONFIG_AUX_DIR([autotools])
AM_INIT_AUTOMAKE([1.15 foreign subdir-objects -Wall])
AM_SILENT_RULES([yes])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([autotools/m4])

AX_IS_RELEASE([git-directory])
AX_CHECK_ENABLE_DEBUG()

AC_PROG_CC
AC_PROG_CXX
AC_USE_SYSTEM_EXTENSIONS

m4_ifndef([PKG_PROG_PKG_CONFIG],
  [m4_fatal([Could not locate the pkg-config autoconf macros])])

PKG_PROG_PKG_CONFIG

AX_CFLAGS_WARN_ALL
AX_CXXFLAGS_WARN_ALL
AX_CXX_COMPILE_STDCXX_11

AC_LANG(C++)
AC_C_INLINE
AC_CHECK_HEADER_STDBOOL
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([atexit bzero memset])
AC_CHECK_FUNCS([strdup strerror strrchr strtol])
AC_CHECK_FUNCS([sqrt sin cos floor pow])

AC_SEARCH_LIBS([sqrt], [m])

AC_ARG_ENABLE([openmp], [AS_HELP_STRING([--enable-openmp], [enable support for OpenMP])])
AS_IF([test "x$enable_openmp" != "xno"],
[
  AX_OPENMP([AC_SUBST([OPENMP_CXXFLAGS])], [AC_MSG_WARN([Cannot enable support for OpenMP on this system])])
],[test "x$enable_openmp" == "xno"],
[
  dnl Avoid warnings about ignored OpenMP flags when it is disabled
  AX_APPEND_COMPILE_FLAGS([-Wno-unknown-pragmas],[OPENMP_CXXFLAGS])
])

AC_ARG_WITH([opengl], AS_HELP_STRING([--with-opengl], [compile with support for OpenGL]))
AS_IF([test "x$with_opengl" != "xno"],
[
  AX_CHECK_GLU([have_glu=yes], [have_glu=no])
  AX_CHECK_GLUT([have_glut=yes], [have_glut=no])
])
AS_IF([test "x$have_glu" == "xyes" && test "x$have_glut" == "xyes"],
[
  AC_DEFINE([HAVE_OPENGL], [1], [Enable OpenGL rendering])
  AC_SUBST([OPENGL_LIBS],["$GLU_LIBS $GLUT_LIBS"])
])

AC_ARG_WITH([gd], AS_HELP_STRING([--with-gd], [compile with support for gd library]))
AS_IF([test "x$with_png" != "xno"],
  [PKG_CHECK_MODULES(GD, [gdlib], [have_gd=yes], [have_gd=no])])
AS_IF([test "x$have_gd" == "xyes"],
  [AC_DEFINE([HAVE_GD], [1], [Enable gd library support])])

AC_ARG_WITH([tiff], AS_HELP_STRING([--with-tiff], [compile with support for tiff library]))
AS_IF([test "x$with_tiff" != "xno"],
  [PKG_CHECK_MODULES(TIFF, [libtiff-4], [have_tiff=yes], [have_tiff=no])])
AS_IF([test "x$have_tiff" == "xyes"],
  [AC_DEFINE([HAVE_TIFF], [1], [Enable TIFF support])])

AC_ARG_WITH([examples], [AS_HELP_STRING([--with-examples], [install examples])])
AM_CONDITIONAL([EXAMPLES], [test "x$with_examples" = "xyes"])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
