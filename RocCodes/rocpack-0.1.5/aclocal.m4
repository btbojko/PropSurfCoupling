# generated automatically by aclocal 1.15 -*- Autoconf -*-

# Copyright (C) 1996-2014 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
m4_ifndef([AC_AUTOCONF_VERSION],
  [m4_copy([m4_PACKAGE_VERSION], [AC_AUTOCONF_VERSION])])dnl
m4_if(m4_defn([AC_AUTOCONF_VERSION]), [2.69],,
[m4_warning([this file was generated for autoconf 2.69.
You have another version of autoconf.  It may work, but is not guaranteed to.
If you have problems, you may need to regenerate the build system entirely.
To do so, use the procedure documented by the package, typically 'autoreconf'.])])

# ===========================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_append_compile_flags.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_APPEND_COMPILE_FLAGS([FLAG1 FLAG2 ...], [FLAGS-VARIABLE], [EXTRA-FLAGS], [INPUT])
#
# DESCRIPTION
#
#   For every FLAG1, FLAG2 it is checked whether the compiler works with the
#   flag.  If it does, the flag is added FLAGS-VARIABLE
#
#   If FLAGS-VARIABLE is not specified, the current language's flags (e.g.
#   CFLAGS) is used.  During the check the flag is always added to the
#   current language's flags.
#
#   If EXTRA-FLAGS is defined, it is added to the current language's default
#   flags (e.g. CFLAGS) when the check is done.  The check is thus made with
#   the flags: "CFLAGS EXTRA-FLAGS FLAG".  This can for example be used to
#   force the compiler to issue an error when a bad flag is given.
#
#   INPUT gives an alternative input source to AC_COMPILE_IFELSE.
#
#   NOTE: This macro depends on the AX_APPEND_FLAG and
#   AX_CHECK_COMPILE_FLAG. Please keep this macro in sync with
#   AX_APPEND_LINK_FLAGS.
#
# LICENSE
#
#   Copyright (c) 2011 Maarten Bosmans <mkbosmans@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 5

AC_DEFUN([AX_APPEND_COMPILE_FLAGS],
[AX_REQUIRE_DEFINED([AX_CHECK_COMPILE_FLAG])
AX_REQUIRE_DEFINED([AX_APPEND_FLAG])
for flag in $1; do
  AX_CHECK_COMPILE_FLAG([$flag], [AX_APPEND_FLAG([$flag], [$2])], [], [$3], [$4])
done
])dnl AX_APPEND_COMPILE_FLAGS

# ===========================================================================
#      http://www.gnu.org/software/autoconf-archive/ax_append_flag.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_APPEND_FLAG(FLAG, [FLAGS-VARIABLE])
#
# DESCRIPTION
#
#   FLAG is appended to the FLAGS-VARIABLE shell variable, with a space
#   added in between.
#
#   If FLAGS-VARIABLE is not specified, the current language's flags (e.g.
#   CFLAGS) is used.  FLAGS-VARIABLE is not changed if it already contains
#   FLAG.  If FLAGS-VARIABLE is unset in the shell, it is set to exactly
#   FLAG.
#
#   NOTE: Implementation based on AX_CFLAGS_GCC_OPTION.
#
# LICENSE
#
#   Copyright (c) 2008 Guido U. Draheim <guidod@gmx.de>
#   Copyright (c) 2011 Maarten Bosmans <mkbosmans@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 6

AC_DEFUN([AX_APPEND_FLAG],
[dnl
AC_PREREQ(2.64)dnl for _AC_LANG_PREFIX and AS_VAR_SET_IF
AS_VAR_PUSHDEF([FLAGS], [m4_default($2,_AC_LANG_PREFIX[FLAGS])])
AS_VAR_SET_IF(FLAGS,[
  AS_CASE([" AS_VAR_GET(FLAGS) "],
    [*" $1 "*], [AC_RUN_LOG([: FLAGS already contains $1])],
    [
     AS_VAR_APPEND(FLAGS,[" $1"])
     AC_RUN_LOG([: FLAGS="$FLAGS"])
    ])
  ],
  [
  AS_VAR_SET(FLAGS,[$1])
  AC_RUN_LOG([: FLAGS="$FLAGS"])
  ])
AS_VAR_POPDEF([FLAGS])dnl
])dnl AX_APPEND_FLAG

# ===========================================================================
#    http://www.gnu.org/software/autoconf-archive/ax_cflags_warn_all.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CFLAGS_WARN_ALL   [(shellvar [,default, [A/NA]])]
#   AX_CXXFLAGS_WARN_ALL [(shellvar [,default, [A/NA]])]
#   AX_FCFLAGS_WARN_ALL  [(shellvar [,default, [A/NA]])]
#
# DESCRIPTION
#
#   Try to find a compiler option that enables most reasonable warnings.
#
#   For the GNU compiler it will be -Wall (and -ansi -pedantic) The result
#   is added to the shellvar being CFLAGS, CXXFLAGS, or FCFLAGS by default.
#
#   Currently this macro knows about the GCC, Solaris, Digital Unix, AIX,
#   HP-UX, IRIX, NEC SX-5 (Super-UX 10), Cray J90 (Unicos 10.0.0.8), and
#   Intel compilers.  For a given compiler, the Fortran flags are much more
#   experimental than their C equivalents.
#
#    - $1 shell-variable-to-add-to : CFLAGS, CXXFLAGS, or FCFLAGS
#    - $2 add-value-if-not-found : nothing
#    - $3 action-if-found : add value to shellvariable
#    - $4 action-if-not-found : nothing
#
#   NOTE: These macros depend on AX_APPEND_FLAG.
#
# LICENSE
#
#   Copyright (c) 2008 Guido U. Draheim <guidod@gmx.de>
#   Copyright (c) 2010 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 15

AC_DEFUN([AX_FLAGS_WARN_ALL],[dnl
AS_VAR_PUSHDEF([FLAGS],[_AC_LANG_PREFIX[]FLAGS])dnl
AS_VAR_PUSHDEF([VAR],[ac_cv_[]_AC_LANG_ABBREV[]flags_warn_all])dnl
AC_CACHE_CHECK([m4_ifval($1,$1,FLAGS) for maximum warnings],
VAR,[VAR="no, unknown"
ac_save_[]FLAGS="$[]FLAGS"
for ac_arg dnl
in "-warn all  % -warn all"   dnl Intel
   "-pedantic  % -Wall"       dnl GCC
   "-xstrconst % -v"          dnl Solaris C
   "-std1      % -verbose -w0 -warnprotos" dnl Digital Unix
   "-qlanglvl=ansi % -qsrcmsg -qinfo=all:noppt:noppc:noobs:nocnd" dnl AIX
   "-ansi -ansiE % -fullwarn" dnl IRIX
   "+ESlit     % +w1"         dnl HP-UX C
   "-Xc        % -pvctl[,]fullmsg" dnl NEC SX-5 (Super-UX 10)
   "-h conform % -h msglevel 2" dnl Cray C (Unicos)
   #
do FLAGS="$ac_save_[]FLAGS "`echo $ac_arg | sed -e 's,%%.*,,' -e 's,%,,'`
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
                     [VAR=`echo $ac_arg | sed -e 's,.*% *,,'` ; break])
done
FLAGS="$ac_save_[]FLAGS"
])
AS_VAR_POPDEF([FLAGS])dnl
AX_REQUIRE_DEFINED([AX_APPEND_FLAG])
case ".$VAR" in
     .ok|.ok,*) m4_ifvaln($3,$3) ;;
   .|.no|.no,*) m4_default($4,[m4_ifval($2,[AX_APPEND_FLAG([$2], [$1])])]) ;;
   *) m4_default($3,[AX_APPEND_FLAG([$VAR], [$1])]) ;;
esac
AS_VAR_POPDEF([VAR])dnl
])dnl AX_FLAGS_WARN_ALL
dnl  implementation tactics:
dnl   the for-argument contains a list of options. The first part of
dnl   these does only exist to detect the compiler - usually it is
dnl   a global option to enable -ansi or -extrawarnings. All other
dnl   compilers will fail about it. That was needed since a lot of
dnl   compilers will give false positives for some option-syntax
dnl   like -Woption or -Xoption as they think of it is a pass-through
dnl   to later compile stages or something. The "%" is used as a
dnl   delimiter. A non-option comment can be given after "%%" marks
dnl   which will be shown but not added to the respective C/CXXFLAGS.

AC_DEFUN([AX_CFLAGS_WARN_ALL],[dnl
AC_LANG_PUSH([C])
AX_FLAGS_WARN_ALL([$1], [$2], [$3], [$4])
AC_LANG_POP([C])
])

AC_DEFUN([AX_CXXFLAGS_WARN_ALL],[dnl
AC_LANG_PUSH([C++])
AX_FLAGS_WARN_ALL([$1], [$2], [$3], [$4])
AC_LANG_POP([C++])
])

AC_DEFUN([AX_FCFLAGS_WARN_ALL],[dnl
AC_LANG_PUSH([Fortran])
AX_FLAGS_WARN_ALL([$1], [$2], [$3], [$4])
AC_LANG_POP([Fortran])
])

# ===========================================================================
#   http://www.gnu.org/software/autoconf-archive/ax_check_compile_flag.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_COMPILE_FLAG(FLAG, [ACTION-SUCCESS], [ACTION-FAILURE], [EXTRA-FLAGS], [INPUT])
#
# DESCRIPTION
#
#   Check whether the given FLAG works with the current language's compiler
#   or gives an error.  (Warnings, however, are ignored)
#
#   ACTION-SUCCESS/ACTION-FAILURE are shell commands to execute on
#   success/failure.
#
#   If EXTRA-FLAGS is defined, it is added to the current language's default
#   flags (e.g. CFLAGS) when the check is done.  The check is thus made with
#   the flags: "CFLAGS EXTRA-FLAGS FLAG".  This can for example be used to
#   force the compiler to issue an error when a bad flag is given.
#
#   INPUT gives an alternative input source to AC_COMPILE_IFELSE.
#
#   NOTE: Implementation based on AX_CFLAGS_GCC_OPTION. Please keep this
#   macro in sync with AX_CHECK_{PREPROC,LINK}_FLAG.
#
# LICENSE
#
#   Copyright (c) 2008 Guido U. Draheim <guidod@gmx.de>
#   Copyright (c) 2011 Maarten Bosmans <mkbosmans@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 4

AC_DEFUN([AX_CHECK_COMPILE_FLAG],
[AC_PREREQ(2.64)dnl for _AC_LANG_PREFIX and AS_VAR_IF
AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether _AC_LANG compiler accepts $1], CACHEVAR, [
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
  AC_COMPILE_IFELSE([m4_default([$5],[AC_LANG_PROGRAM()])],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])
  _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_VAR_IF(CACHEVAR,yes,
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_COMPILE_FLAGS

# ===========================================================================
#   http://www.gnu.org/software/autoconf-archive/ax_check_enable_debug.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_ENABLE_DEBUG([enable by default=yes/info/profile/no], [ENABLE DEBUG VARIABLES ...], [DISABLE DEBUG VARIABLES NDEBUG ...], [IS-RELEASE])
#
# DESCRIPTION
#
#   Check for the presence of an --enable-debug option to configure, with
#   the specified default value used when the option is not present.  Return
#   the value in the variable $ax_enable_debug.
#
#   Specifying 'yes' adds '-g -O0' to the compilation flags for all
#   languages. Specifying 'info' adds '-g' to the compilation flags.
#   Specifying 'profile' adds '-g -pg' to the compilation flags and '-pg' to
#   the linking flags. Otherwise, nothing is added.
#
#   Define the variables listed in the second argument if debug is enabled,
#   defaulting to no variables.  Defines the variables listed in the third
#   argument if debug is disabled, defaulting to NDEBUG.  All lists of
#   variables should be space-separated.
#
#   If debug is not enabled, ensure AC_PROG_* will not add debugging flags.
#   Should be invoked prior to any AC_PROG_* compiler checks.
#
#   IS-RELEASE can be used to change the default to 'no' when making a
#   release.  Set IS-RELEASE to 'yes' or 'no' as appropriate. By default, it
#   uses the value of $ax_is_release, so if you are using the AX_IS_RELEASE
#   macro, there is no need to pass this parameter.
#
#     AX_IS_RELEASE([git-directory])
#     AX_CHECK_ENABLE_DEBUG()
#
# LICENSE
#
#   Copyright (c) 2011 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2014, 2015 Philip Withnall <philip@tecnocode.co.uk>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

#serial 5

AC_DEFUN([AX_CHECK_ENABLE_DEBUG],[
    AC_BEFORE([$0],[AC_PROG_CC])dnl
    AC_BEFORE([$0],[AC_PROG_CXX])dnl
    AC_BEFORE([$0],[AC_PROG_F77])dnl
    AC_BEFORE([$0],[AC_PROG_FC])dnl

    AC_MSG_CHECKING(whether to enable debugging)

    ax_enable_debug_default=m4_tolower(m4_normalize(ifelse([$1],,[no],[$1])))
    ax_enable_debug_is_release=m4_tolower(m4_normalize(ifelse([$4],,
                                                              [$ax_is_release],
                                                              [$4])))

    # If this is a release, override the default.
    AS_IF([test "$ax_enable_debug_is_release" = "yes"],
      [ax_enable_debug_default="no"])

    m4_define(ax_enable_debug_vars,[m4_normalize(ifelse([$2],,,[$2]))])
    m4_define(ax_disable_debug_vars,[m4_normalize(ifelse([$3],,[NDEBUG],[$3]))])

    AC_ARG_ENABLE(debug,
	[AS_HELP_STRING([--enable-debug=]@<:@yes/info/profile/no@:>@,[compile with debugging])],
	[],enable_debug=$ax_enable_debug_default)

    # empty mean debug yes
    AS_IF([test "x$enable_debug" = "x"],
      [enable_debug="yes"])

    # case of debug
    AS_CASE([$enable_debug],
      [yes],[
	AC_MSG_RESULT(yes)
	CFLAGS="${CFLAGS} -g -O0"
	CXXFLAGS="${CXXFLAGS} -g -O0"
	FFLAGS="${FFLAGS} -g -O0"
	FCFLAGS="${FCFLAGS} -g -O0"
	OBJCFLAGS="${OBJCFLAGS} -g -O0"
      ],
      [info],[
	AC_MSG_RESULT(info)
	CFLAGS="${CFLAGS} -g"
	CXXFLAGS="${CXXFLAGS} -g"
	FFLAGS="${FFLAGS} -g"
	FCFLAGS="${FCFLAGS} -g"
	OBJCFLAGS="${OBJCFLAGS} -g"
      ],
      [profile],[
	AC_MSG_RESULT(profile)
	CFLAGS="${CFLAGS} -g -pg"
	CXXFLAGS="${CXXFLAGS} -g -pg"
	FFLAGS="${FFLAGS} -g -pg"
	FCFLAGS="${FCFLAGS} -g -pg"
	OBJCFLAGS="${OBJCFLAGS} -g -pg"
	LDFLAGS="${LDFLAGS} -pg"
      ],
      [
	AC_MSG_RESULT(no)
	dnl Ensure AC_PROG_CC/CXX/F77/FC/OBJC will not enable debug flags
	dnl by setting any unset environment flag variables
	AS_IF([test "x${CFLAGS+set}" != "xset"],
	  [CFLAGS=""])
	AS_IF([test "x${CXXFLAGS+set}" != "xset"],
	  [CXXFLAGS=""])
	AS_IF([test "x${FFLAGS+set}" != "xset"],
	  [FFLAGS=""])
	AS_IF([test "x${FCFLAGS+set}" != "xset"],
	  [FCFLAGS=""])
	AS_IF([test "x${OBJCFLAGS+set}" != "xset"],
	  [OBJCFLAGS=""])
      ])

    dnl Define various variables if debugging is disabled.
    dnl assert.h is a NOP if NDEBUG is defined, so define it by default.
    AS_IF([test "x$enable_debug" = "xyes"],
      [m4_map_args_w(ax_enable_debug_vars, [AC_DEFINE(], [,,[Define if debugging is enabled])])],
      [m4_map_args_w(ax_disable_debug_vars, [AC_DEFINE(], [,,[Define if debugging is disabled])])])
    ax_enable_debug=$enable_debug
])

# ===========================================================================
#        http://www.gnu.org/software/autoconf-archive/ax_check_gl.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_GL([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Checks for an OpenGL implementation. If a valid OpenGL implementation is
#   found, this macro would set C preprocessor symbol HAVE_GL to 1.
#
#   If either a valid OpenGL header or library was not found, by default the
#   configuration would exits on error. This behavior can be overwritten by
#   providing a custom "ACTION-IF-NOT-FOUND" hook.
#
#   If the header, library was found, and been tested for compiling and
#   linking the configuration would export the required compiler flags to
#   "GL_CFLAGS" and "GL_LIBS". These two variables can also be overwritten
#   by user from the command line if they want to link against the library
#   they specified instead of having the configuration script to detect the
#   flags automatically. Note that having "GL_CFLAGS" or "GL_LIBS" set
#   doesn't mean it can compile or link with the flags, since it could be
#   overwritten by user. However the "HAVE_GL" symbol and "ACTION-IF-FOUND"
#   hook is always guaranteed to reflect a valid OpenGL implementation.
#
#   If user didn't specify the "ACTION-IF-FOUND" hook, the configuration
#   would prepend "GL_CFLAGS" and "GL_LIBS" to "CFLAGS" and "LIBS", like
#   many other autoconf macros do.
#
#   OpenGL is one of the libraries that has different header names on
#   different platforms. This macro does the header detection, and will
#   export the following symbol: "HAVE_GL_GL_H" for having "GL/gl.h" or
#   "HAVE_OPENGL_GL_H" for having "OpenGL/gl.h". To write a portable OpenGL
#   code, you should include OpenGL header like so:
#
#     #if defined(HAVE_WINDOWS_H) && defined(_WIN32)
#     # include <windows.h>
#     #endif
#     #ifdef HAVE_GL_GL_H
#     # include <GL/gl.h>
#     #elif defined(HAVE_OPENGL_GL_H)
#     # include <OpenGL/gl.h>
#     #else
#     # error no gl.h
#     #endif
#
#   On the OSX platform, there's two possible OpenGL implementation. One is
#   the OpenGL that ships with OSX, the other comes with X11/XQuartz
#   (http://www.xquartz.org). To use the xquartz variant, user can use the
#   option --with-xquartz-gl[=path to xquartz root]. By default the
#   configuration will check "/opt/X11", which is the default X11 install
#   location on OSX.
#
# LICENSE
#
#   Copyright (c) 2009 Braden McDaniel <braden@endoframe.com>
#   Copyright (c) 2012 Bastien Roucaries <roucaries.bastien+autoconf@gmail.com>
#   Copyright (c) 2016 Felix Chern <idryman@gmail.com>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 19

# example gl program
m4_define([_AX_CHECK_GL_PROGRAM],
          [AC_LANG_PROGRAM([[
# if defined(HAVE_WINDOWS_H) && defined(_WIN32)
#   include <windows.h>
# endif
# ifdef HAVE_GL_GL_H
#   include <GL/gl.h>
# elif defined(HAVE_OPENGL_GL_H)
#   include <OpenGL/gl.h>
# else
#   error no gl.h
# endif
]],[[glBegin(0)]])])

dnl Default include : add windows.h
dnl see http://www.opengl.org/wiki/Platform_specifics:_Windows
dnl (acceded 20120801)
AC_DEFUN([_AX_CHECK_GL_INCLUDES_DEFAULT],dnl
[
  AC_INCLUDES_DEFAULT
  [
  # if defined(HAVE_WINDOWS_H) && defined(_WIN32)
  #   include <windows.h>
  # endif
  ]
])


# _AX_CHECK_GL_SAVE_FLAGS(LIST-OF-FLAGS)
# Use this macro before you modify the flags.
# Restore the flags by _AX_CHECK_GL_RESTORE_FLAGS
#
# Example: _AX_CHECK_GL_SAVE_FLAGS([[CFLAGS],[LIBS]]) expands to
# gl_saved_flag_cflags=$CFLAGS
# gl_saved_flag_libs=$LIBS
# CFLAGS="$GL_CFLAGS $CFLAGS"
# LIBS="$GL_LIBS $LIBS"
AC_DEFUN([_AX_CHECK_GL_SAVE_FLAGS], [
 AX_SAVE_FLAGS_WITH_PREFIX([GL],[$1])
 AC_LANG_PUSH([C])
])

# _AX_CHECK_GL_RESTORE_FLAGS(LIST-OF-FLAGS)
# Use this marcro to restore the flags you saved using
# _AX_CHECK_GL_SAVE_FLAGS
#
# Example: _AX_CHECK_GL_RESTORE_FLAGS([[CFLAGS],[LIBS]]) expands to
# CFLAGS="$gl_saved_flag_cflags"
# LIBS="$gl_saved_flag_libs"
AC_DEFUN([_AX_CHECK_GL_RESTORE_FLAGS], [
 AX_RESTORE_FLAGS_WITH_PREFIX([GL],[$1])
 AC_LANG_POP([C])
])

# Check if the program compiles
AC_DEFUN([_AX_CHECK_GL_COMPILE],
[dnl
 _AX_CHECK_GL_SAVE_FLAGS([CFLAGS])
 AC_COMPILE_IFELSE([_AX_CHECK_GL_PROGRAM],
                   [ax_check_gl_compile_opengl="yes"],
                   [ax_check_gl_compile_opengl="no"])
 _AX_CHECK_GL_RESTORE_FLAGS([CFLAGS])
])

# Compile the example program (cache)
AC_DEFUN([_AX_CHECK_GL_COMPILE_CV],
[dnl
 AC_CACHE_CHECK([for compiling a minimal OpenGL program],[ax_cv_check_gl_compile_opengl],
                [_AX_CHECK_GL_COMPILE()
                 ax_cv_check_gl_compile_opengl="${ax_check_gl_compile_opengl}"])
 ax_check_gl_compile_opengl="${ax_cv_check_gl_compile_opengl}"
])

# Link the example program
AC_DEFUN([_AX_CHECK_GL_LINK],
[dnl
 _AX_CHECK_GL_SAVE_FLAGS([[CFLAGS],[LIBS]])
 AC_LINK_IFELSE([_AX_CHECK_GL_PROGRAM],
                [ax_check_gl_link_opengl="yes"],
                [ax_check_gl_link_opengl="no"])
 _AX_CHECK_GL_RESTORE_FLAGS([[CFLAGS],[LIBS],[LDFLAGS]])
])

# Link the example program (cache)
AC_DEFUN([_AX_CHECK_GL_LINK_CV],
[dnl
 AC_CACHE_CHECK([for linking a minimal OpenGL program],[ax_cv_check_gl_link_opengl],
                [_AX_CHECK_GL_LINK()
                 ax_cv_check_gl_link_opengl="${ax_check_gl_link_opengl}"])
 ax_check_gl_link_opengl="${ax_cv_check_gl_link_opengl}"
])


# _AX_CHECK_GL_MANUAL_LIBS_GENERIC(LIBRARIES-TO-SEARCH)
# Searches library provided in $1, and output the flag
# $ax_check_gl_lib_opengl
AC_DEFUN([_AX_CHECK_GL_MANUAL_LIBS_GENERIC], [
  AS_IF([test -n "$GL_LIBS"],[], [
    ax_check_gl_manual_libs_generic_extra_libs="$1"
    AS_IF([test "X$ax_check_gl_manual_libs_generic_extra_libs" = "X"],
          [AC_MSG_ERROR([AX_CHECK_GL_MANUAL_LIBS_GENERIC argument must no be empty])])

    _AX_CHECK_GL_SAVE_FLAGS([CFLAGS])
    AC_SEARCH_LIBS([glBegin],[$ax_check_gl_manual_libs_generic_extra_libs], [
                   ax_check_gl_lib_opengl="yes"
                   break
                   ])
    AS_IF([test "X$ax_check_gl_lib_opengl"="Xyes"],
          [GL_LIBS="${ac_cv_search_glBegin}"])
    _AX_CHECK_GL_RESTORE_FLAGS([CFLAGS])
 ])
])

# _WITH_XQUARTZ_GL
# ----------------
# Provides an option in command line to specify the XQuartz installation
# path on OSX, so that user can link to it instead of using the default
# OSX OpenGL framework. (Mac OSX only)
AC_DEFUN_ONCE([_WITH_XQUARTZ_GL],[
  AC_ARG_WITH([xquartz-gl],
   [AS_HELP_STRING([--with-xquartz-gl@<:@=DIR@:>@],
                   [On Mac OSX, use opengl provided by X11/XQuartz instead of the built-in framework.
                    If enabled, the default location is @<:@DIR=/opt/X11@:>@.
                    This option is default to false.])],
   [AS_IF([test "X$with_xquartz_gl"="Xyes"],
          [with_xquartz_gl="/opt/X11"])],
   [with_xquartz_gl=no])
  AS_IF([test "X$with_xquartz_gl" != "Xno"],
        [AC_MSG_CHECKING([OSX X11 path])
         AS_IF([test -e "$with_xquartz_gl"],
               [AC_MSG_RESULT(["$with_xquartz_gl"])
                CFLAGS="-I$with_xquartz_gl/include $CFLAGS"
                LIBS="-L$with_xquartz_gl/lib $LIBS"
               ],
               [with_xquartz_gl=no
                AC_MSG_RESULT([no])
                AC_MSG_WARN([--with-xquartz-gl was given, but test for X11 failed. Fallback to system framework])
               ])
        ])
])

# OSX specific setup for OpenGL check
AC_DEFUN([_AX_CHECK_DARWIN_GL], [
 AC_REQUIRE([_WITH_XQUARTZ_GL])
 AS_IF([test "x$with_xquartz_gl" != "xno"],
       [GL_LIBS="${GL_LIBS:--lGL}"],
       [GL_LIBS="${GL_LIBS:--framework OpenGL}"])
])


# AX_CHECK_GL_LIB([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------
# Checks OpenGL headers and library and provides hooks for success and failures.
# When $1 is not set, this macro would modify CFLAGS and LIBS environment variables.
# However, user can override this behavior by providing their own hooks.
# The CFLAGS and LIBS flags required by OpenGL is always exported in
# GL_CFLAGS and GL_LIBS environment variable.
#
# In other words, the default behavior of AX_CHECK_GL_LIB() is equivalent to
# AX_CHECK_GL_LIB(
#   [CFLAGS="$GL_CFLAGS $CFLAGS"
#    LIBS="$GL_LIBS $LIBS"]
# )
AC_DEFUN([AX_CHECK_GL],
[AC_REQUIRE([AC_CANONICAL_HOST])
 AC_REQUIRE([PKG_PROG_PKG_CONFIG])
 AC_ARG_VAR([GL_CFLAGS],[C compiler flags for GL, overriding configure script defaults])
 AC_ARG_VAR([GL_LIBS],[Linker flags for GL, overriding configure script defaults])

 dnl --with-gl or not can be implemented outside of check-gl
 AS_CASE([${host}],
         [*-darwin*],[_AX_CHECK_DARWIN_GL],
         dnl some windows may support X11 opengl, and should be able to linked
         dnl by -lGL. However I have no machine to test it.
         [*-cygwin*|*-mingw*],[
          _AX_CHECK_GL_MANUAL_LIBS_GENERIC([opengl32 GL gl])
          AC_CHECK_HEADERS([windows.h])
          ],
         [PKG_PROG_PKG_CONFIG
          PKG_CHECK_MODULES([GL],[gl],
          [],
          [_AX_CHECK_GL_MANUAL_LIBS_GENERIC([GL gl])])
         ]) dnl host specific checks

 dnl this was cache
 _AX_CHECK_GL_SAVE_FLAGS([CFLAGS])
 AC_CHECK_HEADERS([GL/gl.h OpenGL/gl.h],
   [ax_check_gl_have_headers="yes";break])
 _AX_CHECK_GL_RESTORE_FLAGS([CFLAGS])

 AS_IF([test "X$ax_check_gl_have_headers" = "Xyes"],
       [_AX_CHECK_GL_COMPILE_CV()],
       [no_gl=yes])
 AS_IF([test "X$ax_check_gl_compile_opengl" = "Xyes"],
       [_AX_CHECK_GL_LINK_CV()],
       [no_gl=yes])
 AS_IF([test "X$no_gl" = "X"],
   [AC_DEFINE([HAVE_GL], [1], [Defined if a valid OpenGL implementation is found.])
    m4_ifval([$1],
      [$1],
      [CFLAGS="$GL_CFLAGS $CFLAGS"
       LIBS="$GL_LIBS $LIBS"])
   ],
   [m4_ifval([$2],
     [$2],
     [AC_MSG_ERROR([Could not find a valid OpenGL implementation])])
   ])
])

# ===========================================================================
#       http://www.gnu.org/software/autoconf-archive/ax_check_glu.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_GLU([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Checks for GLUT. If a valid GLU implementation is found, the configure
#   script would export the C preprocessor symbol "HAVE_GLU=1".
#
#   If either a valid GLU header or library was not found, by default the
#   configure script would exit on error. This behavior can be overwritten
#   by providing a custom "ACTION-IF-NOT-FOUND" hook.
#
#   If the header, library was found, and been tested for compiling and
#   linking the configuration would export the required compiler flags to
#   "GLU_CFLAGS" and "GLU_LIBS" environment variables. These two variables
#   can also be overwritten by defining the environment variables before
#   executing the configure program. If it was predefined, configure would
#   not try to overwrite it, but it would still perform the compile and link
#   test. Only when the tests succeeded does the configure script to export
#   "HAVE_GLU=1" and to run "ACTION-IF-FOUND" hook.
#
#   If user didn't specify the "ACTION-IF-FOUND" hook, the configuration
#   would prepend "GLU_CFLAGS" and "GLU_LIBS" to "CFLAGS" and "LIBS", like
#   many other autoconf macros do.
#
#   If the header "GL/glu.h" is found, "HAVE_GL_GLU_H" is defined. If the
#   header "OpenGL/glu.h" is found, HAVE_OPENGL_GLU_H is defined.
#
#   You should use something like this in your headers:
#
#     # if defined(HAVE_WINDOWS_H) && defined(_WIN32)
#     #  include <windows.h>
#     # endif
#     # if defined(HAVE_GL_GLU_H)
#     #  include <GL/glu.h>
#     # elif defined(HAVE_OPENGL_GLU_H)
#     #  include <OpenGL/glu.h>
#     # else
#     #  error no glu.h
#     # endif
#
#   On the OSX platform, you can use the option --with-xquartz-gl to use
#   X11/Xquartz GLU implementation instead of the system built in GLU
#   framework.
#
#   Some implementations (in particular, some versions of Mac OS X) are
#   known to treat the GLU tesselator callback function type as "GLvoid
#   (*)(...)" rather than the standard "GLvoid (*)()". If the former
#   condition is detected, this macro defines "HAVE_VARARGS_GLU_TESSCB".
#
# LICENSE
#
#   Copyright (c) 2009 Braden McDaniel <braden@endoframe.com>
#   Copyright (c) 2013 Bastien Roucaries <roucaries.bastien+autoconf@gmail.com>
#   Copyright (c) 2016 Felix Chern <idryman@gmail.com>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 20

# exemple program
m4_define([_AX_CHECK_GLU_PROGRAM],
          [AC_LANG_PROGRAM([[
# if defined(HAVE_WINDOWS_H) && defined(_WIN32)
#   include <windows.h>
# endif
# ifdef HAVE_GL_GLU_H
#   include <GL/glu.h>
# elif defined(HAVE_OPENGL_GLU_H)
#   include <OpenGL/glu.h>
# else
#   error no glu.h
# endif
]],[[gluBeginCurve(0)]])])


dnl Default include : add windows.h
dnl see http://www.opengl.org/wiki/Platform_specifics:_Windows
dnl (acceded 20120801)
AC_DEFUN([_AX_CHECK_GLU_INCLUDES_DEFAULT],dnl
[
  AC_INCLUDES_DEFAULT
  [
  # if defined(HAVE_WINDOWS_H) && defined(_WIN32)
  #   include <windows.h>
  # endif
  ]
])

# check tesselation callback function signature.
m4_define([_AX_CHECK_GLU_VARARGS_TESSVB_PROGRAM],
[AC_LANG_PROGRAM([[
# if defined(HAVE_WINDOWS_H) && defined(_WIN32)
#   include <windows.h>
# endif
# ifdef HAVE_GL_GLU_H
#   include <GL/glu.h>
# elif defined(HAVE_OPENGL_GLU_H)
#   include <OpenGL/glu.h>
# else
#   error no glu.h
# endif
]],
[[GLvoid (*func)(...); gluTessCallback(0, 0, func)]])
])


# _AX_CHECK_GLU_SAVE_FLAGS(LIST-OF-FLAGS,[LANG])
# ----------------------------------------------
# Save the flags to shell variables.
# Example: _AX_CHECK_GLU_SAVE_FLAGS([[CFLAGS],[LIBS]]) expands to
# AC_LANG_PUSH([C])
# glu_saved_flag_cflags=$CFLAGS
# glu_saved_flag_libs=$LIBS
# CFLAGS="$GLU_CFLAGS $CFLAGS"
# LIBS="$GLU_LIBS $LIBS"
#
# Can optionally support other LANG by specifying $2
AC_DEFUN([_AX_CHECK_GLU_SAVE_FLAGS], [
 m4_ifval([$2],
          [AC_LANG_PUSH([$2])],
          [AC_LANG_PUSH([C])])
 AX_SAVE_FLAGS_WITH_PREFIX([GLU],[$1]) dnl defined in ax_check_gl
])

# _AX_CHECK_GLU_RESTORE_FLAGS(LIST-OF-FLAGS)
# Use this marcro to restore the flags you saved using
# _AX_CHECK_GLU_SAVE_FLAGS
#
# Example: _AX_CHECK_GLU_RESTORE_FLAGS([[CFLAGS],[LIBS]]) expands to
# CFLAGS="$glu_saved_flag_cflags"
# LIBS="$glu_saved_flag_libs"
# AC_LANG_POP([C])
AC_DEFUN([_AX_CHECK_GLU_RESTORE_FLAGS], [
 AX_RESTORE_FLAGS_WITH_PREFIX([GLU],[$1]) dnl defined in ax_check_gl
 m4_ifval([$2],
          [AC_LANG_POP([$2])],
          [AC_LANG_POP([C])])
])


# Search headers and export $ax_check_glu_have_headers
AC_DEFUN([_AX_CHECK_GLU_HEADERS], [
  _AX_CHECK_GLU_SAVE_FLAGS([CFLAGS])
  AC_CHECK_HEADERS([$1],
                   [ax_check_glu_have_headers="yes";],
                   [],
                   [_AX_CHECK_GLU_INCLUDES_DEFAULT()])
  _AX_CHECK_GLU_RESTORE_FLAGS([CFLAGS])
])


# _AX_CHECK_GLU_SEARCH_LIBS(LIBS)
# -------------------------------
# Search for a valid GLU lib from $1 and set
# GLU_LIBS respectively
AC_DEFUN([_AX_CHECK_GLU_SEARCH_LIBS], [
 _AX_CHECK_GLU_SAVE_FLAGS([[CFLAGS],[LIBS]])
 AC_SEARCH_LIBS([gluBeginCurve],[$1],
 	        [GLU_LIBS="${GLU_LIBS:-$ac_cv_search_gluBeginCurve}"])
  _AX_CHECK_GLU_RESTORE_FLAGS([[CFLAGS],[LIBS]])
])

# OSX specific GLU checks
AC_DEFUN([_AX_CHECK_DARWIN_GLU], [
  AC_REQUIRE([_WITH_XQUARTZ_GL])
  AS_IF([test "x$with_xquartz_gl" != "xno"],
        [GLU_LIBS="${GLU_LIBS:--lGLU}"],
        [GLU_LIBS="${GLU_LIBS:--framework OpenGL}"])
])

# AX_CHECK_GLU([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
# -----------------------------------------------------
# Checks GLU and provides hooks for success and failures
AC_DEFUN([AX_CHECK_GLU],[
  AC_REQUIRE([AC_CANONICAL_HOST])
  AC_REQUIRE([_WITH_XQUARTZ_GL])
  AC_REQUIRE([PKG_PROG_PKG_CONFIG])
  AC_ARG_VAR([GLU_CFLAGS],[C compiler flags for GLU, overriding system check])
  AC_ARG_VAR([GLU_LIBS],[Linker flags for GLU, overriding system check])

  dnl Setup GLU_CFLAGS and GLU_LIBS
  AS_CASE([${host}],
          [*-darwin*],[_AX_CHECK_DARWIN_GLU],
          [*-cygwin*],[_AX_CHECK_GLU_SEARCH_LIBS([GLU glu MesaGLU glu32])
                       AC_CHECK_HEADERS([windows.h])],
          # try first native
 	  [*-mingw*],[_AX_CHECK_GLU_SEARCH_LIBS([glu32 GLU glu MesaGLU])
                      AC_CHECK_HEADERS([windows.h])],
          [PKG_PROG_PKG_CONFIG
           PKG_CHECK_MODULES([GLU],[glu],
           [],
           [_AX_CHECK_GLU_SEARCH_LIBS([GLU glu MesaGLU])])
          ])

  AS_CASE([$host],
          [*-darwin*],
            [AS_IF([test "X$with_xquartz_gl" = "Xno"],
                   [_AX_CHECK_GLU_HEADERS([OpenGL/glu.h])],
                   [_AX_CHECK_GLU_HEADERS([GL/glu.h])]
                   )],
          [_AX_CHECK_GLU_HEADERS([GL/glu.h])])

  dnl compile test
  AS_IF([test "X$ax_check_glu_have_headers" = "Xyes"],
        [AC_CACHE_CHECK([for compiling a minimal OpenGL Utility (GLU) program],
                        [ax_cv_check_glu_compile],
                        [_AX_CHECK_GLU_SAVE_FLAGS([CFLAGS])
                         AC_COMPILE_IFELSE([_AX_CHECK_GLU_PROGRAM],
                                           [ax_cv_check_glu_compile="yes"],
                                           [ax_cv_check_glu_compile="no"])
                         _AX_CHECK_GLU_RESTORE_FLAGS([CFLAGS])])
         ])

  dnl link test
  AS_IF([test "X$ax_cv_check_glu_compile" = "Xyes"],
        [AC_CACHE_CHECK([for linking a minimal GLU program],
                        [ax_cv_check_glu_link],
                        [_AX_CHECK_GLU_SAVE_FLAGS([[CFLAGS],[LIBS]])
                         AC_LINK_IFELSE([_AX_CHECK_GLU_PROGRAM],
                                        [ax_cv_check_glu_link="yes"],
                                        [ax_cv_check_glu_link="no"])
                         _AX_CHECK_GLU_RESTORE_FLAGS([[CFLAGS],[LIBS]])])
        ])

#
# Some versions of Mac OS X include a broken interpretation of the GLU
# tesselation callback function signature.
  AS_IF([test "X$ax_cv_check_glu_link" = "Xyes"],
        [AC_CACHE_CHECK([if GLU varargs tesselator is using non-standard form],
                        [ax_cv_varargs_glu_tesscb],
                        [_AX_CHECK_GLU_SAVE_FLAGS([CFLAGS],[C++])
                         AC_COMPILE_IFELSE([_AX_CHECK_GLU_VARARGS_TESSVB_PROGRAM],
                                           [ax_cv_varargs_glu_tesscb="yes"],
                                           [ax_cv_varargs_glu_tesscb="no"])
                         _AX_CHECK_GLU_RESTORE_FLAGS([CFLAGS],[C++])])
        AS_IF([test "X$ax_cv_varargs_glu_tesscb" = "yes"],
              [AC_DEFINE([HAVE_VARARGS_GLU_TESSCB], [1],
                         [Use nonstandard varargs form for the GLU tesselator callback])])
        ])

  dnl hook
  AS_IF([test "X$ax_cv_check_glu_link" = "Xyes"],
        [AC_DEFINE([HAVE_GLU],[1],[Defined if a valid GLU implementation is found.])
         m4_ifval([$1],
                  [$1],
                  [CFLAGS="$GLU_CFLAGS $CFLAGS"
                   LIBS="$GLU_LIBS $LIBS"])],
        [m4_ifval([$2],
                  [$2],
                  [AC_MSG_ERROR([Could not find a valid GLU implementation])])
        ])
])

# ===========================================================================
#       http://www.gnu.org/software/autoconf-archive/ax_check_glut.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_GLUT([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Checks for GLUT. If a valid GLUT implementation is found, the configure
#   script would export the C preprocessor symbol "HAVE_GLUT=1".
#
#   If either a valid GLUT header or library was not found, by default the
#   configure script would exit on error. This behavior can be overwritten
#   by providing a custom "ACTION-IF-NOT-FOUND" hook.
#
#   If the header, library was found, and been tested for compiling and
#   linking the configuration would export the required compiler flags to
#   "GLUT_CFLAGS" and "GLUT_LIBS" environment variables. These two variables
#   can also be overwritten by defining the environment variables before
#   executing the configure program. If it was predefined, configure would
#   not try to overwrite it, but it would still perform the compile and link
#   test. Only when the tests succeeded does the configure script to export
#   "HAVE_GLUT=1" and to run "ACTION-IF-FOUND" hook.
#
#   If user didn't specify the "ACTION-IF-FOUND" hook, the configuration
#   would prepend "GLUT_CFLAGS" and "GLUT_LIBS" to "CFLAGS" and "LIBS", like
#   many other autoconf macros do.
#
#   If the header "GL/glut.h" is found, "HAVE_GL_GLUT_H" is defined. If the
#   header "GLUT/glut.h" is found, HAVE_GLUT_GLUT_H is defined.
#
#   You should use something like this in your headers:
#
#     # if HAVE_WINDOWS_H && defined(_WIN32)
#     #  include <windows.h>
#     # endif
#     # if defined(HAVE_GL_GLUT_H)
#     #  include <GL/glut.h>
#     # elif defined(HAVE_GLUT_GLUT_H)
#     #  include <GLUT/glut.h>
#     # else
#     #  error no glut.h
#     # endif
#
#   On the OSX platform, you can use the option --with-xquartz-gl to use
#   X11/Xquartz GLUT implementation instead of the system built in GLUT
#   framework.
#
# LICENSE
#
#   Copyright (c) 2009 Braden McDaniel <braden@endoframe.com>
#   Copyright (c) 2013 Bastien Roucaries <roucaries.bastien+autoconf@gmail.com>
#   Copyright (c) 2016 Felix Chern <idryman@gmail.com>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 16

AC_DEFUN([_AX_CHECK_GLUT_SAVE_FLAGS], [
 AX_SAVE_FLAGS_WITH_PREFIX([GLUT],[$1]) dnl defined in ax_check_gl
 AC_LANG_PUSH([C])
])

AC_DEFUN([_AX_CHECK_GLUT_RESTORE_FLAGS], [
 AX_RESTORE_FLAGS_WITH_PREFIX([GLUT],[$1]) dnl defined in ax_check_gl
 AC_LANG_POP([C])
])

dnl Default include : add windows.h
dnl see http://www.opengl.org/wiki/Platform_specifics:_Windows
dnl (acceded 20120801)
AC_DEFUN([_AX_CHECK_GLUT_INCLUDES_DEFAULT],dnl
[
  AC_INCLUDES_DEFAULT
  [
  # if defined(HAVE_WINDOWS_H) && defined(_WIN32)
  #   include <windows.h>
  # endif
  ]
])

m4_define([_AX_CHECK_GLUT_PROGRAM],
          [AC_LANG_PROGRAM([[
# if HAVE_WINDOWS_H && defined(_WIN32)
#   include <windows.h>
# endif
# ifdef HAVE_GL_GLUT_H
#   include <GL/glut.h>
# elif defined(HAVE_GLUT_GLUT_H)
#   include <GLUT/glut.h>
# else
#   error no glut.h
# endif]],
[[glutMainLoop()]])])


# _AX_CHECK_GLUT_MANUAL_LIBS_GENERIC(LIST-OF-LIBS)
# ------------------------------------------------
# Searches libraries provided in $1, and export variable
# $ax_check_glut_lib_glut
AC_DEFUN([_AX_CHECK_GLUT_MANUAL_LIBS_GENERIC],
[
 _AX_CHECK_GLUT_SAVE_FLAGS([[CFLAGS],[LIBS]])
 AC_SEARCH_LIBS([glutMainLoop],[$1],
                [GLUT_LIBS="${GLUT_LIBS:-$ac_cv_search_glutMainLoop}"])
 _AX_CHECK_GLUT_RESTORE_FLAGS([[CFLAGS],[LIBS]])
])

# Wrapper macro to check GLUT header
AC_DEFUN([_AX_CHECK_GLUT_HEADER],[
  _AX_CHECK_GLUT_SAVE_FLAGS([CFLAGS])
  AC_CHECK_HEADERS([$1],
                   [ax_check_glut_have_headers=yes])
  _AX_CHECK_GLUT_RESTORE_FLAGS([CFLAGS])
])


# AX_CHECK_GLUT_LIB([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------
# Checks GLUT headers and library and provides hooks for success and failures.
AC_DEFUN([AX_CHECK_GLUT],
[AC_REQUIRE([AC_CANONICAL_HOST])
 AC_REQUIRE([_WITH_XQUARTZ_GL])
 AC_ARG_VAR([GLUT_CFLAGS],[C compiler flags for GLUT, overriding configure script defaults])
 AC_ARG_VAR([GLUT_LIBS],[Linker flags for GLUT, overriding configure script defaults])

 AS_CASE([${host}],
         [*-darwin*],[AS_IF([test "x$with_xquartz_gl" != "xno"],
                            [GLUT_LIBS="${GLUT_LIBS:--lGLUT}"],
                            [GLUT_LIBS="${GLUT_LIBS:--framework GLUT}"])],
         [*-cygwin*|*-mingw*],[
            _AX_CHECK_GLUT_MANUAL_LIBS_GENERIC([glut32 glut])
            AC_CHECK_HEADERS([windows.h])
          ],
         [_AX_CHECK_GLUT_MANUAL_LIBS_GENERIC([glut])
         ]) dnl host specific checks

 dnl checks header
 AS_CASE([${host}],
   [*-darwin*],[AS_IF([test "x$with_xquartz_gl" = "xno"],
                      [_AX_CHECK_GLUT_HEADER([GLUT/glut.h])],
                      [_AX_CHECK_GLUT_HEADER([GL/glut.h])]
                      )],
   [_AX_CHECK_GLUT_HEADER([GL/glut.h])])

 dnl compile
 AS_IF([test "X$ax_check_glut_have_headers" = "Xyes"],
       [AC_CACHE_CHECK([for compiling a minimal GLUT program],
                       [ax_cv_check_glut_compile],
                       [_AX_CHECK_GLUT_SAVE_FLAGS([CFLAGS])
                        AC_COMPILE_IFELSE([_AX_CHECK_GLUT_PROGRAM],
                                          [ax_cv_check_glut_compile="yes"],
                                          [ax_cv_check_glut_compile="no"])
                        _AX_CHECK_GLUT_RESTORE_FLAGS([CFLAGS])
                       ])
      ])

 dnl link
 AS_IF([test "X$ax_cv_check_glut_compile" = "Xyes"],
       [AC_CACHE_CHECK([for linking a minimal GLUT program],
                       [ax_cv_check_glut_link],
                       [_AX_CHECK_GLUT_SAVE_FLAGS([[CFLAGS],[LIBS]])
                        AC_LINK_IFELSE([_AX_CHECK_GLUT_PROGRAM],
                                       [ax_cv_check_glut_link="yes"],
                                       [ax_cv_check_glut_link="no"])
                        _AX_CHECK_GLUT_RESTORE_FLAGS([[CFLAGS],[LIBS]])
                       ])
       ])

 dnl hook
 AS_IF([test "X$ax_cv_check_glut_link" = "Xyes"],
   [AC_DEFINE([HAVE_GLUT], [1], [Defined if a valid GLUT implementation is found])
    m4_ifval([$1],
     [$1],
     [CFLAGS="$GLUT_CFLAGS $CFLAGS"
      LIBS="$GLUT_LIBS $LIBS"])
   ],
   [m4_ifval([$2],
     [$2],
     [AC_MSG_ERROR([Could not find a valid GLUT implementation])]
     )
   ])

])

# ===========================================================================
#   http://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX(VERSION, [ext|noext], [mandatory|optional])
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the specified
#   version of the C++ standard.  If necessary, add switches to CXX and
#   CXXCPP to enable support.  VERSION may be '11' (for the C++11 standard)
#   or '14' (for the C++14 standard).
#
#   The second argument, if specified, indicates whether you insist on an
#   extended mode (e.g. -std=gnu++11) or a strict conformance mode (e.g.
#   -std=c++11).  If neither is specified, you get whatever works, with
#   preference for an extended mode.
#
#   The third argument, if specified 'mandatory' or if left unspecified,
#   indicates that baseline support for the specified C++ standard is
#   required and that the macro should error out if no mode with that
#   support is found.  If specified 'optional', then configuration proceeds
#   regardless, after defining HAVE_CXX${VERSION} if and only if a
#   supporting mode is found.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#   Copyright (c) 2012 Zack Weinberg <zackw@panix.com>
#   Copyright (c) 2013 Roy Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2014, 2015 Google Inc.; contributed by Alexey Sokolov <sokolov@google.com>
#   Copyright (c) 2015 Paul Norman <penorman@mac.com>
#   Copyright (c) 2015 Moritz Klammler <moritz@klammler.eu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.

#serial 4

dnl  This macro is based on the code from the AX_CXX_COMPILE_STDCXX_11 macro
dnl  (serial version number 13).

AC_DEFUN([AX_CXX_COMPILE_STDCXX], [dnl
  m4_if([$1], [11], [],
        [$1], [14], [],
        [$1], [17], [m4_fatal([support for C++17 not yet implemented in AX_CXX_COMPILE_STDCXX])],
        [m4_fatal([invalid first argument `$1' to AX_CXX_COMPILE_STDCXX])])dnl
  m4_if([$2], [], [],
        [$2], [ext], [],
        [$2], [noext], [],
        [m4_fatal([invalid second argument `$2' to AX_CXX_COMPILE_STDCXX])])dnl
  m4_if([$3], [], [ax_cxx_compile_cxx$1_required=true],
        [$3], [mandatory], [ax_cxx_compile_cxx$1_required=true],
        [$3], [optional], [ax_cxx_compile_cxx$1_required=false],
        [m4_fatal([invalid third argument `$3' to AX_CXX_COMPILE_STDCXX])])
  AC_LANG_PUSH([C++])dnl
  ac_success=no
  AC_CACHE_CHECK(whether $CXX supports C++$1 features by default,
  ax_cv_cxx_compile_cxx$1,
  [AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_testbody_$1])],
    [ax_cv_cxx_compile_cxx$1=yes],
    [ax_cv_cxx_compile_cxx$1=no])])
  if test x$ax_cv_cxx_compile_cxx$1 = xyes; then
    ac_success=yes
  fi

  m4_if([$2], [noext], [], [dnl
  if test x$ac_success = xno; then
    for switch in -std=gnu++$1 -std=gnu++0x; do
      cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx$1_$switch])
      AC_CACHE_CHECK(whether $CXX supports C++$1 features with $switch,
                     $cachevar,
        [ac_save_CXX="$CXX"
         CXX="$CXX $switch"
         AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_testbody_$1])],
          [eval $cachevar=yes],
          [eval $cachevar=no])
         CXX="$ac_save_CXX"])
      if eval test x\$$cachevar = xyes; then
        CXX="$CXX $switch"
        if test -n "$CXXCPP" ; then
          CXXCPP="$CXXCPP $switch"
        fi
        ac_success=yes
        break
      fi
    done
  fi])

  m4_if([$2], [ext], [], [dnl
  if test x$ac_success = xno; then
    dnl HP's aCC needs +std=c++11 according to:
    dnl http://h21007.www2.hp.com/portal/download/files/unprot/aCxx/PDF_Release_Notes/769149-001.pdf
    dnl Cray's crayCC needs "-h std=c++11"
    for switch in -std=c++$1 -std=c++0x +std=c++$1 "-h std=c++$1"; do
      cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx$1_$switch])
      AC_CACHE_CHECK(whether $CXX supports C++$1 features with $switch,
                     $cachevar,
        [ac_save_CXX="$CXX"
         CXX="$CXX $switch"
         AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_testbody_$1])],
          [eval $cachevar=yes],
          [eval $cachevar=no])
         CXX="$ac_save_CXX"])
      if eval test x\$$cachevar = xyes; then
        CXX="$CXX $switch"
        if test -n "$CXXCPP" ; then
          CXXCPP="$CXXCPP $switch"
        fi
        ac_success=yes
        break
      fi
    done
  fi])
  AC_LANG_POP([C++])
  if test x$ax_cxx_compile_cxx$1_required = xtrue; then
    if test x$ac_success = xno; then
      AC_MSG_ERROR([*** A compiler with support for C++$1 language features is required.])
    fi
  fi
  if test x$ac_success = xno; then
    HAVE_CXX$1=0
    AC_MSG_NOTICE([No compiler with C++$1 support was found])
  else
    HAVE_CXX$1=1
    AC_DEFINE(HAVE_CXX$1,1,
              [define if the compiler supports basic C++$1 syntax])
  fi
  AC_SUBST(HAVE_CXX$1)
])


dnl  Test body for checking C++11 support

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_11],
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_11
)


dnl  Test body for checking C++14 support

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_14],
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_11
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_14
)


dnl  Tests for new features in C++11

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_new_in_11], [[

// If the compiler admits that it is not ready for C++11, why torture it?
// Hopefully, this will speed up the test.

#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 201103L

#error "This is not a C++11 compiler"

#else

namespace cxx11
{

  namespace test_static_assert
  {

    template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

  }

  namespace test_final_override
  {

    struct Base
    {
      virtual void f() {}
    };

    struct Derived : public Base
    {
      virtual void f() override {}
    };

  }

  namespace test_double_right_angle_brackets
  {

    template < typename T >
    struct check {};

    typedef check<void> single_type;
    typedef check<check<void>> double_type;
    typedef check<check<check<void>>> triple_type;
    typedef check<check<check<check<void>>>> quadruple_type;

  }

  namespace test_decltype
  {

    int
    f()
    {
      int a = 1;
      decltype(a) b = 2;
      return a + b;
    }

  }

  namespace test_type_deduction
  {

    template < typename T1, typename T2 >
    struct is_same
    {
      static const bool value = false;
    };

    template < typename T >
    struct is_same<T, T>
    {
      static const bool value = true;
    };

    template < typename T1, typename T2 >
    auto
    add(T1 a1, T2 a2) -> decltype(a1 + a2)
    {
      return a1 + a2;
    }

    int
    test(const int c, volatile int v)
    {
      static_assert(is_same<int, decltype(0)>::value == true, "");
      static_assert(is_same<int, decltype(c)>::value == false, "");
      static_assert(is_same<int, decltype(v)>::value == false, "");
      auto ac = c;
      auto av = v;
      auto sumi = ac + av + 'x';
      auto sumf = ac + av + 1.0;
      static_assert(is_same<int, decltype(ac)>::value == true, "");
      static_assert(is_same<int, decltype(av)>::value == true, "");
      static_assert(is_same<int, decltype(sumi)>::value == true, "");
      static_assert(is_same<int, decltype(sumf)>::value == false, "");
      static_assert(is_same<int, decltype(add(c, v))>::value == true, "");
      return (sumf > 0.0) ? sumi : add(c, v);
    }

  }

  namespace test_noexcept
  {

    int f() { return 0; }
    int g() noexcept { return 0; }

    static_assert(noexcept(f()) == false, "");
    static_assert(noexcept(g()) == true, "");

  }

  namespace test_constexpr
  {

    template < typename CharT >
    unsigned long constexpr
    strlen_c_r(const CharT *const s, const unsigned long acc) noexcept
    {
      return *s ? strlen_c_r(s + 1, acc + 1) : acc;
    }

    template < typename CharT >
    unsigned long constexpr
    strlen_c(const CharT *const s) noexcept
    {
      return strlen_c_r(s, 0UL);
    }

    static_assert(strlen_c("") == 0UL, "");
    static_assert(strlen_c("1") == 1UL, "");
    static_assert(strlen_c("example") == 7UL, "");
    static_assert(strlen_c("another\0example") == 7UL, "");

  }

  namespace test_rvalue_references
  {

    template < int N >
    struct answer
    {
      static constexpr int value = N;
    };

    answer<1> f(int&)       { return answer<1>(); }
    answer<2> f(const int&) { return answer<2>(); }
    answer<3> f(int&&)      { return answer<3>(); }

    void
    test()
    {
      int i = 0;
      const int c = 0;
      static_assert(decltype(f(i))::value == 1, "");
      static_assert(decltype(f(c))::value == 2, "");
      static_assert(decltype(f(0))::value == 3, "");
    }

  }

  namespace test_uniform_initialization
  {

    struct test
    {
      static const int zero {};
      static const int one {1};
    };

    static_assert(test::zero == 0, "");
    static_assert(test::one == 1, "");

  }

  namespace test_lambdas
  {

    void
    test1()
    {
      auto lambda1 = [](){};
      auto lambda2 = lambda1;
      lambda1();
      lambda2();
    }

    int
    test2()
    {
      auto a = [](int i, int j){ return i + j; }(1, 2);
      auto b = []() -> int { return '0'; }();
      auto c = [=](){ return a + b; }();
      auto d = [&](){ return c; }();
      auto e = [a, &b](int x) mutable {
        const auto identity = [](int y){ return y; };
        for (auto i = 0; i < a; ++i)
          a += b--;
        return x + identity(a + b);
      }(0);
      return a + b + c + d + e;
    }

    int
    test3()
    {
      const auto nullary = [](){ return 0; };
      const auto unary = [](int x){ return x; };
      using nullary_t = decltype(nullary);
      using unary_t = decltype(unary);
      const auto higher1st = [](nullary_t f){ return f(); };
      const auto higher2nd = [unary](nullary_t f1){
        return [unary, f1](unary_t f2){ return f2(unary(f1())); };
      };
      return higher1st(nullary) + higher2nd(nullary)(unary);
    }

  }

  namespace test_variadic_templates
  {

    template <int...>
    struct sum;

    template <int N0, int... N1toN>
    struct sum<N0, N1toN...>
    {
      static constexpr auto value = N0 + sum<N1toN...>::value;
    };

    template <>
    struct sum<>
    {
      static constexpr auto value = 0;
    };

    static_assert(sum<>::value == 0, "");
    static_assert(sum<1>::value == 1, "");
    static_assert(sum<23>::value == 23, "");
    static_assert(sum<1, 2>::value == 3, "");
    static_assert(sum<5, 5, 11>::value == 21, "");
    static_assert(sum<2, 3, 5, 7, 11, 13>::value == 41, "");

  }

  // http://stackoverflow.com/questions/13728184/template-aliases-and-sfinae
  // Clang 3.1 fails with headers of libstd++ 4.8.3 when using std::function
  // because of this.
  namespace test_template_alias_sfinae
  {

    struct foo {};

    template<typename T>
    using member = typename T::member_type;

    template<typename T>
    void func(...) {}

    template<typename T>
    void func(member<T>*) {}

    void test();

    void test() { func<foo>(0); }

  }

}  // namespace cxx11

#endif  // __cplusplus >= 201103L

]])


dnl  Tests for new features in C++14

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_new_in_14], [[

// If the compiler admits that it is not ready for C++14, why torture it?
// Hopefully, this will speed up the test.

#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 201402L

#error "This is not a C++14 compiler"

#else

namespace cxx14
{

  namespace test_polymorphic_lambdas
  {

    int
    test()
    {
      const auto lambda = [](auto&&... args){
        const auto istiny = [](auto x){
          return (sizeof(x) == 1UL) ? 1 : 0;
        };
        const int aretiny[] = { istiny(args)... };
        return aretiny[0];
      };
      return lambda(1, 1L, 1.0f, '1');
    }

  }

  namespace test_binary_literals
  {

    constexpr auto ivii = 0b0000000000101010;
    static_assert(ivii == 42, "wrong value");

  }

  namespace test_generalized_constexpr
  {

    template < typename CharT >
    constexpr unsigned long
    strlen_c(const CharT *const s) noexcept
    {
      auto length = 0UL;
      for (auto p = s; *p; ++p)
        ++length;
      return length;
    }

    static_assert(strlen_c("") == 0UL, "");
    static_assert(strlen_c("x") == 1UL, "");
    static_assert(strlen_c("test") == 4UL, "");
    static_assert(strlen_c("another\0test") == 7UL, "");

  }

  namespace test_lambda_init_capture
  {

    int
    test()
    {
      auto x = 0;
      const auto lambda1 = [a = x](int b){ return a + b; };
      const auto lambda2 = [a = lambda1(x)](){ return a; };
      return lambda2();
    }

  }

  namespace test_digit_seperators
  {

    constexpr auto ten_million = 100'000'000;
    static_assert(ten_million == 100000000, "");

  }

  namespace test_return_type_deduction
  {

    auto f(int& x) { return x; }
    decltype(auto) g(int& x) { return x; }

    template < typename T1, typename T2 >
    struct is_same
    {
      static constexpr auto value = false;
    };

    template < typename T >
    struct is_same<T, T>
    {
      static constexpr auto value = true;
    };

    int
    test()
    {
      auto x = 0;
      static_assert(is_same<int, decltype(f(x))>::value, "");
      static_assert(is_same<int&, decltype(g(x))>::value, "");
      return x;
    }

  }

}  // namespace cxx14

#endif  // __cplusplus >= 201402L

]])

# ============================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx_11.html
# ============================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX_11([ext|noext], [mandatory|optional])
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the C++11
#   standard; if necessary, add switches to CXX and CXXCPP to enable
#   support.
#
#   This macro is a convenience alias for calling the AX_CXX_COMPILE_STDCXX
#   macro with the version set to C++11.  The two optional arguments are
#   forwarded literally as the second and third argument respectively.
#   Please see the documentation for the AX_CXX_COMPILE_STDCXX macro for
#   more information.  If you want to use this macro, you also need to
#   download the ax_cxx_compile_stdcxx.m4 file.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#   Copyright (c) 2012 Zack Weinberg <zackw@panix.com>
#   Copyright (c) 2013 Roy Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2014, 2015 Google Inc.; contributed by Alexey Sokolov <sokolov@google.com>
#   Copyright (c) 2015 Paul Norman <penorman@mac.com>
#   Copyright (c) 2015 Moritz Klammler <moritz@klammler.eu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 17

AX_REQUIRE_DEFINED([AX_CXX_COMPILE_STDCXX])
AC_DEFUN([AX_CXX_COMPILE_STDCXX_11], [AX_CXX_COMPILE_STDCXX([11], [$1], [$2])])

# ===========================================================================
#       http://www.gnu.org/software/autoconf-archive/ax_is_release.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_IS_RELEASE(POLICY)
#
# DESCRIPTION
#
#   Determine whether the code is being configured as a release, or from
#   git. Set the ax_is_release variable to 'yes' or 'no'.
#
#   If building a release version, it is recommended that the configure
#   script disable compiler errors and debug features, by conditionalising
#   them on the ax_is_release variable.  If building from git, these
#   features should be enabled.
#
#   The POLICY parameter specifies how ax_is_release is determined. It can
#   take the following values:
#
#    * git-directory:  ax_is_release will be 'no' if a '.git' directory exists
#    * minor-version:  ax_is_release will be 'no' if the minor version number
#                      in $PACKAGE_VERSION is odd; this assumes
#                      $PACKAGE_VERSION follows the 'major.minor.micro' scheme
#    * micro-version:  ax_is_release will be 'no' if the micro version number
#                      in $PACKAGE_VERSION is odd; this assumes
#                      $PACKAGE_VERSION follows the 'major.minor.micro' scheme
#    * always:         ax_is_release will always be 'yes'
#    * never:          ax_is_release will always be 'no'
#
#   Other policies may be added in future.
#
# LICENSE
#
#   Copyright (c) 2015 Philip Withnall <philip@tecnocode.co.uk>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

#serial 3

AC_DEFUN([AX_IS_RELEASE],[
    AC_BEFORE([AC_INIT],[$0])

    m4_case([$1],
      [git-directory],[
        # $is_release = (.git directory does not exist)
        AS_IF([test -d .git],[ax_is_release=no],[ax_is_release=yes])
      ],
      [minor-version],[
        # $is_release = ($minor_version is even)
        minor_version=`echo "$PACKAGE_VERSION" | sed 's/[[^.]][[^.]]*.\([[^.]][[^.]]*\).*/\1/'`
        AS_IF([test "$(( $minor_version % 2 ))" -ne 0],
              [ax_is_release=no],[ax_is_release=yes])
      ],
      [micro-version],[
        # $is_release = ($micro_version is even)
        micro_version=`echo "$PACKAGE_VERSION" | sed 's/[[^.]]*\.[[^.]]*\.\([[^.]]*\).*/\1/'`
        AS_IF([test "$(( $micro_version % 2 ))" -ne 0],
              [ax_is_release=no],[ax_is_release=yes])
      ],
      [always],[ax_is_release=yes],
      [never],[ax_is_release=no],
      [
        AC_MSG_ERROR([Invalid policy. Valid policies: git-directory, minor-version.])
      ])
])

# ===========================================================================
#         http://www.gnu.org/software/autoconf-archive/ax_openmp.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_OPENMP([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use OpenMP a
#   standard API and set of compiler directives for parallel programming
#   (see http://www-unix.mcs/)
#
#   On success, it sets the OPENMP_CFLAGS/OPENMP_CXXFLAGS/OPENMP_F77FLAGS
#   output variable to the flag (e.g. -omp) used both to compile *and* link
#   OpenMP programs in the current language.
#
#   NOTE: You are assumed to not only compile your program with these flags,
#   but also link it with them as well.
#
#   If you want to compile everything with OpenMP, you should set:
#
#     CFLAGS="$CFLAGS $OPENMP_CFLAGS"
#     #OR#  CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
#     #OR#  FFLAGS="$FFLAGS $OPENMP_FFLAGS"
#
#   (depending on the selected language).
#
#   The user can override the default choice by setting the corresponding
#   environment variable (e.g. OPENMP_CFLAGS).
#
#   ACTION-IF-FOUND is a list of shell commands to run if an OpenMP flag is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_OPENMP.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2015 John W. Peterson <jwpeterson@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 11

AC_DEFUN([AX_OPENMP], [
AC_PREREQ([2.69]) dnl for _AC_LANG_PREFIX

AC_CACHE_CHECK([for OpenMP flag of _AC_LANG compiler], ax_cv_[]_AC_LANG_ABBREV[]_openmp, [save[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
ax_cv_[]_AC_LANG_ABBREV[]_openmp=unknown
# Flags to try:  -fopenmp (gcc), -openmp (icc), -mp (SGI & PGI),
#                -xopenmp (Sun), -omp (Tru64), -qsmp=omp (AIX), none
ax_openmp_flags="-fopenmp -openmp -mp -xopenmp -omp -qsmp=omp none"
if test "x$OPENMP_[]_AC_LANG_PREFIX[]FLAGS" != x; then
  ax_openmp_flags="$OPENMP_[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flags"
fi
for ax_openmp_flag in $ax_openmp_flags; do
  case $ax_openmp_flag in
    none) []_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[] ;;
    *) []_AC_LANG_PREFIX[]FLAGS="$save[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flag" ;;
  esac
  AC_LINK_IFELSE([AC_LANG_SOURCE([[
@%:@include <omp.h>

static void
parallel_fill(int * data, int n)
{
  int i;
@%:@pragma omp parallel for
  for (i = 0; i < n; ++i)
    data[i] = i;
}

int
main()
{
  int arr[100000];
  omp_set_num_threads(2);
  parallel_fill(arr, 100000);
  return 0;
}
]])],[ax_cv_[]_AC_LANG_ABBREV[]_openmp=$ax_openmp_flag; break],[])
done
[]_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[]FLAGS
])
if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" = "xunknown"; then
  m4_default([$2],:)
else
  if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" != "xnone"; then
    OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ax_cv_[]_AC_LANG_ABBREV[]_openmp
  fi
  m4_default([$1], [AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])
fi
])dnl AX_OPENMP

# ===========================================================================
#    http://www.gnu.org/software/autoconf-archive/ax_require_defined.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_REQUIRE_DEFINED(MACRO)
#
# DESCRIPTION
#
#   AX_REQUIRE_DEFINED is a simple helper for making sure other macros have
#   been defined and thus are available for use.  This avoids random issues
#   where a macro isn't expanded.  Instead the configure script emits a
#   non-fatal:
#
#     ./configure: line 1673: AX_CFLAGS_WARN_ALL: command not found
#
#   It's like AC_REQUIRE except it doesn't expand the required macro.
#
#   Here's an example:
#
#     AX_REQUIRE_DEFINED([AX_CHECK_LINK_FLAG])
#
# LICENSE
#
#   Copyright (c) 2014 Mike Frysinger <vapier@gentoo.org>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 1

AC_DEFUN([AX_REQUIRE_DEFINED], [dnl
  m4_ifndef([$1], [m4_fatal([macro ]$1[ is not defined; is a m4 file missing?])])
])dnl AX_REQUIRE_DEFINED

# ================================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_restore_flags_with_prefix.html
# ================================================================================
#
# SYNOPSIS
#
#   AX_RESTORE_FLAGS_WITH_PREFIX(PREFIX, LIST-OF-FLAGS)
#
# DESCRIPTION
#
#   Restore the flags saved by AX_SAVE_FLAGS_WITH_PREFIX.
#
#   Expansion example: AX_RESTORE_FLAGS_WITH_PREFIX([GL], [[CFLAGS],[LIBS]])
#   expands to
#
#     CFLAGS="$gl_saved_flag_cflags"
#     LIBS="$gl_saved_flag_libs"
#
#   One common use case is to define a package specific wrapper macro around
#   this one, and also restore other variables if needed. For example:
#
#     AC_DEFUN([_AX_CHECK_GL_RESTORE_FLAGS], [
#       AX_RESTORE_FLAGS_WITH_PREFIX([GL],[$1])
#       AC_LANG_POP([C])
#     ])
#
#     # Restores CFLAGS, LIBS and language state
#     _AX_CHECK_GL_RESTORE_FLAGS([[CFLAGS],[LIBS]])
#
# LICENSE
#
#   Copyright (c) 2016 Felix Chern <idryman@gmail.com>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 2

AC_DEFUN([AX_RESTORE_FLAGS_WITH_PREFIX],[
m4_ifval([$2], [
m4_car($2)="$_ax_[]m4_tolower($1)_saved_flag_[]m4_tolower(m4_car($2))"
$0($1, m4_cdr($2))])
])

# =============================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_save_flags_with_prefix.html
# =============================================================================
#
# SYNOPSIS
#
#   AX_SAVE_FLAGS_WITH_PREFIX(PREFIX, LIST-OF-FLAGS)
#
# DESCRIPTION
#
#   For each flag in LIST-OF-FLAGS, it expands to lower-cased shell variable
#   with the prefix holding the flag original value.  The saved variables
#   can be restored by AX_RESTORE_FLAGS_WITH_PREFIX
#
#   As an example: AX_SAVE_FLAGS_WITH_PREFIX([GL], [[CFLAGS],[LIBS]])
#   expands to
#
#     gl_saved_flag_cflags="$CFLAGS"
#     gl_saved_flag_libs="$LIBS"
#     CFLAGS="$GL_CFLAGS $CFLAGS"
#     LIBS="$GL_LIBS $LIBS"
#
#   One common use case is to define a package specific wrapper macro around
#   this one, and also setup other variables if needed. For example:
#
#     AC_DEFUN([_AX_CHECK_GL_SAVE_FLAGS], [
#       AX_SAVE_FLAGS_WITH_PREFIX([GL],[$1])
#       AC_LANG_PUSH([C])
#     ])
#
#     # pushes GL_CFLAGS and GL_LIBS to CFLAGS and LIBS
#     # also set the current language to test to C
#     _AX_CHECK_GL_SAVE_FLAGS([[CFLAGS],[LIBS]])
#
# LICENSE
#
#   Copyright (c) 2016 Felix Chern <idryman@gmail.com>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 2

AC_DEFUN([AX_SAVE_FLAGS_WITH_PREFIX],[
m4_ifval([$2], [
_ax_[]m4_tolower($1)_saved_flag_[]m4_tolower(m4_car($2))="$m4_car($2)"
m4_car($2)="$$1_[]m4_car($2) $m4_car($2)"
$0($1, m4_cdr($2))
])])

dnl pkg.m4 - Macros to locate and utilise pkg-config.   -*- Autoconf -*-
dnl serial 11 (pkg-config-0.29.1)
dnl
dnl Copyright © 2004 Scott James Remnant <scott@netsplit.com>.
dnl Copyright © 2012-2015 Dan Nicholson <dbn.lists@gmail.com>
dnl
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl
dnl This program is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
dnl General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception to the GNU General Public License, if you
dnl distribute this file as part of a program that contains a
dnl configuration script generated by Autoconf, you may include it under
dnl the same distribution terms that you use for the rest of that
dnl program.

dnl PKG_PREREQ(MIN-VERSION)
dnl -----------------------
dnl Since: 0.29
dnl
dnl Verify that the version of the pkg-config macros are at least
dnl MIN-VERSION. Unlike PKG_PROG_PKG_CONFIG, which checks the user's
dnl installed version of pkg-config, this checks the developer's version
dnl of pkg.m4 when generating configure.
dnl
dnl To ensure that this macro is defined, also add:
dnl m4_ifndef([PKG_PREREQ],
dnl     [m4_fatal([must install pkg-config 0.29 or later before running autoconf/autogen])])
dnl
dnl See the "Since" comment for each macro you use to see what version
dnl of the macros you require.
m4_defun([PKG_PREREQ],
[m4_define([PKG_MACROS_VERSION], [0.29.1])
m4_if(m4_version_compare(PKG_MACROS_VERSION, [$1]), -1,
    [m4_fatal([pkg.m4 version $1 or higher is required but ]PKG_MACROS_VERSION[ found])])
])dnl PKG_PREREQ

dnl PKG_PROG_PKG_CONFIG([MIN-VERSION])
dnl ----------------------------------
dnl Since: 0.16
dnl
dnl Search for the pkg-config tool and set the PKG_CONFIG variable to
dnl first found in the path. Checks that the version of pkg-config found
dnl is at least MIN-VERSION. If MIN-VERSION is not specified, 0.9.0 is
dnl used since that's the first version where most current features of
dnl pkg-config existed.
AC_DEFUN([PKG_PROG_PKG_CONFIG],
[m4_pattern_forbid([^_?PKG_[A-Z_]+$])
m4_pattern_allow([^PKG_CONFIG(_(PATH|LIBDIR|SYSROOT_DIR|ALLOW_SYSTEM_(CFLAGS|LIBS)))?$])
m4_pattern_allow([^PKG_CONFIG_(DISABLE_UNINSTALLED|TOP_BUILD_DIR|DEBUG_SPEW)$])
AC_ARG_VAR([PKG_CONFIG], [path to pkg-config utility])
AC_ARG_VAR([PKG_CONFIG_PATH], [directories to add to pkg-config's search path])
AC_ARG_VAR([PKG_CONFIG_LIBDIR], [path overriding pkg-config's built-in search path])

if test "x$ac_cv_env_PKG_CONFIG_set" != "xset"; then
	AC_PATH_TOOL([PKG_CONFIG], [pkg-config])
fi
if test -n "$PKG_CONFIG"; then
	_pkg_min_version=m4_default([$1], [0.9.0])
	AC_MSG_CHECKING([pkg-config is at least version $_pkg_min_version])
	if $PKG_CONFIG --atleast-pkgconfig-version $_pkg_min_version; then
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
		PKG_CONFIG=""
	fi
fi[]dnl
])dnl PKG_PROG_PKG_CONFIG

dnl PKG_CHECK_EXISTS(MODULES, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl -------------------------------------------------------------------
dnl Since: 0.18
dnl
dnl Check to see whether a particular set of modules exists. Similar to
dnl PKG_CHECK_MODULES(), but does not set variables or print errors.
dnl
dnl Please remember that m4 expands AC_REQUIRE([PKG_PROG_PKG_CONFIG])
dnl only at the first occurence in configure.ac, so if the first place
dnl it's called might be skipped (such as if it is within an "if", you
dnl have to call PKG_CHECK_EXISTS manually
AC_DEFUN([PKG_CHECK_EXISTS],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
if test -n "$PKG_CONFIG" && \
    AC_RUN_LOG([$PKG_CONFIG --exists --print-errors "$1"]); then
  m4_default([$2], [:])
m4_ifvaln([$3], [else
  $3])dnl
fi])

dnl _PKG_CONFIG([VARIABLE], [COMMAND], [MODULES])
dnl ---------------------------------------------
dnl Internal wrapper calling pkg-config via PKG_CONFIG and setting
dnl pkg_failed based on the result.
m4_define([_PKG_CONFIG],
[if test -n "$$1"; then
    pkg_cv_[]$1="$$1"
 elif test -n "$PKG_CONFIG"; then
    PKG_CHECK_EXISTS([$3],
                     [pkg_cv_[]$1=`$PKG_CONFIG --[]$2 "$3" 2>/dev/null`
		      test "x$?" != "x0" && pkg_failed=yes ],
		     [pkg_failed=yes])
 else
    pkg_failed=untried
fi[]dnl
])dnl _PKG_CONFIG

dnl _PKG_SHORT_ERRORS_SUPPORTED
dnl ---------------------------
dnl Internal check to see if pkg-config supports short errors.
AC_DEFUN([_PKG_SHORT_ERRORS_SUPPORTED],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])
if $PKG_CONFIG --atleast-pkgconfig-version 0.20; then
        _pkg_short_errors_supported=yes
else
        _pkg_short_errors_supported=no
fi[]dnl
])dnl _PKG_SHORT_ERRORS_SUPPORTED


dnl PKG_CHECK_MODULES(VARIABLE-PREFIX, MODULES, [ACTION-IF-FOUND],
dnl   [ACTION-IF-NOT-FOUND])
dnl --------------------------------------------------------------
dnl Since: 0.4.0
dnl
dnl Note that if there is a possibility the first call to
dnl PKG_CHECK_MODULES might not happen, you should be sure to include an
dnl explicit call to PKG_PROG_PKG_CONFIG in your configure.ac
AC_DEFUN([PKG_CHECK_MODULES],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
AC_ARG_VAR([$1][_CFLAGS], [C compiler flags for $1, overriding pkg-config])dnl
AC_ARG_VAR([$1][_LIBS], [linker flags for $1, overriding pkg-config])dnl

pkg_failed=no
AC_MSG_CHECKING([for $1])

_PKG_CONFIG([$1][_CFLAGS], [cflags], [$2])
_PKG_CONFIG([$1][_LIBS], [libs], [$2])

m4_define([_PKG_TEXT], [Alternatively, you may set the environment variables $1[]_CFLAGS
and $1[]_LIBS to avoid the need to call pkg-config.
See the pkg-config man page for more details.])

if test $pkg_failed = yes; then
   	AC_MSG_RESULT([no])
        _PKG_SHORT_ERRORS_SUPPORTED
        if test $_pkg_short_errors_supported = yes; then
	        $1[]_PKG_ERRORS=`$PKG_CONFIG --short-errors --print-errors --cflags --libs "$2" 2>&1`
        else 
	        $1[]_PKG_ERRORS=`$PKG_CONFIG --print-errors --cflags --libs "$2" 2>&1`
        fi
	# Put the nasty error message in config.log where it belongs
	echo "$$1[]_PKG_ERRORS" >&AS_MESSAGE_LOG_FD

	m4_default([$4], [AC_MSG_ERROR(
[Package requirements ($2) were not met:

$$1_PKG_ERRORS

Consider adjusting the PKG_CONFIG_PATH environment variable if you
installed software in a non-standard prefix.

_PKG_TEXT])[]dnl
        ])
elif test $pkg_failed = untried; then
     	AC_MSG_RESULT([no])
	m4_default([$4], [AC_MSG_FAILURE(
[The pkg-config script could not be found or is too old.  Make sure it
is in your PATH or set the PKG_CONFIG environment variable to the full
path to pkg-config.

_PKG_TEXT

To get pkg-config, see <http://pkg-config.freedesktop.org/>.])[]dnl
        ])
else
	$1[]_CFLAGS=$pkg_cv_[]$1[]_CFLAGS
	$1[]_LIBS=$pkg_cv_[]$1[]_LIBS
        AC_MSG_RESULT([yes])
	$3
fi[]dnl
])dnl PKG_CHECK_MODULES


dnl PKG_CHECK_MODULES_STATIC(VARIABLE-PREFIX, MODULES, [ACTION-IF-FOUND],
dnl   [ACTION-IF-NOT-FOUND])
dnl ---------------------------------------------------------------------
dnl Since: 0.29
dnl
dnl Checks for existence of MODULES and gathers its build flags with
dnl static libraries enabled. Sets VARIABLE-PREFIX_CFLAGS from --cflags
dnl and VARIABLE-PREFIX_LIBS from --libs.
dnl
dnl Note that if there is a possibility the first call to
dnl PKG_CHECK_MODULES_STATIC might not happen, you should be sure to
dnl include an explicit call to PKG_PROG_PKG_CONFIG in your
dnl configure.ac.
AC_DEFUN([PKG_CHECK_MODULES_STATIC],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
_save_PKG_CONFIG=$PKG_CONFIG
PKG_CONFIG="$PKG_CONFIG --static"
PKG_CHECK_MODULES($@)
PKG_CONFIG=$_save_PKG_CONFIG[]dnl
])dnl PKG_CHECK_MODULES_STATIC


dnl PKG_INSTALLDIR([DIRECTORY])
dnl -------------------------
dnl Since: 0.27
dnl
dnl Substitutes the variable pkgconfigdir as the location where a module
dnl should install pkg-config .pc files. By default the directory is
dnl $libdir/pkgconfig, but the default can be changed by passing
dnl DIRECTORY. The user can override through the --with-pkgconfigdir
dnl parameter.
AC_DEFUN([PKG_INSTALLDIR],
[m4_pushdef([pkg_default], [m4_default([$1], ['${libdir}/pkgconfig'])])
m4_pushdef([pkg_description],
    [pkg-config installation directory @<:@]pkg_default[@:>@])
AC_ARG_WITH([pkgconfigdir],
    [AS_HELP_STRING([--with-pkgconfigdir], pkg_description)],,
    [with_pkgconfigdir=]pkg_default)
AC_SUBST([pkgconfigdir], [$with_pkgconfigdir])
m4_popdef([pkg_default])
m4_popdef([pkg_description])
])dnl PKG_INSTALLDIR


dnl PKG_NOARCH_INSTALLDIR([DIRECTORY])
dnl --------------------------------
dnl Since: 0.27
dnl
dnl Substitutes the variable noarch_pkgconfigdir as the location where a
dnl module should install arch-independent pkg-config .pc files. By
dnl default the directory is $datadir/pkgconfig, but the default can be
dnl changed by passing DIRECTORY. The user can override through the
dnl --with-noarch-pkgconfigdir parameter.
AC_DEFUN([PKG_NOARCH_INSTALLDIR],
[m4_pushdef([pkg_default], [m4_default([$1], ['${datadir}/pkgconfig'])])
m4_pushdef([pkg_description],
    [pkg-config arch-independent installation directory @<:@]pkg_default[@:>@])
AC_ARG_WITH([noarch-pkgconfigdir],
    [AS_HELP_STRING([--with-noarch-pkgconfigdir], pkg_description)],,
    [with_noarch_pkgconfigdir=]pkg_default)
AC_SUBST([noarch_pkgconfigdir], [$with_noarch_pkgconfigdir])
m4_popdef([pkg_default])
m4_popdef([pkg_description])
])dnl PKG_NOARCH_INSTALLDIR


dnl PKG_CHECK_VAR(VARIABLE, MODULE, CONFIG-VARIABLE,
dnl [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl -------------------------------------------
dnl Since: 0.28
dnl
dnl Retrieves the value of the pkg-config variable for the given module.
AC_DEFUN([PKG_CHECK_VAR],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
AC_ARG_VAR([$1], [value of $3 for $2, overriding pkg-config])dnl

_PKG_CONFIG([$1], [variable="][$3]["], [$2])
AS_VAR_COPY([$1], [pkg_cv_][$1])

AS_VAR_IF([$1], [""], [$5], [$4])dnl
])dnl PKG_CHECK_VAR

# Copyright (C) 2002-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_AUTOMAKE_VERSION(VERSION)
# ----------------------------
# Automake X.Y traces this macro to ensure aclocal.m4 has been
# generated from the m4 files accompanying Automake X.Y.
# (This private macro should not be called outside this file.)
AC_DEFUN([AM_AUTOMAKE_VERSION],
[am__api_version='1.15'
dnl Some users find AM_AUTOMAKE_VERSION and mistake it for a way to
dnl require some minimum version.  Point them to the right macro.
m4_if([$1], [1.15], [],
      [AC_FATAL([Do not call $0, use AM_INIT_AUTOMAKE([$1]).])])dnl
])

# _AM_AUTOCONF_VERSION(VERSION)
# -----------------------------
# aclocal traces this macro to find the Autoconf version.
# This is a private macro too.  Using m4_define simplifies
# the logic in aclocal, which can simply ignore this definition.
m4_define([_AM_AUTOCONF_VERSION], [])

# AM_SET_CURRENT_AUTOMAKE_VERSION
# -------------------------------
# Call AM_AUTOMAKE_VERSION and AM_AUTOMAKE_VERSION so they can be traced.
# This function is AC_REQUIREd by AM_INIT_AUTOMAKE.
AC_DEFUN([AM_SET_CURRENT_AUTOMAKE_VERSION],
[AM_AUTOMAKE_VERSION([1.15])dnl
m4_ifndef([AC_AUTOCONF_VERSION],
  [m4_copy([m4_PACKAGE_VERSION], [AC_AUTOCONF_VERSION])])dnl
_AM_AUTOCONF_VERSION(m4_defn([AC_AUTOCONF_VERSION]))])

# AM_AUX_DIR_EXPAND                                         -*- Autoconf -*-

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# For projects using AC_CONFIG_AUX_DIR([foo]), Autoconf sets
# $ac_aux_dir to '$srcdir/foo'.  In other projects, it is set to
# '$srcdir', '$srcdir/..', or '$srcdir/../..'.
#
# Of course, Automake must honor this variable whenever it calls a
# tool from the auxiliary directory.  The problem is that $srcdir (and
# therefore $ac_aux_dir as well) can be either absolute or relative,
# depending on how configure is run.  This is pretty annoying, since
# it makes $ac_aux_dir quite unusable in subdirectories: in the top
# source directory, any form will work fine, but in subdirectories a
# relative path needs to be adjusted first.
#
# $ac_aux_dir/missing
#    fails when called from a subdirectory if $ac_aux_dir is relative
# $top_srcdir/$ac_aux_dir/missing
#    fails if $ac_aux_dir is absolute,
#    fails when called from a subdirectory in a VPATH build with
#          a relative $ac_aux_dir
#
# The reason of the latter failure is that $top_srcdir and $ac_aux_dir
# are both prefixed by $srcdir.  In an in-source build this is usually
# harmless because $srcdir is '.', but things will broke when you
# start a VPATH build or use an absolute $srcdir.
#
# So we could use something similar to $top_srcdir/$ac_aux_dir/missing,
# iff we strip the leading $srcdir from $ac_aux_dir.  That would be:
#   am_aux_dir='\$(top_srcdir)/'`expr "$ac_aux_dir" : "$srcdir//*\(.*\)"`
# and then we would define $MISSING as
#   MISSING="\${SHELL} $am_aux_dir/missing"
# This will work as long as MISSING is not called from configure, because
# unfortunately $(top_srcdir) has no meaning in configure.
# However there are other variables, like CC, which are often used in
# configure, and could therefore not use this "fixed" $ac_aux_dir.
#
# Another solution, used here, is to always expand $ac_aux_dir to an
# absolute PATH.  The drawback is that using absolute paths prevent a
# configured tree to be moved without reconfiguration.

AC_DEFUN([AM_AUX_DIR_EXPAND],
[AC_REQUIRE([AC_CONFIG_AUX_DIR_DEFAULT])dnl
# Expand $ac_aux_dir to an absolute path.
am_aux_dir=`cd "$ac_aux_dir" && pwd`
])

# AM_CONDITIONAL                                            -*- Autoconf -*-

# Copyright (C) 1997-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_CONDITIONAL(NAME, SHELL-CONDITION)
# -------------------------------------
# Define a conditional.
AC_DEFUN([AM_CONDITIONAL],
[AC_PREREQ([2.52])dnl
 m4_if([$1], [TRUE],  [AC_FATAL([$0: invalid condition: $1])],
       [$1], [FALSE], [AC_FATAL([$0: invalid condition: $1])])dnl
AC_SUBST([$1_TRUE])dnl
AC_SUBST([$1_FALSE])dnl
_AM_SUBST_NOTMAKE([$1_TRUE])dnl
_AM_SUBST_NOTMAKE([$1_FALSE])dnl
m4_define([_AM_COND_VALUE_$1], [$2])dnl
if $2; then
  $1_TRUE=
  $1_FALSE='#'
else
  $1_TRUE='#'
  $1_FALSE=
fi
AC_CONFIG_COMMANDS_PRE(
[if test -z "${$1_TRUE}" && test -z "${$1_FALSE}"; then
  AC_MSG_ERROR([[conditional "$1" was never defined.
Usually this means the macro was only invoked conditionally.]])
fi])])

# Copyright (C) 1999-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.


# There are a few dirty hacks below to avoid letting 'AC_PROG_CC' be
# written in clear, in which case automake, when reading aclocal.m4,
# will think it sees a *use*, and therefore will trigger all it's
# C support machinery.  Also note that it means that autoscan, seeing
# CC etc. in the Makefile, will ask for an AC_PROG_CC use...


# _AM_DEPENDENCIES(NAME)
# ----------------------
# See how the compiler implements dependency checking.
# NAME is "CC", "CXX", "OBJC", "OBJCXX", "UPC", or "GJC".
# We try a few techniques and use that to set a single cache variable.
#
# We don't AC_REQUIRE the corresponding AC_PROG_CC since the latter was
# modified to invoke _AM_DEPENDENCIES(CC); we would have a circular
# dependency, and given that the user is not expected to run this macro,
# just rely on AC_PROG_CC.
AC_DEFUN([_AM_DEPENDENCIES],
[AC_REQUIRE([AM_SET_DEPDIR])dnl
AC_REQUIRE([AM_OUTPUT_DEPENDENCY_COMMANDS])dnl
AC_REQUIRE([AM_MAKE_INCLUDE])dnl
AC_REQUIRE([AM_DEP_TRACK])dnl

m4_if([$1], [CC],   [depcc="$CC"   am_compiler_list=],
      [$1], [CXX],  [depcc="$CXX"  am_compiler_list=],
      [$1], [OBJC], [depcc="$OBJC" am_compiler_list='gcc3 gcc'],
      [$1], [OBJCXX], [depcc="$OBJCXX" am_compiler_list='gcc3 gcc'],
      [$1], [UPC],  [depcc="$UPC"  am_compiler_list=],
      [$1], [GCJ],  [depcc="$GCJ"  am_compiler_list='gcc3 gcc'],
                    [depcc="$$1"   am_compiler_list=])

AC_CACHE_CHECK([dependency style of $depcc],
               [am_cv_$1_dependencies_compiler_type],
[if test -z "$AMDEP_TRUE" && test -f "$am_depcomp"; then
  # We make a subdir and do the tests there.  Otherwise we can end up
  # making bogus files that we don't know about and never remove.  For
  # instance it was reported that on HP-UX the gcc test will end up
  # making a dummy file named 'D' -- because '-MD' means "put the output
  # in D".
  rm -rf conftest.dir
  mkdir conftest.dir
  # Copy depcomp to subdir because otherwise we won't find it if we're
  # using a relative directory.
  cp "$am_depcomp" conftest.dir
  cd conftest.dir
  # We will build objects and dependencies in a subdirectory because
  # it helps to detect inapplicable dependency modes.  For instance
  # both Tru64's cc and ICC support -MD to output dependencies as a
  # side effect of compilation, but ICC will put the dependencies in
  # the current directory while Tru64 will put them in the object
  # directory.
  mkdir sub

  am_cv_$1_dependencies_compiler_type=none
  if test "$am_compiler_list" = ""; then
     am_compiler_list=`sed -n ['s/^#*\([a-zA-Z0-9]*\))$/\1/p'] < ./depcomp`
  fi
  am__universal=false
  m4_case([$1], [CC],
    [case " $depcc " in #(
     *\ -arch\ *\ -arch\ *) am__universal=true ;;
     esac],
    [CXX],
    [case " $depcc " in #(
     *\ -arch\ *\ -arch\ *) am__universal=true ;;
     esac])

  for depmode in $am_compiler_list; do
    # Setup a source with many dependencies, because some compilers
    # like to wrap large dependency lists on column 80 (with \), and
    # we should not choose a depcomp mode which is confused by this.
    #
    # We need to recreate these files for each test, as the compiler may
    # overwrite some of them when testing with obscure command lines.
    # This happens at least with the AIX C compiler.
    : > sub/conftest.c
    for i in 1 2 3 4 5 6; do
      echo '#include "conftst'$i'.h"' >> sub/conftest.c
      # Using ": > sub/conftst$i.h" creates only sub/conftst1.h with
      # Solaris 10 /bin/sh.
      echo '/* dummy */' > sub/conftst$i.h
    done
    echo "${am__include} ${am__quote}sub/conftest.Po${am__quote}" > confmf

    # We check with '-c' and '-o' for the sake of the "dashmstdout"
    # mode.  It turns out that the SunPro C++ compiler does not properly
    # handle '-M -o', and we need to detect this.  Also, some Intel
    # versions had trouble with output in subdirs.
    am__obj=sub/conftest.${OBJEXT-o}
    am__minus_obj="-o $am__obj"
    case $depmode in
    gcc)
      # This depmode causes a compiler race in universal mode.
      test "$am__universal" = false || continue
      ;;
    nosideeffect)
      # After this tag, mechanisms are not by side-effect, so they'll
      # only be used when explicitly requested.
      if test "x$enable_dependency_tracking" = xyes; then
	continue
      else
	break
      fi
      ;;
    msvc7 | msvc7msys | msvisualcpp | msvcmsys)
      # This compiler won't grok '-c -o', but also, the minuso test has
      # not run yet.  These depmodes are late enough in the game, and
      # so weak that their functioning should not be impacted.
      am__obj=conftest.${OBJEXT-o}
      am__minus_obj=
      ;;
    none) break ;;
    esac
    if depmode=$depmode \
       source=sub/conftest.c object=$am__obj \
       depfile=sub/conftest.Po tmpdepfile=sub/conftest.TPo \
       $SHELL ./depcomp $depcc -c $am__minus_obj sub/conftest.c \
         >/dev/null 2>conftest.err &&
       grep sub/conftst1.h sub/conftest.Po > /dev/null 2>&1 &&
       grep sub/conftst6.h sub/conftest.Po > /dev/null 2>&1 &&
       grep $am__obj sub/conftest.Po > /dev/null 2>&1 &&
       ${MAKE-make} -s -f confmf > /dev/null 2>&1; then
      # icc doesn't choke on unknown options, it will just issue warnings
      # or remarks (even with -Werror).  So we grep stderr for any message
      # that says an option was ignored or not supported.
      # When given -MP, icc 7.0 and 7.1 complain thusly:
      #   icc: Command line warning: ignoring option '-M'; no argument required
      # The diagnosis changed in icc 8.0:
      #   icc: Command line remark: option '-MP' not supported
      if (grep 'ignoring option' conftest.err ||
          grep 'not supported' conftest.err) >/dev/null 2>&1; then :; else
        am_cv_$1_dependencies_compiler_type=$depmode
        break
      fi
    fi
  done

  cd ..
  rm -rf conftest.dir
else
  am_cv_$1_dependencies_compiler_type=none
fi
])
AC_SUBST([$1DEPMODE], [depmode=$am_cv_$1_dependencies_compiler_type])
AM_CONDITIONAL([am__fastdep$1], [
  test "x$enable_dependency_tracking" != xno \
  && test "$am_cv_$1_dependencies_compiler_type" = gcc3])
])


# AM_SET_DEPDIR
# -------------
# Choose a directory name for dependency files.
# This macro is AC_REQUIREd in _AM_DEPENDENCIES.
AC_DEFUN([AM_SET_DEPDIR],
[AC_REQUIRE([AM_SET_LEADING_DOT])dnl
AC_SUBST([DEPDIR], ["${am__leading_dot}deps"])dnl
])


# AM_DEP_TRACK
# ------------
AC_DEFUN([AM_DEP_TRACK],
[AC_ARG_ENABLE([dependency-tracking], [dnl
AS_HELP_STRING(
  [--enable-dependency-tracking],
  [do not reject slow dependency extractors])
AS_HELP_STRING(
  [--disable-dependency-tracking],
  [speeds up one-time build])])
if test "x$enable_dependency_tracking" != xno; then
  am_depcomp="$ac_aux_dir/depcomp"
  AMDEPBACKSLASH='\'
  am__nodep='_no'
fi
AM_CONDITIONAL([AMDEP], [test "x$enable_dependency_tracking" != xno])
AC_SUBST([AMDEPBACKSLASH])dnl
_AM_SUBST_NOTMAKE([AMDEPBACKSLASH])dnl
AC_SUBST([am__nodep])dnl
_AM_SUBST_NOTMAKE([am__nodep])dnl
])

# Generate code to set up dependency tracking.              -*- Autoconf -*-

# Copyright (C) 1999-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.


# _AM_OUTPUT_DEPENDENCY_COMMANDS
# ------------------------------
AC_DEFUN([_AM_OUTPUT_DEPENDENCY_COMMANDS],
[{
  # Older Autoconf quotes --file arguments for eval, but not when files
  # are listed without --file.  Let's play safe and only enable the eval
  # if we detect the quoting.
  case $CONFIG_FILES in
  *\'*) eval set x "$CONFIG_FILES" ;;
  *)   set x $CONFIG_FILES ;;
  esac
  shift
  for mf
  do
    # Strip MF so we end up with the name of the file.
    mf=`echo "$mf" | sed -e 's/:.*$//'`
    # Check whether this is an Automake generated Makefile or not.
    # We used to match only the files named 'Makefile.in', but
    # some people rename them; so instead we look at the file content.
    # Grep'ing the first line is not enough: some people post-process
    # each Makefile.in and add a new line on top of each file to say so.
    # Grep'ing the whole file is not good either: AIX grep has a line
    # limit of 2048, but all sed's we know have understand at least 4000.
    if sed -n 's,^#.*generated by automake.*,X,p' "$mf" | grep X >/dev/null 2>&1; then
      dirpart=`AS_DIRNAME("$mf")`
    else
      continue
    fi
    # Extract the definition of DEPDIR, am__include, and am__quote
    # from the Makefile without running 'make'.
    DEPDIR=`sed -n 's/^DEPDIR = //p' < "$mf"`
    test -z "$DEPDIR" && continue
    am__include=`sed -n 's/^am__include = //p' < "$mf"`
    test -z "$am__include" && continue
    am__quote=`sed -n 's/^am__quote = //p' < "$mf"`
    # Find all dependency output files, they are included files with
    # $(DEPDIR) in their names.  We invoke sed twice because it is the
    # simplest approach to changing $(DEPDIR) to its actual value in the
    # expansion.
    for file in `sed -n "
      s/^$am__include $am__quote\(.*(DEPDIR).*\)$am__quote"'$/\1/p' <"$mf" | \
	 sed -e 's/\$(DEPDIR)/'"$DEPDIR"'/g'`; do
      # Make sure the directory exists.
      test -f "$dirpart/$file" && continue
      fdir=`AS_DIRNAME(["$file"])`
      AS_MKDIR_P([$dirpart/$fdir])
      # echo "creating $dirpart/$file"
      echo '# dummy' > "$dirpart/$file"
    done
  done
}
])# _AM_OUTPUT_DEPENDENCY_COMMANDS


# AM_OUTPUT_DEPENDENCY_COMMANDS
# -----------------------------
# This macro should only be invoked once -- use via AC_REQUIRE.
#
# This code is only required when automatic dependency tracking
# is enabled.  FIXME.  This creates each '.P' file that we will
# need in order to bootstrap the dependency handling code.
AC_DEFUN([AM_OUTPUT_DEPENDENCY_COMMANDS],
[AC_CONFIG_COMMANDS([depfiles],
     [test x"$AMDEP_TRUE" != x"" || _AM_OUTPUT_DEPENDENCY_COMMANDS],
     [AMDEP_TRUE="$AMDEP_TRUE" ac_aux_dir="$ac_aux_dir"])
])

# Do all the work for Automake.                             -*- Autoconf -*-

# Copyright (C) 1996-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This macro actually does too much.  Some checks are only needed if
# your package does certain things.  But this isn't really a big deal.

dnl Redefine AC_PROG_CC to automatically invoke _AM_PROG_CC_C_O.
m4_define([AC_PROG_CC],
m4_defn([AC_PROG_CC])
[_AM_PROG_CC_C_O
])

# AM_INIT_AUTOMAKE(PACKAGE, VERSION, [NO-DEFINE])
# AM_INIT_AUTOMAKE([OPTIONS])
# -----------------------------------------------
# The call with PACKAGE and VERSION arguments is the old style
# call (pre autoconf-2.50), which is being phased out.  PACKAGE
# and VERSION should now be passed to AC_INIT and removed from
# the call to AM_INIT_AUTOMAKE.
# We support both call styles for the transition.  After
# the next Automake release, Autoconf can make the AC_INIT
# arguments mandatory, and then we can depend on a new Autoconf
# release and drop the old call support.
AC_DEFUN([AM_INIT_AUTOMAKE],
[AC_PREREQ([2.65])dnl
dnl Autoconf wants to disallow AM_ names.  We explicitly allow
dnl the ones we care about.
m4_pattern_allow([^AM_[A-Z]+FLAGS$])dnl
AC_REQUIRE([AM_SET_CURRENT_AUTOMAKE_VERSION])dnl
AC_REQUIRE([AC_PROG_INSTALL])dnl
if test "`cd $srcdir && pwd`" != "`pwd`"; then
  # Use -I$(srcdir) only when $(srcdir) != ., so that make's output
  # is not polluted with repeated "-I."
  AC_SUBST([am__isrc], [' -I$(srcdir)'])_AM_SUBST_NOTMAKE([am__isrc])dnl
  # test to see if srcdir already configured
  if test -f $srcdir/config.status; then
    AC_MSG_ERROR([source directory already configured; run "make distclean" there first])
  fi
fi

# test whether we have cygpath
if test -z "$CYGPATH_W"; then
  if (cygpath --version) >/dev/null 2>/dev/null; then
    CYGPATH_W='cygpath -w'
  else
    CYGPATH_W=echo
  fi
fi
AC_SUBST([CYGPATH_W])

# Define the identity of the package.
dnl Distinguish between old-style and new-style calls.
m4_ifval([$2],
[AC_DIAGNOSE([obsolete],
             [$0: two- and three-arguments forms are deprecated.])
m4_ifval([$3], [_AM_SET_OPTION([no-define])])dnl
 AC_SUBST([PACKAGE], [$1])dnl
 AC_SUBST([VERSION], [$2])],
[_AM_SET_OPTIONS([$1])dnl
dnl Diagnose old-style AC_INIT with new-style AM_AUTOMAKE_INIT.
m4_if(
  m4_ifdef([AC_PACKAGE_NAME], [ok]):m4_ifdef([AC_PACKAGE_VERSION], [ok]),
  [ok:ok],,
  [m4_fatal([AC_INIT should be called with package and version arguments])])dnl
 AC_SUBST([PACKAGE], ['AC_PACKAGE_TARNAME'])dnl
 AC_SUBST([VERSION], ['AC_PACKAGE_VERSION'])])dnl

_AM_IF_OPTION([no-define],,
[AC_DEFINE_UNQUOTED([PACKAGE], ["$PACKAGE"], [Name of package])
 AC_DEFINE_UNQUOTED([VERSION], ["$VERSION"], [Version number of package])])dnl

# Some tools Automake needs.
AC_REQUIRE([AM_SANITY_CHECK])dnl
AC_REQUIRE([AC_ARG_PROGRAM])dnl
AM_MISSING_PROG([ACLOCAL], [aclocal-${am__api_version}])
AM_MISSING_PROG([AUTOCONF], [autoconf])
AM_MISSING_PROG([AUTOMAKE], [automake-${am__api_version}])
AM_MISSING_PROG([AUTOHEADER], [autoheader])
AM_MISSING_PROG([MAKEINFO], [makeinfo])
AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
AC_REQUIRE([AM_PROG_INSTALL_STRIP])dnl
AC_REQUIRE([AC_PROG_MKDIR_P])dnl
# For better backward compatibility.  To be removed once Automake 1.9.x
# dies out for good.  For more background, see:
# <http://lists.gnu.org/archive/html/automake/2012-07/msg00001.html>
# <http://lists.gnu.org/archive/html/automake/2012-07/msg00014.html>
AC_SUBST([mkdir_p], ['$(MKDIR_P)'])
# We need awk for the "check" target (and possibly the TAP driver).  The
# system "awk" is bad on some platforms.
AC_REQUIRE([AC_PROG_AWK])dnl
AC_REQUIRE([AC_PROG_MAKE_SET])dnl
AC_REQUIRE([AM_SET_LEADING_DOT])dnl
_AM_IF_OPTION([tar-ustar], [_AM_PROG_TAR([ustar])],
	      [_AM_IF_OPTION([tar-pax], [_AM_PROG_TAR([pax])],
			     [_AM_PROG_TAR([v7])])])
_AM_IF_OPTION([no-dependencies],,
[AC_PROVIDE_IFELSE([AC_PROG_CC],
		  [_AM_DEPENDENCIES([CC])],
		  [m4_define([AC_PROG_CC],
			     m4_defn([AC_PROG_CC])[_AM_DEPENDENCIES([CC])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_CXX],
		  [_AM_DEPENDENCIES([CXX])],
		  [m4_define([AC_PROG_CXX],
			     m4_defn([AC_PROG_CXX])[_AM_DEPENDENCIES([CXX])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_OBJC],
		  [_AM_DEPENDENCIES([OBJC])],
		  [m4_define([AC_PROG_OBJC],
			     m4_defn([AC_PROG_OBJC])[_AM_DEPENDENCIES([OBJC])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_OBJCXX],
		  [_AM_DEPENDENCIES([OBJCXX])],
		  [m4_define([AC_PROG_OBJCXX],
			     m4_defn([AC_PROG_OBJCXX])[_AM_DEPENDENCIES([OBJCXX])])])dnl
])
AC_REQUIRE([AM_SILENT_RULES])dnl
dnl The testsuite driver may need to know about EXEEXT, so add the
dnl 'am__EXEEXT' conditional if _AM_COMPILER_EXEEXT was seen.  This
dnl macro is hooked onto _AC_COMPILER_EXEEXT early, see below.
AC_CONFIG_COMMANDS_PRE(dnl
[m4_provide_if([_AM_COMPILER_EXEEXT],
  [AM_CONDITIONAL([am__EXEEXT], [test -n "$EXEEXT"])])])dnl

# POSIX will say in a future version that running "rm -f" with no argument
# is OK; and we want to be able to make that assumption in our Makefile
# recipes.  So use an aggressive probe to check that the usage we want is
# actually supported "in the wild" to an acceptable degree.
# See automake bug#10828.
# To make any issue more visible, cause the running configure to be aborted
# by default if the 'rm' program in use doesn't match our expectations; the
# user can still override this though.
if rm -f && rm -fr && rm -rf; then : OK; else
  cat >&2 <<'END'
Oops!

Your 'rm' program seems unable to run without file operands specified
on the command line, even when the '-f' option is present.  This is contrary
to the behaviour of most rm programs out there, and not conforming with
the upcoming POSIX standard: <http://austingroupbugs.net/view.php?id=542>

Please tell bug-automake@gnu.org about your system, including the value
of your $PATH and any error possibly output before this message.  This
can help us improve future automake versions.

END
  if test x"$ACCEPT_INFERIOR_RM_PROGRAM" = x"yes"; then
    echo 'Configuration will proceed anyway, since you have set the' >&2
    echo 'ACCEPT_INFERIOR_RM_PROGRAM variable to "yes"' >&2
    echo >&2
  else
    cat >&2 <<'END'
Aborting the configuration process, to ensure you take notice of the issue.

You can download and install GNU coreutils to get an 'rm' implementation
that behaves properly: <http://www.gnu.org/software/coreutils/>.

If you want to complete the configuration process using your problematic
'rm' anyway, export the environment variable ACCEPT_INFERIOR_RM_PROGRAM
to "yes", and re-run configure.

END
    AC_MSG_ERROR([Your 'rm' program is bad, sorry.])
  fi
fi
dnl The trailing newline in this macro's definition is deliberate, for
dnl backward compatibility and to allow trailing 'dnl'-style comments
dnl after the AM_INIT_AUTOMAKE invocation. See automake bug#16841.
])

dnl Hook into '_AC_COMPILER_EXEEXT' early to learn its expansion.  Do not
dnl add the conditional right here, as _AC_COMPILER_EXEEXT may be further
dnl mangled by Autoconf and run in a shell conditional statement.
m4_define([_AC_COMPILER_EXEEXT],
m4_defn([_AC_COMPILER_EXEEXT])[m4_provide([_AM_COMPILER_EXEEXT])])

# When config.status generates a header, we must update the stamp-h file.
# This file resides in the same directory as the config header
# that is generated.  The stamp files are numbered to have different names.

# Autoconf calls _AC_AM_CONFIG_HEADER_HOOK (when defined) in the
# loop where config.status creates the headers, so we can generate
# our stamp files there.
AC_DEFUN([_AC_AM_CONFIG_HEADER_HOOK],
[# Compute $1's index in $config_headers.
_am_arg=$1
_am_stamp_count=1
for _am_header in $config_headers :; do
  case $_am_header in
    $_am_arg | $_am_arg:* )
      break ;;
    * )
      _am_stamp_count=`expr $_am_stamp_count + 1` ;;
  esac
done
echo "timestamp for $_am_arg" >`AS_DIRNAME(["$_am_arg"])`/stamp-h[]$_am_stamp_count])

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_PROG_INSTALL_SH
# ------------------
# Define $install_sh.
AC_DEFUN([AM_PROG_INSTALL_SH],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
if test x"${install_sh+set}" != xset; then
  case $am_aux_dir in
  *\ * | *\	*)
    install_sh="\${SHELL} '$am_aux_dir/install-sh'" ;;
  *)
    install_sh="\${SHELL} $am_aux_dir/install-sh"
  esac
fi
AC_SUBST([install_sh])])

# Copyright (C) 2003-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# Check whether the underlying file-system supports filenames
# with a leading dot.  For instance MS-DOS doesn't.
AC_DEFUN([AM_SET_LEADING_DOT],
[rm -rf .tst 2>/dev/null
mkdir .tst 2>/dev/null
if test -d .tst; then
  am__leading_dot=.
else
  am__leading_dot=_
fi
rmdir .tst 2>/dev/null
AC_SUBST([am__leading_dot])])

# Check to see how 'make' treats includes.	            -*- Autoconf -*-

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_MAKE_INCLUDE()
# -----------------
# Check to see how make treats includes.
AC_DEFUN([AM_MAKE_INCLUDE],
[am_make=${MAKE-make}
cat > confinc << 'END'
am__doit:
	@echo this is the am__doit target
.PHONY: am__doit
END
# If we don't find an include directive, just comment out the code.
AC_MSG_CHECKING([for style of include used by $am_make])
am__include="#"
am__quote=
_am_result=none
# First try GNU make style include.
echo "include confinc" > confmf
# Ignore all kinds of additional output from 'make'.
case `$am_make -s -f confmf 2> /dev/null` in #(
*the\ am__doit\ target*)
  am__include=include
  am__quote=
  _am_result=GNU
  ;;
esac
# Now try BSD make style include.
if test "$am__include" = "#"; then
   echo '.include "confinc"' > confmf
   case `$am_make -s -f confmf 2> /dev/null` in #(
   *the\ am__doit\ target*)
     am__include=.include
     am__quote="\""
     _am_result=BSD
     ;;
   esac
fi
AC_SUBST([am__include])
AC_SUBST([am__quote])
AC_MSG_RESULT([$_am_result])
rm -f confinc confmf
])

# Fake the existence of programs that GNU maintainers use.  -*- Autoconf -*-

# Copyright (C) 1997-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_MISSING_PROG(NAME, PROGRAM)
# ------------------------------
AC_DEFUN([AM_MISSING_PROG],
[AC_REQUIRE([AM_MISSING_HAS_RUN])
$1=${$1-"${am_missing_run}$2"}
AC_SUBST($1)])

# AM_MISSING_HAS_RUN
# ------------------
# Define MISSING if not defined so far and test if it is modern enough.
# If it is, set am_missing_run to use it, otherwise, to nothing.
AC_DEFUN([AM_MISSING_HAS_RUN],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
AC_REQUIRE_AUX_FILE([missing])dnl
if test x"${MISSING+set}" != xset; then
  case $am_aux_dir in
  *\ * | *\	*)
    MISSING="\${SHELL} \"$am_aux_dir/missing\"" ;;
  *)
    MISSING="\${SHELL} $am_aux_dir/missing" ;;
  esac
fi
# Use eval to expand $SHELL
if eval "$MISSING --is-lightweight"; then
  am_missing_run="$MISSING "
else
  am_missing_run=
  AC_MSG_WARN(['missing' script is too old or missing])
fi
])

# Helper functions for option handling.                     -*- Autoconf -*-

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_MANGLE_OPTION(NAME)
# -----------------------
AC_DEFUN([_AM_MANGLE_OPTION],
[[_AM_OPTION_]m4_bpatsubst($1, [[^a-zA-Z0-9_]], [_])])

# _AM_SET_OPTION(NAME)
# --------------------
# Set option NAME.  Presently that only means defining a flag for this option.
AC_DEFUN([_AM_SET_OPTION],
[m4_define(_AM_MANGLE_OPTION([$1]), [1])])

# _AM_SET_OPTIONS(OPTIONS)
# ------------------------
# OPTIONS is a space-separated list of Automake options.
AC_DEFUN([_AM_SET_OPTIONS],
[m4_foreach_w([_AM_Option], [$1], [_AM_SET_OPTION(_AM_Option)])])

# _AM_IF_OPTION(OPTION, IF-SET, [IF-NOT-SET])
# -------------------------------------------
# Execute IF-SET if OPTION is set, IF-NOT-SET otherwise.
AC_DEFUN([_AM_IF_OPTION],
[m4_ifset(_AM_MANGLE_OPTION([$1]), [$2], [$3])])

# Copyright (C) 1999-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_PROG_CC_C_O
# ---------------
# Like AC_PROG_CC_C_O, but changed for automake.  We rewrite AC_PROG_CC
# to automatically call this.
AC_DEFUN([_AM_PROG_CC_C_O],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
AC_REQUIRE_AUX_FILE([compile])dnl
AC_LANG_PUSH([C])dnl
AC_CACHE_CHECK(
  [whether $CC understands -c and -o together],
  [am_cv_prog_cc_c_o],
  [AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])
  # Make sure it works both with $CC and with simple cc.
  # Following AC_PROG_CC_C_O, we do the test twice because some
  # compilers refuse to overwrite an existing .o file with -o,
  # though they will create one.
  am_cv_prog_cc_c_o=yes
  for am_i in 1 2; do
    if AM_RUN_LOG([$CC -c conftest.$ac_ext -o conftest2.$ac_objext]) \
         && test -f conftest2.$ac_objext; then
      : OK
    else
      am_cv_prog_cc_c_o=no
      break
    fi
  done
  rm -f core conftest*
  unset am_i])
if test "$am_cv_prog_cc_c_o" != yes; then
   # Losing compiler, so override with the script.
   # FIXME: It is wrong to rewrite CC.
   # But if we don't then we get into trouble of one sort or another.
   # A longer-term fix would be to have automake use am__CC in this case,
   # and then we could set am__CC="\$(top_srcdir)/compile \$(CC)"
   CC="$am_aux_dir/compile $CC"
fi
AC_LANG_POP([C])])

# For backward compatibility.
AC_DEFUN_ONCE([AM_PROG_CC_C_O], [AC_REQUIRE([AC_PROG_CC])])

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_RUN_LOG(COMMAND)
# -------------------
# Run COMMAND, save the exit status in ac_status, and log it.
# (This has been adapted from Autoconf's _AC_RUN_LOG macro.)
AC_DEFUN([AM_RUN_LOG],
[{ echo "$as_me:$LINENO: $1" >&AS_MESSAGE_LOG_FD
   ($1) >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD
   ac_status=$?
   echo "$as_me:$LINENO: \$? = $ac_status" >&AS_MESSAGE_LOG_FD
   (exit $ac_status); }])

# Check to make sure that the build environment is sane.    -*- Autoconf -*-

# Copyright (C) 1996-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_SANITY_CHECK
# ---------------
AC_DEFUN([AM_SANITY_CHECK],
[AC_MSG_CHECKING([whether build environment is sane])
# Reject unsafe characters in $srcdir or the absolute working directory
# name.  Accept space and tab only in the latter.
am_lf='
'
case `pwd` in
  *[[\\\"\#\$\&\'\`$am_lf]]*)
    AC_MSG_ERROR([unsafe absolute working directory name]);;
esac
case $srcdir in
  *[[\\\"\#\$\&\'\`$am_lf\ \	]]*)
    AC_MSG_ERROR([unsafe srcdir value: '$srcdir']);;
esac

# Do 'set' in a subshell so we don't clobber the current shell's
# arguments.  Must try -L first in case configure is actually a
# symlink; some systems play weird games with the mod time of symlinks
# (eg FreeBSD returns the mod time of the symlink's containing
# directory).
if (
   am_has_slept=no
   for am_try in 1 2; do
     echo "timestamp, slept: $am_has_slept" > conftest.file
     set X `ls -Lt "$srcdir/configure" conftest.file 2> /dev/null`
     if test "$[*]" = "X"; then
	# -L didn't work.
	set X `ls -t "$srcdir/configure" conftest.file`
     fi
     if test "$[*]" != "X $srcdir/configure conftest.file" \
	&& test "$[*]" != "X conftest.file $srcdir/configure"; then

	# If neither matched, then we have a broken ls.  This can happen
	# if, for instance, CONFIG_SHELL is bash and it inherits a
	# broken ls alias from the environment.  This has actually
	# happened.  Such a system could not be considered "sane".
	AC_MSG_ERROR([ls -t appears to fail.  Make sure there is not a broken
  alias in your environment])
     fi
     if test "$[2]" = conftest.file || test $am_try -eq 2; then
       break
     fi
     # Just in case.
     sleep 1
     am_has_slept=yes
   done
   test "$[2]" = conftest.file
   )
then
   # Ok.
   :
else
   AC_MSG_ERROR([newly created file is older than distributed files!
Check your system clock])
fi
AC_MSG_RESULT([yes])
# If we didn't sleep, we still need to ensure time stamps of config.status and
# generated files are strictly newer.
am_sleep_pid=
if grep 'slept: no' conftest.file >/dev/null 2>&1; then
  ( sleep 1 ) &
  am_sleep_pid=$!
fi
AC_CONFIG_COMMANDS_PRE(
  [AC_MSG_CHECKING([that generated files are newer than configure])
   if test -n "$am_sleep_pid"; then
     # Hide warnings about reused PIDs.
     wait $am_sleep_pid 2>/dev/null
   fi
   AC_MSG_RESULT([done])])
rm -f conftest.file
])

# Copyright (C) 2009-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_SILENT_RULES([DEFAULT])
# --------------------------
# Enable less verbose build rules; with the default set to DEFAULT
# ("yes" being less verbose, "no" or empty being verbose).
AC_DEFUN([AM_SILENT_RULES],
[AC_ARG_ENABLE([silent-rules], [dnl
AS_HELP_STRING(
  [--enable-silent-rules],
  [less verbose build output (undo: "make V=1")])
AS_HELP_STRING(
  [--disable-silent-rules],
  [verbose build output (undo: "make V=0")])dnl
])
case $enable_silent_rules in @%:@ (((
  yes) AM_DEFAULT_VERBOSITY=0;;
   no) AM_DEFAULT_VERBOSITY=1;;
    *) AM_DEFAULT_VERBOSITY=m4_if([$1], [yes], [0], [1]);;
esac
dnl
dnl A few 'make' implementations (e.g., NonStop OS and NextStep)
dnl do not support nested variable expansions.
dnl See automake bug#9928 and bug#10237.
am_make=${MAKE-make}
AC_CACHE_CHECK([whether $am_make supports nested variables],
   [am_cv_make_support_nested_variables],
   [if AS_ECHO([['TRUE=$(BAR$(V))
BAR0=false
BAR1=true
V=1
am__doit:
	@$(TRUE)
.PHONY: am__doit']]) | $am_make -f - >/dev/null 2>&1; then
  am_cv_make_support_nested_variables=yes
else
  am_cv_make_support_nested_variables=no
fi])
if test $am_cv_make_support_nested_variables = yes; then
  dnl Using '$V' instead of '$(V)' breaks IRIX make.
  AM_V='$(V)'
  AM_DEFAULT_V='$(AM_DEFAULT_VERBOSITY)'
else
  AM_V=$AM_DEFAULT_VERBOSITY
  AM_DEFAULT_V=$AM_DEFAULT_VERBOSITY
fi
AC_SUBST([AM_V])dnl
AM_SUBST_NOTMAKE([AM_V])dnl
AC_SUBST([AM_DEFAULT_V])dnl
AM_SUBST_NOTMAKE([AM_DEFAULT_V])dnl
AC_SUBST([AM_DEFAULT_VERBOSITY])dnl
AM_BACKSLASH='\'
AC_SUBST([AM_BACKSLASH])dnl
_AM_SUBST_NOTMAKE([AM_BACKSLASH])dnl
])

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_PROG_INSTALL_STRIP
# ---------------------
# One issue with vendor 'install' (even GNU) is that you can't
# specify the program used to strip binaries.  This is especially
# annoying in cross-compiling environments, where the build's strip
# is unlikely to handle the host's binaries.
# Fortunately install-sh will honor a STRIPPROG variable, so we
# always use install-sh in "make install-strip", and initialize
# STRIPPROG with the value of the STRIP variable (set by the user).
AC_DEFUN([AM_PROG_INSTALL_STRIP],
[AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
# Installed binaries are usually stripped using 'strip' when the user
# run "make install-strip".  However 'strip' might not be the right
# tool to use in cross-compilation environments, therefore Automake
# will honor the 'STRIP' environment variable to overrule this program.
dnl Don't test for $cross_compiling = yes, because it might be 'maybe'.
if test "$cross_compiling" != no; then
  AC_CHECK_TOOL([STRIP], [strip], :)
fi
INSTALL_STRIP_PROGRAM="\$(install_sh) -c -s"
AC_SUBST([INSTALL_STRIP_PROGRAM])])

# Copyright (C) 2006-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_SUBST_NOTMAKE(VARIABLE)
# ---------------------------
# Prevent Automake from outputting VARIABLE = @VARIABLE@ in Makefile.in.
# This macro is traced by Automake.
AC_DEFUN([_AM_SUBST_NOTMAKE])

# AM_SUBST_NOTMAKE(VARIABLE)
# --------------------------
# Public sister of _AM_SUBST_NOTMAKE.
AC_DEFUN([AM_SUBST_NOTMAKE], [_AM_SUBST_NOTMAKE($@)])

# Check how to create a tarball.                            -*- Autoconf -*-

# Copyright (C) 2004-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_PROG_TAR(FORMAT)
# --------------------
# Check how to create a tarball in format FORMAT.
# FORMAT should be one of 'v7', 'ustar', or 'pax'.
#
# Substitute a variable $(am__tar) that is a command
# writing to stdout a FORMAT-tarball containing the directory
# $tardir.
#     tardir=directory && $(am__tar) > result.tar
#
# Substitute a variable $(am__untar) that extract such
# a tarball read from stdin.
#     $(am__untar) < result.tar
#
AC_DEFUN([_AM_PROG_TAR],
[# Always define AMTAR for backward compatibility.  Yes, it's still used
# in the wild :-(  We should find a proper way to deprecate it ...
AC_SUBST([AMTAR], ['$${TAR-tar}'])

# We'll loop over all known methods to create a tar archive until one works.
_am_tools='gnutar m4_if([$1], [ustar], [plaintar]) pax cpio none'

m4_if([$1], [v7],
  [am__tar='$${TAR-tar} chof - "$$tardir"' am__untar='$${TAR-tar} xf -'],

  [m4_case([$1],
    [ustar],
     [# The POSIX 1988 'ustar' format is defined with fixed-size fields.
      # There is notably a 21 bits limit for the UID and the GID.  In fact,
      # the 'pax' utility can hang on bigger UID/GID (see automake bug#8343
      # and bug#13588).
      am_max_uid=2097151 # 2^21 - 1
      am_max_gid=$am_max_uid
      # The $UID and $GID variables are not portable, so we need to resort
      # to the POSIX-mandated id(1) utility.  Errors in the 'id' calls
      # below are definitely unexpected, so allow the users to see them
      # (that is, avoid stderr redirection).
      am_uid=`id -u || echo unknown`
      am_gid=`id -g || echo unknown`
      AC_MSG_CHECKING([whether UID '$am_uid' is supported by ustar format])
      if test $am_uid -le $am_max_uid; then
         AC_MSG_RESULT([yes])
      else
         AC_MSG_RESULT([no])
         _am_tools=none
      fi
      AC_MSG_CHECKING([whether GID '$am_gid' is supported by ustar format])
      if test $am_gid -le $am_max_gid; then
         AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no])
        _am_tools=none
      fi],

  [pax],
    [],

  [m4_fatal([Unknown tar format])])

  AC_MSG_CHECKING([how to create a $1 tar archive])

  # Go ahead even if we have the value already cached.  We do so because we
  # need to set the values for the 'am__tar' and 'am__untar' variables.
  _am_tools=${am_cv_prog_tar_$1-$_am_tools}

  for _am_tool in $_am_tools; do
    case $_am_tool in
    gnutar)
      for _am_tar in tar gnutar gtar; do
        AM_RUN_LOG([$_am_tar --version]) && break
      done
      am__tar="$_am_tar --format=m4_if([$1], [pax], [posix], [$1]) -chf - "'"$$tardir"'
      am__tar_="$_am_tar --format=m4_if([$1], [pax], [posix], [$1]) -chf - "'"$tardir"'
      am__untar="$_am_tar -xf -"
      ;;
    plaintar)
      # Must skip GNU tar: if it does not support --format= it doesn't create
      # ustar tarball either.
      (tar --version) >/dev/null 2>&1 && continue
      am__tar='tar chf - "$$tardir"'
      am__tar_='tar chf - "$tardir"'
      am__untar='tar xf -'
      ;;
    pax)
      am__tar='pax -L -x $1 -w "$$tardir"'
      am__tar_='pax -L -x $1 -w "$tardir"'
      am__untar='pax -r'
      ;;
    cpio)
      am__tar='find "$$tardir" -print | cpio -o -H $1 -L'
      am__tar_='find "$tardir" -print | cpio -o -H $1 -L'
      am__untar='cpio -i -H $1 -d'
      ;;
    none)
      am__tar=false
      am__tar_=false
      am__untar=false
      ;;
    esac

    # If the value was cached, stop now.  We just wanted to have am__tar
    # and am__untar set.
    test -n "${am_cv_prog_tar_$1}" && break

    # tar/untar a dummy directory, and stop if the command works.
    rm -rf conftest.dir
    mkdir conftest.dir
    echo GrepMe > conftest.dir/file
    AM_RUN_LOG([tardir=conftest.dir && eval $am__tar_ >conftest.tar])
    rm -rf conftest.dir
    if test -s conftest.tar; then
      AM_RUN_LOG([$am__untar <conftest.tar])
      AM_RUN_LOG([cat conftest.dir/file])
      grep GrepMe conftest.dir/file >/dev/null 2>&1 && break
    fi
  done
  rm -rf conftest.dir

  AC_CACHE_VAL([am_cv_prog_tar_$1], [am_cv_prog_tar_$1=$_am_tool])
  AC_MSG_RESULT([$am_cv_prog_tar_$1])])

AC_SUBST([am__tar])
AC_SUBST([am__untar])
]) # _AM_PROG_TAR
