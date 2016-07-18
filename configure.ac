
# Check for type traits in <type_traits> and std namespace

AC_DEFUN([ACX_CHECK_TYPE_TRAITS_CXX0X], [
  AC_LANG_PUSH([C++])

  AC_MSG_CHECKING([for type traits <tr1/type_traits> or <type_trait>, for std or std:tr1, and if compiler option -std=c++0x is needed])
  acx_type_traits=no

# (1) =========================== no compiler option ===========================

  if test "$acx_type_traits" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
      	  [[#include <type_traits>]],
          [[typedef std::is_arithmetic<double> isArithmeticT;]]
      	)	  
      ],		    
      [
        AC_DEFINE([USE_TYPE_TRAITS],[1],[define if using <type_traits>.])
        AC_DEFINE([HAS_STD_TYPE_TRAITS],[1],[define if std type traits are available.])
        acx_type_traits=yes
      ]
    )
  fi

  if test "$acx_type_traits" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <type_traits>]],
          [[typedef std::tr1::is_arithmetic<double> isArithmeticT;]]
        )
      ],
      [
        AC_DEFINE([USE_TYPE_TRAITS],[1],[define if using <type_traits>.])
        AC_DEFINE([HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
        acx_type_traits=yes
      ]
    )
  fi
  
  if test "$acx_type_traits" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/type_traits>]],
          [[typedef std::is_arithmetic<double> isArithmeticT;]]
        )
      ],
      [
        AC_DEFINE([USE_TR1_TYPE_TRAITS],[1],[define if using <tr1/type_traits>.])
        AC_DEFINE([HAS_STD_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
        acx_type_traits=yes
      ]
    )
  fi
  
  if test "$acx_type_traits" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/type_traits>]],
          [[typedef std::tr1::is_arithmetic<double> isArithmeticT;]]
        )
      ],
      [
        AC_DEFINE([USE_TR1_TYPE_TRAITS],[1],[define if using <tr1/type_traits>.])
        AC_DEFINE([HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
        acx_type_traits=yes
      ]
    )
  fi

# (2) ========================= c++0x compiler option ===========================
  
  if test "$acx_type_traits" = no; then
    ac_save_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS -std=c++0x"
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
      	  [[#include <type_traits>]],
          [[typedef std::is_arithmetic<double> isArithmeticT;]]
      	)	  
      ],		    
      [
        AC_DEFINE([USE_TYPE_TRAITS],[1],[define if using <type_traits>.])
        AC_DEFINE([HAS_STD_TYPE_TRAITS],[1],[define if std type traits are available.])
        acx_type_traits=yes
      ]
    )
    if test "$acx_type_traits" = no; then
       CXXFLAGS="$ac_save_CXXFLAGS"
    fi
  fi



  if test "$acx_type_traits" = no; then
    ac_save_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS -std=c++0x"
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <type_traits>]],
          [[typedef std::tr1::is_arithmetic<double> isArithmeticT;]]
        )
      ],
      [
        AC_DEFINE([USE_TYPE_TRAITS],[1],[define if using <type_traits>.])
        AC_DEFINE([HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
        acx_type_traits=yes
      ]
    )
    if test "$acx_type_traits" = no; then
       CXXFLAGS="$ac_save_CXXFLAGS"
    fi
  fi


  
  if test "$acx_type_traits" = no; then
    ac_save_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS -std=c++0x"
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/type_traits>]],
          [[typedef std::is_arithmetic<double> isArithmeticT;]]
        )
      ],
      [
        AC_DEFINE([USE_TR1_TYPE_TRAITS],[1],[define if using <tr1/type_traits>.])
        AC_DEFINE([HAS_STD_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
        acx_type_traits=yes
      ]
    )
    if test "$acx_type_traits" = no; then
       CXXFLAGS="$ac_save_CXXFLAGS"
    fi
  fi


  
  if test "$acx_type_traits" = no; then
    ac_save_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS -std=c++0x"
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/type_traits>]],
          [[typedef std::tr1::is_arithmetic<double> isArithmeticT;]]
        )
      ],
      [
        AC_DEFINE([USE_TR1_TYPE_TRAITS],[1],[define if using <tr1/type_traits>.])
        AC_DEFINE([HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
        acx_type_traits=yes
      ]
    )
    if test "$acx_type_traits" = no; then
       CXXFLAGS="$ac_save_CXXFLAGS"
    fi
  fi


  # post type traits results
  AC_MSG_RESULT([$acx_type_traits])
  if test "$acx_type_traits" = no; then
    AC_MSG_ERROR([std and std::tr1 <type_traits> or <tr1/type_traits> are not supported.])
  fi  
  AC_LANG_POP([])  
])



# ===========================================================================





AC_INIT([iedera],[1.06],[laurent DT noe AT univ-lille1 DT fr],[iedera],[http://bioinfo.lifl.fr/yass/iedera])
AC_CONFIG_SRCDIR([main.cc])
AC_LANG([C++])
AM_INIT_AUTOMAKE
AC_PROG_CXX
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h])
AC_CHECK_FUNCS([gettimeofday sqrt strdup strtol])
ACX_CHECK_TYPE_TRAITS_CXX0X
CXXFLAGS="$CXXFLAGS -O3 -pipe -funroll-loops -Wall"
LIBS="$LIBS -lm"
AC_SUBST(CXXFLAGS)
AC_SUBST(LIBS)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
