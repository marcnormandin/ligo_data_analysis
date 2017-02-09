# LIGO Data Analysis
# by Marc Normandin.

# ****************************************************************************************************
AC_INIT([lda],[0.1],[normandin.utb@gmail.com])
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_LIBTOOL
AC_PROG_CC([cc])
AM_INIT_AUTOMAKE()

# ****************************************************************************************************
# Google test library
# We can't use the standard AC_CHECK_LIB([gtest],[some_fun]) because "some_fun" will have C++
# name-mangling, and we need the mangled name. Instead we check for the main()
# in the accompanied gtest_main library. If we have gtest_main, we should also have gtest.
# Note: AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
#AC_CHECK_LIB([gtest_main], [main], [HAVE_GTEST=1] [TEST_LIBS="$TEST_LIBS -lgtest_main"], AC_MSG_WARN([libgtest is not installed.]))
#AM_CONDITIONAL([HAVE_GTEST], [test "$HAVE_GTEST" = 1])
#AC_SUBST([HAVE_GTEST])

# ****************************************************************************************************
# C/C++ math library
AC_CHECK_LIB([m],[cos])

# ****************************************************************************************************
# GNU Scientific Library
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

# ****************************************************************************************************
# FFTW3 library
#AC_CHECK_LIB([fftw3],[fftw_plan_dft_1d],[], [AC_MSG_ERROR([fftw_plan_dft_1d was not found in libfftw3])])

AC_TYPE_SIZE_T
AC_SUBST([TEST_LIBS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST(ac_top_srcdir)



# ****************************************************************************************************
# List of Makefile.am's
AC_CONFIG_SRCDIR([src_c/main.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src_c/Makefile])

# ****************************************************************************************************
# Done!
AC_OUTPUT
