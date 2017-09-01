# LIGO Data Analysis
# by Marc Normandin.

# ****************************************************************************************************
AC_INIT([lda],[2.5-production],[normandin.utb@gmail.com])
AC_CONFIG_MACRO_DIR([m4])
#AC_PROG_LIBTOOL
AC_PROG_CC
AC_PROG_CXX

#AC_PROG_CXX([mpicxx])
AC_PROG_CC([gcc-6])
AC_PROG_CXX([g++-6])
#AC_PROG_CC([icc])
#AC_PROG_CXX([icc])
#AM_PROG_CC_C_O
AM_INIT_AUTOMAKE()
AC_PROG_LIBTOOL

AC_OPENMP()
if test "$enable_openmp" = yes; then
	AC_OPENMP()
	if test "x$ac_cv_prog_cc_openmp" = "xunsupported"; then
		AC_MSG_NOTICE(OpenMP requested, but is not supported. Multithreaded operation is disabled.)
	else
		CFLAGS="$CFLAGS $OPENMP_CFLAGS"
		AC_DEFINE(HAVE_OPENMP, 1, Define to 1 if you have OpenMP)
	fi
else
	AC_MSG_NOTICE(OpenMP disabled.)
fi

#AX_OPENMP()
#AC_SUBST([AM_CPPFLAGS], ["$CFLAGS $OPENMP_CFLAGS"])
AC_SUBST([OPENMP_CFLAGS])
AM_CONDITIONAL([HAVE_OPENMP], [test "x$ac_cv_prog_cc_openmp" != "xunsupported"])
#AC_DEFINE([HAVE_OPENMP], ["$HAVE_OPENMP"], [Have OpenMP])
#AC_SUBST([HAVE_OPENMP])

#AC_DEFINE([HAVE_OMP], [], [Have OpenMP])
#AX_MPI()
#AC_DEFINE([HAVE_MPI], [1], [Enable MPI usage.])
 # If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
  # If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
  # If --with-mpi=no is used, use a standard C compiler instead.

#AC_DEFINE([HAVE_MPI], [1], [Enable MPI])

#AC_SUBST([HAVE_MPI])
#AM_CONDITIONAL([HAVE_MPI], [test x"$use_mpi" = xyes], [MPI Enabled])
#AM_CONDITIONAL([HAVE_MPI], [test ])

#AX_LIB_HDF5([serial])
#if test "$with_hdf5" = "no"; then
#	AC_MSG_ERROR([Unable to find HDF5, which is required.])
#fi

# ****************************************************************************************************
# Google test library
# We can't use the standard AC_CHECK_LIB([gtest],[some_fun]) because "some_fun" will have C++
# name-mangling, and we need the mangled name. Instead we check for the main()
# in the accompanied gtest_main library. If we have gtest_main, we should also have gtest.
# Note: AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
AC_CHECK_LIB([gtest_main], [main], [HAVE_GTEST=yes] [TEST_LIBS="$TEST_LIBS -lgtest_main"], AC_MSG_WARN([libgtest is not installed.]))
AM_CONDITIONAL([HAVE_GTEST], [test "x$HAVE_GTEST" = xyes])
AC_SUBST([HAVE_GTEST])
AC_DEFINE([HAVE_GTEST], [], [Have GoogleTest Framework])

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
AC_CONFIG_SRCDIR([libcore/detector_antenna_patterns.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
#AC_CONFIG_FILES([programs/histogram/Makefile])
AC_CONFIG_FILES([settings/Makefile])
AC_CONFIG_FILES([libcore/Makefile])
AC_CONFIG_FILES([libpso/Makefile])
AC_CONFIG_FILES([programs/Makefile])
#AC_CONFIG_FILES([programs/matlab_data_mpi/Makefile])
AC_CONFIG_FILES([programs/matlab_data_serial/Makefile])
AC_CONFIG_FILES([programs/simulate_data/Makefile])
#AC_CONFIG_FILES([programs/simulate_matlab_data/Makefile])
#AC_CONFIG_FILES([programs/diagnostics/Makefile])
AC_CONFIG_FILES([tests/Makefile])

# ****************************************************************************************************
# Done!
AC_OUTPUT
