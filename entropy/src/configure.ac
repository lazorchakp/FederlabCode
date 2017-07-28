dnl $ autoreconf --install to compile configure
dnl To configure and make entropy
dnl  ./configure
dnl  make
dnl  make install


  AC_INIT([entropy], [1.1])

  dnl beings automake
  AM_INIT_AUTOMAKE([-Wall foreign])

  AC_CONFIG_HEADERS([config.h])

dnl according to Appendix C in the GSl Reference Manual all three of these lines are needed in this order to properly check for GSL
echo "Testing  for gsl"

 sinclude(gsl.m4)
 AX_PATH_GSL
 LIBS="$GSL_LIBS $LIBS"
 CFLAGS="$GSL_CFLAGS $CFLAGS"

dnl     AC_CHECK_LIB([m],[cos])
dnl     AC_CHECK_LIB([gslcblas],[cblas_dgemm])
dnl     AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

dnl   AC_SEARCH_LIBS([cos],[m], [], [AC_MSG_ERROR([unable to find the cos() function])])
  AC_SEARCH_LIBS([cblas_dgemm],[gslcblas m], [], [AC_MSG_ERROR([cblas_dgemm not found: exiting configure])])
  AC_SEARCH_LIBS([gsl_blas_dgemm], [gsl m], [], [AC_MSG_ERROR([gsl_blas_dgemm not found: exiting configure])])

echo "Testing  for a C++ compiler"
  AC_PROG_CXX([h5c++])
  AC_PROG_CXX_C_O
  AC_PROG_CXXCPP

echo "Setting  language to C++"
  AC_LANG(C++)

echo "Setting  CXXFLAGS -DGSL_RANGE_CHECK=0"
  CXXFLAGS="$CXXFLAGS -DGSL_RANGE_CHECK=0"

echo "Testing  for iostream"
  AC_CHECK_HEADERS(iostream)
echo "Testing  for iomanip"
  AC_CHECK_HEADERS(iomanip)
echo "Testing  for cmath"
  AC_CHECK_HEADERS(cmath)
echo "Testing  for fstream"
  AC_CHECK_HEADERS(fstream)
echo "Testing  for sstream"
  AC_CHECK_HEADERS(sstream)
echo "Testing  for functions"
  AC_CHECK_FUNCS([pow])
  AC_CHECK_FUNCS([sqrt])
echo "Testing  for float.h, sting.h, & unistd.hrm "
  AC_CHECK_HEADERS([float.h string.h unistd.h])
echo "Testing  for STDBOOL"
  AC_HEADER_STDBOOL
echo "Testing  for TYPE_SIZE_T"
  AC_TYPE_SIZE_T

  AC_PREREQ


echo "Creating Makefile"
  AC_CONFIG_FILES([Makefile])

echo "Creating config.h, config.h.in, and config.status"
  dnl AC_OUTPUT creates files config.h and config.h.in, and config.status 
  dnl which are used by automake
  AC_OUTPUT

