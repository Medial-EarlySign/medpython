cmake_minimum_required(VERSION 3.5.0)
set(SRC_FILES
	blas/f2c.hpp
	blas.hpp
	earth.hpp
     	blas/daxpy.cpp
	blas/dcopy.cpp
	blas/ddot.cpp
	blas/dnrm2.cpp
	blas/dqrls.cpp
	blas/dqrsl.cpp
	blas/dscal.cpp
	blas/d_sign.cpp
	blas/dtrsl.cpp
	earth.cpp
	R/dqrdc2.cpp

)
add_library(Mars STATIC ${SRC_FILES})
