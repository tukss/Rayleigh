#!/bin/csh
###  RAYLEIGH BUILD SCRIPT
###  Usage:  ./build_rayleigh machine_name
###  The Makefile employed is machine dependent
###  See the Makefiles directory for example Makefiles
###  The value of machine_name should correspond to one of the Makefiles
###  in that directory.  Machine name for Makefile is indicated following underscore.

###  For example, to use Makefile_LCD
###  ./build_rayleigh LCD

###  To use Makefile_Pleiades
###  ./build_rayleigh Pleiades

###  A debugging flag may be passed to makefiles that support it.
### ./build_rayleigh Pleiades debug  (for example)

setenv MACHINE $1
setenv RAYLEIGH_OPT1 $2
setenv RAYLEIGH_OPT2 $3
setenv RAYLEIGH_OPT3 $4
mkdir build
cd build
rm *.F90
rm Makefile
cp ../parallel_framework/*.F90 .
cp ../data_structures/*.F90 .
cp ../math_layer/*.F90 .
cp ../IO/*.F90 .
cp ../test_suite/*.F90 .
cp ../physics/*.F90 .
cp ../Utility/*.F90 .
cp ../Utility/*.c .
cp ../Makefiles/Makefile_$MACHINE Makefile
cp ../Makefiles/object_list .
if($MACHINE == "Mira") then
cp ../Utility/MakeDir.F90_IBM MakeDir.F90
if ($RAYLEIGH_OPT1 == "omp") then
	echo Builing OpenMP-version-Mira
	cp ../math_layer/Fourier_Transform.mira_omp Fourier_Transform.F90
	cp ../math_layer/Chebyshev_Polynomials.mira_omp Chebyshev_Polynomials.F90
	cp ../parallel_framework/Parallel_Framework.mira_omp Parallel_Framework.F90
endif
endif
if($MACHINE == "Fermi") then
cp ../Utility/MakeDir.F90_IBM MakeDir.F90
if ($RAYLEIGH_OPT1 == "omp") then
	echo Builing OpenMP-version-Mira
	cp ../math_layer/Fourier_Transform.mira_omp Fourier_Transform.F90
	cp ../math_layer/Chebyshev_Polynomials.mira_omp Chebyshev_Polynomials.F90
	cp ../parallel_framework/Parallel_Framework.mira_omp Parallel_Framework.F90
endif
endif

make clean
if ($RAYLEIGH_OPT1 == "mic") then
make rayleigh.mic
else
make all
endif
