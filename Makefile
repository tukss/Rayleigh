FFTWINC = /home3/feathern/FFTW3/include
FFTWLIB = /home3/feathern/FFTW3/lib
F90  = /usr/local/mpich2-1.0.7/bin/mpif90
F90FLAGS = -r8 -O3 -fpp -FR -fp-model strict -I$(FFTWINC)

MKLPATH = /central/intel/mkl_10.1.0.015
INTERFACE = ${MKLPATH}/lib/em64t/libmkl_intel_lp64.a
COMP =  ${MKLPATH}/lib/em64t/libmkl_lapack.a
THREADING = ${MKLPATH}/lib/em64t/libmkl_intel_thread.a
RUNTIME = ${MKLPATH}/lib/em64t/libguide.a -lpthread
MKL_START = -L${MKLPATH}/lib/em64t -I${MKLPATH}/include -Wl, --start-group
MKL_END = -l -l -l -Wl, --end-group

#   Flags for the MKL libraries
LIBFLAGS = $(MKL_START) $(INTERFACE) $(THREADING) $(COMP) $(MKL_END) $(RUNTIME) -L/home3/feathern/FFTW3/lib/ -lfftw3 -lm -lc



$LIBFLAGS = -L/home3/feathern/FFTW3/lib/ -lfftw3 -lm -lc 
F90FLAGS = -r8 -O0 -CB -g -traceback -fpp -FR 
##LIBFLAGS = -lmpi -lmpio
PROG =	pseudo

.SUFFIXES: .o .F90 .f .F
OBJ = Finite_Difference.o Linear_Solve.o Legendre_Polynomials.o Legendre_Transforms.o Fourier_Transform.o Fourier_Derivatives.o MPI_Base.o \
	All_to_All.o MPI_LAYER.o ProblemSize.o Input.o Load_Balance.o Parallel_Framework.o \
	PseudoPhysics.o Main.o \

$(PROG) :$(OBJ)
	$(F90) $(F90FLAGS) -o  $(PROG) $(OBJ) $(LIBFLAGS)
.F90.o :
	$(F90) $(F90FLAGS) -c $<
clean : 
	rm -f *.o *.mod pseudo
ready :
	ln -s parallel_framework/*.mod .
	ln -s parallel_framework/*.o .
	ln -s math_layer/fft/*.o .
	ln -s math_layer/fft/*.mod .
	ln -s math_layer/linear_solve/*.mod .
	ln -s math_layer/linear_solve/*.o .
	ln -s math_layer/*.o .
	ln -s math_layer/*.mod .
Main.o : Parallel_Framework.o PseudoPhysics.o

