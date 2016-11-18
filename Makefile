#CC = xlc_r
CPP = g++
CFLAGS = -O2 -Wno-deprecated  #-qsmp #-O2 #-qarch=auto -qtune=auto  

all: plshuffle
	
plshuffle: Par-Lab-Shuffle-Com.cpp Par-Lab-Shuffle-Com.hpp Timer.hpp Serial_Multinomial.hpp utility.hpp utility2.hpp TSet.hpp parmergesort.hpp
	module load mpi/mvapich2/gcc/4.7.2/1.9-psm;	mpicxx -o plshufflec Par-Lab-Shuffle-Com.cpp

clean:
	rm -f plshufflec
	

