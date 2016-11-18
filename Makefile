
#CC = xlc_r
CPP = g++
CFLAGS = -O2 -Wno-deprecated  #-qsmp #-O2 #-qarch=auto -qtune=auto
BOOST = -I /vbi/packages/boost/include/boost-1_38

all: ashuffle

ashuffle: ashuffle.cpp utility.hpp Timer.hpp Serial_Multinomial.hpp
	$(CPP) $(CFLAGS) ashuffle.cpp -o as

clean:
	rm -f as
	