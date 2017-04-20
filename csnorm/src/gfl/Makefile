CC=gcc
NVCC=nvcc
CFLAGS=-O3 -Wall -Wextra -ansi -std=c99 -pedantic -fopenmp
OBJ=obj/*.o
IDIR=../include/

all:
	#mkdir -p lib
	#mkdir -p obj
	cd src; ${CC} ${CFLAGS}  -c -fPIC *.c -I${IDIR}
	#${CC} -shared -o lib/libgraphfl.so ${OBJ} -lgsl -lgslcblas
	
docs:
	doxygen Doxyfile

clean:
	#rm -f obj/*.o
	rm -f src/*.o

cleanall:
	#rm -f obj/*.o
	rm -f src/*.o
	#rm lib/*.so
