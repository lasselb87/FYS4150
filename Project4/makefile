CC = gcc
CXX = mpic++
CPPflags = -g -Wall -std=c++11 -O3

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
INCL := -I /usr/local/Cellar/armadillo/9.100.5/include
LIB := -DARMA_DONT_USE_WRAPPER -framework Accelerate
endif

ifeq ($(UNAME), Linux)
INCL :=
LIB := -DARMA_DONT_USE_WRAPPER -lblas -llapack
endif

simulation.x: simulation.o ising.o metropolis.o
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

simulation.o: simulation.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

ising.o: Classes/ising.cpp
	g++ ${CPPflags} -c $^ ${INCL} ${LIB}

metropolis.o: Classes/metropolis.cpp
	g++ ${CPPflags} -c $^ ${INCL} ${LIB}

test.x: test.o ising.o metropolis.o
	g++ ${CPPflags} -o $@ $^ ${INCL} ${LIB}
	./test.x

test.o: Classes/test.cpp
	g++ ${CPPflags} -c $^ ${INCL} ${LIB}

clean:
	rm -f *.o *.txt *.x
