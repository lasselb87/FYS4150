CC = gcc
CXX = g++
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

main.o: main.cpp
	${CXX} ${CPPflags} -c $^

main.x: main.o VMC.o
	${CXX} ${CPPflags} -o $@ $^

VMC.o: classes/VMC.cpp
	${CXX} ${CPPflags} -c $^

test.x: test.o VMC.o
	g++ ${CPPflags} -o $@ $^ ${INCL} ${LIB}
	./test.x

test.o: classes/test.cpp
	g++ ${CPPflags} -c $^ ${INCL} ${LIB}

clean:
	rm -f *.o *.txt *.x
