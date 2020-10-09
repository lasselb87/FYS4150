CPPflags = g++ -std=c++11 -O3

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
INCLUDES := -I /usr/local/Cellar/armadillo/9.100.5/include
LIB := -DARMA_DONT_USE_WRAPPER -framework Accelerate
endif

ifeq ($(UNAME), Linux)
LIB := -DARMA_DONT_USE_WRAPPER -lblas -llapack
endif

TARGETS = bucklingBeam quantumDot waveFunctions benchmark test

all: $(TARGETS)

bucklingBeam: bucklingBeam.o jacobi.o
	${CPPflags} -o bucklingBeam.x bucklingBeam.cpp jacobi.cpp ${LIB}
	rm -f *.o

quantumDot: quantumDot.o jacobi.o
	${CPPflags} -o quantumDot.x quantumDot.cpp jacobi.cpp ${LIB}
	rm -f *.o

waveFunctions: waveFunctions.o jacobi.o
	${CPPflags} -o waveFunctions.x waveFunctions.cpp jacobi.cpp ${LIB}
	rm -f *.o

benchmark: benchmark.o jacobi.o
	${CPPflags} -o benchmark.x benchmark.cpp jacobi.cpp ${LIB}
	rm -f *.o
	./benchmark.x
	python plot.py benchmark

test: test.o jacobi.o
	${CPPflags} -o test.x test.cpp jacobi.cpp ${LIB}
	rm -f *.o
	./test.x
	rm -f ./test.x

clean:
	rm -f *.o *.txt *.x

cleanResults:
	rm -f ./Results/*.txt ./Results/*.png
