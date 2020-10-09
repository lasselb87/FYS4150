CC = gcc
CXX = g++
CPPflags = -g -Wall -std=c++11 -O3

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
INCL := -I /usr/local/Cellar/armadillo/9.100.5/include
LIB += -DARMA_DONT_USE_WRAPPER -framework Accelerate
endif

ifeq ($(UNAME), Linux)
INCL :=
LIB := -DARMA_DONT_USE_WRAPPER -lblas -llapack
endif

OBJS = planet.o solver.o
TARGET1 = twoBody.x
TARGET2 = allPlanets.x
TARGET3 = benchmark.x
TARGET4 = escapeVelocity.x
TARGET5 = threeBody.x
TARGET6 = relativistic.x
TARGET7 = test.x

TARGETS = $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7)

all: $(TARGETS)

twoBody.x: twoBody.o $(OBJS)
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

twoBody.o: twoBody.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

benchmark.x: benchmark.o
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

benchmark.o: benchmark.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

allPlanets.x: allPlanets.o $(OBJS)
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

allPlanets.o: allPlanets.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

escapeVelocity.x: escapeVelocity.o $(OBJS)
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

escapeVelocity.o: escapeVelocity.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

threeBody.x: threeBody.o $(OBJS)
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

threeBody.o: threeBody.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

relativistic.x: relativistic.o $(OBJS)
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

relativistic.o: relativistic.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

solver.o:
	${CXX} ${CPPflags} -c ./Classes/solver.cpp ${INCL} ${LIB}

planet.o:
	${CXX} ${CPPflags} -c ./Classes/planet.cpp ${INCL} ${LIB}

test.x: test.o $(OBJS)
		${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}
		./test.x
		rm -f test.o test.x

test.o: test.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

clean:
	rm -f *.o *.txt *.x

cleanRaw:
	rm -f ./Raw_data/*.txt

cleanResults:
	rm -f ./Results/*.txt ./Results/*.png ./Results/*.pdf
