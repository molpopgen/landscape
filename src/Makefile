CXX=c++
CXXFLAGS=-std=c++11 -O2 -Wall -W

all: rtree_example.o wflandscape.o
	$(CXX) $(CXXFLAGS) -o rtree_example rtree_example.o -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS) -o wflandscape wflandscape.o -lgsl -lgslcblas

clean:
	rm -f *.o

wflandscape.o: simtypes.hpp wfrules.hpp
