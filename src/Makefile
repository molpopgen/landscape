CXXFLAGS=-std=c++11 -O2 -Wall -W -DNDEBUG

all: rtree_example.o rtree_wtf.o wflandscape.o wflandscape_timing.o
	$(CXX) $(CXXFLAGS) -o rtree_example rtree_example.o -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS) -o rtree_wtf rtree_wtf.o -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS) -o wflandscape wflandscape.o -lgsl -lgslcblas -lsequence -lpthread
	$(CXX) $(CXXFLAGS) -o wflandscape_timing wflandscape_timing.o -lgsl -lgslcblas -lsequence -lpthread

clean:
	rm -f *.o

wflandscape.o: simtypes.hpp wfrules.hpp
