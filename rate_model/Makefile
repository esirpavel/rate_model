# CXXFLAGS = -O0 -g
CXXFLAGS = -O3

all: cycover_ring.so

cycover_ring.so: ring.o cycover_ring.o
	g++ $(CXXFLAGS) -shared cycover_ring.o ring.o -lblas -o $@

cycover_ring.o: cycover_ring.cpp
	g++ -c -Wno-cpp -fPIC $(CXXFLAGS) `python3-config --includes` -I`python3-config --prefix`/lib/python3.7/site-packages/numpy/core/include/ cycover_ring.cpp -o $@

cycover_ring.cpp: cycover_ring.pyx
	cython -3 --cplus cycover_ring.pyx -o $@

ring.o: ring.cpp ring.h
	g++ ring.cpp $(CXXFLAGS) -c -fPIC -o $@

clean:
	rm cycover_ring.cpp ring.o cycover_ring.o cycover_ring.so
