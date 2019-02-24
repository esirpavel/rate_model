all: cycover_ring.so

cycover_ring.so: ring.o cycover_ring.o
	g++ -O3 -shared cycover_ring.o ring.o -lblas -o $@

cycover_ring.o: cycover_ring.cpp
	g++ -c -Wno-cpp -fPIC -O3 `python3-config --includes` cycover_ring.cpp -o $@

cycover_ring.cpp: cycover_ring.pyx
	cython -3 --cplus cycover_ring.pyx -o $@

ring.o: ring.cpp ring.h
	g++ ring.cpp -O3 -c -fPIC -o $@

clean:
	rm cycover_ring.cpp ring.o cycover_ring.o cycover_ring.so
