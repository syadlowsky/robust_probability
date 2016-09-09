CXX=g++
CPPFLAGS:= -O3 -msse3 -fPIC -fopenmp -march=native -shared -std=c++0x -pedantic -Wall -Wcast-qual -mavx -g3 -I /usr/include/eigen3
RM= rm -f
LDFLAGS= -shared
.PHONY: all clean


all: robust_probability.so

robust_probability.so: robust_probability.o
	$(CXX) -shared $< -o $@

clean:
	rm -rf *.so
	rm -rf *.o



