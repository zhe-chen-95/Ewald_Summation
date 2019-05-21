CXX = nvcc
CXXFLAGS = -std=c++11
CXXFLAGSDEBUG = -std=c++11 -O0 -g



TARGETS = $(basename $(wildcard *.cu *.cpp))

all : main

main:main.cpp *.h *.cu
	$(CXX) $(CXXFLAGS) -Xcompiler "-O3 -march=native -fopenmp" $< $(LIBS) ewald.cu global.cpp -lcufft -o $@

clean:
	-$(RM) $(TARGETS) *~

.PHONY: all, clean
