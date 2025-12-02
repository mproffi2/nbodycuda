# Makefile for N-Body Simulation (auto-detect CUDA/CPU)

TARGET = nbody

# Detect source file
ifeq ($(wildcard nbody.cu), nbody.cu)
    SRC = nbody.cu
    NVCC = nvcc
    NVCCFLAGS = -O3 -std=c++17 -arch=sm_61
    BUILD_CMD = $(NVCC) $(NVCCFLAGS) -o $(TARGET) $(SRC)
else ifeq ($(wildcard nbody.cpp), nbody.cpp)
    SRC = nbody.cpp
    CXX = g++
    CXXFLAGS = -O3 -std=c++17
    BUILD_CMD = $(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)
else
    $(error "No nbody source file found (nbody.cu or nbody.cpp)")
endif

all:
	$(BUILD_CMD)

clean:
	rm -f $(TARGET) *.o