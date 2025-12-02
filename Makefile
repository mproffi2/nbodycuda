# Makefile for N-Body Simulation (CPU + optional GPU)

# Detect if nbody.cu exists; if not, use nbody.cpp
ifeq ($(wildcard nbody.cu),)
  SRC = nbody.cpp
  USE_NVCC = false
else
  SRC = nbody.cu
  USE_NVCC = true
endif

CXX      = g++
CXXFLAGS = -O3 -std=c++17
NVCC     = nvcc
NVCCFLAGS= -O3 -std=c++17 -arch=sm_61

TARGET   = nbody

all: $(TARGET)

$(TARGET):
ifeq ($(USE_NVCC),true)
	$(NVCC) $(NVCCFLAGS) -o $(TARGET) $(SRC)
else
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)
endif

clean:
	rm -f $(TARGET) *.o