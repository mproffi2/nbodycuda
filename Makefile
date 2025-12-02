Makefile for N-Body Simulation (auto-detect CUDA/CPU)

TARGET = nbody

Auto-detect source file

ifneq ("$(wildcard nbody.cu)","")
SRC = nbody.cu
NVCC = nvcc
NVCCFLAGS = -O3 -std=c++17 -arch=sm_61
BUILD_CMD = $(NVCC) $(NVCCFLAGS) -o $(TARGET) $(SRC)
else ifneq ("$(wildcard nbody.cpp)","")
SRC = nbody.cpp
CXX = g++
CXXFLAGS = -O3 -std=c++17
BUILD_CMD = $(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)
else
$(error "No nbody source file found (nbody.cu or nbody.cpp)")
endif

all: $(TARGET)

$(TARGET):
$(BUILD_CMD)

clean:
rm -f $(TARGET) *.o *.log