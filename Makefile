# Makefile for N-Body Simulation (CPU + GPU)

NVCC        = nvcc
CXX         = g++
CXXFLAGS    = -O3 -std=c++17
NVCCFLAGS   = -O3 -std=c++17 -arch=sm_61

TARGET      = nbody
SRC         = nbody.cu

all: $(TARGET)

$(TARGET): $(SRC)
	$(NVCC) $(NVCCFLAGS) -o $(TARGET) $(SRC)

cpu: $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET)_cpu $(SRC)

clean:
	rm -f $(TARGET) $(TARGET)_cpu *.o