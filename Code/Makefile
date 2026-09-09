CXX = g++
CXXFLAGS = -O2 -std=c++11 -Wall -Wextra -pedantic
INCLUDES = -isystem /usr/local/include/CCfits
LIBS = -lcfitsio -lCCfits -L/usr/local/lib 

COMMON_SOURCES = \
    Voxelise.cpp \
    SmoothingKernel.cpp \
    CubicSplineSmoothingKernel.cpp \
    QuarticSplineSmoothingKernel.cpp \
    ScaledGaussianSmoothingKernel.cpp

COMMON_OBJECTS = $(COMMON_SOURCES:.cpp=.o)

all: VoxeliseSimulation CalculateCumulativeKernel

VoxeliseSimulation: VoxeliseSimulation.o $(COMMON_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

CalculateCumulativeKernel: CalculateCumulativeKernel.o $(COMMON_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f *.o VoxeliseSimulation CalculateCumulativeKernel
