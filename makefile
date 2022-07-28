CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -Werror -pedantic -O3 -ffast-math

SOURCES = $(wildcard *.cpp)

OBJECTS = $(SOURCES: %.cpp = %.o)

EXECUTABLE = mlp

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)

%.o:
	$(CXX) $(CXXFLAGS) -c $*.cpp
