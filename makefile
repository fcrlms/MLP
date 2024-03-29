CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Werror -pedantic -O3 -ffast-math

SOURCES = $(wildcard src/*.cc)
OBJECTS = $(patsubst src/%.cc, obj/%.o, $(SOURCES))
HEADERS = $(wildcard src/*.h)

EXECUTABLE = mlp

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)

$(OBJECTS): obj/%.o : src/%.cc | obj
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj:
	mkdir -p $@
