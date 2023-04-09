CXX        = g++
CXXFLAGS   = -O3 -mavx2 -msse3
SOURCES    = mandelbrot.cpp
SFML	   = -lsfml-graphics -lsfml-window -lsfml-system
OBJECTS    = $(SOURCES:.cpp=.o)
EXECUTABLE = mandelbrot

.PHONY: compile
compile: $(SOURCES) $(EXECUTABLE)
				./mandelbrot

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(SFML)

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

.PHONY: run
run:
	./mandelbrot

.PHONY: clear
clear:
	rm -rf $(OBJECTS) $(EXECUTABLE)

clrObj:
	rm -rf $(OBJECTS)
