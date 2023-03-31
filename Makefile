all:
	g++ -c mandelbrot.cpp -mavx2 -o main.o
	g++ main.o -o mandelbrot -lsfml-graphics -lsfml-window -lsfml-system -lsfml-audio
	./mandelbrot

clear:
	rm -f *.o