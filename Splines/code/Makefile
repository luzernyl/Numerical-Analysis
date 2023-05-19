output : B.o 
	g++ -o B B.o -std=c++14

B.o : B.cpp Splines.h Splines.o Vector.h
	g++ -c B.cpp -llapacke -llapack -std=c++14

Splines.o : Splines.cpp Splines.h Vector.h
	g++ -c Splines.cpp -llapacke -llapack -std=c++14

clean :
	rm *.o output


