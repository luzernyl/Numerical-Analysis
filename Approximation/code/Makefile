all : A B

A : A.o Approx.o
	g++ -o A A.o -std=c++11 -llapacke -llapack -lcblas -lblas

B : B.o Approx.o
	g++ -o B B.o -std=c++11 -llapacke -llapack -lcblas -lblas

A.o : A.cpp Approx.cpp Approx.h
	g++ -c A.cpp -std=c++11 -llapacke -llapack -lcblas -lblas

B.o : B.cpp Approx.cpp Approx.h
	g++ -c B.cpp -std=c++11 -llapacke -llapack -lcblas -lblas

Approx.o : Approx.cpp Approx.h
	g++ -c Approx.cpp -std=c++11 -llapacke -llapack -lcblas -lblas

clean :
	rm A.o B.o A B


