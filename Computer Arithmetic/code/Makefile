output : A.o B.o 
	g++ -o A A.o -std=c++11
	g++ -o B B.o -std=c++11

A.o : A.cpp
	g++ -c A.cpp -std=c++11

B.o : B.cpp
	g++ -c B.cpp -std=c++11
clean :
	rm *.o output


