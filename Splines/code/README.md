Structure of the programs:
Splines.h contains the declaration of classes and functions, while Splines.cpp is the implementation. Vector.h is the header file for vectors.
Compile using "g++ -c Splines.h -llapacke -llapack -std=c++14"

testpoly.cpp is the testing program for the class Polynomial
Compile using "g++ -o testpoly testpoly.cpp -std=c++14"

testvector.cpp is the testing program for the class Vec in Vector.h
Compile using "g++ -o testvector testvector.cpp -std=c++14"

B.cpp is the program for Problem B
Compile using "g++ -o B B.cpp -llapacke -llapack -std=c++14"

