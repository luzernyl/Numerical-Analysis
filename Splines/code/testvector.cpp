#include "Vector.h"
#include <iostream>
#include <vector>
#include <initializer_list>

int main(int argc, char* argv[])
{
  std::initializer_list<double> l{1.5,2.7,3};
  std::initializer_list<double> l1{2, 2, 2};
  Vec<double,3> vec(l);
  Vec<double,3> vec1(l1);
  std::cout << vec << std::endl;
  std::cout << vec1 << std::endl;
  //Vec<int,2> vec1(vec);
  //std::cout << vec1 << std::endl;
  //std::cout << vec.unit(1) << std::endl;
  std::cout << vec-vec1 << std::endl;
  std::cout << vec*vec1 << std::endl;

  Vec<double,1> minus(-1);
  std::cout << minus << std::endl;
  std::cout << vec*(-1) << std::endl;
  return 0;
}
