#include "HermiteVec.h"
#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>

int main(int argc, char* argv[])
{
  std::vector<double> day{0, 6, 10, 13, 17, 20, 28};
  std::vector<double> sp1{6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7};
  std::vector<double> sp2{6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89};
  std::vector<double> zero(7);

  InterpConditions data1(day, sp1, zero);
  InterpConditions data2(day, sp2, zero);
  
  NewtonInterp result1(data1), result2(data2);
  Polynomial poly1(result1.NewtonFormula());
  Polynomial poly2(result2.NewtonFormula());
  poly1.print();
  std::cout << std::endl << std::endl;
  poly2.print();
  std::cout << std::endl << std::endl;
  for(int i = 0; i < 45; i++)
    {
      std::cout << poly1.eval(i) << " ";
    }
  std::cout << std::endl << std::endl;
  for(int i = 0; i < 45; i++)
    {
      std::cout << poly2.eval(i) << " ";
    }
  std::cout << std::endl << std::endl;
  std::cout << "After another 15 days : " << std::endl;
  std::cout << "Sample 1 : " << poly1.eval(43) << std::endl;
  std::cout << "Sample 2 : " << poly2.eval(43) << std::endl;
  
  return 0;
}
