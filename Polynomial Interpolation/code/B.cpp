#include "HermiteVec.h"
#include <iostream>
#include <vector>
#include <math.h>

double f1(double x)
{
  return 1 / (1+x*x);
}

int main(int argc, char* argv[])
{
  InterpConditions data2, data4, data6, data8;
  
  double n = 2;
  for(int i = 0; i <= n; i++)
  {
    double x = -5 + (10*i/n);
    data2.input_site(x);
    data2.input_value(f1(x));
    data2.input_der(0);
  }
  
  n=4;
  for(int i = 0; i <= n; i++)
  {
    double x = -5 + (10*i/n);
    data4.input_site(x);
    data4.input_value(f1(x));
    data4.input_der(0);
  }

  n = 6;
  for(int i = 0; i <= n; i++)
  {
    double x = -5 + (10*i/n);
    data6.input_site(x);
    data6.input_value(f1(x));
    data6.input_der(0);
  }

  n = 8;
  for(int i = 0; i <= n; i++)
  {
    double x = -5 + (10*i/n);
    data8.input_site(x);
    data8.input_value(f1(x));
    data8.input_der(0);
  }
  
  NewtonInterp result2(data2);
  Polynomial poly(result2.NewtonFormula());
  poly.print();
  std::cout << std::endl;
  for(int j = 0; j <= 2; j++)
  {
    double x = -5 + (10*j/2);
    std::cout << poly.eval(x) << std::endl;
  }

  NewtonInterp result4(data4);
  Polynomial poly1(result4.NewtonFormula());
  poly1.print();
  std::cout << std::endl;
  for(int j = 0; j <= 4; j++)
  {
    double x = -5 + (10*j/4);
    std::cout << poly1.eval(x) << std::endl;
  }
  
  
  NewtonInterp result6(data6);
  Polynomial poly2(result6.NewtonFormula());
  poly2.print();
  std::cout << std::endl;
  for(int j = 0; j <= 6; j++)
  {
    double x = -5 + (10*j/6);
    std::cout << poly2.eval(x) << std::endl;
  }
  
  
  NewtonInterp result8(data8);
  Polynomial poly3(result8.NewtonFormula());
  poly3.print();
  std::cout << std::endl;
  for(int j = 0; j <= 8; j++)
  {
    double x = -5 + (10*j/8);
    std::cout << poly3.eval(x) << std::endl;
  }
return 0;
}
