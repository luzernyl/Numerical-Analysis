#include "HermiteVec.h"
#include <iostream>
#include <math.h>
#include <cmath>

double f(double x)
{
  return 1 / (double)(1 + (25*x*x));
}

int main(int argc, char*argv[])
{
  InterpConditions data5, data10, data15, data20;
  for(int i = 1; i <= 5; i++)
  {
    double x = cos((double)(2*i - 1)*M_PI / (double)(2*5));
    data5.input_site(x);
    data5.input_value(f(x));
    data5.input_der(0);
  }
  for(int i = 1; i <= 10; i++)
  {
    double x = cos((double)(2*i - 1)*M_PI / (double)(2*10));
    data10.input_site(x);
    data10.input_value(f(x));
    data10.input_der(0);
  }
  for(int i = 1; i <= 15; i++)
  {
    double x =  cos((double)(2*i - 1)*M_PI / (double)(2*15));
    data15.input_site(x);
    data15.input_value(f(x));
    data15.input_der(0);
  }
  for(int i = 1; i <= 20; i++)
  {
    double x =  cos((double)(2*i - 1)*M_PI / (double)(2*20));
    data20.input_site(x);
    data20.input_value(f(x));
    data20.input_der(0);
  }

  NewtonInterp result5(data5);
  Polynomial poly5(result5.NewtonFormula());
  poly5.print();
  std::cout << std::endl;
  for(int i = 1; i <= 5; i++)
  {
    double x = cos((double)(2*i-1)*M_PI / (double)(2*5));
    std::cout << poly5.eval(x) << std::endl;
  }
  std::cout << std::endl;

  NewtonInterp result10(data10);
  Polynomial poly10(result10.NewtonFormula());
  poly10.print();
  std::cout << std::endl;
  for(int i = 1; i <= 10; i++)
  {
    double x = cos((double)(2*i-1)*M_PI / (double)(2*10));
    std::cout << poly5.eval(x) << std::endl;
  }
  std::cout << std::endl;

  NewtonInterp result15(data15);
  Polynomial poly15(result15.NewtonFormula());
  poly15.print();
  std::cout << std::endl;
  for(int i = 1; i <= 15; i++)
  {
    double x = cos((double)(2*i-1)*M_PI / (double)(2*15));
    std::cout << poly5.eval(x) << std::endl;
  }
  std::cout << std::endl;

  NewtonInterp result20(data20);
  Polynomial poly20(result20.NewtonFormula());
  poly20.print();
  std::cout << std::endl;
  for(int i = 1; i <= 20; i++)
  {
    double x = cos((double)(2*i-1)*M_PI / (double)(2*20));
    std::cout << poly5.eval(x) << std::endl;
  }
  
  return 0;
}
