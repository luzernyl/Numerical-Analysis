#include "Splines.h"
#include "Vector.h"
#include <math.h>
#include <iostream>

double f1(double x)
{
  return 1 / (double)(1 + (25*x*x));
}

double derf1(double x)
{
  return (-50*x) / pow((25*x*x)+1, 2);
}

InterpConditions construct(double N)
{
  InterpConditions data;
  std::vector<double> value;
  for(double x = -1; x <= 1; x += 2.0/(N-1))
  {
    data.input_site(x);
    value.clear();
    value.push_back(f1(x));
    if(x == -1 or (x + 2.0/(N-1) >= 1))
    {
      value.push_back(derf1(x));
    }
    data.input_value(value);
    data.input_der(derf1(x));
  }

  return data;
}

int main(int argc, char* argv[])
{
  InterpConditions data6 = construct(6);
  InterpConditions data11 = construct(11);
  InterpConditions data21 = construct(21);
  
  Spline<1,3,SplineType::ppForm> spl6(interpolate<SplineType::ppForm>(data6, BCType::complete));
  std::cout << "N = 6 : " << std::endl;
  for(int i = 0; i < 5; i++)
  {
    std::cout << spl6.get_poly(i) << std::endl;
  }
  std::cout << std::endl;

  Spline<1,3,SplineType::ppForm> spl11(interpolate<SplineType::ppForm>(data11, BCType::complete));
  std::cout << "N = 11 : " << std::endl;
  for(int i = 0; i < 10; i++)
  {
    std::cout << spl11.get_poly(i) << std::endl;
  }
  std::cout << std::endl;

  Spline<1,3,SplineType::ppForm> spl21(interpolate<SplineType::ppForm>(data21, BCType::complete));
  std::cout << "N = 21 : " << std::endl;
  for(int i = 0; i < 20; i++)
  {
    std::cout << spl21.get_poly(i) << std::endl;
  }
  return 0;
}
