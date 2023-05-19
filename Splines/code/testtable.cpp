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


int main(int argc, char* argv[])
{
  /*
  double N = 6;
  for(double i = -1; i <= 1; i+= 2.0/(N-1))
  {
    data.input_site(i);
    data.input_value(f1(i));
    std::vector<double> df{derf1(i)};
    data.input_der(df);
  }
  */
  std::vector<double> site{1,2,3,4,6};
  std::vector<std::vector<double> > value;
  std::vector<double> der;
  value.clear();
  der.clear();
  for(double i = 0; i < site.size(); i++)
  {
    std::vector<double > t {log(site[i])};
    if(i == 0 || i == site.size()-1)
      t.push_back((double)(1/site[i]));
    value.push_back(t);
    der.push_back((double)(1/site[i]));
  }
  InterpConditions data(site,value,der);
  Spline<1,3,SplineType::ppForm> sp;
  sp.make_table(data);
  sp.print_table();
  std::cout << std::endl;
  Spline<1,3,SplineType::ppForm> spl(interpolate<SplineType::ppForm>(data, BCType::complete));
  for(int i = 0; i < site.size() - 1; i++)
  {
    std::cout << spl.get_poly(i) << std::endl;
  }
  std::cout << spl.get_poly(3).eval(5) << std::endl;
  return 0;
}
