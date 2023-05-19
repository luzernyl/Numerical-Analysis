#include "HermiteVec.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>

double f1(double x)
{
  return 1 / (1+x*x);
}

double f(double x)
{
  return 1 / (double)(1 + (25*x*x));
}

int main(int argc, char* argv[])
{
  std::cout << "Problem B :" << std::endl;
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
  std::cout << "n = 2 : " << std::endl;
  poly.print();
  std::cout << std::endl << std::endl;
  /*
  for(int j = 0; j <= 2; j++)
  {
    double x = -5 + (10*j/2);
    std::cout << poly.eval(x) << std::endl;
  }
  */

  NewtonInterp result4(data4);
  Polynomial poly1(result4.NewtonFormula());
  std::cout << "n = 4 : " << std::endl;
  poly1.print();
  std::cout << std::endl << std::endl;
  /*
  for(int j = 0; j <= 4; j++)
  {
    double x = -5 + (10*j/4);
    std::cout << poly1.eval(x) << std::endl;
  }
  */
  
  NewtonInterp result6(data6);
  Polynomial poly2(result6.NewtonFormula());
  std::cout << "n = 6 : " << std::endl;
  poly2.print();
  std::cout << std::endl << std::endl;
  /*
  for(int j = 0; j <= 6; j++)
  {
    double x = -5 + (10*j/6);
    std::cout << poly2.eval(x) << std::endl;
  }
  */
  
  NewtonInterp result8(data8);
  Polynomial poly3(result8.NewtonFormula());
  std::cout << "n = 8 : " << std::endl;
  poly3.print();
  std::cout << std::endl << std::endl;
  /*
  for(int j = 0; j <= 8; j++)
  {
    double x = -5 + (10*j/8);
    std::cout << poly3.eval(x) << std::endl;
  }
  */

  std::cout << "Problem C :" << std::endl;
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
  std::cout << "n = 5 : " << std::endl;
  poly5.print();
  std::cout << std::endl << std::endl;
  /*
  for(int i = 1; i <= 5; i++)
  {
    double x = cos((double)(2*i-1)*M_PI / (double)(2*5));
    std::cout << poly5.eval(x) << std::endl;
  }
  std::cout << std::endl;
  */

  NewtonInterp result10(data10);
  Polynomial poly10(result10.NewtonFormula());
  std::cout << "n = 10 : " << std::endl;
  poly10.print();
  std::cout << std::endl << std::endl;
  /*
  for(int i = 1; i <= 10; i++)
  {
    double x = cos((double)(2*i-1)*M_PI / (double)(2*10));
    std::cout << poly5.eval(x) << std::endl;
  }
  std::cout << std::endl;
  */

  NewtonInterp result15(data15);
  Polynomial poly15(result15.NewtonFormula());
  std::cout << "n = 15 : " << std::endl;
  poly15.print();
  std::cout << std::endl << std::endl;
  /*
  for(int i = 1; i <= 15; i++)
  {
    double x = cos((double)(2*i-1)*M_PI / (double)(2*15));
    std::cout << poly5.eval(x) << std::endl;
  }
  std::cout << std::endl;
  */

  NewtonInterp result20(data20);
  Polynomial poly20(result20.NewtonFormula());
  std::cout << "n = 20 : " << std::endl;
  poly20.print();
  std::cout << std::endl << std::endl;
  /*
  for(int i = 1; i <= 20; i++)
  {
    double x = cos((double)(2*i-1)*M_PI / (double)(2*20));
    std::cout << poly5.eval(x) << std::endl;
  }
  */

  std::cout << "Problem D :" << std::endl;
  std::vector<double> time{0, 3, 5, 8, 13};
  std::vector<double> dist{0, 225, 383, 623, 993};
  std::vector<double> speed{75, 77, 80, 74, 72};
  InterpConditions data(time, dist, speed);

  NewtonInterp result(data); 
  Polynomial positionpoly(result.NewtonFormula());
  Polynomial speedpoly(positionpoly.derivative());
  std::cout << "Position polynomial : " << std::endl;
  positionpoly.print();
  std::cout << std::endl << std::endl;
  std::cout << "Speed polynomial : " << std::endl;
  speedpoly.print();
  std::cout << std::endl << std::endl;
  
  std::cout << "Position when t = 10s : " << positionpoly.eval(10) << std::endl;
  std::cout << "Speed when t = 10s : " << speedpoly.eval(10) << std::endl;
  std::cout << std::endl;
  for(int i = 0; i<=13; i++)
  {
    std::cout << speedpoly.eval(i) << " ";
  }
  std::cout << std::endl << std::endl;

  std::cout << "Problem E :" << std::endl;
  std::vector<double> day{0, 6, 10, 13, 17, 20, 28};
  std::vector<double> sp1{6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7};
  std::vector<double> sp2{6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89};
  std::vector<double> zero(7);

  InterpConditions datasp1(day, sp1, zero);
  InterpConditions datasp2(day, sp2, zero);
  
  NewtonInterp resultsp1(datasp1), resultsp2(datasp2);
  Polynomial polysp1(resultsp1.NewtonFormula());
  Polynomial polysp2(resultsp2.NewtonFormula());
  std::cout << "Sample 1 : " << std::endl;
  polysp1.print();
  std::cout << std::endl << std::endl;
  std::cout << "Sample 2 : " << std::endl;
  polysp2.print();
  std::cout << std::endl << std::endl;
  /*
  for(int i = 0; i < 45; i++)
    {
      std::cout << polysp1.eval(i) << " ";
    }
  std::cout << std::endl << std::endl;
  for(int i = 0; i < 45; i++)
    {
      std::cout << polysp2.eval(i) << " ";
    }
  std::cout << std::endl << std::endl;
  */
  std::cout << "After another 15 days : " << std::endl;
  std::cout << "Sample 1 : " << polysp1.eval(43) << std::endl;
  std::cout << "Sample 2 : " << polysp2.eval(43) << std::endl;
  return 0;
}
