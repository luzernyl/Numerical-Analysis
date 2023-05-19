/**
 * @file   NLE_F.cpp
 * @author Luzern Yuven Luis <student@student>
 * @date   Mon Oct 11 12:38:34 2021
 * 
 * @brief  Program for Problem F
 * 
 * 
 */

#include "NLE.h"
#include <iostream>
#include <math.h>
#include <cmath>

/// Function for F (a)
double f1(double alpha)
{
  return 89*sin(11.5)*sin(alpha)*cos(alpha) + 89*cos(11.5)*pow(sin(alpha),2) -
    ( (49+0.5*55)*sin(11.5)*cos(alpha) ) + (0.5*55*tan(11.5)*cos(alpha)) -
     ( (49+0.5*55)*cos(11.5)-(0.5*55) )*sin(alpha);
}

/// Derivative of function for F (a)
double df1(double alpha)
{
  return 89*sin(11.5)*2*cos(2*alpha) + 89*cos(11.5)*sin(2*alpha) +
    ( (49+0.5*55)*sin(11.5)*sin(alpha) ) - (0.5*55*tan(11.5)*sin(alpha)) -
     ( (49+0.5*55)*cos(11.5)-(0.5*55) )*cos(alpha);
}

/// Function for F (b)
double f2(double alpha)
{
  return 89*sin(11.5)*sin(alpha)*cos(alpha) + 89*cos(11.5)*pow(sin(alpha),2) -
    ( (49+0.5*30)*sin(11.5)*cos(alpha) ) + (0.5*30*tan(11.5)*cos(alpha)) -
     ( (49+0.5*30)*cos(11.5)-(0.5*30) )*sin(alpha);
}

/// Derivative of function for F (b)
double df2(double alpha)
{
  return 89*sin(11.5)*2*cos(2*alpha) + 89*cos(11.5)*sin(2*alpha) +
    (49+0.5*30)*sin(11.5)*sin(alpha)-(0.5*30*tan(11.5)*sin(alpha)) -
     ( (49+0.5*30)*cos(11.5)-(0.5*30) )*cos(alpha);
}

int main(int args, char* argv[])
{
  EquationSolver<double> *nptr1, *nptr2;
  Newton<double> sol1,sol2;
  sol1.input(f1, df1, 33, 100, pow(10,-16));
  sol2.input(f2, df2, 33, 100, pow(10,-16));
  nptr1 = &sol1;
  nptr2 = &sol2;
  nptr1->solve();
  nptr2->solve();
  std::cout << "Using Newton's Method :" << std::endl;
  std::cout << "D = 55in., Alpha = " << nptr1->get_root() << std::endl;
  std::cout << "D = 30in., Alpha = " << nptr2->get_root() << std::endl;
  std::cout << std::endl;

  EquationSolver<double> *sptr3, *sptr4, *sptr5, *sptr6;
  Secant<double> sol3,sol4,sol5,sol6;
  sol3.input(f1, 32, 34, 100, pow(10,-8), pow(10,-16));
  sol4.input(f2, 32, 34, 100, pow(10,-8), pow(10,-16));
  sol5.input(f1, 0, 100, 100, pow(10,-8), pow(10,-16));
  sol6.input(f2, 0, 100, 100, pow(10,-8), pow(10,-16));
  sptr3 = &sol3;
  sptr4 = &sol4;
  sptr5 = &sol5;
  sptr6 = &sol6;
  sptr3->solve();
  sptr4->solve();
  sptr5->solve();
  sptr6->solve();
  std::cout << "Using Secant Method :" << std::endl;
  std::cout << "D = 55in., Alpha = " << sptr3->get_root() << std::endl;
  std::cout << "D = 30in., Alpha = " << sptr4->get_root() << std::endl;
  std::cout << "When the initial values are far away from 33 degrees : " << std::endl;
  std::cout << "D = 55in., Alpha = " << sptr5->get_root() << std::endl;
  std::cout << "D = 30in., Alpha = " << sptr6->get_root() << std::endl;
  return 0;
}
