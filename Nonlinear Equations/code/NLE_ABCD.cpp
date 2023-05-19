/**
 * @file   NLE_ABCD.cpp
 * @author Luzern Yuven Luis <student@student>
 * @date   Mon Oct 11 11:58:39 2021
 * 
 * @brief  Program for Problems A,B,C,D
 * 
 * 
 */

#include "NLE.h"
#include <iostream>
#include <math.h>
#include <cmath>

/// function for B.1
double bisfunc1(double x)
{
  return 1/x - tan(x);
}

/// function for B.2
double bisfunc2(double x)
{
  return 1/x - pow(2,x);
}

/// function for B.3
double bisfunc3(double x)
{
  return 1/(pow(2,x)) + exp(x) + 2*cos(x) - 6;
}

/// function for B.4
double bisfunc4(double x)
{
  return (pow(x,3)+4*pow(x,2)+3*x+5) / (2*pow(x,3)-9*pow(x,2)+18*x-2);
}

/// function for C
double newfunc(double x)
{
  return x - tan(x);
}

/// derivative of the function for C
double dernewfunc(double x)
{
  return 1 - (1/pow(cos(x),2));
}

/// function for D.1
double secfunc1(double x)
{
  return sin(x/2) - 1;
}

/// function for D.2
double secfunc2(double x)
{
  return exp(x) - tan(x);
}

/// function for D.3
double secfunc3(double x)
{
  return pow(x,3) - 12*pow(x,2) + 3*x + 1;
}

int main(int argc, char* argv[])
{
  EquationSolver<double> *bptr1, *bptr2, *bptr3, *bptr4;
  Bisection<double> bissol1,bissol2,bissol3,bissol4;
  bissol1.input(bisfunc1, 0, M_PI/2, 100, pow(10,-8), pow(10,-16));
  bissol2.input(bisfunc2, 0, 1, 100, pow(10,-8), pow(10,-16));
  bissol3.input(bisfunc3, 1, 3, 100, pow(10,-8), pow(10,-16));
  bissol4.input(bisfunc4, 0, 4, 100, pow(10,-8), pow(10,-16));
  bptr1 = &bissol1;
  bptr2 = &bissol2;
  bptr3 = &bissol3;
  bptr4 = &bissol4;
  bptr1->solve();
  bptr2->solve();
  bptr3->solve();
  bptr4->solve();
  std::cout << "Bisection Method" << std::endl;
  std::cout << "Solution 1 = " << bptr1->get_root() << std::endl;
  std::cout << "Solution 2 = " << bptr2->get_root() << std::endl;
  std::cout << "Solution 3 = " << bptr3->get_root() << std::endl;
  std::cout << "Solution 4 = " << bptr4->get_root() << std::endl;
  std::cout << std::endl;

  EquationSolver<double> *nptr1, *nptr2;
  Newton<double> newsol1, newsol2;
  newsol1.input(newfunc, dernewfunc, 4.5, 100, pow(10,-16));
  newsol2.input(newfunc, dernewfunc, 7.7, 100, pow(10,-16));
  nptr1 = &newsol1;
  nptr2 = &newsol2;
  nptr1->solve();
  nptr2->solve();
  std::cout << "Newton Method" << std::endl;
  std::cout << "Root near 4.5 = " << nptr1->get_root() << std::endl;
  std::cout << "Root near 7.7 = " << nptr2->get_root() << std::endl;
  std::cout << std::endl;

  EquationSolver<double> *sptr1, *sptr1a, *sptr2, *sptr2a, *sptr3, *sptr3a;
  Secant<double> secsol1,secsol1a,secsol2,secsol2a,secsol3,secsol3a;
  secsol1.input(secfunc1, 0, M_PI/2, 100, pow(10,-8), pow(10,-16));
  secsol1a.input(secfunc1, -M_PI/2, M_PI/2, 100, pow(10,-8), pow(10,-16));
  secsol2.input(secfunc2, 1, 1.4, 100, pow(10,-8), pow(10,-16));
  secsol2a.input(secfunc2, 0, 5, 100, pow(10,-8), pow(10,-16));
  secsol3.input(secfunc3, 0, -0.5, 100, pow(10,-8), pow(10,-16));
  secsol3a.input(secfunc3, 2, 3, 100, pow(10,-8), pow(10,-16));
  sptr1 = &secsol1;
  sptr1a = &secsol1a;
  sptr2 = &secsol2;
  sptr2a = &secsol2a;
  sptr3 = &secsol3;
  sptr3a = &secsol3a;
  sptr1->solve();
  sptr1a->solve();
  sptr2->solve();
  sptr2a->solve();
  sptr3->solve();
  sptr3a->solve();
  std::cout << "Secant Method" << std::endl;
  std::cout << "Solution 1 = " << sptr1->get_root() << std::endl;
  std::cout << "Solution 1 (x0 = 0, x1 = PI/4) = " << sptr1a->get_root() << std::endl;
  std::cout << "Solution 2 = " << sptr2->get_root() << std::endl;
  std::cout << "Solution 2a (x0 = -5, x1 = 5) = " << sptr2a->get_root() << std::endl;
  std::cout << "Solution 3 = " << sptr3->get_root() << std::endl;
  std::cout << "Solution 3a (x0 = 2, x1 = 3) = " << sptr3a->get_root() << std::endl;
 
return 0;
}
