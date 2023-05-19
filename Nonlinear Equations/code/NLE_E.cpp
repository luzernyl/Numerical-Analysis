/**
 * @file   NLE_E.cpp
 * @author Luzern Yuven Luis <student@student>
 * @date   Mon Oct 11 12:30:01 2021
 * 
 * @brief  Program for Problem E
 * 
 * 
 */

#include <iostream>
#include <iomanip>
#include "NLE.h"
#include <cmath>
#include <math.h>

double func(double h)
{
  return 12.4 - 10*(0.5*M_PI*1 - 1*asin(h) - h*sqrt((1-pow(h,2))));
  //return 12.4 - 10*(0.5*M_PI*1 - 1*asin(1-h) - (1-h)*sqrt((1-pow(1-h,2))));
}

double derfunc(double h)
{
  return 10*(1/sqrt(1-pow(h,2))+sqrt(1-pow(h,2))+(h*0.5*(1/pow(1-pow(h,2),1/2)*2*h)));
}

int main(int args, char* argv[])
{
  EquationSolver<double> *bptr, *nptr, *sptr;
  Bisection<double> bissol;
  Newton<double> newsol;
  Secant<double> secsol;
  bissol.input(func, 0, 1, 100, pow(10,-8), pow(10,-16));
  bptr = &bissol;
  bptr->solve();
  std::cout << std::fixed << std::setprecision(2);
  std::cout << "Using Bisection Method : Depth = " << 1 - bptr->get_root() << " ft"<< std::endl;
  newsol.input(func, derfunc, 0.8, 100, pow(10,-16));
  nptr = &newsol;
  nptr->solve();
  std::cout << "Using Newton's Method : Depth = " << 1 - nptr->get_root() << " ft" << std::endl;
  secsol.input(func, 0.8, 0.85, 100, pow(10,-8), pow(10,-16));
  sptr = &secsol;
  sptr->solve();
  std::cout << "Using Secant Method : Depth = " << 1 - sptr->get_root() << " ft" << std::endl;
  return 0;
}
