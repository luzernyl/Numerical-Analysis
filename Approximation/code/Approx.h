/**
 * @file   Approx.h
 * @author Luzern Yuven Luis <student@student>
 * @date   Sun Dec 19 23:25:32 2021
 * 
 * @brief  Header file for Approx.cpp
 * 
 * 
 */


#ifndef _APPROX_H_
#define _APPROX_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>

class Approx
{
 private:
  // input of x
  std::vector<double> xdata;

  // input of y
  std::vector<double> ydata;

  // Gram matrix
  double Gram[3*3];

  // RHS of normal equations
  double c[3];

  // Result of QR Factorization
  double R1[3*3];

  /** 
   * Function f(x) = 1
   * 
   * @param x input
   * 
   * @return 1
   */
  static double one(double x);

  /** 
   * Function f(x) = x
   * 
   * @param _x input
   * 
   * @return x
   */
  static double x(double _x);

  /** 
   * Function f(x) = x^2
   * 
   * @param _x input
   * 
   * @return x^2
   */
  static double x2(double _x);

  /** 
   * Discrete inner product of two vectors
   * 
   * @param u first vector
   * @param v second vector
   * @param p weight function
   * 
   * @return the inner product <u,v>
   */
  double innerProd(double (*u)(double), double (*v)(double), double (*p)(double));

  /** 
   * Discrete inner product of a vector and y
   * 
   * @param v the vector
   * @param p weight function
   * 
   * @return the inner product <u,y>
   */
  double innerProd(double (*v)(double), double (*p)(double));

 public:
  // Default constructor
  Approx();

  /** 
   * Constructor to create xdata and ydata
   * 
   * @param _xdata xdata
   * @param _ydata ydata
   */
  Approx(std::vector<double> _xdata, std::vector<double> _ydata);

  /** 
   * Solve DLS using Normal Equations
   * 
   */
  void solveNormal();

  /** 
   * Solve DLS using QR Factorization
   * 
   * 
   * @return 
   */
  int solveQR();

  /** 
   * Condition number of Gram matrix
   * 
   * 
   * @return condG
   */
  double condG();

  /** 
   * Condition number of matrix R1
   * 
   * 
   * @return condR1
   */
  double condR1();
};

#include "Approx.cpp"

#endif
