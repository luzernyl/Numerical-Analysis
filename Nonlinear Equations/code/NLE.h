/**
 * @file   NLE.h
 * @author Luzern Yuven Luis <student@student>
 * @date   Mon Oct 11 11:38:50 2021
 * 
 * @brief  An implementation of nonlinear equation solvers.
 *         Bisection Method, Newton's Method, and the Secant Method.
 * 
 * 
 */

#ifndef NLE_H_
#define NLE_H_

#include <iostream>
#include <math.h>
#include <cmath>

/** 
 * The sign function
 * 
 * @param val value
 * 
 * @return sign of the value
 */
template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

/** 
 * An abstract base class with a pure virtual method solve.
 * 
 */
template <typename T>
class EquationSolver
{
 public:
  virtual void solve() = 0;
  virtual T get_root() = 0;
};

/** 
 * A derived class of EquationSolver for the bisection method
 * 
 */
template <typename T>
class Bisection : public EquationSolver<T>
{
 private:
  double c,h,k;
  T (*f)(double);
  double a,b,delta,eps;
  int M;

 public :
  /// Constructor
  Bisection()
  {
    this->c = 0;
    this->h = 0;
    this->k = 1;
  };

  /** 
   * The inputs of the bisection method
   * 
   * @param _f function
   * @param _a left endpoint of interval
   * @param _b right endpoint of interval
   * @param _M maximum number of iterations
   * @param _delta desired length of interval
   * @param _eps desired accuracy
   */
  void input(T (*_f)(double), double _a, double _b, int _M, double _delta, double _eps)
  {
    this->f = _f;
    this->a = _a;
    this->b = _b;
    this->M = _M;
    this->delta = _delta;
    this->eps = _eps;
  };
  
  /// returns the root of the function
  T get_root()
  {
    return this->c;
  };

  /// returns the length of the interval from the last iteration
  double get_h()
  {
    return this->h;
  };
  
  /// returns the number of iterations needed to obtain the root 
  double get_k()
  {
    return this->k;
  };

  /// Implementation of the bisection method
  void solve()
  {
    this->h = this->b - this->a;
    double u = f(this->a);
    this->c = 0;
    double w = 0.0;
    for(this->k = 1; this->k <= this->M; this->k++)
      {
	this->h = this->h/2;
	this->c = this->a + this->h;
	w = f(this->c);
	if (fabs(this->h) < this->delta || fabs(w) < this->eps) break;
	else if (sgn(w) == sgn(u)) this->a = this->c;
      }
  };
};

/**
 *  A derived class of EquationSolver for Newton's method
 * 
 */
template <typename T>
class Newton : public EquationSolver<T>
{
 private:
  double x,k;
  T (*f)(double), (*df)(double);
  double x0, eps;
  int M;
  
 public :
  /// Constructor
  Newton()
  {
    this->x = 0;
    this->k = 0;
  };

  /** 
   * Inputs of Newton's Method
   * 
   * @param _f function
   * @param _df derivative of function
   * @param _x0 initial guess
   * @param _M maximum number of iterations
   * @param _eps desired accuracy
   */
  void input(T (*_f)(double), T (*_df)(double), double _x0, int _M, double _eps)
  {
    this->f = _f;
    this->df = _df;
    this->x0 = _x0;
    this->M = _M;
    this->eps = _eps;
  };
  
  /// returns the root of the function 
  T get_root()
  {
    return this->x;
  };
  /// returns the number of iterations needed to obtain the root
  double get_k()
  {
    return this->k;
  };

  /// An implementation of Newton's Method
  void solve()
  {
    this->x = x0;
    double u = 0;
    for(this->k = 0; k <= M; k++)
      {
	u = f(this->x);
	if (fabs(u) < eps) break;
	this->x = this->x - (u / df(x));
      }
  };
};

/**
 * A derived class of EquationSolver for secant method
 * 
 */
template <typename T>
class Secant : public EquationSolver<T>
{
 private:
  double xn, xn_1, k;
  T (*f)(double);
  double x0, x1, delta, eps;
  int M;
  
 public:
  /// Constructor
  Secant()
  {
    this->xn = 0;
    this->xn_1 = 0;
    this->k = 2;
  };

  /** 
   * Inputs of the secant method
   * 
   * @param _f function
   * @param _x0 initial guess 1
   * @param _x1 initial guess 2
   * @param _M maximum number of iterations
   * @param _delta desired difference of iterates
   * @param _eps desired accuracy
   */
  void input(T (*_f)(double), double _x0, double _x1, int _M, double _delta, double _eps)
  {
    this->f = _f;
    this->x0 = _x0;
    this->x1 = _x1;
    this->M = _M;
    this->delta = _delta;
    this->eps = _eps;
  };

  /// returns the root of the function
  T get_root()
  {
    return this->xn;
  };
  
  /// returns the iterate of the (n-1)th iteration
  double get_xn_1()
  {
    return this->xn_1;
  };

  /// returns the number of iterations needed to obtain the root
  double get_k()
  {
    return k;
  };

  /// An implementation of the secant method
  void solve()
  {
    this->xn = this->x1;
    this->xn_1 = this->x0;
    double u = f(this->xn);
    double v = f(this->xn_1);
    double s = 0;
    for(this->k = 2; k <= this->M; this->k++)
      {
	if(fabs(u) > fabs(v))
	  {
	    double tmp1 = this->xn;
	    this->xn = this->xn_1;
	    this->xn_1 = tmp1;
	    double tmp2 = u;
	    u = v;
	    v = tmp2;
	  }
	s = (this->xn - this->xn_1) / (u-v);
	this->xn_1 = this->xn;
	v = u;
	this->xn = this->xn-u*s;
	u = f(this->xn);
	if(fabs(this->xn - this->xn_1) < this->delta || fabs(u) < this->eps) break;
      }
  };
};

#endif
