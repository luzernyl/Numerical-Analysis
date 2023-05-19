#ifndef _SPLINES_H_
#define _SPLINES_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <ostream>
#include <lapacke.h>
#include <algorithm>
#include "Vector.h"

struct InterpConditions
{
  std::vector<double> site;
  std::vector<std::vector<double> > value;
  std::vector<double> der;
  
  //bool distinct(std::vector<double> _site);
  
  InterpConditions();
  
  InterpConditions(std::vector<double> _site, std::vector<std::vector<double> > _value, std::vector<double> der);

  void set_site_size(int n);

  void set_value_size(int n);

  int get_site_size();
  
  std::vector<double> get_value_size();

  void input_site(double x);

  void input_value(std::vector<double> v);

  void input_der(double x);

  int get_N() const;

  int find_site(double x) const;
  
};

template<int Order, class CoefType>
class Polynomial
{
 private:
  std::vector<CoefType> degCoeff;

 public:
  Polynomial();
  
  Polynomial(const std::vector<CoefType> &_degCoeff);

  //template<int _Order, class _CoefType>
  Polynomial(Polynomial<Order, CoefType> const &_p);
  
  void set_order(int _order);
  
  int get_order();
  
  void set_coeff(int k, CoefType _coeff);
  
  CoefType get_coeff(int k) const;

  template<int _Order, class _CoefType>
  Polynomial<std::max(Order,_Order),CoefType> operator+(Polynomial<_Order,_CoefType> const &p2);

  template<int _Order, class _CoefType>
  friend Polynomial<std::max(Order,_Order),CoefType> operator+(Polynomial<Order,CoefType> const &p1, Polynomial<_Order,_CoefType> const &p2);

  template<int _Order, class _CoefType>
  Polynomial<std::max(Order,_Order),CoefType> operator-(Polynomial<_Order, _CoefType> const &p2);

  template<int _Order, class _CoefType>
  friend Polynomial<std::max(Order,_Order),CoefType> operator-(Polynomial<Order,CoefType> const &p1, Polynomial<_Order, _CoefType> const &p2);

  template<int _Order, class _CoefType>
  Polynomial<Order+_Order,CoefType> operator*(Polynomial<_Order,_CoefType> const &p2);

  template<int _Order, class _CoefType>
  friend Polynomial<Order+_Order,CoefType> operator*(Polynomial<Order,CoefType> const &p1, Polynomial<_Order,_CoefType> const &p2);

  //template<int _Order, class _CoefType>
  void operator=(Polynomial<Order,CoefType> const &p);
  
  CoefType eval(double x);
  
  Polynomial<Order-1,CoefType> derivative();

  template<int _Order, class _CoefType>
  friend std::ostream& operator<<(std::ostream &os, const Polynomial<_Order,_CoefType> &p);
  
};


enum class SplineType {ppForm, cardinalB};
enum BCType {complete, notAknot, periodic};

template<int Dim, int Order, SplineType t>
class Spline
{
 public:
  Spline();
};

// template specialization for Dim = 1
template<SplineType t>
class Spline<1,3,t>
{
 private:
  std::vector<std::vector<double> > tableOfDividedDiffs;
  int N;
  std::vector<Polynomial<3,double> > spline;
  
 public:
  Spline<1,3,t>(){};
  Spline<1,3,t>(const Spline<1,3,t> &spl);
  int make_table(const InterpConditions &data);
  void print_table();

  template<SplineType tt> friend Spline<1,3,tt>
  interpolate(const InterpConditions &, BCType);

  Polynomial<3,double> get_poly(int i);
};


/*
// template specialization for SplineType of ppForm
template<int Dim, int Order>
class Spline<Dim, Order, SplineType::ppForm>
{
 public:
  template<int Ord> friend Spline<2,Ord,SplineType::ppForm>
  fitcurve(const std::vector<Vec<double,2> > &, BCType);
};
*/

#include "Splines.cpp"

#endif
