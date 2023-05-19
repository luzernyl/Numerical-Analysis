#ifndef _SPLINES_CPP_
#define _SPLINES_CPP_

#include "Splines.h"

InterpConditions::InterpConditions(){};

InterpConditions::InterpConditions(std::vector<double> _site, std::vector<std::vector<double> > _value, std::vector<double> _der)
{
  site.resize(_site.size());
  value.resize(_value.size());
  //if(distinct(_site) == true)
  //{
    site = _site;
    value = _value;
    der = _der;
    //}
    //else
    // std::cerr << "Interpolation sites must be distinct !" << std::endl;
}
/*
bool InterpConditions::distinct(std::vector<double> _site)
{
  for(int i = 0; i < _site.size(); i++)
    {
      if(std::cout(_site.begin(), _site.end(), _site[i]) > 1)
	return false;
    }
  return true;
}
*/

int InterpConditions::get_site_size()
{
  return site.size();
}

std::vector<double> InterpConditions::get_value_size()
{
  std::vector<double> sz;
  for(int i = 0; i < value.size(); i++)
    sz.push_back(value[i].size());
  return sz;
}

void InterpConditions::input_site(double x)
{
  site.push_back(x);
}

void InterpConditions::input_value(std::vector<double> v)
{
  value.push_back(v);
}

void InterpConditions::input_der(double x)
{
  der.push_back(x);
}

int InterpConditions::get_N() const
{
  int s = 0;
  for(int i = 0; i < value.size(); i++)
    s += value[i].size();
  
  return s;
}

int InterpConditions::find_site(double x) const
{
  for(int i = 0; i < site.size(); i++)
  {
    if(site[i] == x)
      return i;
  }
  return -1;
}

template<int Order, class CoefType>
Polynomial<Order,CoefType>::Polynomial()
{
  degCoeff.resize(Order+1);
}


template<int Order, class CoefType>
Polynomial<Order,CoefType>::Polynomial(const std::vector<CoefType> &_degCoeff)
{
  degCoeff.resize(_degCoeff.size());
  degCoeff = _degCoeff;
}

template<int Order, class CoefType>
//template<int _Order, class _CoefType>
Polynomial<Order,CoefType>::Polynomial(Polynomial<Order, CoefType> const &_p)
{
  std::vector<CoefType> newdeg(Order+1, 0);
  for(int i = 0; i <= Order+1; i++)
    newdeg[i] = _p.get_coeff(i);

  degCoeff = newdeg;
  //Order = _Order+1;
}

template<int Order, class CoefType>
void Polynomial<Order,CoefType>::set_order(int _order)
{
  Order = _order;
}

template<int Order, class CoefType>
int Polynomial<Order,CoefType>::get_order()
{
  return Order;
}

template<int Order, class CoefType>
void Polynomial<Order,CoefType>::set_coeff(int k, CoefType _coeff)
{
  degCoeff[k] = _coeff;
}

template<int Order, class CoefType>
CoefType Polynomial<Order,CoefType>::get_coeff(int k) const
{
  return degCoeff[k];
}

template<int Order, class CoefType>
template<int _Order, class _CoefType>
Polynomial<std::max(Order,_Order),CoefType> Polynomial<Order,CoefType>::operator+(Polynomial<_Order,_CoefType> const &p2)
{
  const int NewOrder = std::max(Order,_Order);
  Polynomial<NewOrder,CoefType> p3;

  for(int i = 0; i <= NewOrder; i++)
  {
    if(i <= Order && i <= _Order)
      p3.degCoeff[i] = degCoeff[i] + p2.get_coeff(i);
    else if(i > _Order)
      p3.degCoeff[i] = degCoeff[i];
    else
      p3.degCoeff[i] = p2.get_coeff(i);
  }

  return p3;
}
// friend operator
template<int Order, class CoefType, int _Order, class _CoefType>
Polynomial<std::max(Order,_Order),CoefType> operator+(Polynomial<Order,CoefType> const &p1, Polynomial<_Order,_CoefType> const &p2)
{
  const int NewOrder = std::max(Order,_Order);
  Polynomial<NewOrder,CoefType> p3;

  for(int i = 0; i <= NewOrder; i++)
  {
    if(i <= Order && i <= _Order)
      p3.degCoeff[i] = p1.degCoeff[i] + p2.degCoeff[i];
    else if(i > _Order)
      p3.degCoeff[i] = p1.degCoeff[i];
    else
      p3.degCoeff[i] = p2.degCoeff[i];
  }

  return p3;
}

template<int Order, class CoefType>
template<int _Order, class _CoefType>
Polynomial<std::max(Order,_Order),CoefType> Polynomial<Order,CoefType>::operator-(Polynomial<_Order, _CoefType> const &p2)
{
  const int NewOrder = std::max(Order,_Order);
  Polynomial<NewOrder,CoefType> p3;

  for(int i = 0; i <= NewOrder; i++)
  {
    if(i <= Order && i <= _Order)
      p3.degCoeff[i] = degCoeff[i] - p2.degCoeff[i];
    else if (i > _Order)
      p3.degCoeff[i] = degCoeff[i];
    else
      p3.degCoeff[i] = (p2.degCoeff[i])*(-1);
  }

  return p3;
}

// friend operator
template<int Order, class CoefType, int _Order, class _CoefType>
Polynomial<std::max(Order,_Order),CoefType> operator-(Polynomial<Order,CoefType> const &p1, Polynomial<_Order, _CoefType> const &p2)
{
  const int NewOrder = std::max(Order,_Order);
  Polynomial<NewOrder,CoefType> p3;

  for(int i = 0; i <= NewOrder; i++)
  {
    if(i <= Order && i <= _Order)
      p3.degCoeff[i] = p1.degCoeff[i] - p2.degCoeff[i];
    else if(i > _Order)
      p3.degCoeff[i] = p1.degCoeff[i];
    else
      p3.degCoeff[i] = (p2.degCoeff[i])*(-1);
  }

  return p3;
}

template<int Order, class CoefType>
template<int _Order, class _CoefType>
Polynomial<(Order+_Order),CoefType> Polynomial<Order,CoefType>::operator*(Polynomial<_Order,_CoefType> const &p2)
{
  const int NewOrder = Order + _Order;
  Polynomial<NewOrder,CoefType> p3;

  for(int i = 0; i <= Order; i++)
  {
    for(int j = 0; j <= _Order; j++)
      p3.set_coeff(i+j, p3.get_coeff(i+j) + degCoeff[i] * (p2.get_coeff(j)));
  }

  return p3;
}

// friend operator
template<int Order, class CoefType, int _Order, class _CoefType>
Polynomial<(Order+_Order),CoefType> operator*(Polynomial<Order,CoefType> const &p1, Polynomial<_Order,_CoefType> const &p2)
{
  const int NewOrder = Order + _Order;
  Polynomial<NewOrder,CoefType> p3;

  for(int i = 0; i <= Order; i++)
  {
    for(int j = 0; j <= _Order; j++)
      p3.set_coeff(i+j, p3.get_coeff(i+j) + p1.get_coeff(i)*p2.degCoeff[j]);
  }

  return p3;
}

template<int Order, class CoefType>
//template<int _Order, class _CoefType>
void Polynomial<Order,CoefType>::operator=(Polynomial<Order,CoefType> const &p)
{
  std::vector<CoefType> newdeg(Order+1,0);

  for(int i = 0; i <= Order; i++)
    newdeg[i] = p.get_coeff(i);

  degCoeff = newdeg;
}

template<int Order, class CoefType>
CoefType Polynomial<Order,CoefType>::eval(double x)
{
  CoefType result = 0;
  for(int i = 0; i <= Order; i++)
  {
    result = result + degCoeff[i]*pow(x,i);
  }
  return result;
}

template<int Order, class CoefType>
Polynomial<Order-1,CoefType> Polynomial<Order,CoefType>::derivative()
{
  Polynomial<Order-1,CoefType> der;
  for(int i = 0; i < Order; i++)
  {
    der.set_coeff(i, degCoeff[i+1]*(i+1));
    //der.degCoeff[i] = (degCoeff[i+1])*(i+1);
  }
  return der;
}

template<int _Order, class _CoefType>
std::ostream& operator<<(std::ostream &os, Polynomial<_Order,_CoefType> const &p)
{
  for(int i = _Order; i >= 1; i--)
  {
    os << p.degCoeff[i] << "*x^" << i << " + ";
  }
  os << p.degCoeff[0];
  return os;
}

double factorial(int n)
{
    double sum = 1;
    for(int i = 1 ; i <= n ; i ++)
        sum *= i;
    return sum;
}

template<SplineType t>
Spline<1,3,t>::Spline(const Spline<1,3,t> &spl)
{
  tableOfDividedDiffs = spl.tableOfDividedDiffs;
  N = spl.N;
  spline = spl.spline;
}

template<SplineType t>
int Spline<1,3,t>::make_table(const InterpConditions &data)
{
  tableOfDividedDiffs.clear();
  for(int i = 0; i < data.value.size(); i++)
  {
    std::vector<double> v{data.site[i], data.value[i][0]};
    for(int j = 1; j <= data.value[i].size(); j++)
      tableOfDividedDiffs.push_back(v);
  }

  for(int i = 2; i <= data.get_N(); i++)
  {
    for(int j = i - 1; j <= data.get_N() - 1; j++)
    {
      if(tableOfDividedDiffs[j-i+1][0] == tableOfDividedDiffs[j][0])
	tableOfDividedDiffs[j].push_back(data.value[data.find_site(tableOfDividedDiffs[j][0])][i-1] / factorial(i-1));
      else
      {
	double tmp = (tableOfDividedDiffs[j][i-1]- tableOfDividedDiffs[j-1][i-1]) / (tableOfDividedDiffs[j][0] - tableOfDividedDiffs[j-i+1][0]);
	tableOfDividedDiffs[j].push_back(tmp);
      }
    }
  }

  return 0;
}

template<SplineType t>
void Spline<1,3,t>::print_table()
{
  for(int i = 0; i < tableOfDividedDiffs.size(); i++)
  {
    for(int j = 0; j < tableOfDividedDiffs[i].size(); j++)
      std::cout << tableOfDividedDiffs[i][j] << " ";
    std::cout << std::endl;
  }
}

template<SplineType tt>
Spline<1,3,tt> interpolate(const InterpConditions &data, BCType type)
{
  const std::vector<double> &site = data.site;
  const std::vector<std::vector<double> > &value = data.value;
  const std::vector<double> &der = data.der;
  Spline<1,3,tt> spl;
  spl.make_table(data);
  std::vector<std::vector<double> > &table = spl.tableOfDividedDiffs;
  int N = site.size();
  int n = N-1;
  double DL[n], D[N], DU[n], b[N];

  for(int i = 0; i < N; i++)
  {
    D[i] = 2;
    b[i] = table[i+2][3]*6;
    //int miu = (site[i+1] - site[i]) / (site[i+2] - site[i]);
    //int lambda = (site[i+2] - site[i+1]) / (site[i+2] - site[i]);
    //b[i] = 3*miu*table[i+2][2] - 3*lambda*table[i+1][2];
    //std::cout << b[i] << " ";
  }
  //std::cout << std::endl;
  
  if(type == BCType::complete)
  {
    for(int i = 0; i < n; i++)
    {
	DL[i] = (site[i+1] - site[i]) / (site[i+2] - site[i]);
	DU[i+1] = (site[i+2] - site[i+1]) / (site[i+2] - site[i]);
    }
    DU[0] = 1;
    DL[n-1] = 1;
  }

  lapack_int info = LAPACKE_dgtsv(LAPACK_COL_MAJOR, N, 1, DL, D, DU, b, N);

  spl.spline.resize(N-1);
  spl.spline.clear();
  for(int i = 0; i < N-1; i++)
  {
    Polynomial<3,double> result;
    Polynomial<0,double> p;
    result.set_coeff(0, table[i+1][1]);
    std::vector<double> factor{-site[i], 1};
    Polynomial<1,double> pfac(factor);
    Polynomial<0,double> x;
    Polynomial<1,double> p1;
    Polynomial<2,double> p2;
    Polynomial<3,double> p3;
    for(int j = 1; j <= 3; j++)
    {
      //p = p * p1;
      if(j == 1)
      {
	p1 = pfac;
	x.set_coeff(0, der[i]);
	result = result + x * p1;
      }
      else if(j == 2)
      {
	p2 = p1 * pfac;
	x.set_coeff(0, b[i] / 2.0);
	result = result + x*p2;
      }
      else
      {
	p3 = p2 * pfac;
	x.set_coeff(0, (b[i+1] - b[i]) / (6 * (site[i+1] - site[i])));
	result = result + x*p3;
      }
      //result = result + x*p;
    }
    spl.spline.push_back(result);
  }

  return spl;
}

template<SplineType t>
Polynomial<3,double> Spline<1,3,t>::get_poly(int i)
{
  return spline[i];
}

#endif
