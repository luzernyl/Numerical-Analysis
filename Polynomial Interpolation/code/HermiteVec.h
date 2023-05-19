/**
 * @file   HermiteVec.h
 * @author Luzern Yuven Luis<student@student>
 * @date   Sat Oct 30 23:40:27 2021
 * 
 * @brief  A C++ package to solve Hermite interpolation problem.
 *         Consists of three classes :
 *         struct InterpConditions, class Polynomial, class NewtonInterp
 * 
 * 
 */

#ifndef HERMITE
#define HERMITE

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <ostream>

/**
 * Stores the interpolation sites and corresponding function values
 * and derivatives
 * 
 */
struct InterpConditions
{
  std::vector<double> site, value, der;

  /** 
   * Determines whether interpolation sites are distinct from each other.
   * 
   * @param site_ interpolation sites
   * 
   * @return true if sites are distinct
   */
  bool distinct(std::vector<double> site_)
  {
    for(int i = 0; i < site_.size(); i++)
      {
	if(std::count(site_.begin(), site_.end(), site_[i]) > 1)
	   return false;
      }
    return true;
  }

  // Default Constructor
  InterpConditions(){};

  /** 
   * Constructor to store interpolation sites, corresponding function values and derivatives
   * 
   * @param site_ interpolation sites
   * @param value_ function values
   * @param der_ function derivatives
   * 
   */
  InterpConditions(std::vector<double> site_, std::vector<double> value_, std::vector<double> der_)
  {
    this->site.resize(site_.size());
    this->value.resize(value_.size());
    this->der.resize(der_.size());
    if(distinct(site_) == true)
    {
	this->site = site_;
	this->value = value_;
	this->der = der_;
    }
    else
      std::cerr << "Interpolation sites must be distinct !" << std::endl;
  }

  // returns the number of sites
  int get_site_size()
  {
    return this->site.size();
  }

  // returns the number of function values
  int get_value_size()
  {
    return this->value.size();
  }

  // returns the number of derivative values
  int get_der_size()
  {
    return this->der.size();
  }

  // input site manually
  void input_site(double x)
  {
    this->site.push_back(x);
  }

  // input function value manually
  void input_value(double x)
  {
    this->value.push_back(x);
  }

  // input derivative value manually
  void input_der(double x)
  {
    this->der.push_back(x);
  }

  /** 
   * Determines whether an interpolation problem is 
   * a Hermite interpolation problem
   * 
   * 
   * @return true if it is a Hermite problem
   */
  bool isHermite()
  {
    bool flag = false;
    for(int i = 0; i < this->der.size(); i++)
      {
	if(der[i] != 0)
	{
	  flag = true;
	  break;
	}
      }
    return flag;
  }
  
};

/**
 * An implementation of Polynomial. Consists of 
 * addition, substraction, multiplication, evaluation
 * at some x, output to an ostream, returning the 
 * degree and coefficients of a polynomial
 * 
 */
class Polynomial
{
 private:
  std::vector<double> degCoeff;
  int capacity;

 public:
  /** 
   * Default constructor.
   * Initializes a zero-degree polynomial.
   * 
   */
  Polynomial()
  {
    this->capacity = 0;
    this->degCoeff.resize(1);
  }

  /** 
   * Constructor to initialize a polynomial with
   * specified capacity
   * 
   * @param capacity_ capacity of polynomial
   */
  Polynomial(int capacity_)
  {
    this->capacity = capacity_;
    this->degCoeff.resize(capacity+1);
  }

  /** 
   * Constructor to initialize a polynomial with a
   * given polynomial
   * 
   * @param p given polynomial
   */
  Polynomial (Polynomial const &p)
  {
    std::vector<double> newdeg(p.capacity+1,0);

    for(int i = 0; i <= p.capacity; i++)
      newdeg[i] = p.degCoeff[i];
    
    this->degCoeff=newdeg;
        
    this->capacity=p.capacity;
  }

  // set all coefficients to zero
  void set_to_zero()
  {
    for(int i = 0; i <= this->capacity; i++)
    {
      this->degCoeff[i] = 0;
    }
  }

  /** 
   * Set a coefficient of the k-th degree
   * 
   * @param deg_ k-th degree
   * @param coeff_ new coefficient
   */
  void set_coeff(int deg_, double coeff_)
  {
    if(deg_ > this->capacity)
    {
      std::cerr<< "Over Capacity !";
    }
    else
    {
      this->degCoeff[deg_] = coeff_;
    }
  }

  // set capacity of polynomial to cap
  void set_capacity(int cap)
  {
    this->capacity = cap;
    this->degCoeff.resize(cap+1);
  }

  // returns the capacity of polynomial
  int get_capacity()
  {
    return this->capacity;
  }

  // returns the coefficient of the i-th degree
  double get_coeff(int i)
  {
    return this->degCoeff[i];
  }

  /** 
   * The addition operator for polynomials
   * 
   * @param p2 the second polynomial
   * 
   * @return polynomial after addition
   */
  Polynomial operator+(Polynomial const &p2)
  {
    int newcap = std::max(this->capacity, p2.capacity);
    Polynomial p3(newcap);

    for(int i = 0; i <= newcap; i++)
    {
      if(i <= this->capacity && i <= p2.capacity)
	p3.degCoeff[i] = this->degCoeff[i] + p2.degCoeff[i];
      else if(i > p2.capacity)
	p3.degCoeff[i] = this->degCoeff[i];
      else
	p3.degCoeff[i] = p2.degCoeff[i];
    }

    return p3;
  }

  friend Polynomial operator+(Polynomial const &p1, Polynomial const &p2);

  /** 
   * The substraction operator for polynomials
   * 
   * @param p2 the second polynomial
   * 
   * @return the polynomial after subtraction
   */
  Polynomial operator-(Polynomial const &p2)
  {
    int newcap = std::max(this->capacity, p2.capacity);
    Polynomial p3(newcap);

    for(int i = 0; i <= newcap; i++)
    {
      if(i <= this->capacity && i <= p2.capacity)
	p3.degCoeff[i] = this->degCoeff[i] - p2.degCoeff[i];
      else if(i > p2.capacity)
	p3.degCoeff[i] = this->degCoeff[i];
      else
	p3.degCoeff[i] = (-1)*(p2.degCoeff[i]);
    }

    return p3;
  }

  friend Polynomial operator-(Polynomial const &p1, Polynomial const &p2);

  /** 
   * The multiplication operator for polynomials
   * 
   * @param p2 the second polynomial
   * 
   * @return the polynomial after multiplication
   */
  Polynomial operator*(Polynomial const &p2)
  {
    int newcap = this->capacity + p2.capacity;
    Polynomial p3(newcap);

    for(int i = 0; i <= this->capacity; i++)
    {
	for(int j = 0; j <= p2.capacity; j++)
	{
	  p3.degCoeff[i+j] += this->degCoeff[i]*p2.degCoeff[j];
	}
    }

    return p3;
  }
  
  friend Polynomial operator*(Polynomial const &p1, Polynomial const &p2);

  /** 
   * The assignment operator for polynomials
   * 
   * @param p the second polynomial
   */
  void operator=(Polynomial const &p)
  {
    std::vector<double> newdeg(p.capacity+1,0);
    //Copy the contents
    for(int i = 0; i <= p.capacity;i++)
        newdeg[i]=p.degCoeff[i];
                
    this->degCoeff=newdeg;
        
    this->capacity=p.capacity;
  }

  /** 
   * Evalution of a polynomial at some x
   * 
   * @param x a point
   * 
   * @return the value of the polynomial at x
   */
  double eval(double x)
  {
    double result = 0;
    for(int i = 0; i <= this->capacity; i++)
    {
      result += this->degCoeff[i]*pow(x,i);
    }
    return result;
  }

  // returns the derivative of a polynomial
  Polynomial derivative()
  {
    Polynomial der(this->capacity - 1);
    for(int i = 0; i < this->capacity; i++)
    {
      der.degCoeff[i] = (this->degCoeff[i+1])*(i+1);
    }
    return der;
  }

  // output to an ostream for polynomials
  std::ostream& operator<<(std::ostream& os)
  {
    for(int i = this->capacity; i >= 0; i--)
    {
      if (i != 0)
	os << this->degCoeff[i] << "x^" << i << " + ";
      else
	os << this->degCoeff[0]; 
    }
    return os;
  }

  // prints a polynomial
  void print()
  {      
    for(int i=0;i<=this->capacity;i++)
    {
       if(degCoeff[i]!=0)
       std::cout<<degCoeff[i]<<"*x.^"<<i<<"  + ";
    }
  }
};

// friend addition operator
Polynomial operator+(Polynomial const &p1, Polynomial const &p2)
{
    int newcap = std::max(p1.capacity, p2.capacity);
    Polynomial p3(newcap);

    for(int i = 0; i <= newcap; i++)
    {
      if(i <= p1.capacity && i <= p2.capacity)
	p3.degCoeff[i] = p1.degCoeff[i] + p2.degCoeff[i];
      else if(i > p2.capacity)
	p3.degCoeff[i] = p1.degCoeff[i];
      else
	p3.degCoeff[i] = p2.degCoeff[i];
    }

    return p3;
}

// friend substraction operator
Polynomial operator-(Polynomial const &p1, Polynomial const &p2)
{
    int newcap = std::max(p1.capacity, p2.capacity);
    Polynomial p3(newcap);

    for(int i = 0; i <= newcap; i++)
    {
      if(i <= p1.capacity && i <= p2.capacity)
	p3.degCoeff[i] = p1.degCoeff[i] - p2.degCoeff[i];
      else if(i > p2.capacity)
	p3.degCoeff[i] = p1.degCoeff[i];
      else
	p3.degCoeff[i] = (-1)*(p2.degCoeff[i]);
    }

    return p3;
}

// friend multiplication operator
Polynomial operator*(Polynomial const &p1, Polynomial const &p2)
{
    int newcap = p1.capacity + p2.capacity;
    Polynomial p3(newcap);

    for(int i = 0; i <= p1.capacity; i++)
    {
	for(int j = 0; j <= p2.capacity; j++)
	{
	  p3.degCoeff[i+j] += p1.degCoeff[i]*p2.degCoeff[j];
	}
    }

    return p3;
}

/**
 * Class for implementing Newton's formula and 
 * Neville-Aitken algorithm
 * 
 */
class NewtonInterp : public Polynomial
{
 private:
  Polynomial interPoly_; // the interpolating polynomial
  static Polynomial NAPoly_; // interpolating polynomial for Neville-Aitken
  std::vector< std::vector<double> > tableOfDividedDiffs_;
  std::vector<double> newsite;
  Polynomial symmPoly; // symmetric polynomial

 public:
  
  /** 
   * Default Constructor. Initializes the interpolating
   * polynomial and the table of divided differences.
   * 
   */
  NewtonInterp()
  {
    Polynomial interPoly_();
    std::vector< std::vector<double> > tableOfDividedDiffs_;
    }

  /** 
   * Constructor to initialize table of divided differences
   * based on interpolation data
   * 
   * @param data interpolation data
   *  
   */
  NewtonInterp(InterpConditions data)
  {
    // Determines whether the problem is Hermite
    if(data.isHermite() == true)
    {
      this->newsite.resize(2*data.site.size());
      this->tableOfDividedDiffs_.resize(newsite.size());
      for(int i = 0; i < newsite.size(); i++)
      {
	this->tableOfDividedDiffs_[i].resize(newsite.size());
      }
      this->interPoly_.set_capacity(newsite.size());
    
      for(int i = 0; i < data.site.size(); i++)
      {
	this->newsite[2*i] = data.site[i];
        this->newsite[(2*i)+1] = data.site[i];
	this->tableOfDividedDiffs_[2*i][0] = data.value[i];
        this->tableOfDividedDiffs_[(2*i)+1][0] = data.value[i];
	this->tableOfDividedDiffs_[(2*i)+1][1] = data.der[i];
      
	if(i != 0)
	{
	  this->tableOfDividedDiffs_[2*i][1] = (this->tableOfDividedDiffs_[2*i][0] - this->tableOfDividedDiffs_[2*i-1][0]) / (this->newsite[2*i] - this->newsite[2*i-1]);
	}
      }

      for(int i = 2; i < this->newsite.size(); i++)
      {
	 for(int j = 2; j <= i; j++)
	 {
	   this->tableOfDividedDiffs_[i][j] = (this->tableOfDividedDiffs_[i][j-1] - this->tableOfDividedDiffs_[i-1][j-1]) / (this->newsite[i] - this->newsite[i-j]);
	 }
      }
    }
    // Problem is Newton interpolation
    else
    {
      this->newsite = data.site;
      this->tableOfDividedDiffs_.resize(newsite.size());
      for(int i = 0; i < newsite.size(); i++)
      {
	this->tableOfDividedDiffs_[i].resize(newsite.size());
      }
      this->interPoly_.set_capacity(newsite.size());

      for(int i = 0; i < this->newsite.size(); i++)
      {
	this->tableOfDividedDiffs_[i][0] = data.value[i];
	for(int j = 1; j <= i; j++)
	{
	  this->tableOfDividedDiffs_[i][j] = (this->tableOfDividedDiffs_[i][j-1] - this->tableOfDividedDiffs_[i-1][j-1]) / (data.site[i] - data.site[i-j]);
	}
      }
    }
    
  }

  /** 
   * Implementation of Newton Formula
   * Interpolating from scratch to overwrite current data members
   * 
   * 
   * @return interpolating polynomial
   */
  Polynomial NewtonFormula()
  {
    this->symmPoly.set_capacity(0);
    Polynomial poly(1);
    Polynomial p1(0);
    Polynomial coef(0);

    this->interPoly_.set_to_zero();
    
    for(int i = 0; i < this->tableOfDividedDiffs_.size(); i++)
    {
      // Set polynomial to be x-0
      poly.set_coeff(0,0);
      poly.set_coeff(1,0);
      if(i == 0)
      {
	poly.set_coeff(0,1);
	//std::cout << "poly : ";
	//poly.print();
	//std::cout << std::endl;
	symmPoly.set_coeff(0,1);
      }
      else
      {
	// Set polynomial to be x-x_(i-1)
	poly.set_coeff(1,1);
	poly.set_coeff(0, (-1)*this->newsite[i-1]);
	//std::cout << "poly : ";
	//poly.print();
	//std::cout << std::endl;
	symmPoly.set_capacity(i-1);
	symmPoly = symmPoly*poly;
      }
      //std::cout << i << " : tmpPoly = ";
      //tmpPoly.print();
      //std::cout << std::endl;
      p1.set_capacity(i-1);
      p1.set_to_zero();
      p1 = symmPoly;
      //std::cout << "p1 : ";
      //p1.print();
      //std::cout << std::endl;

      coef.set_coeff(0,0);
      coef.set_coeff(0, this->tableOfDividedDiffs_[i][i]);
      p1 = coef*p1;
      interPoly_ = interPoly_ + p1;
      //std::cout << "interPoly_ : ";
      //interPoly_.print();
      //std::cout << std::endl << std::endl;
    }
    return interPoly_;
  }

  /** 
   * Incrementally build upon a polynomial obtained from previous
   * interpolation conditions
   * 
   * @param previous the previous interpolation conditions
   * @param data new interpolation data
   * 
   * @return new interpolating polynomial
   */
  Polynomial NewtonFormula(NewtonInterp previous, InterpConditions data)
  {
    this->tableOfDividedDiffs_ = previous.tableOfDividedDiffs_;
    this->interPoly_ = previous.interPoly_;
    int lastline = this->tableOfDividedDiffs_.size() - 1;
    if(data.isHermite() == true)
    {
      this->newsite.resize(this->newsite.size() + (2*data.site.size()));
      this->tableOfDividedDiffs_.resize(this->newsite.size());
      for(int i = 0; i < newsite.size(); i++)
      {
	this->tableOfDividedDiffs_[i].resize(newsite.size());
      }
      this->interPoly_.set_capacity(this->newsite.size());

      for(int i = lastline; i < lastline + data.site.size(); i++)
      {
	this->newsite[2*i] = data.site[i];
	this->newsite[(2*i)+1] = data.site[i];
	this->tableOfDividedDiffs_[2*i][0] = data.value[i];
	this->tableOfDividedDiffs_[(2*i)+1][0] = data.value[i];
	this->tableOfDividedDiffs_[(2*i)+1][1] = data.der[i];

	if(i != lastline)
	{
	  this->tableOfDividedDiffs_[2*i][1] = (this->tableOfDividedDiffs_[2*i][0] - this->tableOfDividedDiffs_[2*i-1][0]) / (this->newsite[2*i] - this->newsite[2*i-1]);
	}
      }

      for(int i = lastline + 2; i < this->newsite.size(); i++)
      {
	for(int j = 2; j <= i; j++)
	{
	  this->tableOfDividedDiffs_[i][j] = (this->tableOfDividedDiffs_[i][j-1] - this->tableOfDividedDiffs_[i-1][j-1]) / (this->newsite[i] - this->newsite[i-j]);
	}
      }
    }
    else
    {
      int lastline = this->newsite.size();
      this->newsite.resize(lastline + data.site.size());
      for(int i = lastline; i < this->newsite.size(); i++)
      {
	this->newsite[i] = data.site[i-lastline];
      }
      this->tableOfDividedDiffs_.resize(this->newsite.size());
      for(int i = lastline; i < this->newsite.size(); i++)
      {
	this->tableOfDividedDiffs_[i][0] = data.value[i];
	for(int j = 1; j <= 1; j++)
	{
	  this->tableOfDividedDiffs_[i][j] = (this->tableOfDividedDiffs_[i][j-1] - this->tableOfDividedDiffs_[i-1][j-1]) / (this->newsite[i] - this->newsite[i-j]);
	}
      }
    }

    this->symmPoly = previous.symmPoly;
    Polynomial poly(0);
    Polynomial coef(0);
    Polynomial p1(lastline);
    for(int i = lastline; i < this->tableOfDividedDiffs_.size(); i++)
    {
      poly.set_coeff(0,0);
      poly.set_coeff(1,0);
      if(i == 0)
      {
	poly.set_coeff(0,1);
	//std::cout << "poly : ";
	//poly.print();
	//std::cout << std::endl;
	symmPoly.set_coeff(0,1);
      }
      else
      {
	poly.set_coeff(1,1);
	poly.set_coeff(0, (-1)*this->newsite[i-1]);
	//std::cout << "poly : ";
	//poly.print();
	//std::cout << std::endl;
	symmPoly.set_capacity(i-1);
	symmPoly = symmPoly*poly;
      }
      //std::cout << i << " : tmpPoly = ";
      //tmpPoly.print();
      //std::cout << std::endl;
      p1.set_capacity(i-1);
      p1.set_to_zero();
      p1 = symmPoly;
      //std::cout << "p1 : ";
      //p1.print();
      //std::cout << std::endl;

      coef.set_coeff(0,0);
      coef.set_coeff(0, this->tableOfDividedDiffs_[i][i]);
      p1 = coef*p1;
      interPoly_ = interPoly_ + p1;
      //std::cout << "interPoly_ : ";
      //interPoly_.print();
      //std::cout << std::endl << std::endl;
    }

    return interPoly_;
  }

  /** 
   * Implementation of Neville-Aitken algorithm
   * 
   * @param data_ interpolation data
   * @param k interpolate until x_(i+k)
   * @param i interpolate from x_i
   * 
   * @return 
   */
  static Polynomial NAFormula(InterpConditions data_, int k, int i)
  {
    Polynomial p1(1), p2(1);
    p1.set_coeff(0, (-1)*data_.site[i]);
    p1.set_coeff(1, 1);
    p2.set_coeff(0, (-1)*data_.site[i+k]);
    p2.set_coeff(1, 1);

    if(k+1 == 0)
      return data_.value[i];
    else
    {
       Polynomial poly(k+1);
       poly = (p1*NAFormula(data_,k-1,i+1) - p2*NAFormula(data_,k-1,i));
       double deno = data_.site[i+k] - data_.site[i];
       for(int i = 0; i<=k+1; i++)
       {
	 poly.set_coeff(i, poly.get_coeff(i)/deno);
       }
       return poly;
    }
  }

  /** 
   * Neville-Aitken algorithm for i = 0 and k = n, 
   * which returns the interpolating polynomial of degree n
   * for the function f at the points x_0,...,x_n
   * 
   * @param data interpolation data
   * 
   * @return interpolating polynomial
   */
  static Polynomial NevilleAitken(InterpConditions data)
  {
    NAPoly_ = NAFormula(data, data.site.size(), 0);
    return NAPoly_;
  }
  
};


#endif
