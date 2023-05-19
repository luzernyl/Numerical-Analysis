#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <utility>
#include <initializer_list>
#include <cmath>
#include <cassert>
#include <ostream>
#include <vector>

template <class T, int _Dim>
class Vec
{
 protected:
  //std::vector<T> coord(_Dim);
  T coord[_Dim];

 public:
  enum { Dim = _Dim};

  Vec(const T &t = T())
  {
    for(int i = 0; i < Dim; coord[i++] = t);
  }

  Vec(std::initializer_list<T> l)
  {
    auto j = l.begin();
    for(int d = 0; d < Dim; ++d)
      coord[d] =* j++;
  }

  //Conversion Constructor
  template <class T2>
  explicit Vec(const Vec<T2, Dim> &rhs)
  {
    for(int d  = 0; d < Dim; d++)
      coord[d] = static_cast<T>(rhs[d]);
  }

  static Vec<T,Dim> unit(int D)
  {
    Vec<T,Dim> r;
    r[D] = static_cast<T>(1);
    return r;
  }

  // accessors
 public:
  T &operator[](int _d) { return coord[_d]; }

  const T &operator[](int _d) const { return coord[_d]; }

  const T *data() const { return &coord[0]; }

 public :

#define ELMWISE_BINARY_OP(OpName, Op)          \
  template <class T2>                          \
  auto OpName(const Vec<T2,Dim> &rhs) const    \
  {                                            \
    using Tx = decltype(coord[0] Op rhs[0]);   \
    Vec<Tx,Dim> res;                           \
    for(int i = 0; i < Dim; i++)               \
      res[i] = coord[i] Op rhs[i];             \
    return res;                                \
  }

  ELMWISE_BINARY_OP(operator+, +)
  ELMWISE_BINARY_OP(operator-, -)
  ELMWISE_BINARY_OP(operator*, *)
  ELMWISE_BINARY_OP(operator/, /)
#undef ELMWISE_BINARY_OP
    
#define RIGHT_BROADCAST(OpName, Op)            \
    template <class T2>                        \
    auto OpName (const T2 &rhs) const          \
    {                                          \
      using Tx = decltype(coord[0] Op rhs);    \
      Vec<Tx,Dim> res;                         \
      for(int d = 0; d < Dim; ++d)             \
	res[d] = coord[d] Op rhs;              \
      return res;                              \
    }

  RIGHT_BROADCAST(operator+, +)
  RIGHT_BROADCAST(operator-, -)
  RIGHT_BROADCAST(operator*, *)
  RIGHT_BROADCAST(operator/, /)
#undef RIGHT_BROADCAST

    friend std::ostream &operator<<(std::ostream &os, const Vec<T,Dim> &p)
    {
      os << "(" << p[0];
      for(int d = 1; d < Dim; d++)
	os << "," << p[d];
      os << ")";
      return os;
    }
  
  template<class T2>
  void operator=(const Vec<T2,_Dim> &rhs)
  {
    for(int i = 0; i < Dim; i++)
      coord[i] = rhs.coord[i];
  }
};

#endif
