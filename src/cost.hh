#ifndef __COST_HH__
#define __COST_HH__

/** @defgroup cost cost class template
 *  @brief costs defined over the templated C elements (int, long int, long long int, or other ...),
 *  and in the @f$ (\oplus = min, \otimes = plus) @f$ tropical semi-ring
 *
 *
 *  @see matrix, automaton
 */

// @{

//STL
#include <functional>
#include <algorithm>
#include <vector>
//STD
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
//STR
using namespace std;

/**
 * @class cost
 * @tparam C
 * @brief costs defined in the @f$ (\oplus = min, \otimes = plus) @f$
 *   tropical semi-ring
 *
 *  This class defines common operations in the @f$ (\oplus = min, \otimes = plus) @f$ tropical semi-ring to get
 *  compatibility issues with the classical  @f$ (\oplus = +, \otimes = \times) @f$ semi-ring used for probabilistic
 *  computations. As such, any algorithm programmed with  @f$ (\oplus, \otimes) @f$
 *  (C++ operators  @f$+@f$ and @f$\times@f$) can be used in both cases : either on costs, otherwise on more
 *  traditional floating point values...
 *
 *  @see matrix, automaton
 */

template<typename C> class cost {
 public:

  /// Build a cost
 cost(const C c): _c(c) {};

  /// Erase a cost
  ~cost() {};

  // @{
  /// Operator + (min) for two costs
  template<typename U> friend cost<U> operator+ (const cost<U> l,const cost<U> r);
  /// Operator @f$ \times @f$ (add) for two costs
  template<typename U> friend cost<U> operator* (const cost<U> l,const cost<U> r);
  /// Operator @f$ / @f$ (sub) for two costs
  template<typename U> friend cost<U> operator/ (const cost<U> l,const cost<U> r);
  /// Operator @f$ != @f$ for two costs
  template<typename U> friend bool    operator!= (const cost<U> l,const cost<U> r);
  /// Operator @f$ == @f$ for two costs
  template<typename U> friend bool    operator== (const cost<U> l,const cost<U> r);
  /// Operator @f$ < @f$ for two costs
  template<typename U> friend bool    operator< (const cost<U> l,const cost<U> r);
  /// Operator @f$ > @f$ for two costs
  template<typename U> friend bool    operator> (const cost<U> l,const cost<U> r);
  // @}

  // @{
  /// Print a cost
  template<typename U> friend ostream& operator<< (ostream& os, const cost<U>& c);
  /// Load a cost
  template<typename U> friend istream& operator>> (istream& is, cost<U>& c);
  // @}



protected:
  /// memorized cost
  C _c;
};

/// Operator + (min) for two costs
template<typename C> inline cost<C> operator+ (const cost<C> l, const cost<C> r) {
  cost<C> x(MIN(l._c,r._c));
  return x;
}

/// Operator @f$ \times @f$ (add) for two costs
template<typename C> inline cost<C> operator* (const cost<C> l, const cost<C> r) {
  cost<C> x(l._c + r._c);
  return x;
}

/// Operator @f$ / @f$ (sub) for two costs
template<typename C> inline cost<C> operator* (const cost<C> l, const cost<C> r) {
  cost<C> x(l._c - r._c);
  return x;
}

/// Operator @f$ == @f$ for two costs
template<typename C> inline bool    operator== (const cost<C> l, const cost<C> r) {
  return l._c == r._c;
}

/// Operator @f$ != @f$ for two costs
template<typename C> inline bool    operator!= (const cost<C> l, const cost<C> r) {
  return l._c != r._c;
}

/// Operator @f$ < @f$ for two costs
template<typename C> inline bool    operator< (const cost<C> l, const cost<C> r) {
  return l._c < r._c;
}

/// Operator @f$ > @f$ for two costs
template<typename C> inline bool    operator> (const cost<C> l, const cost<C> r) {
  return l._c > r._c;
}

/// Print a cost
template<typename C> ostream& operator<< (ostream& os, const cost<C>& c) {
  os << c._c;
  return os;
}

/// Load a cost
template<typename C> istream& operator>> (istream& is, cost<C>& c) {
  is >> c._c;
  return is;
}

// @}

#endif


