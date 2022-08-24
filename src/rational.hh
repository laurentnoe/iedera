#ifndef __RATIONAL_HH__
#define __RATIONAL_HH__

/** @page rational Rational
 *  @brief Rational description and functions
 *  @tableofcontents
 *
 *  @section rational-description Description
 *
 *  This part describes a rational\<C\> : each @ref rational<C> is mainly represented by a pair of @e C.
 *
 *  This pair is stored in the @ref rational::_num and rational::_den.
 *
 *  Several methods are also proposed to build (@ref rational-construction), or operate (@ref rational-operators) theses rational,
 *
 *  @section rational-construction Construction
 *
 *  Three methods are proposed to build or read rationals :
 *
 *    @li @ref rational::rational(const C & num = C(0), const C & den = C(1))     gives a constructor with one or two arguments
 *    @li @ref rational::rational(const rational<C> & r) gives a separate copy constructor
 *    @li @ref rational::operator= (const rational<C> & other ) gives an assignemnt operator to an already built rational
 *    @li @ref rational::operator>>(istream& is, rational<U>& c) can also be used to read a rational
 *
 *  @section rational-operators Operators
 *
 *  Three methods are proposed to compute new rationals from already built ones :
 *
 *  @li @ref rational::operator+()
 *  @li @ref rational::operator-()
 *  @li @ref rational::operator*()
 *  @li @ref rational::operator/()
 *
 * all are depending on the operators +, -,* or / that are defined for \<C\> (this could be a @e min, @e plus semi-ring for example)
 *
 * @todo FIXME : to be continued
 *
 */

/** @defgroup rational rational class template
 *  @brief rational defined over the templated C coefficients (int, long long int, polynomial<C> or other ...)
 *
 *
 *  @see polynomial
 */

// @{

//STL
#include <algorithm>
#include <vector>
#include <functional>
#include <string>
#include <limits>
//STD
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <set>
//STR
#include "macro.h"
#include "infint.hh"
using namespace std;

/**
 * @class rational
 * @tparam C
 * @brief rational are defined over the templated C coefficients (int, long int, long long int,  polynomial<C> or other ...) ; this template can handle multivariate polynomials if needed
 *        and is provided for positive or negative coefficients if C enables them
 */

template<typename C> class rational {
 public:

  /// Normalize a rational with a "gcd" algorithm (adapted to "polynomials" and "integers" : they have the "%" and "/" defined)
  void normalize() {
    C a = _num;
    C b = _den;
    if (C(b) == C(0))
      return; // NOT CORRECT ANYWAY...
    if (C(a) == C(0)) {
      _den = C(1);
      return;
    }

    set<C> already_seen;
    while (C(b) != C(0)) {
      C tmp = a % b;
      // for polynomials loops when no reduction occurs ...
      if (already_seen.find(tmp) == already_seen.end()) {
        already_seen.insert(tmp);
      } else {
        already_seen.clear();
        return;
      }
      a = b;
      b = tmp;
    }
    _num = _num / a;
    _den = _den / a;
    return;
  }


  /// Build a constant rational
  rational(const C & num = C(0), const C & den = C(1)) {
    _num = C(num);
    _den = C(den);
    normalize();
  };

  /// Copy Constructor
  rational(const rational<C> &other) {
    _num = C(other._num);
    _den = C(other._den);
    normalize();
  };

  /// Assignment Operator
  rational<C> & operator= (const rational<C> & other ) {
    _num = C(other._num);
    _den = C(other._den);
    normalize();
    return *this;
  }

  // @{
  /// Operator + for two rationals
  template<typename U> friend rational<U> operator+ (const rational<U> & l, const rational<U> & r);
  /// Operator -  for two rationals
  template<typename U> friend rational<U> operator- (const rational<U> & l, const rational<U> & r);
  /// Operator @f$ \times @f$ for two rationals
  template<typename U> friend rational<U> operator* (const rational<U> & l, const rational<U> & r);
  /// Operator /  for two rationals
  template<typename U> friend rational<U> operator/ (const rational<U> & l, const rational<U> & r);
  /// Operator @f$ \neq  @f$ for two rationals
  template<typename U> friend bool  operator!= (const rational<U> & l, const rational<U> & r);
  /// Operator @f$  == @f$ for two rationals
  template<typename U> friend bool  operator== (const rational<U> & l, const rational<U> & r);
  /// Operator @f$ < @f$ for two rationals, given to provide an order "without evaluation", only here to classify rationals
  template<typename U> friend bool  operator< (const rational<U> & l, const rational<U> & r);
  /// Operator @f$ > @f$ for two rationals, given to provide an order "without evaluation", only here to classify rationals
  template<typename U> friend bool  operator> (const rational<U> & l, const rational<U> & r);
  // @}

  // @{
  /// Print a rational with the string " { <num> } / { <den> } "
  template<typename U> friend ostream& operator<< (ostream& os, const rational<U>& c);
  /// Load a rational with the string " { <num> } / { <den> } "
  template<typename U> friend istream& operator>> (istream& is, rational<U>& c);
  // @}

protected:
  C _num,_den;
};

/// Operator @f$ == @f$ for two rationals
template<typename C> inline bool    operator== (const rational<C> & l, const rational<C> & r) {
  return l._num * r._den == r._num * l._den;
}

/// Operator @f$ != @f$ for two rationals
template<typename C> inline bool    operator!= (const rational<C> & l, const rational<C> & r) {
  return !(l == r);
}

/// Operator @f$ < @f$ for two positive rationals
template<typename C> inline bool    operator< (const rational<C> & l, const rational<C> & r) {
  return l._num * r._den < r._num * l._den;
}

/// Operator @f$ > @f$ for two positive rationals
template<typename C> inline bool    operator> (const rational<C> l, const rational<C> r) {
  return l._num * r._den > r._num * l._den;
}


/// Operator + for two rationals
template<typename C> inline rational<C> operator+ (const rational<C> & l, const rational<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] + [" << r << "]"););
  rational<C> result = rational<C>((l._num * r._den) + (r._num * l._den), l._den * r._den);
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("\t = [" << result << "]"););
  return result;
}

/// Operator - for two rationals
template<typename C> inline rational<C> operator- (const rational<C> & l, const rational<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] - [" << r << "]"););
  rational<C> result = rational<C>((l._num * r._den) - (r._num * l._den), l._den * r._den);
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("\t = [" << result << "]"););
  return result;
}

/// Operator * for two rationals
template<typename C> inline rational<C> operator* (const rational<C> & l, const rational<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] * [" << r << "]"););
  rational<C> result = rational<C>(l._num * r._num, l._den * r._den);
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("\t = [" << result << "]"););
  return result;
}

/// Operator / for two rationals
template<typename C> inline rational<C> operator/ (const rational<C> & l, const rational<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] * [" << r << "]"););
  rational<C> result = rational<C>(l._num * r._den, l._den * r._num);
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("\t = [" << result << "]"););
  return result;
}

/// Print a rational<C>
template<typename C> ostream& operator<< (ostream& os, const rational<C> & r) {
  os << "{" << (r._num) << "} / {" << (r._den) << "}";
  return os;
}

/// Load a rational ( { \<C\>(num) } / \<C\>(den) } )
template<typename C> istream& operator>> (istream& is, rational<C> & r) {
  char c;

  r._num = C(0);
  r._den = C(0);

  // read spaces until a symbol '}' occurs
  c = 0xff;
  while (is.get(c) && (c == ' ' || c == '\t'));
  if (is.eof() || c != '{') {
    _ERROR("operator>>"," invalid rational numerator not starting with '{'"<< endl << "\t rational format : { C(num) } / C(den) }" << endl);
  }

  // read numerator
  C num;
  is >> num;
  if (is.fail()) {
    _ERROR("operator>>"," invalid rational numerator coefficient" << endl << "\t rational format : { C(num) } / C(den) }" << endl);
  }
  r._num = C(num);
  cerr << "num:" << num << endl;

  // read spaces until a symbol '}' occurs
  c = 0xff;
  while (is.get(c) && (c == ' ' || c == '\t'));
  if (is.eof() || c != '}') {
    _ERROR("operator>>"," invalid rational numerator not starting with '}'" << endl << "\t rational format : { C(num) } / C(den) }" << endl);
  }

  // read spaces until a symbol '/' occurs
  c = 0xff;
  while (is.get(c) && (c == ' ' || c == '\t'));
  if (is.eof() || c != '/') {
    _ERROR("operator>>"," invalid rational separator in between '/'" << endl << "\t rational format : { C(num) } / C(den) }" << endl);
  }

  // read spaces until a symbol '{' occurs
  c = 0xff;
  while (is.get(c) && (c == ' ' || c == '\t'));
  if (is.eof() || c != '{') {
    _ERROR("operator>>"," invalid rational denominator not starting with '{'" << endl << "\t rational format : { C(num) } / C(den) }" << endl);
  }

  // read denominator
  C den;
  is >> den;
  if (is.fail()) {
    _ERROR("operator>>"," invalid rational denominator coefficient" << endl << "\t rational format : { C(num) } / { C(den) }" << endl);
  }
  r._den = C(den);
  cerr << "den:" << den << endl;

  // read spaces until a symbol '{' occurs
  c = 0xff;
  while (is.get(c) && (c == ' ' || c == '\t'));
  if (is.eof() || c != '}') {
    _ERROR("operator>>"," invalid rational denominator not ending with '}'" << endl << "\t rational format : { C(num) } / C(den) }" << endl);
  }

  r.normalize();
  return is;
}

// @}

#endif

