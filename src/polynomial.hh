#ifndef __POLYNOMIAL_HH__
#define __POLYNOMIAL_HH__

/** @page polynomial Polynomial
 *  @brief Polynomial description and functions
 *  @tableofcontents
 *
 *  @section polynomial-description Description
 *
 *  This part describes a polynomial\<C\> : each @ref polynomial<C> is mainly represented by an inner vector of pairs, where each pair represents a @e monomial.
 *
 *  Each monomial is composed (as a pair) by those two elements :
 *  @li a vector of int, to store the different degrees of variables @e x, @e y @e z (variables names are stored in the static @ref polynomial::_var_names list).
 *  @li a templated coefficient \<C\> (that could be, for example, a \<long long\> or an @ref infint<long long>) : this last one could be @b negative @b if C @b allows @b them.
 *
 *  The set of monomials is stored in the @ref polynomial::_coefs vector. The set of variables used is stored in the static variable @ref polynomial::_var_names list.
 *
 *  Several methods are also proposed to build (@ref polynomial-construction), or operate (@ref polynomial-operators) theses polynomials,
 *
 *  @section polynomial-construction Construction
 *
 *  Three methods are proposed to build or read polynomials :
 *
 *    @li @ref polynomial::polynomial(int u)   gives a constant polynomial : degree is zero for  @e x, @e y @e z ... so the monomial coefficient is the integer C(u)
 *    @li @ref polynomial::polynomial(C c, std::vector<int> &pows)  gives a monomial polynomial : degree is set for @e x, @e y @e z with pows, ... and the monomial coefficient is the integer C(u)
 *    @li @ref polynomial::polynomial(const polynomial<C> & other)  gives a separate copy constructor
 *    @li @ref polynomial::operator= (const polynomial<C> & other ) gives an assignemnt operator to an already built polynomial
 *    @li @ref polynomial::operator>>(istream& is, polynomial<U>& c) can also be used to read a polynomial : variables will be added to the @ref polynomial::_var_names list if not already present
 *
 *  @section polynomial-operators Operators
 *
 *  Three methods are proposed to compute new polynomials from already built ones :
 *
 *  @li @ref polynomial::operator+()
 *  @li @ref polynomial::operator-()
 *  @li @ref polynomial::operator*()
 *  @li @ref polynomial::div()
 *  @li @ref polynomial::operator/()
 *  @li @ref polynomial::operator%()
 *
 * all are depending on the operators +, - or * that are defined for \<C\> (this could be a @e min, @e plus semi-ring for example)
 *
 * @todo FIXME : to be continued
 *
 */

/** @defgroup polynomial polynomial class template
 *  @brief polynomials defined over the templated C coefficients (int, long int, long long int, or other ...) with several variables
 *
 *
 *  @see matrix
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
 * @class polynomial
 * @tparam C
 * @brief polynomials are defined over the templated C coefficients (int, long int, long long int, or other ...) ; this template can handle multivariate polynomials if needed
 *        and is provided for positive or negative coefficients if C enables them
 */

template<typename C> bool sortbypowervec(const std::pair<std::vector<int>, C > & coef1, const std::pair<std::vector<int>, C > &coef2) {
  typename std::vector<int>::const_iterator i_l = coef1.first.begin(),  i_r = coef2.first.begin();
  while (i_l != coef1.first.end() && i_r != coef2.first.end()) {
    if ((*i_l) < (*i_r)) return true;
    if ((*i_l) > (*i_r)) return false;
    i_l++, i_r++;
  }
  while (i_l != coef1.first.end()) {
    if ((*i_l) < 0) return false;
    if ((*i_l) > 0) return true;
    i_l++;
  }
  while (i_r != coef2.first.end()) {
    if ((*i_r) < 0) return true;
    if ((*i_r) > 0) return false;
    i_r++;
  }
  return false;
}

template<typename C> class polynomial {


 public:
  /// Normalize a polynomial
  void normalize() {
    /* extends all the coef vectors to size "polynomial<C>::_var_names.size()" by adding extra zeros, sort them */
    for (typename std::vector<std::pair<std::vector<int>, C > >::iterator i_p = _coefs.begin(); i_p != _coefs.end(); i_p++)
      i_p->first.resize(polynomial<C>::_var_names.size(),0);
    sort(_coefs.begin(), _coefs.end(),sortbypowervec<C>);
  }


  /// Build a constant polynomial
  polynomial(int u = 0) {
    if (u != 0)
      _coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(polynomial<C>::_var_names.size(),0), C(u)));
    normalize();
  }

  /// Build a constant polynomial
  polynomial(C c, std::vector<int> &pows) {
    if (c != C(0)) {
      std::pair<std::vector<int>, C> coef;
      coef.second = C(c);
      for (unsigned u = 0; u < pows.size(); u++) {
        coef.first.push_back(pows[u]);
      }
      coef.first.resize(polynomial<C>::_var_names.size(),0);
      _coefs.push_back(coef);
    }
    normalize();
  }


  /// Copy Constructor
  polynomial(const polynomial<C> &other) {
    // Manual Deep Copy
    for (unsigned i = 0; i < other._coefs.size(); i++) {
      std::pair<std::vector<int>, C> coef;
      coef.second = C(other._coefs[i].second);
      for (unsigned u = 0; u < other._coefs[i].first.size(); u++) {
        coef.first.push_back(other._coefs[i].first[u]);
      }
      coef.first.resize(polynomial<C>::_var_names.size(),0);
      _coefs.push_back(coef);
    }
    normalize();
  }

  /// Assignment Operator
  polynomial<C> & operator= (const polynomial<C> & other ) {
    // Deep cleaning of previous elements cumulated in "this"
    for(unsigned i = 0; i < _coefs.size(); i++)
      _coefs[i].first.clear();
    _coefs.clear();

    // Manual Deep Copy
    for (unsigned i = 0; i < other._coefs.size(); i++) {
      std::pair<std::vector<int>, C> coef;
      coef.second = C(other._coefs[i].second);
      for (unsigned u = 0; u < other._coefs[i].first.size(); u++) {
        coef.first.push_back(other._coefs[i].first[u]);
      }
      coef.first.resize(polynomial<C>::_var_names.size(),0);
      _coefs.push_back(coef);
    }

    normalize();
    return *this;
  }


  /// Destructor
  ~polynomial() {
    for(unsigned i = 0; i < _coefs.size(); i++)
      _coefs[i].first.clear();
    _coefs.clear();
  }



  // @{
  /// Operator @f$ + @f$ for two polynomials
  template<typename U> friend polynomial<U> operator+ (const polynomial<U> & l, const polynomial<U> & r);
  /// Operator @f$ - @f$ for two polynomials
  template<typename U> friend polynomial<U> operator- (const polynomial<U> & l, const polynomial<U> & r);
  /// Operator @f$ \times @f$ for two polynomials
  template<typename U> friend polynomial<U> operator* (const polynomial<U> & l, const polynomial<U> & r);
  /// Operator @f$ div @f$ for two polynomials
  template<typename U> friend std::pair<polynomial<U>,polynomial<U> > div(const polynomial<U> & l, const polynomial<U> & r);
  /// Operator @f$ / @f$ for two polynomials
  template<typename U> friend polynomial<U> operator/ (const polynomial<U> & l, const polynomial<U> & r);
  /// Operator @f$ % @f$ for two polynomials
  template<typename U> friend polynomial<U> operator% (const polynomial<U> & l, const polynomial<U> & r);
  /// Operator @f$ gcd @f$ for two polynomials
  template<typename U> friend polynomial<U> gcd(const polynomial<U> & l, const polynomial<U> & r);

  /// Operator @f$ != @f$ for two polynomials
  template<typename U> friend bool    operator!= (const polynomial<U> & l, const polynomial<U> & r);
  /// Operator @f$ == @f$ for two polynomials
  template<typename U> friend bool    operator== (const polynomial<U> & l, const polynomial<U> & r);
  /// Operator @f$ < @f$ for two polynomials, given to provide an order "without evaluation", only here to classify polynomials
  template<typename U> friend bool    operator< (const polynomial<U> & l, const polynomial<U> & r);
  /// Operator @f$ > @f$ for two polynomials, given to provide an order "without evaluation", only here to classify polynomials
  template<typename U> friend bool    operator> (const polynomial<U> & l, const polynomial<U> & r);
  // @}

  // @{
  /// Print a polynomial
  template<typename U> friend ostream& operator<< (ostream& os, const polynomial<U>& c);
  /// Load a polynomial
  template<typename U> friend istream& operator>> (istream& is, polynomial<U>& c);
  // @}

protected:
  /// give a list of the coefficients C for a polynomial, each vector associated with C represents the list of variables (each variable being addressed by its proper index in the vector) and its degree (integer given by the vector at the index given by the variable)
  std::vector<std::pair<std::vector<int>, C> > _coefs;
  /// give the index of each variable name, or the name of the variable from its index
  static std::vector<string> _var_names;
};


/// Operator @f$ == @f$ for two Monomials POWS ONLY
template<typename C> inline bool equal_pows (const std::pair<std::vector<int>, C> & coef1, const std::pair<std::vector<int>, C> & coef2) {
  typename std::vector<int>::const_iterator i_l = coef1.first.begin(),  i_r = coef2.first.begin();
  while (i_l != coef1.first.end() && i_r != coef2.first.end()) {
    if ((*i_l) != (*i_r)) return false;
    i_l++, i_r++;
  }
  while (i_l != coef1.first.end()) {
    if ((*i_l) != 0) return false;
    i_l++;
  }
  while (i_r != coef2.first.end()) {
    if ((*i_r) != 0) return false;
    i_r++;
  }
  return true;
}

/// Operator @f$ diff_pows @f$ for two Monomial POWS ONLY
template<typename C> inline std::vector<int> diff_pows (const std::pair<std::vector<int>, C> & coef1, const std::pair<std::vector<int>, C> & coef2) {
  std::vector<int> result;
  typename std::vector<int>::const_iterator i_l = coef1.first.begin(),  i_r = coef2.first.begin();
  while (i_l != coef1.first.end() && i_r != coef2.first.end()) {
    result.push_back((*i_l) - (*i_r));
    i_l++, i_r++;
  }
  while (i_l != coef1.first.end()) {
    result.push_back(*i_l);
    i_l++;
  }
  while (i_r != coef2.first.end()) {
    result.push_back(-(*i_r));
    i_r++;
  }
  return result;
}

/// Operator @f$ computable_diff_pows @f$ for two Monomial POWS ONLY
template<typename C> inline bool computable_diff_pows (const std::pair<std::vector<int>, C> & coef1, const std::pair<std::vector<int>, C> & coef2) {
  typename std::vector<int>::const_iterator i_l = coef1.first.begin(),  i_r = coef2.first.begin();
  while (i_l != coef1.first.end() && i_r != coef2.first.end()) {
    if ((*i_l) < (*i_r)) return false;
    i_l++, i_r++;
  }
  while (i_l != coef1.first.end()) {
    if ((*i_l) < 0) return false;
    i_l++;
  }
  while (i_r != coef2.first.end()) {
    if ((*i_r) > 0) return false;
    i_r++;
  }
  return true;
}


/// Operator @f$ cumul_pows @f$ for two Monomial POWS ONLY
template<typename C> inline std::vector<int> cumul_pows (const std::pair<std::vector<int>, C> & coef1, const std::pair<std::vector<int>, C> & coef2) {
  std::vector<int> result;
  typename std::vector<int>::const_iterator i_l = coef1.first.begin(),  i_r = coef2.first.begin();
  while (i_l != coef1.first.end() && i_r != coef2.first.end()) {
    result.push_back((*i_l) + (*i_r));
    i_l++, i_r++;
  }
  while (i_l != coef1.first.end()) {
    result.push_back(*i_l);
    i_l++;
  }
  while (i_r != coef2.first.end()) {
    result.push_back(-(*i_r));
    i_r++;
  }
  return result;
}


/// Operator @f$ == @f$ for two polynomials
template<typename C> inline bool equal_coefs (const std::pair<std::vector<int>, C> & coef1, const std::pair<std::vector<int>, C> & coef2) {
  return coef1.second == coef2.second && equal_pows(coef1, coef2);
}




/// Operator @f$ == @f$ for two polynomials
template<typename C> inline bool    operator== (const polynomial<C> & l, const polynomial<C> & r) {
  return l._coefs.size() == r._coefs.size() && std::equal(l._coefs.begin(), l._coefs.end(), r._coefs.begin(), equal_coefs<C>);
}

/// Operator @f$ != @f$ for two polynomials
template<typename C> inline bool    operator!= (const polynomial<C> & l, const polynomial<C> & r) {
  return !(l._coefs.size() == r._coefs.size() && std::equal(l._coefs.begin(), l._coefs.end(), r._coefs.begin(), equal_coefs<C>));
}

/// Operator @f$ < @f$ for two polynomials
template<typename C> inline bool    operator< (const polynomial<C> & l, const polynomial<C> & r) {
  return std::lexicographical_compare(l._coefs.begin(), l._coefs.end(), r._coefs.begin(), r._coefs.end());
}

/// Operator @f$ > @f$ for two polynomials
template<typename C> inline bool    operator> (const polynomial<C> l, const polynomial<C> r) {
  return !(std::lexicographical_compare(l._coefs.begin(), l._coefs.end(), r._coefs.begin(), r._coefs.end()));
}


/// Operator + for two polynomials
template<typename C> inline polynomial<C> operator+ (const polynomial<C> & l, const polynomial<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] + [" << r << "]"););
  polynomial<C> result = polynomial<C>();
  typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_l = l._coefs.begin();
  typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_r = r._coefs.begin();
  while (i_l != l._coefs.end() && i_r != r._coefs.end()) {
    if (equal_pows(*i_l,*i_r)) {
      C val = C(i_l->second) + C(i_r->second);
      std::vector<int> degree_vector = std::vector<int>(i_l->first);
#ifndef USEINFINT
      if ( (i_l->second > 0) && (i_r->second > 0) && (val < 0) ) {
        _ERROR("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","polynomial<C> operator+ DOES OVERFLOW ...");
      }
#endif
      if (val != C(0))
        result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(degree_vector), C(val)));
      i_l++;
      i_r++;
    } else {
      if (i_l->first < i_r->first) {
        C val = C(i_l->second);
        std::vector<int> degree_vector = std::vector<int>(i_l->first);
        if (val != C(0))
          result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(degree_vector), C(val)));
        i_l++;
      } else {
        C val = C(i_r->second);
        std::vector<int> degree_vector = std::vector<int>(i_r->first);
        if (val != C(0))
          result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(degree_vector), C(val)));
        i_r++;
      }
    }
  }
  while (i_l != l._coefs.end()) {
    C val = C(i_l->second);
    std::vector<int> degree_vector = std::vector<int>(i_l->first);
    if (val != C(0))
      result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(degree_vector), C(val)));
    i_l++;
  }
  while (i_r != r._coefs.end()) {
    C val = C(i_r->second);
    std::vector<int> degree_vector = std::vector<int>(i_r->first);
    if (val != C(0))
      result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(degree_vector), C(val)));
    i_r++;
  }
  result.normalize();
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("\t = [" << result << "]"););
  return result;
}



/// Operator - for two polynomials
template<typename C> inline polynomial<C> operator- (const polynomial<C> & l, const polynomial<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] - [" << r << "]"););
  polynomial<C> result = polynomial<C>();
  typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_l = l._coefs.begin();
  typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_r = r._coefs.begin();
  while (i_l != l._coefs.end() && i_r != r._coefs.end()) {
    if (equal_pows(*i_l,*i_r)) {
      C val = C(i_l->second) - C(i_r->second);
      std::vector<int> degree_vector = std::vector<int>(i_l->first);
      if (val != C(0))
        result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(degree_vector), C(val)));
      i_l++;
      i_r++;
    } else {
      if (i_l->first < i_r->first) {
        C val = C(i_l->second);
        std::vector<int> degree_vector = std::vector<int>(i_l->first);
        if (val != C(0))
          result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(degree_vector), C(val)));
        i_l++;
      } else {
        C val = C(0) - C(i_r->second);
        std::vector<int> degree_vector = std::vector<int>(i_r->first);
        if (val != C(0))
          result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(degree_vector), C(val)));
        i_r++;
      }
    }
  }
  while (i_l != l._coefs.end()) {
    C val = C(i_l->second);
    std::vector<int> degree_vector = std::vector<int>(i_l->first);
    if (val != C(0))
      result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(degree_vector), C(val)));
    i_l++;
  }
  while (i_r != r._coefs.end()) {
    C val = C(0) - C(i_r->second);
    std::vector<int> degree_vector = std::vector<int>(i_r->first);
    if (val != C(0))
      result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(degree_vector), C(val)));
    i_r++;
  }
  result.normalize();
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("\t = [" << result << "]"););
  return result;
}


/** Operator @f$ \times @f$ for two polynomials
 *  @todo{FIXME : this must be updated with Karatsuba, Toom-k or FFT/NNT methods (Strassen)}
 */
template<typename C> inline polynomial<C> operator* (const polynomial<C> & l, const polynomial<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] * [" << r << "]"););
  map<std::vector<int>, C> tmp_coefs;
  for (typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_l = l._coefs.begin(); i_l != l._coefs.end(); i_l++) {
    for (typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_r = r._coefs.begin(); i_r != r._coefs.end(); i_r++) {
      // build the degree_vector
      std::vector<int> degree_vector = cumul_pows(*i_l,*i_r);

      // check already existance in result (log search)
      typename map<std::vector<int>, C>::iterator low = tmp_coefs.lower_bound(degree_vector);
      if (low != tmp_coefs.end()) {
        if (low->first.size() == degree_vector.size() && std::equal(low->first.begin(),low->first.end(),degree_vector.begin())) {
#ifndef USEINFINT
          if ( (i_l->second > 0) && (i_r->second > 0) && ((i_l->second * i_r->second) < 0) ) {
            _ERROR("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","polynomial<C> operator* DOES OVERFLOW ...");
          }
          if ( (low->second > 0) && ((i_l->second * i_r->second) > 0) && (low->second + (i_l->second * i_r->second) < 0) ) {
            _ERROR("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","polynomial<C> operator* DOES OVERFLOW ...");
          }
#endif
          low->second = C(low->second) + C(i_l->second * i_r->second);
        } else {
#ifndef USEINFINT
          if ( (i_l->second > 0) && (i_r->second > 0) && ((i_l->second * i_r->second) < 0) ) {
            _ERROR("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","polynomial<C> operator* DOES OVERFLOW ...");
          }
#endif
          tmp_coefs.insert(low, std::pair<std::vector<int>,C >(std::vector<int>(degree_vector), C(i_l->second * i_r->second)));
        }
      } else {
#ifndef USEINFINT
        if ( (i_l->second > 0) && (i_r->second > 0) && ((i_l->second * i_r->second) < 0) ) {
          _ERROR("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","polynomial<C> operator* DOES OVERFLOW ...");
        }
#endif
        tmp_coefs.insert(std::pair<std::vector<int>, C >(std::vector<int>(degree_vector), C(i_l->second * i_r->second)));
      }
    }
  }

  // simplify zeros
  polynomial<C> result = polynomial<C>();
  for(typename map<std::vector<int>, C>::iterator it = tmp_coefs.begin(); it != tmp_coefs.end(); ++it)
    if (it->second != C(0))
      result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(it->first),C(it->second)));

  tmp_coefs.clear();
  result.normalize();

  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("\t = [" << result << "]"););
  return result;
}


/** Operator @f$ \div @f$ "with remainder" for two multivariate polynomials
 */
template<typename C> inline std::pair<polynomial<C>,polynomial<C> > div(const polynomial<C> & l, const polynomial<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] div [" << r << "]"););
  polynomial<C> p = polynomial<C>(l);
  polynomial<C> d = polynomial<C>(r);
  if (d._coefs.size() == 0)
    return std::pair<polynomial<C>,polynomial<C> >(polynomial<C>(),polynomial<C>(l)); // To avoid div by zero ... this is an artefact
  if (p._coefs.size() == 0)
    return std::pair<polynomial<C>,polynomial<C> >(polynomial<C>(),polynomial<C>());  // "zero" / whatever" = "zero" ...
  polynomial<C> res = polynomial<C>(); // res for result (and not rest) is set to "zero"

  std::pair<std::vector<int>, C> d_high = d._coefs.back();
  std::pair<std::vector<int>, C> p_high = p._coefs.back();

  /// Compute the difference between degrees and check if any < number
  while (computable_diff_pows(p_high,d_high)) {

    /// If OK, build C_p_over_d  * x^diff_x * y^diff_y
    std::vector<int> p_minus_d_pows = diff_pows(p_high,d_high);
    C                p_minus_d_coef = C(p_high.second) / C(d_high.second);
    // Check if the coefficient is not zero ... otherwise, we must stop
    if (p_minus_d_coef == C(0))
      break;

    // Build the monomial
    polynomial<C> q = polynomial<C>(p_minus_d_coef,p_minus_d_pows);
    // Add it to the result
    res = res + q;
    // And substract it from "p" (after multiplying it by "d")
    p   = p -  q * d;
    //cerr << "\t\t" <<l << " =  (" << r << ") * (" << res << ") + (" << p << ")" << endl;
    p_minus_d_pows.clear();

    // Check if "p" cannot be reduced anymore
    if (p._coefs.size() == 0)
      break;
    else
      p_high = p._coefs.back();
  }
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("\res = [" << res << "], mod = [" << p << "]"););
  return std::pair<polynomial<C>,polynomial<C> >(polynomial<C>(res),polynomial<C>(p));
}


/** Partial operator @f$ / @f$ (div without remainer) for two multivariate polynomials
 */
template<typename C> inline polynomial<C> operator/ (const polynomial<C> & l, const polynomial<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] / [" << r << "]"););
  return (div(l,r)).first;
}

/** Partial operator @f$ % @f$ (modulus or mod) for two multivariate polynomials
 */
template<typename C> inline polynomial<C> operator% (const polynomial<C> & l, const polynomial<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] / [" << r << "]"););
  return (div(l,r)).second;
}

/** Partial operator @f$ % @f$ (modulus or mod) for two multivariate polynomials
 */
template<typename C> inline polynomial<C> gcd (const polynomial<C> & l, const polynomial<C> & r) {
  // https://mathsci2.appstate.edu/~cookwj/sage/algebra/Euclidean_algorithm-poly.html#:~:text=The%20Polynomial%20Euclidean%20Algorithm%20computes,gcd%20(%20b%20%2C%20r%20)%20.
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] gcd [" << r << "]"););
  cerr << endl;
  polynomial<C> p = polynomial<C>(l);
  polynomial<C> d = polynomial<C>(r);
  if (polynomial<C>(d) == polynomial<C>(0) || polynomial<C>(p) == polynomial<C>(0))
    return polynomial<C>(1);
  set<polynomial<C> > already_seen;
  while (polynomial<C>(d) != polynomial<C>(0)) {
    polynomial<C> tmp = p % d;
    // for polynomials loops when no reduction occurs ...
    if (already_seen.find(tmp) == already_seen.end()) {
      already_seen.insert(tmp);
    } else {
      already_seen.clear();
      return polynomial<C>(1);
    }
    p = d;
    d = tmp;
  }
  return polynomial<C>(p);
}




/// Print a polynomial
template<typename C> ostream& operator<< (ostream& os, const polynomial<C> & p) {
  if (p._coefs.size() == 0)
    os << C(0);

  for (unsigned i = 0 ; i < p._coefs.size() ; i++) {
    C sign = C(1);
    if (p._coefs[i].second < C(0))
      sign = C(-1);

    if (sign < C(0))
      os << " - ";
    else
      if (i > 0)
        os << " + ";

    bool first_prod = true;
    // dont print '1 *'  "something ..."
    C abs_coef = sign * (p._coefs[i].second);
    if (abs_coef != C(1)) {
      os << abs_coef;
      first_prod = false;
    }
    for (unsigned j = 0 ; j < p._coefs[i].first.size(); j++) {
      if (p._coefs[i].first[j] != 0) {
        if (first_prod)
          first_prod = false;
        else
          os << " * ";
        os << polynomial<C>::_var_names[j];
        if (p._coefs[i].first[j] != 1)
          os << " ^ " << p._coefs[i].first[j];
      }
    }
    if (first_prod) {
      os << C(1);
    }
  }
  return os;
}

/// Load a polynomial ( \<C\>"coef" * variable1 [^power1]  [* variable2 [^power2] ... ] + ...)
//  @todo{FIXME : this kind of Ugly and Non secure code will not be used from command lined iedera !! unless one uses already generated polynomials from this code ...}
template<typename C> istream& operator>> (istream& is, polynomial<C> & p) {
  char c;
  /// for parsing and keeping
  C sign = C(1);
  C coef = C(1);
  std::vector<int> var_degree;

  p._coefs          = std::vector<std::pair<std::vector<int>, C> > ();

  // read sign
 read_sign:
  c = 0xff;
  while (is.get(c) && (c == ' ' || c == '\t'));
  if (is.eof()) {
    _ERROR("operator>>"," invalid start for a sign/variable/coefficient" << endl << "\t polynom format : <C>\"coef\" * variable1 [^power1]  [* variable2 [^power2] ... ] +  <C>\"coef\" * ..." << endl);
  }
  if (c >= '0' && c <= '9') {
    is.unget();
    goto read_coef;
  } else {
    if (c == '+' || c == '-') {
      if (c == '-') {
        sign = sign * C(-1);
      } else {
        sign = sign * C(1);
      }
      goto read_sign;
    } else {
      is.unget();
      coef = C(1);
      goto read_variable;
    }
  }


 read_coef:
  long long int coef_long_long;
  is >> coef_long_long;
  coef = C(coef_long_long);
  if (is.fail()) {
    _ERROR("operator>>"," invalid coefficient" << endl << "\t polynom format : <C>\"coef\" * variable1 [^power1]  [* variable2 [^power2] ... ] +  <C>\"coef\" * ..." << endl);
  }
  //cerr << "coef:" << coef << endl;

 read_product:
  c = 0xff;
  while (is.get(c) && (c == ' ' || c == '\t'));
  if (is.eof())
    goto end_of_polynomial;
  else if (c == '*')
    goto read_variable;
  else if (c == '+' || c == '-') {
    is.unget();
    p._coefs.push_back(std::pair<std::vector<int>, C>(var_degree, coef * sign));
    var_degree.clear();
    sign = C(1);
    coef = C(1);
    goto read_sign;
  } else {
    is.unget();
    goto end_of_polynomial;
  }

 read_variable:
    {
      string var_symbol;
      is >> var_symbol;
      for(string::iterator it = var_symbol.begin(); it != var_symbol.end(); it++) {
        if ( ! ((*it >= 'a' && *it <= 'z') || (*it >= 'A' && *it <= 'Z') || *it == '_')) {
          _ERROR("operator>>"," variable name \""<< var_symbol<<"\" include non alpha symbols \'"<< (*it) <<"\'" << endl << "\t polynom format : <C>\"coef\" * variable1 [^power1]  [* variable2 [^power2] ... ] +  <C>\"coef\" * ..." << endl);
        }
      }



    //cerr << "variable:" << var_symbol << endl;
    int i_var = -1;

    // extend var_degree up to this variable if it is already defined
    bool defined_var = false;
    for (unsigned i = 0; i < polynomial<C>::_var_names.size(); i++) {
      if (!polynomial<C>::_var_names[i].compare(var_symbol)) {
        i_var = (int) i;
        defined_var = true;
        // already in the monomial list, so no need to extend the monomial list
        if (i_var < (int) var_degree.size()) {
          if (var_degree[i_var] != 0) { // check if possibly not at zero, because this is a classical error is to set it twice
            _ERROR("operator>>"," variable \""<< var_symbol<<"\" occuring twice in a monomial" << endl << "\t polynom format : <C>\"coef\" * variable1 [^power1]  [* variable2 [^power2] ... ] +  <C>\"coef\" * ..." << endl);
          } else {
            var_degree[i_var] = 1;
          }
        } else {
          // or need to extend it up to "i_var"
        for (int j = (int) var_degree.size(); j < i_var; j++)
          var_degree.push_back(0);
        var_degree.push_back(1);
        }
        break;
      }
    }

    // if not defined, create its name and extend "var_degree"
    if (!defined_var) {
      unsigned i = polynomial<C>::_var_names.size();
      i_var = (int) i;
      polynomial<C>::_var_names.push_back(var_symbol);
      //cerr << "new symbol " << var_symbol <<  " at index "  << (polynomial<C>::_var_names.size()) << endl;
      for (unsigned j = var_degree.size(); j < i; j++)
        var_degree.push_back(0);
      var_degree.push_back(1);
    }


    // read_product_or_power_or_sign
    c = 0xff;
    while (is.get(c) && (c == ' ' || c == '\t'));
    if (is.eof())
      goto end_of_polynomial;
    else if (c == '*')
      goto read_variable;
    else if (c == '+' || c == '-') {
      is.unget();
      p._coefs.push_back(std::pair<std::vector<int>, C>(var_degree, coef * sign));
      var_degree.clear();
      sign = C(1);
      coef = C(1);
      goto read_sign;
    } else if (c == '^') {
      int power_value = 0;
      is >> power_value;
      if (is.fail()) {
        _ERROR("operator>>"," invalid power coefficient" << endl << "\t polynom format : <C>\"coef\" * variable1 [^ power1]  [* variable2 [^ power2] ... ] +  <C>\"coef\" * ..." << endl);
      }
      //cerr << "power:" << power_value << endl;
      var_degree[i_var] = power_value;
      goto read_product;
    } else {
      is.unget();
      goto end_of_polynomial;
    }
  }

 end_of_polynomial:
  p._coefs.push_back(std::pair<std::vector<int>, C>(var_degree, coef * sign));
  p.normalize();
  return is;
}

template<typename C> std::vector<string> polynomial<C>::_var_names = std::vector<string>(0);
// @}

#endif
