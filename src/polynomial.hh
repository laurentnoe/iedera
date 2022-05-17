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
 *  @li a templated coefficient \<C\> (that could be, for example, a \<long long\> or an @ref infint<long long>) : this last one could be @b positive @b is C @b allows @b them.
 *
 *  The set of monomials is stored in the @ref polynomial::_coefs vector. The set of variables used is stored in the static variable @ref polynomial::_var_names list.
 *
 *  Several methods are also proposed to build (@ref polynomial-construction), or operate (@ref polynomial-operators) theses polynomials,
 *
 *  @section polynomial-construction Construction
 *
 *  Three methods are proposed to build or read polynomials :
 *
 *    @li @ref polynomial::polynomial(int u)   gives a constant polynomial : degree is zero for  @e x, @e y @e z ...
 *    @li @ref polynomial::polynomial(const polynomial<C> & other)  gives a separate copy constructor,
 *    @li @ref polynomial::operator>>(istream& is, polynomial<U>& c) can also be used to read polynomials : variables will be added to the @ref polynomial::_var_names list if not already present.
 *
 *  @section polynomial-operators Operators
 *
 *  Three methods are proposed to compute new polynomials from already built ones :
 *
 *  @li @ref polynomial::operator+()
 *  @li @ref polynomial::operator-()
 *  @li @ref polynomial::operator*()
 *
 * all are depending on the operators that are defined for \<C\> (could be a @e min, @e plus semi-ring for example)
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

template<typename C> class polynomial {
 public:

  /// Build an empty polynomial
  polynomial() {
    for(unsigned i = 0; i < _coefs.size(); i++)
      _coefs[i].first.clear();
    _coefs.clear();
  };

  /// Build a constant polynomial
  polynomial(int u) {
    for(unsigned i = 0; i < _coefs.size(); i++)
      _coefs[i].first.clear();
    _coefs.clear();
    if (u != 0)
      _coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(0), C(u)));

    /* extends all the coef vectors to size "_var_names.size()" by adding extra zeros, sort them */
    for (typename std::vector<std::pair<std::vector<int>, C > >::iterator i_p = _coefs.begin(); i_p != _coefs.end(); i_p++)
      i_p->first.resize(_var_names.size(),0);
    sort(_coefs.begin(), _coefs.end());
  };

  /// Copy Constructor
  polynomial(const polynomial<C> &other) {
    for(unsigned i = 0; i < _coefs.size(); i++)
      _coefs[i].first.clear();
    _coefs.clear();

    /* copy */
    for (unsigned i = 0; i < other._coefs.size(); i++)
      _coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(other._coefs[i].first), C(other._coefs[i].second)));

    /* extends all the coef vectors to size "p._var_names.size()" by adding extra zeros, sort them */
    for (typename std::vector<std::pair<std::vector<int>, C > >::iterator i_p = _coefs.begin(); i_p != _coefs.end(); i_p++)
      i_p->first.resize(_var_names.size(),0);
    sort(_coefs.begin(), _coefs.end());
  };



  /// Erase a polynomial
  ~polynomial() {for(unsigned i = 0; i < _coefs.size(); i++)
      _coefs[i].first.clear();
    _coefs.clear();};

  // @{
  /// Operator @f$ + @f$ for two polynomials
  template<typename U> friend polynomial<U> operator+ (const polynomial<U> & l, const polynomial<U> & r);
  /// Operator @f$ - @f$ for two polynomials
  template<typename U> friend polynomial<U> operator- (const polynomial<U> & l, const polynomial<U> & r);
  /// Operator @f$ \times @f$ for two polynomials
  template<typename U> friend polynomial<U> operator* (const polynomial<U> & l, const polynomial<U> & r);
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


/// Operator @f$ == @f$ for two polynomials
template<typename C> inline bool equal_coefs (const std::pair<std::vector<int>, C> & coef1, const std::pair<std::vector<int>, C> & coef2) {
  return
    coef1.second == coef2.second
    &&
    coef1.first.size() == coef2.first.size()
    &&
    std::equal(coef1.first.begin(), coef1.first.end(), coef2.first.begin());
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
  // Sorting is supposed to be done before; since polynomials l and r are const here, must not be done here
  //sort(l._coefs.begin(),l._coefs.end());
  //sort(r._coefs.begin(),r._coefs.end());
  polynomial<C> result = polynomial<C>();
  typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_l = l._coefs.begin();
  typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_r = r._coefs.begin();
  while (i_l != l._coefs.end() && i_r != r._coefs.end()) {
    if (i_l->first.size() == i_r->first.size() && std::equal(i_l->first.begin(),i_l->first.end(),i_r->first.begin())) {
      C val = C(i_l->second) + C(i_r->second);
#ifndef USEINFINT
      if ( (i_l->second > 0) && (i_r->second > 0) && (val < 0) ) {
        _ERROR("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","polynomial<C> operator+ DOES OVERFLOW ...");
      }
#endif
      if (val != C(0))
        result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(i_l->first), val));
      i_l++;
      i_r++;
    } else {
      if (i_l->first < i_r->first) {
        C val = C(i_l->second);
        if (val != C(0))
          result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(i_l->first), val));
        i_l++;
      } else {
        C val = C(i_r->second);
        if (val != C(0))
          result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(i_r->first), val));
        i_r++;
      }
    }
  }
  while (i_l != l._coefs.end()) {
    C val = C(i_l->second);
    if (val != C(0))
      result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(i_l->first), val));
    i_l++;
  }
  while (i_r != r._coefs.end()) {
    C val = C(i_r->second);
    if (val != C(0))
      result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(i_r->first), val));
    i_r++;
  }
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("\t = [" << result << "]"););
  return result;
}



/// Operator - for two polynomials
template<typename C> inline polynomial<C> operator- (const polynomial<C> & l, const polynomial<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] - [" << r << "]"););
  // Sorting is supposed to be done before; since polynomials l and r are const here, must not be done here
  //sort(l._coefs.begin(),l._coefs.end());
  //sort(r._coefs.begin(),r._coefs.end());
  polynomial<C> result = polynomial<C>();
  typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_l = l._coefs.begin();
  typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_r = r._coefs.begin();
  while (i_l != l._coefs.end() && i_r != r._coefs.end()) {
    if (i_l->first.size() == i_r->first.size() && std::equal(i_l->first.begin(),i_l->first.end(),i_r->first.begin())) {
      C val = C(i_l->second) - C(i_r->second);
      if (val != C(0))
        result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(i_l->first), val));
      i_l++;
      i_r++;
    } else {
      if (i_l->first < i_r->first) {
        C val = C(i_l->second);
        if (val != C(0))
          result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(i_l->first), val));
        i_l++;
      } else {
        C val = C(0) - C(i_r->second);
        if (val != C(0))
          result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(i_r->first), val));
        i_r++;
      }
    }
  }
  while (i_l != l._coefs.end()) {
    C val = C(i_l->second);
    if (val != C(0))
      result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(i_l->first), val));
    i_l++;
  }
  while (i_r != r._coefs.end()) {
    C val = C(0) - C(i_r->second);
    if (val != C(0))
      result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(i_r->first), val));
    i_r++;
  }
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("\t = [" << result << "]"););
  return result;
}


/** Operator @f$ \times @f$ for two polynomials
 *  @todo{FIXME : this must be updated with Karatsuba, Toom-k or FFT/NNT methods (Strassen)}
 */
template<typename C> inline polynomial<C> operator* (const polynomial<C> & l, const polynomial<C> & r) {
  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("[" << l << "] * [" << r << "]"););
  // Sorting is supposed to be done before; since polynomials l and r are const here, must not be done here
  //sort(l._coefs.begin(),l._coefs.end());
  //sort(r._coefs.begin(),r._coefs.end());
  map<std::vector<int>, C> tmp_coefs;
  for (typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_l = l._coefs.begin(); i_l != l._coefs.end(); i_l++) {
    for (typename std::vector<std::pair<std::vector<int>, C > >::const_iterator i_r = r._coefs.begin(); i_r != r._coefs.end(); i_r++) {
      // build the degree_vector
      std::vector<int> degree_vector;
      if (i_l->first.size() > i_r->first.size()) {
        degree_vector = std::vector<int>(i_l->first);
        transform(i_r->first.begin(), i_r->first.end(), degree_vector.begin(), degree_vector.begin(), plus<int>());
      } else {
        degree_vector = std::vector<int>(i_r->first);
        transform(i_l->first.begin(), i_l->first.end(), degree_vector.begin(), degree_vector.begin(), plus<int>());
      }
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
          low->second = low->second + (i_l->second * i_r->second);
        } else {
#ifndef USEINFINT
          if ( (i_l->second > 0) && (i_r->second > 0) && ((i_l->second * i_r->second) < 0) ) {
            _ERROR("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","polynomial<C> operator* DOES OVERFLOW ...");
          }
#endif
          tmp_coefs.insert(low, std::pair<std::vector<int>,C >(std::vector<int>(degree_vector), (i_l->second * i_r->second)));
        }
      } else {
#ifndef USEINFINT
        if ( (i_l->second > 0) && (i_r->second > 0) && ((i_l->second * i_r->second) < 0) ) {
          _ERROR("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","polynomial<C> operator* DOES OVERFLOW ...");
        }
#endif
        tmp_coefs.insert(std::pair<std::vector<int>, C >(std::vector<int>(degree_vector), (i_l->second * i_r->second)));
      }
    }
  }

  // simplify zeros
  polynomial<C> result = polynomial<C>();
  for(typename map<std::vector<int>, C>::iterator it = tmp_coefs.begin(); it != tmp_coefs.end(); ++it)
    if (it->second != C(0))
      result._coefs.push_back(std::pair<std::vector<int>, C> (std::vector<int>(it->first),it->second));

  tmp_coefs.clear();

  VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("\t = [" << result << "]"););
  return result;
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
        os << p._var_names[j];
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
  p._coefs          = std::vector<std::pair<std::vector<int>, C> > ();
  string line;
  //cerr << "\"" << line << "\"" << endl;
  C sign = C(1);
  // must be formatted only on one single line
  if (getline(is, line)) {

    istringstream row(line);
    while (!row.eof()) {
      char c_times_power_plus_symbol = ' ';
      C coef;
      std::vector<int> var_degree;

      // read spaces until a symbol or a number occurs
      while (!row.eof()) {
        char c = row.peek();
        if (c == ' ' || c == '\t') {
          row.ignore();
          continue;
        } else {
          if (c >= '0' && c <= '9') {
            break;
          } else {
            if (c == '+' || c == '-') {
              if (c == '-') {
                sign = C(-1);
              } else {
                sign = C(1);
              }
              row.ignore();
              continue;
            } else {
              coef = C(1);
              goto next_variable;
            }
          }
        }
      }
      long long int cof;
      row >> cof;
      coef = C(cof);

      if (row.fail()) {
        _ERROR("operator>>"," invalid coefficient (found inside : \"" << line << "\")" << endl << "\t polynom format : <C>\"coef\" * variable1 [^power1]  [* variable2 [^power2] ... ] +  <C>\"coef\" * ..." << endl);
      }
      //cerr << "coef:" << coef << endl;

    next_symbol:
      row >> c_times_power_plus_symbol;
      //cerr << "symbol:" << c_times_power_plus_symbol << endl;
      if (!row.eof() && c_times_power_plus_symbol != '+' && c_times_power_plus_symbol != '-' && c_times_power_plus_symbol != '*') {
        _ERROR("operator>>"," missing * or + symbol (found : \""<< c_times_power_plus_symbol << "\")" << endl << "\t polynom format : <C>\"coef\" * variable1 [^ power1]  [* variable2 [^ power2] ... ] +  <C>\"coef\" * ..." << endl);
      }
      if (row.eof() || c_times_power_plus_symbol == '+' || c_times_power_plus_symbol == '-') {
        p._coefs.push_back(std::pair<std::vector<int>, C>(var_degree, coef * sign));
        if (c_times_power_plus_symbol == '-')
          sign = C(-1);
        else
          sign = C(1);
        continue;
      }

    next_variable:
      string var_symbol;
      row >> var_symbol;
      for(string::iterator it = var_symbol.begin(); it != var_symbol.end(); it++) {
        if ( ! ((*it >= 'a' && *it <= 'z') || (*it >= 'A' && *it <= 'Z') || *it == '_')) {
          _ERROR("operator>>"," variable name \""<< var_symbol<<"\" include non alpha symbols \'"<< (*it) <<"\' (found inside : \"" << line << "\")" << endl << "\t polynom format : <C>\"coef\" * variable1 [^power1]  [* variable2 [^power2] ... ] +  <C>\"coef\" * ..." << endl);
        }
      }
      //cerr << "variable:" << var_symbol << endl;
      int i_var = -1;

      // extend var_degree up to this variable if it is already defined
      bool defined_var = false;
      for (unsigned i = 0; i < p._var_names.size(); i++) {
        if (!p._var_names[i].compare(var_symbol)) {
          i_var = (int) i;
          defined_var = true;
          // already in the monomial list, so no need to extend the monomial list
          if (i_var < (int) var_degree.size()) {
            if (var_degree[i_var] != 0) { // check if possibly not at zero, because this is a classical error is to set it twice
              _ERROR("operator>>"," variable \""<< var_symbol<<"\" occuring twice in a monomial (found inside : \"" << line << "\")" << endl << "\t polynom format : <C>\"coef\" * variable1 [^power1]  [* variable2 [^power2] ... ] +  <C>\"coef\" * ..." << endl);
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
        unsigned i = p._var_names.size();
        i_var = (int) i;
        p._var_names.push_back(var_symbol);
        //cerr << "new symbol " << var_symbol <<  " at index "  << (p._var_names.size()) << endl;
        for (unsigned j = var_degree.size(); j < i; j++)
          var_degree.push_back(0);
        var_degree.push_back(1);
      }

      // possible next symbol if + or - (continue)
      row >> c_times_power_plus_symbol;
      //cerr << "symbol2:" << c_times_power_plus_symbol << endl;
      if (row.eof() || c_times_power_plus_symbol == '+' || c_times_power_plus_symbol == '-') {
        p._coefs.push_back(std::pair<std::vector<int>, C>(var_degree, coef * sign));
        if (c_times_power_plus_symbol == '-')
          sign = C(-1);
        else
          sign = C(1);
        continue;
      }
      if (c_times_power_plus_symbol == '*') {
        goto next_variable;
      }
      if (c_times_power_plus_symbol != '^') {
        _ERROR("operator>>"," missing * or + or ^ symbol  (found : \"" << c_times_power_plus_symbol << "\")" << endl << "\t polynom format : <C>\"coef\" * variable1 [^ power1]  [* variable2 [^ power2] ... ] +  <C>\"coef\" * ..." << endl);
      }
      int power_value = 0;
      row >> power_value;
      if (row.fail()) {
      _ERROR("operator>>"," invalid power coefficient (found inside : \"" << line << "\")" << endl << "\t polynom format : <C>\"coef\" * variable1 [^ power1]  [* variable2 [^ power2] ... ] +  <C>\"coef\" * ..." << endl);
      }
      //cerr << "power:" << power_value << endl;
      var_degree[i_var] = power_value;
      goto next_symbol;
    }
  }
  /* extends all the coef vectors to size "p._var_names.size()" by adding extra zeros, sort them */
  for (typename std::vector<std::pair<std::vector<int>, C > >::iterator i_p = p._coefs.begin(); i_p != p._coefs.end(); i_p++)
    i_p->first.resize(p._var_names.size(),0);
  sort(p._coefs.begin(), p._coefs.end());
  return is;
}

template<typename C> std::vector<string> polynomial<C>::_var_names = std::vector<string>(0);
// @}

#endif
