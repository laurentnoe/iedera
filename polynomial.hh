#ifndef __POLYNOMIAL_HH__
#define __POLYNOMIAL_HH__

/** @defgroup polynomial polynomial class template
 *  @brief polynomials defined over the templated C coefficients (int, long int, long long int, or other ...) with several variables
 *
 *
 *  @see matrix
 */

// @{

//STL
#include <assert.h>
#include <algorithm>
#include <vector>
#include <functional>
#include <string>
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
using namespace std;


/**
 * @class polynomial
 * @brief polynomials are defined over the templated C coefficients (int, long int, long long int, or other ...) ; this template can handle multivariate polynomials if needed
 *        but is only for positive coefficients (no minus symbol allowed)
 */

template<typename C> class polynomial {
 public:

  /// Set the list of variable to work on for all the polynomials
  static void setvars(vector<string> & var_names) {_var_names = var_names;};


  /// Build an empty polynomial
  polynomial() {_coefs.clear();};
  /// Build a constant polynomial
  polynomial(int u) {_coefs.clear(); _coefs.push_back(pair<vector<int>, C> (vector<int>(_var_names.size(),0), C(u)));};

  /// Erase a polynomial
  ~polynomial() {};

  // @{
  /// Operator + (min) for two polynomials
  template<typename U> friend polynomial<U> operator+ (polynomial<U> & l, polynomial<U> & r);
  /// Operator @f$ \times @f$ (add) for two polynomials
  template<typename U> friend polynomial<U> operator* (polynomial<U> & l, polynomial<U> & r);
  /// Operator @f$ != @f$ for two polynomials
  template<typename U> friend bool    operator!= (polynomial<U> & l, polynomial<U> & r);
  /// Operator @f$ == @f$ for two polynomials
  template<typename U> friend bool    operator== (polynomial<U> & l, polynomial<U> & r);
  /// Operator @f$ < @f$ for two polynomials
  template<typename U> friend bool    operator< (polynomial<U> & l, polynomial<U> & r);
  /// Operator @f$ > @f$ for two polynomials
  template<typename U> friend bool    operator> (polynomial<U> & l, polynomial<U> & r);
  // @}

  // @{
  /// Print a polynomial
  template<typename U> friend ostream& operator<< (ostream& os, const polynomial<U>& c);
  /// Load a polynomial
  template<typename U> friend istream& operator>> (istream& is, polynomial<U>& c);
  // @}



protected:
  /// give a list of the coefficients C for a polynomial, each vector associated with C represents the list of variables (each variable being addressed by its proper index in the vector) and its degree (integer given by the vector at the index given by the variable)
  vector<pair<vector<int>, C> > _coefs;
  /// give the index of each variable name, or the name of the variable from its index
  static vector<string> _var_names;
};

/// Operator + for two polynomials
template<typename C> inline polynomial<C> operator+ (polynomial<C> & l, polynomial<C> & r) {
  sort(l._coefs.begin(),l._coefs.end());
  sort(r._coefs.begin(),r._coefs.end());
  assert(l._var_names == r._var_names);
  polynomial<C> result = polynomial<C>();
  typename vector<pair<vector<int>, C > >::const_iterator i_l = l._coefs.begin();
  typename vector<pair<vector<int>, C > >::const_iterator i_r = r._coefs.begin();
  while ( i_l != l._coefs.end() && i_r != r._coefs.end()) {
    if (i_l->first == i_r->first) {
      result._coefs.push_back(pair<vector<int>, C> (vector<int>(i_l->first), C(i_l->second) + C(i_r->second)));
      i_l++;
      i_r++;
    } else {
      if (i_l->first < i_r->first) {
        result._coefs.push_back(pair<vector<int>, C> (vector<int>(i_l->first), C(i_l->second)));
        i_l++;
      } else {
        result._coefs.push_back(pair<vector<int>, C> (vector<int>(i_r->first), C(i_r->second)));
        i_r++;
      }
    }
  }
  while (i_l != l._coefs.end()) {
    result._coefs.push_back(pair<vector<int>, C> (vector<int>(i_l->first), C(i_l->second)));
    i_l++;
  }
  while (i_r != r._coefs.end()) {
    result._coefs.push_back(pair<vector<int>, C> (vector<int>(i_r->first), C(i_r->second)));
    i_r++;
  }
  return result;
}

/// Operator @f$ \times @f$ for two polynomials
template<typename C> inline polynomial<C> operator* (polynomial<C> & l, polynomial<C> & r) {
  sort(l._coefs.begin(),l._coefs.end());
  sort(r._coefs.begin(),r._coefs.end());
  //assert(l._var_names == r._var_names);
  polynomial<C> result = polynomial<C>();
  for (typename vector<pair<vector<int>, C > >::const_iterator i_l = l._coefs.begin(); i_l != l._coefs.end(); i_l++) {
    for (typename vector<pair<vector<int>, C > >::const_iterator i_r = r._coefs.begin();  i_r != r._coefs.end(); i_r++) {
      // build the degree_vector
      assert(i_l->first.size() == i_r->first.size());
      vector<int> degree_vector;
      degree_vector.reserve(i_l->first.size());
      transform(i_l->first.begin(), i_l->first.end(), i_r->first.begin(), back_inserter(degree_vector), plus<int>());
      // check already existance in result (log search)
      typename vector<pair<vector<int>, C > >::iterator low = lower_bound(result._coefs.begin(), result._coefs.end(), pair<vector<int>, C > (degree_vector,C(0)));
      if (low != result._coefs.end() && low->first == degree_vector) {
        low->second = low->second + (i_l->second * i_r->second);
      } else {
        result._coefs.insert(low, pair<vector<int>, C > (degree_vector,C(i_l->second * i_r->second)));
      }
    }
  }
  return result;
}

/// Print a polynomial
template<typename C> ostream& operator<< (ostream& os, const polynomial<C> & p) {
  for (int i = 0 ; i < p._coefs.size() ; i++) {
    if (i > 0)
      os << " + ";
    bool first_prod = true;
    // dont print '1 *'  "something ..."
    if (p._coefs[i].second != C(1)) {
      os << p._coefs[i].second;
      first_prod = false;
    }
    for (int j = 0 ; j <  p._var_names.size() ; j++) {
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
  os << endl;
  return os;
}

/// Load a polynomial ( \<C\>"coef" * variable1 [^power1]  [* variable2 [^power2] ... ] + ...)
//  @todo{FIXME : this kind of Ugly and Non secure code will not be used from command lined iedera !! unless one uses already generated polynomials from this code ...}
template<typename C> istream& operator>> (istream& is, polynomial<C>& p) {
  p._coefs          = vector<pair<vector<int>, C> > ();
  string line;
  // must be formatted only on one single line
  if (getline(is, line)) {

    istringstream row(line);
    while (!row.eof()) {
      char c_times_power_plus_symbol = ' ';
      C coef;
      vector<int> var_degree(p._var_names.size(),0);

      // read spaces until a symbol or a number occurs
      while (!row.eof()) {
      char c = row.peek();
      if (c == ' ' || c == '\t') {
        row.ignore();
        continue;
      } else {
        if ( c >= '0' && c <= '9' ) {
          break;
        } else {
          coef = C(1);
          goto next_variable;
        }
      }
      }

      row >> coef;
      if ( row.fail() ) {
        _ERROR("operator>>"," invalid coefficient (found inside : \"" << line << "\")" << endl << "\t polynom format : <C>\"coef\" * variable1 [^power1]  [* variable2 [^power2] ... ] +  <C>\"coef\" * ..." << endl);
      }
      //cerr << "coef:" << coef << endl;

    next_symbol:
      row >> c_times_power_plus_symbol;
      //cerr << "symbol:" << c_times_power_plus_symbol << endl;
      if (!row.eof() && c_times_power_plus_symbol != '+' && c_times_power_plus_symbol != '*') {
        _ERROR("operator>>"," missing * or + symbol (found : \""<< c_times_power_plus_symbol << "\")" << endl << "\t polynom format : <C>\"coef\" * variable1 [^ power1]  [* variable2 [^ power2] ... ] +  <C>\"coef\" * ..." << endl);
      }
      if (row.eof() || c_times_power_plus_symbol == '+') {
        p._coefs.push_back( pair<vector<int>, C>(var_degree, coef));
        continue;
      }

    next_variable:
      string var_symbol;
      row >> var_symbol;
      //cerr << "variable:" << var_symbol << endl;
      int i_var = -1;
      for (int i=0; i<p._var_names.size(); i++) {
        if (!p._var_names[i].compare(var_symbol)) {
          i_var = i;
          var_degree[i_var] = 1;
        }
      }
      if (i_var == -1) {
        _ERROR("operator>>"," missing correct variable name (found : \"" << var_symbol << "\")" << endl << "\t polynom format : <C>\"coef\" * variable1 [^ power1]  [* variable2 [^ power2] ... ] +  <C>\"coef\" * ..." << endl);
      }

      row >> c_times_power_plus_symbol;
      //cerr << "symbol2:" << c_times_power_plus_symbol << endl;
      if (row.eof() || c_times_power_plus_symbol == '+') {
        p._coefs.push_back( pair<vector<int>, C>(var_degree, coef));
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
      if ( row.fail() ) {
      _ERROR("operator>>"," invalid power coefficient (found inside : \"" << line << "\")" << endl << "\t polynom format : <C>\"coef\" * variable1 [^ power1]  [* variable2 [^ power2] ... ] +  <C>\"coef\" * ..." << endl);
      }
      //cerr << "power:" << power_value << endl;
      var_degree[i_var] = power_value;
      goto next_symbol;
    }
  }

  sort(p._coefs.begin(), p._coefs.end());
  return is;
}


template<typename C> vector<string> polynomial<C>::_var_names = vector<string>(0);
// @}

#endif


