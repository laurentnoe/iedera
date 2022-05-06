#ifndef __INFINT_HH__
#define __INFINT_HH__

/** @defgroup infint infint class template
 *  @brief large integers defined over the templated C elements (int, long int, long long int, or other ...)
 *
 *
 *  @see matrix
 */

// @{

/*
 * infint<T> - Arbitrary-Precision Integer Arithmetic Library
 * Copyright (C) 2013 Sercan Tutar
 * [** with quick&dirty templating by L.Noe : 2015 **]
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 *
 * USAGE:
 *   It is pretty straight forward to use the library. Just create an instance of
 *   infint<T> class and start using it.
 *
 *   Useful methods:
 *      intSqrt:        integer square root operation
 *      digitAt:        returns digit at index
 *      numberOfDigits: returns number of digits
 *      size:           returns size in bytes
 *      toString:       converts it to a string
 *
 *   There are also conversion methods which allow conversion to primitive types:
 *   toInt, toLong, toLongLong, toUnsignedInt, toUnsignedLong, toUnsignedLongLong.
 *
 *   See ReadMe.txt for more info.
 *
 *
 * No overflows, happy programmers!
 *
 */

#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>

#include <limits.h>
#include <stdlib.h>

#ifndef LONG_LONG_MIN
#define LONG_LONG_MIN LLONG_MIN
#endif

#ifndef LONG_LONG_MAX
#define LONG_LONG_MAX LLONG_MAX
#endif

#ifndef ULONG_LONG_MIN
#define ULONG_LONG_MIN ULLONG_MIN
#endif

#ifndef ULONG_LONG_MAX
#define ULONG_LONG_MAX ULLONG_MAX
#endif



/**
 * infint Element storage type and its storage size + some other constants (defined only for \<int\> and \<long long int\> types)
 */
//@{

/** Elem template in general */
template<typename T> struct infint_elem;

/** Elem template for \<int\> parameter */
template<> struct infint_elem<int>           {
  /** type used for storage */
  typedef short int type;
  /** decimal base used for storage */
  const static short int base =            10000;
  /** maximal decimal value stored */
  const static short int upper_bound =      9999;
  /** number of decimal digits stored */
  const static short int digits =              4;
};

/** Elem template for \<long long int\> parameter */
template<> struct infint_elem<long long int> {
  /** type used for storage */
  typedef       int type;
  /** decimal base used for storage */
  const static       int base =        1000000000;
  /** maximal decimal value stored */
  const static       int upper_bound =  999999999;
  /** number of decimal digits stored */
  const static       int digits =               9;
};
//@}

/**
 * infint class to store large integers
 */
template<typename T> class infint {

  template<typename U> friend std::ostream& operator<<(std::ostream &s, const infint<U> &v);
  template<typename U> friend std::istream& operator>>(std::istream &s,       infint<U> &v);
public:

  /**
   * Some constants
   */
  //@{
  static const infint<T> zero;
  static const infint<T> one;
  static const infint<T> two;
  //@}

  /**
   * Constructors
   */
  //@{
  infint();
  infint(const char* c);
  infint(const std::string& s);
  infint(int l);
  infint(long l);
  infint(long long l);
  infint(unsigned int l);
  infint(unsigned long l);
  infint(unsigned long long l);
  //@}

  /**
   * Assignment operators
   */
  //@{
  /**
   * @param s is an input string value to assign
   * @return a reference on the assigned infint\<T\>
   */
  const infint<T>& operator=(const char* s);
  const infint<T>& operator=(const std::string& s);
  /**
   * @param l is an input integer value to assign
   * @return a reference on the assigned infint\<T\>
   */
  const infint<T>& operator=(int l);
  const infint<T>& operator=(long l);
  const infint<T>& operator=(long long l);
  const infint<T>& operator=(unsigned int l);
  const infint<T>& operator=(unsigned long l);
  const infint<T>& operator=(unsigned long long l);
  //@}

  /**
   * Unary increment/decrement operators
   */
  //@{
  /**
   * @return a reference on the current infint\<T\> "*this", once incremeted or decremented
   */
  const infint<T>& operator++();
  const infint<T>& operator--();
  /**
   * @param * is an input value to modify "*this"
   * @return a reference on the current infint\<T\> "*this", once incremeted or decremented
   */
  infint<T> operator++(int);
  infint<T> operator--(int);
  //@}

  /**
   * Operational assignments
   */
  //@{
  /**
   * @param rhs is an input value to modify "*this" (rhs is not modified)
   * @return a reference on the current infint\<T\> "*this", once modified by the operator
   */
  const infint<T>& operator+=(const infint<T>& rhs);
  const infint<T>& operator-=(const infint<T>& rhs);
  const infint<T>& operator*=(const infint<T>& rhs);
  const infint<T>& operator/=(const infint<T>& rhs);
  const infint<T>& operator%=(const infint<T>& rhs);
  const infint<T>& operator*=(typename infint_elem<T>::type rhs);
  //@}

  /**
   * Operations
   */
  //@{
  /**
   * @return a new infint\<T\>, negation of "*this"
   */
  infint<T> operator-() const;
  /**
   * @param rhs is a second right member for the operator with the left member "*this" (both remain unmodified)
   * @return a new assigned infint\<T\>
   */
  infint<T> operator+(const infint<T>& rhs) const;
  infint<T> operator-(const infint<T>& rhs) const;
  infint<T> operator*(const infint<T>& rhs) const;
  infint<T> operator/(const infint<T>& rhs) const;
  infint<T> operator%(const infint<T>& rhs) const;
  infint<T> operator*(typename infint_elem<T>::type rhs) const;
  //@}

  /**
   * Relational operations
   */
  //@{
  /**
   * @param rhs is a second right member for the boolean operator with the left member "*this" (both remain unmodified)
   * @return the boolean value of the relation being tested
   */
  bool operator==(const infint<T>& rhs) const;
  bool operator!=(const infint<T>& rhs) const;
  bool operator<(const infint<T>& rhs) const;
  bool operator<=(const infint<T>& rhs) const;
  bool operator>(const infint<T>& rhs) const;
  bool operator>=(const infint<T>& rhs) const;
  //@}


  /** Integer square root
   * @return the square root of the current infint\<T\>
   */
  infint<T> intSqrt() const;

  /**
   * Digit operations
   */
  //@{
  /**
   * @param i digit index
   * @return the char at this index
   */
  char digitAt(size_t i) const;
  size_t numberOfDigits() const;
  //@}

  /**
   * Size in bytes
   * @return the size in bytes
   */
  size_t size() const;

  /**
   * String conversion
   * @return the decimal representation of the current infint\<T\>
   */
  std::string toString() const;

  /**
   * Conversion to primitive types
   */
  //@{
  /**
   * @return the integer representation of the current infint\<T\>
   */
  int toInt() const;
  long toLong() const;
  long long toLongLong() const;
  unsigned int toUnsignedInt() const;
  unsigned long toUnsignedLong() const;
  unsigned long long toUnsignedLongLong() const;
  //@}

private:
  static typename infint_elem<T>::type dInR(const infint<T>& R, const infint<T>& D);
  static void multiplyByDigit(typename infint_elem<T>::type factor, std::vector<typename infint_elem<T>::type>& val);

  void correct(bool justCheckLeadingZeros = false, bool hasValidSign = false);
  void fromString(const std::string& s);
  void optimizeSqrtSearchBounds(infint<T>& lo, infint<T>& hi) const;
  void truncateToBase();
  bool equalizeSigns();
  void removeLeadingZeros();

  std::vector<typename infint_elem<T>::type> val; // number with base FACTOR
  bool pos; // true if number is positive

};


template<typename T> const infint<T> infint<T>::zero = 0;
template<typename T> const infint<T> infint<T>::one = 1;
template<typename T> const infint<T> infint<T>::two = 2;

template<typename T> inline infint<T>::infint() : pos(true)
{
  val.push_back((typename infint_elem<T>::type) 0);
}

template<typename T> inline infint<T>::infint(const char* s)
{
  fromString(s);
}

template<typename T> inline infint<T>::infint(const std::string& s)
{
  fromString(s);
}

template<typename T> inline infint<T>::infint(int l) : pos(l >= 0)
{
  if (!pos)
    {
      l = -l;
    }
  do
    {
      div_t dt = div(l, infint_elem<T>::base);
      val.push_back((typename infint_elem<T>::type) dt.rem);
      l = dt.quot;
    } while (l > 0);
}

template<typename T> inline infint<T>::infint(long l) : pos(l >= 0)
{
  if (!pos)
    {
      l = -l;
    }
  do
    {
      ldiv_t dt = ldiv(l, infint_elem<T>::base);
      val.push_back((typename infint_elem<T>::type) dt.rem);
      l = dt.quot;
    } while (l > 0);
}

template<typename T> inline infint<T>::infint(long long l) : pos(l >= 0)
{
  if (!pos)
    {
      l = -l;
    }
  do
    {
#ifndef _WIN32
      lldiv_t dt = lldiv(l, infint_elem<T>::base);
      val.push_back((typename infint_elem<T>::type) dt.rem);
      l = dt.quot;
#else
      val.push_back((typename infint_elem<T>::type) (l % infint_elem<T>::base));
      l = l / infint_elem<T>::base;
#endif
    } while (l > 0);
}

template<typename T> inline infint<T>::infint(unsigned int l) : pos(true)
{
  do
    {
      val.push_back((typename infint_elem<T>::type) (l % infint_elem<T>::base));
      l = l / infint_elem<T>::base;
    } while (l > 0);
}

template<typename T> inline infint<T>::infint(unsigned long l) : pos(true)
{
  do
    {
      val.push_back((typename infint_elem<T>::type) (l % infint_elem<T>::base));
      l = l / infint_elem<T>::base;
    } while (l > 0);
}

template<typename T> inline infint<T>::infint(unsigned long long l) : pos(true)
{
  do
    {
      val.push_back((typename infint_elem<T>::type) (l % infint_elem<T>::base));
      l = l / infint_elem<T>::base;
    } while (l > 0);
}

template<typename T> inline const infint<T>& infint<T>::operator=(const char* c)
{
  fromString(c);
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator=(const std::string& s)
{
  fromString(s);
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator=(int l)
{
  pos = l >= 0;
  val.clear();
  if (!pos)
    {
      l = -l;
    }
  do
    {
      div_t dt = div(l, infint_elem<T>::base);
      val.push_back((typename infint_elem<T>::type) dt.rem);
      l = dt.quot;
    } while (l > 0);
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator=(long l)
{
  pos = l >= 0;
  val.clear();
  if (!pos)
    {
      l = -l;
    }
  do
    {
      ldiv_t dt = ldiv(l, infint_elem<T>::base);
      val.push_back((typename infint_elem<T>::type) dt.rem);
      l = dt.quot;
    } while (l > 0);
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator=(long long l)
{
  pos = l >= 0;
  val.clear();
  if (!pos)
    {
      l = -l;
    }
  do
    {
#ifndef _WIN32
      lldiv_t dt = lldiv(l, infint_elem<T>::base);
      val.push_back((typename infint_elem<T>::type) dt.rem);
      l = dt.quot;
#else
      val.push_back((typename infint_elem<T>::type) (l % infint_elem<T>::base));
      l = l / infint_elem<T>::base;
#endif
    } while (l > 0);
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator=(unsigned int l)
{
  pos = true;
  val.clear();
  do
    {
      val.push_back((typename infint_elem<T>::type) (l % infint_elem<T>::base));
      l = l / infint_elem<T>::base;
    } while (l > 0);
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator=(unsigned long l)
{
  pos = true;
  val.clear();
  do
    {
      val.push_back((typename infint_elem<T>::type) (l % infint_elem<T>::base));
      l = l / infint_elem<T>::base;
    } while (l > 0);
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator=(unsigned long long l)
{
  pos = true;
  val.clear();
  do
    {
      val.push_back((typename infint_elem<T>::type) (l % infint_elem<T>::base));
      l = l / infint_elem<T>::base;
    } while (l > 0);
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator++()
{
  val[0] += (pos ? 1 : -1);
  this->correct(false, true);
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator--()
{
  val[0] -= (pos ? 1 : -1);
  this->correct(false, true);
  return *this;
}

template<typename T> inline infint<T> infint<T>::operator++(int)
{
  infint<T> result = *this;
  val[0] += (pos ? 1 : -1);
  this->correct(false, true);
  return result;
}

template<typename T> inline infint<T> infint<T>::operator--(int)
{
  infint<T> result = *this;
  val[0] -= (pos ? 1 : -1);
  this->correct(false, true);
  return result;
}

template<typename T> inline const infint<T>& infint<T>::operator+=(const infint<T>& rhs)
{
  if (rhs.val.size() > val.size())
    {
      val.resize(rhs.val.size(), 0);
    }
  for (size_t i = 0; i < val.size(); ++i)
    {
      val[i] = (pos ? val[i] : -val[i]) + (i < rhs.val.size() ? (rhs.pos ? rhs.val[i] : -rhs.val[i]) : 0);
    }
  correct();
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator-=(const infint<T>& rhs)
{
  if (rhs.val.size() > val.size())
    {
      val.resize(rhs.val.size(), 0);
    }
  for (size_t i = 0; i < val.size(); ++i)
    {
      val[i] = (pos ? val[i] : -val[i]) - (i < rhs.val.size() ? (rhs.pos ? rhs.val[i] : -rhs.val[i]) : 0);
    }
  correct();
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator*=(const infint<T>& rhs)
{
  // TODO: optimize (do not use operator*)
  *this = *this * rhs;
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator/=(const infint<T>& rhs)
{
  if (rhs == zero)
    {
      std::cerr << "Division by zero!" << std::endl;
      return *this;
    }
  infint<T> R, D = (rhs.pos ? rhs : -rhs), N = (pos ? *this : -*this);
  bool oldpos = pos;
  val.clear();
  val.resize(N.val.size(), 0);
  for (int i = (int) N.val.size() - 1; i >= 0; --i)
    {
      R.val.insert(R.val.begin(), (typename infint_elem<T>::type) 0);
      R.val[0] = N.val[i];
      R.correct(true);
      typename infint_elem<T>::type cnt = dInR(R, D);
      R -= D * cnt;
      val[i] += cnt;
    }
  correct();
  pos = (val.size() == 1 && val[0] == 0) ? true : (oldpos == rhs.pos);
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator%=(const infint<T>& rhs)
{
  if (rhs == zero)
    {
      std::cerr << "Division by zero!" << std::endl;
      return zero;
    }
  infint<T> D = (rhs.pos ? rhs : -rhs), N = (pos ? *this : -*this);
  bool oldpos = pos;
  val.clear();
  for (int i = (int) N.val.size() - 1; i >= 0; --i)
    {
      val.insert(val.begin(), (typename infint_elem<T>::type) 0);
      val[0] = N.val[i];
      correct(true);
      *this -= D * dInR(*this, D);
    }
  correct();
  pos = (val.size() == 1 && val[0] == 0) ? true : oldpos;
  return *this;
}

template<typename T> inline const infint<T>& infint<T>::operator*=(typename infint_elem<T>::type rhs)
{
  typename infint_elem<T>::type factor = rhs < 0 ? -rhs : rhs;
  bool oldpos = pos;
  multiplyByDigit(factor, val);
  correct();
  pos = (val.size() == 1 && val[0] == 0) ? true : (oldpos == (rhs >= 0));
  return *this;
}

template<typename T> inline infint<T> infint<T>::operator-() const
{//PROFILED_SCOPE
  infint<T> result = *this;
  result.pos = !pos;
  return result;
}

template<typename T> inline infint<T> infint<T>::operator+(const infint<T>& rhs) const
{//PROFILED_SCOPE
  infint<T> result;
  result.val.resize(val.size() > rhs.val.size() ? val.size() : rhs.val.size(), 0);
  for (size_t i = 0; i < val.size() || i < rhs.val.size(); ++i)
    {
      result.val[i] = (i < val.size() ? (pos ? val[i] : -val[i]) : 0) + (i < rhs.val.size() ? (rhs.pos ? rhs.val[i] : -rhs.val[i]) : 0);
    }
  result.correct();
  return result;
}

template<typename T> inline infint<T> infint<T>::operator-(const infint<T>& rhs) const
{//PROFILED_SCOPE
  infint<T> result;
  result.val.resize(val.size() > rhs.val.size() ? val.size() : rhs.val.size(), 0);
  for (size_t i = 0; i < val.size() || i < rhs.val.size(); ++i)
    {
      result.val[i] = (i < val.size() ? (pos ? val[i] : -val[i]) : 0) - (i < rhs.val.size() ? (rhs.pos ? rhs.val[i] : -rhs.val[i]) : 0);
    }
  result.correct();
  return result;
}

template<typename T> inline infint<T> infint<T>::operator*(const infint<T>& rhs) const
{//PROFILED_SCOPE
  infint<T> result;
  result.val.resize(val.size() + rhs.val.size(), 0);
  T carry = 0;
  size_t digit = 0;
  for (;; ++digit)
    {//PROFILED_SCOPE
      //result.val[digit] = (typename infint_elem<T>::type) (carry % infint_elem<T>::base);
      //carry /= infint_elem<T>::base;

      T oldcarry = carry;
      carry /= infint_elem<T>::base;
      result.val[digit] = (typename infint_elem<T>::type) (oldcarry - carry * infint_elem<T>::base);

      bool found = false;
      for (size_t i = digit < rhs.val.size() ? 0 : digit - rhs.val.size() + 1; i < val.size() && i <= digit; ++i)
        {//PROFILED_SCOPE
          T pval = result.val[digit] + val[i] * (T) rhs.val[digit - i];
          if (pval >= infint_elem<T>::base || pval <= -infint_elem<T>::base)
            {//PROFILED_SCOPE
              //carry += pval / infint_elem<T>::base;
              //pval %= infint_elem<T>::base;

              T quot = pval / infint_elem<T>::base;
              carry += quot;
              pval -= quot * infint_elem<T>::base;
            }
          result.val[digit] = (typename infint_elem<T>::type) pval;
          found = true;
        }
      if (!found)
        {//PROFILED_SCOPE
          break;
        }
    }
  for (; carry > 0; ++digit)
    {//PROFILED_SCOPE
      result.val[digit] = (typename infint_elem<T>::type) (carry % infint_elem<T>::base);
      carry /= infint_elem<T>::base;
    }
  result.correct();
  result.pos = (result.val.size() == 1 && result.val[0] == 0) ? true : (pos == rhs.pos);
  return result;
}

template<typename T> inline infint<T> infint<T>::operator/(const infint<T>& rhs) const
{//PROFILED_SCOPE
  if (rhs == zero)
    {
      std::cerr << "Division by zero!" << std::endl;
      return zero;
    }
  infint<T> Q, R, D = (rhs.pos ? rhs : -rhs), N = (pos ? *this : -*this);
  Q.val.resize(N.val.size(), 0);
  for (int i = (int) N.val.size() - 1; i >= 0; --i)
    {//PROFILED_SCOPE
      R.val.insert(R.val.begin(), (typename infint_elem<T>::type) 0);
      R.val[0] = N.val[i];
      R.correct(true);
      typename infint_elem<T>::type cnt = dInR(R, D);
      R -= D * cnt;
      Q.val[i] += cnt;
    }
  Q.correct();
  Q.pos = (Q.val.size() == 1 && Q.val[0] == 0) ? true : (pos == rhs.pos);
  return Q;
}

template<typename T> inline infint<T> infint<T>::operator%(const infint<T>& rhs) const
{//PROFILED_SCOPE
  if (rhs == zero)
    {
      std::cerr << "Division by zero!" << std::endl;
      return zero;
    }
  infint<T> R, D = (rhs.pos ? rhs : -rhs), N = (pos ? *this : -*this);
  for (int i = (int) N.val.size() - 1; i >= 0; --i)
    {
      R.val.insert(R.val.begin(), (typename infint_elem<T>::type) 0);
      R.val[0] = N.val[i];
      R.correct(true);
      R -= D * dInR(R, D);
    }
  R.correct();
  R.pos = (R.val.size() == 1 && R.val[0] == 0) ? true : pos;
  return R;
}

template<typename T> inline infint<T> infint<T>::operator*(typename infint_elem<T>::type rhs) const
{//PROFILED_SCOPE
  infint<T> result = *this;
  typename infint_elem<T>::type factor = rhs < 0 ? -rhs : rhs;
  multiplyByDigit(factor, result.val);
  result.correct();
  result.pos = (result.val.size() == 1 && result.val[0] == 0) ? true : (pos == (rhs >= 0));
  return result;
}

template<typename T> inline bool infint<T>::operator==(const infint<T>& rhs) const
{//PROFILED_SCOPE
  if (pos != rhs.pos || val.size() != rhs.val.size())
    {
      return false;
    }
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      if (val[i] != rhs.val[i])
        {
          return false;
        }
    }
  return true;
}

template<typename T> inline bool infint<T>::operator!=(const infint<T>& rhs) const
{//PROFILED_SCOPE
  if (pos != rhs.pos || val.size() != rhs.val.size())
    {
      return true;
    }
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      if (val[i] != rhs.val[i])
        {
          return true;
        }
    }
  return false;
}

template<typename T> inline bool infint<T>::operator<(const infint<T>& rhs) const
{//PROFILED_SCOPE
  if (pos && !rhs.pos)
    {
      return false;
    }
  if (!pos && rhs.pos)
    {
      return true;
    }
  if (val.size() > rhs.val.size())
    {
      return pos ? false : true;
    }
  if (val.size() < rhs.val.size())
    {
      return pos ? true : false;
    }
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      if (val[i] < rhs.val[i])
        {
          return pos ? true : false;
        }
      if (val[i] > rhs.val[i])
        {
          return pos ? false : true;
        }
    }
  return false;
}

template<typename T> inline bool infint<T>::operator<=(const infint<T>& rhs) const
{//PROFILED_SCOPE
  if (pos && !rhs.pos)
    {
      return false;
    }
  if (!pos && rhs.pos)
    {
      return true;
    }
  if (val.size() > rhs.val.size())
    {
      return pos ? false : true;
    }
  if (val.size() < rhs.val.size())
    {
      return pos ? true : false;
    }
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      if (val[i] < rhs.val[i])
        {
          return pos ? true : false;
        }
      if (val[i] > rhs.val[i])
        {
          return pos ? false : true;
        }
    }
  return true;
}

template<typename T> inline bool infint<T>::operator>(const infint<T>& rhs) const
{//PROFILED_SCOPE
  if (pos && !rhs.pos)
    {
      return true;
    }
  if (!pos && rhs.pos)
    {
      return false;
    }
  if (val.size() > rhs.val.size())
    {
      return pos ? true : false;
    }
  if (val.size() < rhs.val.size())
    {
      return pos ? false : true;
    }
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      if (val[i] < rhs.val[i])
        {
          return pos ? false : true;
        }
      if (val[i] > rhs.val[i])
        {
          return pos ? true : false;
        }
    }
  return false;
}

template<typename T> inline bool infint<T>::operator>=(const infint<T>& rhs) const
{//PROFILED_SCOPE
  if (pos && !rhs.pos)
    {
      return true;
    }
  if (!pos && rhs.pos)
    {
      return false;
    }
  if (val.size() > rhs.val.size())
    {
      return pos ? true : false;
    }
  if (val.size() < rhs.val.size())
    {
      return pos ? false : true;
    }
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      if (val[i] < rhs.val[i])
        {
          return pos ? false : true;
        }
      if (val[i] > rhs.val[i])
        {
          return pos ? true : false;
        }
    }
  return true;
}

template<typename T> inline void infint<T>::optimizeSqrtSearchBounds(infint<T>& lo, infint<T>& hi) const
{//PROFILED_SCOPE
  infint<T> hdn = one;
  for (int i = (int) this->numberOfDigits() / 2; i >= 2; --i)
    {
      hdn *= 10;
    }
  if (lo < hdn)
    {
      lo = hdn;
    }
  hdn *= 100;
  if (hi > hdn)
    {
      hi = hdn;
    }
}

template<typename T> inline infint<T> infint<T>::intSqrt() const
{//PROFILED_SCOPE
  if (*this <= zero)
    {
      std::cerr << "intSqrt called for non-positive integer: " << *this << std::endl;
      return zero;
    }
  infint<T> hi = *this / two + one, lo = zero, mid, mid2;
  optimizeSqrtSearchBounds(lo, hi);
  do
    {
      mid = (hi + lo) / two; // 8 factor
      mid2 = mid * mid; // 1 factor
      if (mid2 == *this)
        {
          lo = mid;
          break;
        }
      else if (mid2 < *this)
        {
          lo = mid;
        }
      else
        {
          hi = mid;
        }
    } while (lo < hi - one && mid2 != *this);
  return lo;
}

template<typename T> inline char infint<T>::digitAt(size_t i) const
{//PROFILED_SCOPE
  /*
    #ifdef INFINT_USE_SHORT infint_elem<T>::base // uses 10^4 (short) as the base
    static const int powersOfTen[] = { 1, 10, 100, 1000};
    #else
    static const int powersOfTen[] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};
    #endif

    if (numberOfDigits() <= i)
    {
    std::cerr << "Invalid digit index: " << i << std::endl;
    return -1;
    }
    return (val[i / infint_elem<T>::digits] / powersOfTen[i % infint_elem<T>::digits]) % 10;
  */
  typename infint_elem<T>::type powersOfTen = 1;
  size_t j = i % infint_elem<T>::digits;
  while (j > 0) {
    j--;
    powersOfTen *= 10;
  }
  return (val[i / infint_elem<T>::digits] / powersOfTen) % 10;
}

template<typename T> inline size_t infint<T>::numberOfDigits() const
{//PROFILED_SCOPE

  /*
    return (val.size() - 1) * infint_elem<T>::digits
    #ifdef INFINT_USE_SHORT
    + (val.back() > 999 ? 4 : (val.back() > 99 ? 3 : (val.back() > 9 ? 2 : 1)));
    #else
    + (val.back() > 99999999 ? 9 : (val.back() > 9999999 ? 8 : (val.back() > 999999 ? 7 : (val.back() > 99999 ? 6 :
    (val.back() > 9999 ? 5 : (val.back() > 999 ? 4 : (val.back() > 99 ? 3 : (val.back() > 9 ? 2 : 1))))))));
    #endif
  */
  size_t u = (val.size() - 1) * infint_elem<T>::digits;
  typename infint_elem<T>::type  b = infint_elem<T>::base / 10;
  size_t d = infint_elem<T>::digits;
  while (val.back() < b) {
    b = b / 10;
    d = d - 1;
  }
  return u + d;
}

template<typename T> inline std::string infint<T>::toString() const
{//PROFILED_SCOPE
  std::ostringstream oss;
  oss << *this;
  return oss.str();
}

template<typename T> inline size_t infint<T>::size() const
{//PROFILED_SCOPE
  return val.size() * sizeof(typename infint_elem<T>::type) + sizeof(bool);
}

template<typename T> inline int infint<T>::toInt() const
{//PROFILED_SCOPE
  if (*this > INT_MAX || *this < INT_MIN)
    std::cerr << "Out of INT bounds: " << *this << std::endl;
  int result = 0;
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      result = result * infint_elem<T>::base + val[i];
    }
  return pos ? result : -result;
}

template<typename T> inline long infint<T>::toLong() const
{//PROFILED_SCOPE
  if (*this > LONG_MAX || *this < LONG_MIN)
    std::cerr << "Out of LONG bounds: " << *this << std::endl;
  long result = 0;
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      result = result * infint_elem<T>::base + val[i];
    }
  return pos ? result : -result;
}

template<typename T> inline long long infint<T>::toLongLong() const
{//PROFILED_SCOPE
  if (*this > LONG_LONG_MAX || *this < LONG_LONG_MIN)
#ifdef INFINT_USE_EXCEPTIONS
    throw infintException("out of bounds");
#else
  std::cerr << "Out of LLONG bounds: " << *this << std::endl;
#endif
  long long result = 0;
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      result = result * infint_elem<T>::base + val[i];
    }
  return pos ? result : -result;
}

template<typename T> inline unsigned int infint<T>::toUnsignedInt() const
{//PROFILED_SCOPE
  if (!pos || *this > UINT_MAX)
#ifdef INFINT_USE_EXCEPTIONS
    throw infintException("out of bounds");
#else
  std::cerr << "Out of UINT bounds: " << *this << std::endl;
#endif
  unsigned int result = 0;
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      result = result * infint_elem<T>::base + val[i];
    }
  return result;
}

template<typename T> inline unsigned long infint<T>::toUnsignedLong() const
{//PROFILED_SCOPE
  if (!pos || *this > ULONG_MAX)
#ifdef INFINT_USE_EXCEPTIONS
    throw infintException("out of bounds");
#else
  std::cerr << "Out of ULONG bounds: " << *this << std::endl;
#endif
  unsigned long result = 0;
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      result = result * infint_elem<T>::base + val[i];
    }
  return result;
}

template<typename T> inline unsigned long long infint<T>::toUnsignedLongLong() const
{//PROFILED_SCOPE
  if (!pos || *this > ULONG_LONG_MAX)
#ifdef INFINT_USE_EXCEPTIONS
    throw infintException("out of bounds");
#else
  std::cerr << "Out of ULLONG bounds: " << *this << std::endl;
#endif
  unsigned long long result = 0;
  for (int i = (int) val.size() - 1; i >= 0; --i)
    {
      result = result * infint_elem<T>::base + val[i];
    }
  return result;
}

template<typename T> inline void infint<T>::truncateToBase()
{//PROFILED_SCOPE
  for (size_t i = 0; i < val.size(); ++i) // truncate each
    {
      if (val[i] >= infint_elem<T>::base || val[i] <= -infint_elem<T>::base)
        {//PROFILED_SCOPE
          div_t dt = div(val[i], infint_elem<T>::base);
          val[i] = dt.rem;
          if (i + 1 >= val.size())
            {//PROFILED_SCOPE
              val.push_back(dt.quot);
            }
          else
            {//PROFILED_SCOPE
              val[i + 1] += dt.quot;
            }
        }
    }
}

template<typename T> inline bool infint<T>::equalizeSigns()
{//PROFILED_SCOPE
  bool isPositive = true;
  int i = (int) ((val.size())) - 1;
  for (; i >= 0; --i)
    {
      if (val[i] != 0)
        {
          isPositive = val[i--] > 0;
          break;
        }
    }

  if (isPositive)
    {
      for (; i >= 0; --i)
        {
          if (val[i] < 0)
            {
              int k = 0, index = i + 1;
              for (; (size_t)(index) < val.size() && val[index] == 0; ++k) ++index; // count adjacent zeros on left
              //if ((size_t)(index) < val.size() && val[index] > 0)
              { // number on the left is positive
                val[index] -= 1;
                val[i] += infint_elem<T>::base;
                for (; k > 0; --k)
                  {
                    val[i + k] = infint_elem<T>::upper_bound;
                  }
              }
            }
        }
    }
  else
    {
      for (; i >= 0; --i)
        {
          if (val[i] > 0)
            {
              int k = 0, index = i + 1;
              for (; (size_t)(index) < val.size() && val[index] == 0; ++k) ++index; // count adjacent zeros on right
              //if ((size_t)(index) < val.size() && val[index] < 0)
              { // number on the left is negative
                val[index] += 1;
                val[i] -= infint_elem<T>::base;
                for (; k > 0; --k)
                  {
                    val[i + k] = -infint_elem<T>::upper_bound;
                  }
              }
            }
        }
    }

  return isPositive;
}

template<typename T> inline void infint<T>::removeLeadingZeros()
{//PROFILED_SCOPE
  for (int i = (int) (val.size()) - 1; i > 0; --i) // remove leading 0's
    {
      if (val[i] != 0)
        {
          return;
        }
      else
        {
          val.erase(val.begin() + i);
        }
    }
}

template<typename T> inline void infint<T>::correct(bool justCheckLeadingZeros, bool hasValidSign)
{//PROFILED_SCOPE
  if (!justCheckLeadingZeros)
    {
      truncateToBase();

      if (equalizeSigns())
        {
          pos = ((val.size() == 1 && val[0] == 0) || !hasValidSign) ? true : pos;
        }
      else
        {
          pos = hasValidSign ? !pos : false;
          for (size_t i = 0; i < val.size(); ++i)
            {
              val[i] = abs(val[i]);
            }
        }
    }

  removeLeadingZeros();
}

template<typename T> inline void infint<T>::fromString(const std::string& s)
{//PROFILED_SCOPE
  pos = true;
  val.clear();
  // TODO use resize
  val.reserve(s.size() / infint_elem<T>::digits + 1);
  int i = (int) s.size() - infint_elem<T>::digits;
  for (; i >= 0; i -= infint_elem<T>::digits)
    {
      val.push_back(atoi(s.substr(i, infint_elem<T>::digits).c_str()));
    }
  if (i > -infint_elem<T>::digits)
    {
      std::string ss = s.substr(0, i + infint_elem<T>::digits);
      if (ss.size() == 1 && ss[0] == '-')
        {
          pos = false;
        }
      else
        {
          val.push_back(atoi(ss.c_str()));
        }
    }
  if (val.back() < 0)
    {
      val.back() = -val.back();
      pos = false;
    }
  correct(true);
}

template<typename T> inline typename infint_elem<T>::type infint<T>::dInR(const infint<T>& R, const infint<T>& D)
{//PROFILED_SCOPE
  typename infint_elem<T>::type min = 0, max = infint_elem<T>::upper_bound;
  while (max - min > 0)
    {
      typename infint_elem<T>::type avg = max + min;
      //div_t dt = div(avg, 2);
      //avg = dt.rem ? (dt.quot + 1) : dt.quot;
      typename infint_elem<T>::type havg = avg / 2;
      avg = (avg - havg * 2) ? (havg + 1) : havg;
      infint<T> prod = D * avg;
      if (R == prod)
        {//PROFILED_SCOPE
          return avg;
        }
      else if (R > prod)
        {//PROFILED_SCOPE
          min = avg;
        }
      else
        {//PROFILED_SCOPE
          max = avg - 1;
        }
    }
  return min;
}

template<typename T> inline void infint<T>::multiplyByDigit(typename infint_elem<T>::type factor, std::vector<typename infint_elem<T>::type>& val)
{//PROFILED_SCOPE
  typename infint_elem<T>::type carry = 0;
  for (size_t i = 0; i < val.size(); ++i)
    {
      T pval = val[i] * (T) factor + carry;
      if (pval >= infint_elem<T>::base || pval <= -infint_elem<T>::base)
        {
          //carry = (typename infint_elem<T>::type) (pval / infint_elem<T>::base);
          //pval %= infint_elem<T>::base;

          carry = (typename infint_elem<T>::type) (pval / infint_elem<T>::base);
          pval -= carry * infint_elem<T>::base;
        }
      else
        {
          carry = 0;
        }
      val[i] = (typename infint_elem<T>::type) pval;
    }
  if (carry > 0)
    {
      val.push_back(carry);
    }
}

/**************************************************************/
/******************** NON-MEMBER OPERATORS ********************/
/**************************************************************/


//@{
/**
 * @brief Input/output
 * @param s is an istream or ostream
 * @param v is the infint\<U\> to input or output
 * @return a reference on the current istream or ostream
 */
template<typename T> inline std::istream& operator>>(std::istream &s, infint<T> &v)
{//PROFILED_SCOPE
  std::string str;
  char c = ' ';
  bool n = s.get(c);
  if (n && c == '-') {
    str.push_back(c);
    n = s.get(c);
  }

  while (n && (c >= '0') && (c <= '9')){
    str.push_back(c);
    n = s.get(c);
  }
  s.unget();
  v.fromString(str);
  return s;
}
template<typename T> inline std::ostream& operator<<(std::ostream &s, const infint<T> &v)
{//PROFILED_SCOPE
  if (!v.pos)
    {
      s << '-';
    }
  bool first = true;
  for (int i = (int) v.val.size() - 1; i >= 0; --i)
    {
      if (first)
        {
          s << v.val[i];
          first = false;
        }
      else
        {
          s << std::setfill('0') << std::setw(infint_elem<T>::digits) << v.val[i];
        }
    }
  return s;
}
//@}
//@}

#endif
