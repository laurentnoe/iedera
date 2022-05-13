#ifndef __MATRIX_HH__
#define __MATRIX_HH__

/** @page matrix Matrix
 *  @brief Matrix description and functions
 *  @tableofcontents
 *
 *  @section matrix-description Description
 *
 *  This part describes a @ref matrix<T> : each @ref matrix<T> is mainly represented by a set of row\<T\>, each @ref row<T> can be represented as a @b sparse or @b non-sparse set of \<T\> cells.
 *  @li each  @ref row<T> has an additional integer attribute @ref row<T>::_final, to match its equivalent @ref transition<T> attributes (see @ref automaton<T>).
 *  @li each  @ref row<T> can be stored in a @e sparse or @e non-sparse form ( @ref row<T>::_sparse) with a way to revert the storage selection (@ref row<T>::setsparse(const bool sparse)).
 *  @li each  @ref matrix<T> may bring probabilities (T = double, T = polynomial\<long long int\>), costs (T = cost\<int\>), counts (T = unsigned long long).
 *
 *  @note @ref matrix<T> are (just) a more compact way to store @ref automata<T> attributes, once letters are not needed anymore...
 *
 *   Default @ref matrix<T> constructor is almost empty, but several methods from @ref automaton<T> are proposed to produce matrices. They must be used first!
 *
 *  Several methods are also proposed to manipulate theses matrices (@ref local-matrix-manipulation or @ref global-matrices-manipulation), compute properties (@ref matrix-computed-properties).
 *
 *  @section local-matrix-manipulation Local matrix manipulation
 *
 *  Two methods are proposed to manipulate matrices localy :
 *  @li @ref matrix::addNewRow(const int final, const bool sparse) to append a new row at the end of the @ref matrix<T>
 *  @li @ref matrix::addNewCell(const int i, const int j, const T v) to a add a cell on row @e i , @b provided @b that the coordinate @e i for the row\<T\> @b is @b correct.
 *
 *  @section global-matrices-manipulation Global matrices manipulation
 *
 *  Two methods are proposed to manipulate matrices globaly :
 *
 *  @li @ref matrix::Transpose() to reverse lines/columns,
 *  @li @ref matrix::Compose() for the product of two, compatible in size, matrices.
 *
 *  @section matrix-computed-properties Matrices computed properties
 *
 *  @li @ref matrix::Pr() is the most classical computation to reach @e final or @e non-final states after @e nbSteps transitions
 *  @li @ref matrix::Pr_transitive_final() is more complex, it computes "the transitive sum" of the @e final states values (@e final could be 1, but also more) that are crossed during the walk
 *  @li @ref matrix::Pr_one_step_from_one() is doing one single compuation step, but use the @ final values from a second matrix @m_final that is passed as a parameter
 *
 *  @ref matrices_slicer is a class provided when several differents matrices have to be multiplied with a sliding windows moving along them ...
 *
 *      Spaced seed design on profile HMMs for precise HTS read-mapping
 *        efficient sliding window product on the matrix semi-group
 *
 * You can also find three "stepwise equivalent" methods in the @ref automaton<T> class :
 * @see automaton::matrices_step_pr_product,  @see automaton::matrices_step_cost_product, and @see automaton::matrices_step_count_product
 * these three methods give the "breadth first" product as an ordered set of matrices  @f$M_1,M_2,M_3\ldots,M_l@f$, thus enabling any computation @f$M_i,M_{i+1}\ldots,M_{j}@f$ @f$\forall 0 \leq i < j \leq l@f$.
 *
 *  @todo FIXME : to be continued
 *
 *
 *
 */

/** @defgroup matrix matrix and row class templates
 *  @brief sparse / dense matrix template, with "automaton like" attributes ("final" integer attribute per row), and with semi-ring templates for each cell
 */

// @{

/// A try with classical TYPE_TRAITS
#if !defined(USE_TYPE_TRAITS) && !defined(USE_TR1_TYPE_TRAITS)
#define USE_TYPE_TRAITS
#endif

/// A try with classical STD
#if !defined(HAS_STD_TYPE_TRAITS) && !defined(HAS_STD_TR1_TYPE_TRAITS)
#define HAS_STD_TYPE_TRAITS
#endif


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
//STR1
#ifdef USE_TYPE_TRAITS
  #include <type_traits>
#else
  #ifdef USE_TR1_TYPE_TRAITS
    #include <tr1/type_traits>
  #else
    #error  <type_traits> or <tr1/type_traits> not selected in "matrix.hh" : try g++ ... -DUSE_TR1_TYPE_TRAITS or -DUSE_TYPE_TRAITS, and use g++ version >= 4.3 with "-std=c++0x" or "-std=gnu++0x"
  #endif
#endif


#ifdef HAS_STD_TYPE_TRAITS
   using namespace std;
#else
   #ifdef HAS_STD_TR1_TYPE_TRAITS
      using namespace std::tr1;
   #else
      #error  std or std::tr1 not selected in "matrix.hh" : try g++ ... -DHAS_STD_TYPE_TRAITS or -DHAS_STD_TR1_TYPE_TRAITS, and use g++ version >= 4.3 with "-std=c++0x" or "-std=gnu++0x"
   #endif
#endif


//STR
#include "macro.h"
#include "infint.hh"
#include "cost.hh"
#include "polynomial.hh"

/**@name  Enable/Disable if
 * @brief only for definition and selection of good functions
 */
// @{
/** @brief struct template is enabled for bool = true
 */
template <bool B, typename T = void> struct enable_if_ca { /** type enabled */ typedef T type; };
/** @brief struct template is disabled for bool = false
 */
template <typename T> struct enable_if_ca <false, T> {};

/** @brief struct template is enabled for bool = false
 */
template <bool B, typename T = void>  struct disable_if_ca { /** type enabled */  typedef T type; };
/** @brief struct template is disabled for bool = true
 */
template <typename T> struct disable_if_ca <true, T> {};
// @}



/** Zero and One are defined as "copies" for "arithmetic" types (but not for "Cost" types that are <min,plus> defined ... ),
 *  Transition gives the correct weight for a transition for "arithmetic types" (but for "Cost" type, gives the cost of the letter and not the letter)
 */
// @{
#ifdef HAS_STD_TYPE_TRAITS

/// arithmetic templates (for double probabilities or for integer count)
/** @brief Zero constant in the (+,x) semi-ring
 *  @return 0
 */
/** @addtogroup polynomial
 *  @addtogroup infint
 */
template<typename T> typename enable_if_ca < std::is_arithmetic<T>::value || std::is_same<T, infint<long long> >::value || std::is_same<T, polynomial<long long> >::value || std::is_same<T, polynomial<infint<long long> > >::value, T >::type                                    Zero()                   {return T(0);}
/** @brief One constant in the (+,x) semi-ring
 *  @return 1
 */
template<typename T> typename enable_if_ca < std::is_arithmetic<T>::value || std::is_same<T, infint<long long> >::value || std::is_same<T, polynomial<long long> >::value || std::is_same<T, polynomial<infint<long long> > >::value, T >::type                                    One()                    {return T(1);}
/// cost template (cost is not defined "arithmetic")
/** @addtogroup cost
 */
// @{
/** @brief Zero constant in the (min,+) semi-ring
 *  @return +infinity
 */
template<typename T> typename disable_if_ca < std::is_arithmetic<T>::value || std::is_same<T, infint<long long> >::value || std::is_same<T, polynomial<long long> >::value || std::is_same<T, polynomial<infint<long long> > >::value, T>::type                                     Zero()                   {return T(0x7fffffff);}
/** @brief One constant in the (min,+) semi-ring
 *  @return 0
 */
template<typename T> typename disable_if_ca < std::is_arithmetic<T>::value || std::is_same<T, infint<long long> >::value || std::is_same<T, polynomial<long long> >::value || std::is_same<T, polynomial<infint<long long> > >::value, T>::type                                     One()                    {return T(0x00000000);}
// @}

/** @brief IsProb() check if \<T\> can be interpreted as a probability
 *  @return true if \<T\> is a probability, false otherwise
 */
/// floating point and polynomials
template<typename T> typename enable_if_ca< (std::is_arithmetic<T>::value && std::is_floating_point<T>::value) || std::is_same<T, polynomial<long long int> >::value || std::is_same<T, polynomial<infint<long long int> > >::value, bool>::type                                    IsProb()                 {return true;}
/// integer count and costs
template<typename T> typename disable_if_ca< (std::is_arithmetic<T>::value && std::is_floating_point<T>::value) || std::is_same<T, polynomial<long long int> >::value || std::is_same<T, polynomial<infint<long long int> > >::value, bool>::type                                   IsProb()                 {return false;}

/** @brief IsVoid() check if \<T\> is a void type
 *  @return true if \<T\> is a void, false otherwise
 */
/// floating point and polynomial
template<typename T> typename  enable_if_ca< std::is_void<T>::value && true, bool>::type                                         IsVoid()                 {return true;}
/// integer count and costs
template<typename T> typename disable_if_ca< std::is_void<T>::value && true, bool>::type                                         IsVoid()                 {return false;}

#else

/// arithmetic templates (for double)
/** @brief Zero constant in the (+,x) semi-ring
 *  @return 0
 */
/** @addtogroup polynomial
 *  @addtogroup infint
 */
template<typename T> typename enable_if_ca < std::tr1::is_arithmetic<T>::value || std::tr1::is_same<T, infint<long long> >::value || std::tr1::is_same<T, polynomial<long long> >::value || std::tr1::is_same<T, polynomial<infint<long long> > >::value, T >::type                                         Zero()                    {return T(0);}
/** @brief One constant in the (+,x) semi-ring
 *  @return 1
 */
template<typename T> typename enable_if_ca < std::tr1::is_arithmetic<T>::value || std::tr1::is_same<T, infint<long long> >::value || std::tr1::is_same<T, polynomial<long long> >::value || std::tr1::is_same<T, polynomial<infint<long long> > >::value, T >::type                                         One()                     {return T(1);}
/// cost template (cost is not defined "arithmetic")
/** @addtogroup cost
 */
// @{
/** @brief Zero constant in the (min,+) semi-ring
 *  @return +infinity
 */
template<typename T> typename disable_if_ca < std::tr1::is_arithmetic<T>::value || std::tr1::is_same<T, infint<long long> >::value || std::tr1::is_same<T, polynomial<long long> >::value || std::tr1::is_same<T, polynomial<infint<long long> > >::value, T>::type                                          Zero()                   {return T(0x7fffffff);}
/** @brief One constant in the (min,+) semi-ring
 *  @return 0
 */
template<typename T> typename disable_if_ca < std::tr1::is_arithmetic<T>::value || std::tr1::is_same<T, infint<long long> >::value || std::tr1::is_same<T, polynomial<long long> >::value || std::tr1::is_same<T, polynomial<infint<long long> > >::value, T>::type                                          One()                    {return T(0x00000000);}
// @}

/** @brief IsProb() check if \<T\> can be interpreted as a probability
 *  @return true if \<T\> is a probability, false otherwise
 */
/// floating point and polynomials
template<typename T> typename  enable_if_ca< (std::tr1::is_arithmetic<T>::value && std::tr1::is_floating_point<T>::value) || std::tr1::is_same<T, polynomial<long long int> >::value || std::tr1::is_same<T, polynomial<infint<long long int> > >::value, bool>::type                                         IsProb()                 {return true;}
/// integer count and costs
template<typename T> typename disable_if_ca< (std::tr1::is_arithmetic<T>::value && std::tr1::is_floating_point<T>::value) || std::tr1::is_same<T, polynomial<long long int> >::value || std::tr1::is_same<T, polynomial<infint<long long int> > >::value, bool>::type                                        IsProb()                 {return false;}

/** @brief IsVoid() check if \<T\> is a void type
 *  @return true if \<T\> is a void, false otherwise
 */
/// floating point and polynomials
template<typename T> typename  enable_if_ca< std::tr1::is_void<T>::value && true, bool>::type                                         IsVoid()                 {return true;}
/// integer count and costs
template<typename T> typename disable_if_ca< std::tr1::is_void<T>::value && true, bool>::type                                         IsVoid()                 {return false;}

#endif
// @}



/**
 * @class row
 * @brief sparse / dense vector, represented by an ordered vector of cells, each cell containing either a couple @f$ ( integer \times T ) @f$ or a single element @f$ T @f$
 * T may represents here
 * - a probability
 * - a cost
 *  each row has a also a flag to together with a "final" state for this row.
 */

template<typename T> class row {

public:

  /** @brief Constructor for an empty row (empty cell lists)
   *  @param final is the row state (final or not)
   *  @param sparse is the implementation choosen for the row being created
   */
  inline row(const int final = 0, const bool sparse = true): _final(final), _sparse(sparse), _cells_sparse(), _cells_dense() {};

  /** @brief Copy constructor
   *  @param r is the row to be cloned (vector reallocation)
   */
  inline row(const row<T> & r): _final(r._final), _sparse(r._sparse), _cells_sparse(r._cells_sparse), _cells_dense(r._cells_dense) {};

  /** @brief Erase a row (clear cell lists first)
   */
  inline ~row() { _cells_sparse.clear(); _cells_dense.clear(); _final = 0; _sparse = true;};

  /** @brief Clear a row (clear the vector _cells and set _final to 0)
   */
  inline void clear() { _cells_sparse.clear(); _cells_dense.clear(); _final = 0; _sparse = true;};


  /** @brief add a new cell (j) if does not exists, otherwise "add" (according to T) the "v" value to the current cell (j)
   *  @param j is the column number
   *  @param v is the element to be inserted (or added if it already exists)
   */
  void insert(const int j, const T v);

  /** @brief final state
   *  @return final state of the current row
   */
  inline int final() const { return _final; }

  /** @brief set final state
   *  @param final is the value that will be set to the current row
   */
  inline void setfinal(const int final) { _final = final; }

  /** @brief get the implementation of the current row
   *  @return true if the row is sparse (dense otherwise)
   */
  inline bool sparse() const { return _sparse; }

  /** @brief change the implementation of the current row (sparse to dense and dense to sparse)
   *  @param sparse is the implementation choosen for the row being modified
   */
  inline void setsparse(const bool sparse) {
    if (!_sparse && sparse) {
      for (unsigned i = 0; i < _cells_dense.size(); i++) {
        if (_cells_dense[i] != Zero<T>()) {
          std::pair<int,T> p = make_pair(i,_cells_dense[i]);
          _cells_sparse.push_back(p);
        }
      }
      _cells_dense.clear();
      _sparse = true;
    } else if (_sparse && !sparse) {
      for (typename std::vector< std::pair<int,T> >::const_iterator i = _cells_sparse.begin(); i != _cells_sparse.end(); i++) {
        int j = i->first;
        T   v = i->second;
        if (_cells_dense.size() <= (unsigned)j)
          _cells_dense.resize(j+1, Zero<T>());
        _cells_dense[j] = _cells_dense[j] + v;
      }
      _cells_sparse.clear();
      _sparse = false;
    }
  }

  /** @brief return the size as the number of cells
   *  @return the number of cells
   */
  inline int size() const { if (_sparse) return _cells_sparse.size(); else return _cells_dense.size(); }

  /** @brief return the size as the number of cells
   *  @return the number of cells
   */
  inline int max_index() const {
    if (_sparse) {
      int max_column = 0;
      for (typename std::vector< std::pair<int,T> >::const_iterator i = _cells_sparse.begin(); i != _cells_sparse.end(); i++) {
        max_column = MAX(max_column,i->first);
      }
      return max_column;
    } else {
      return _cells_dense.size() - 1;
    }
  }

  /** @brief return the density as the number of cells over the max index (virtual size)
   *  @return the density
   */
  double density() const {
    if (_sparse) {
      return (double) _cells_sparse.size() / (max_index()+1);
    } else {
      int n = 0;
      for (unsigned i = 0; i < _cells_dense.size(); i++)
        if (_cells_dense[i] != Zero<T>())  n++;
      return (double) n / _cells_dense.size(); //FIXME if (_cell_dense.size() == 0)
    }
  }


  /** @brief return the reference to the vector of sparse cells
   *  @return a reference to the vector
   */
  const std::vector< std::pair<int,T> > & cells_sparse() const { return _cells_sparse; }

  /** @brief return the reference to the vector of dense cells
   *  @return a reference to the vector
   */
  const std::vector< T > & cells_dense() const { return _cells_dense; }

  // @{
  /// print row information
  template<typename U> friend ostream& operator<< (ostream& os, const row<U>& m);
  /** @brief print a row in "maple" format
   *  @param os is the outputstream
   *  @param max_size is used to pad with Zeros<T>() on a sparse row
   *  @param ignore_indices_sorted  is a sorted list of elements indices that must not be in the output (because not needed for example)
   *  @param v_end_cell is a cumulative cell that keep the values of elements not present in the output between two outputed elements
   */
  void maple(ostream& os, int max_size, std::vector<int> & ignore_indices_sorted, T & v_end_cell) const;

  /** @brief print a row in "maple" recursive format
   *  @param os is the outputstream
   */
  void maple_recursive(ostream& os) const;

  /// load row information
  template<typename U> friend istream& operator>> (istream& is, row<U>& m);
  // @}

 protected :
  /// final state
  int _final;
  /// is the row sparse or dense represented ?
  bool _sparse;
  /// cells list in sparse mode
  std::vector< std::pair<int,T> > _cells_sparse;
  /// cells vector in dense mode
  std::vector< T > _cells_dense;

  /// matrix is a friend class to ease access
  template<typename U> friend class matrix;
  /// matrices_slicer is a friend class to ease access
  template<typename U> friend class matrices_slicer;

};




/**
 * @class matrix
 *
 * @brief sparse matrix, represented by a (dense) vector of (sparse / dense) rows.
 *        This class has strong link with automaton as it usualy represents such objects (final/non-final states) after
 *        being processed with a HMM model.
 * @see row
 * @see automaton
 */

template<typename T> class matrix {

public:


  /** @brief Constructor for an empty matrix
   */

  matrix() : _rows() {};

  /** @brief Copy constructor
   *  @param m is the matrix to be cloned (vector reallocation)
   */

  matrix(const matrix<T> & m) {
    _rows = std::vector<row<T> >(m._rows.size());
    for (unsigned i = 0; i < m._rows.size(); i++)
      _rows[i] = row<T>(m._rows[i]);
  };


  /** @brief Erase a matrix (clear _cells_sparse and _cells_dense for each row first)
   */
  ~matrix() {
    for (unsigned i = 0; i < _rows.size(); i++)
      _rows[i].clear();
    _rows.clear();
  };

  /** @brief Clear a matrix (clear _cells_sparse and _cells_dense for each row first)
   */
  void clear() {
    for (unsigned i = 0; i < _rows.size(); i++)
      _rows[i].clear();
    _rows.clear();
  };

  /** @brief Resize the number of line of the matrix (clear the rows)
   *  @param s is the new number of rows.
   */
  void resize(int s) {
    _rows.resize(s);
  };

  /** @brief add a new cell (i,j) if does not exists, otherwise "add" (according to T) the "v" value to the current cell (i,j)
   *  @param i is the row number
   *  @param j is the column number
   *  @param v is the element to be inserted (or added if it already exists)
   */
  inline void insert(const int i, const int j, const T v) { _rows[i].insert(j,v);}

  /** @brief add a new row at the end of the matrix
   *  @param final is the row state (final or not)
   *  @param sparse is the implementation choosen for the row being added
   *  @return the "new" index of the row being added
   */
  inline int addNewRow(const int final = 0, const bool sparse = true) { _rows.push_back(row<T>(final,sparse)); return _rows.size() - 1;}

  /** @brief add a new cell (i,j) if does not exists, otherwise "add" (according to T) the "v" value to the current cell (i,j)
   *  @param i is the row number
   *  @param j is the column number
   *  @param v is the element to be inserted (or added if it already exists)
   */
  inline void addNewCell(const int i, const int j, const T v) { _rows[i].insert(j,v); }

  /** @brief return the size as the number of rows
   *  @return the number of rows
   */
  inline int size() const { return _rows.size(); }

  /** @brief return the fullsize as the full number of cells
   *  @return the numrber of cells
   */
  inline int fullsize() const {
    int fsize = 0;
    for (unsigned i = 0; i < _rows.size(); i++)
      fsize += _rows[i].size();
    return fsize;
  }

  /** @brief selfLoop on row i
   *  @param i the row and also the column number (since this is a selfloop)
   *  @param p is the value to be inserted (by default it is a One\<T\>(), i.e. the neutral element of the *\<T\>)
   */
  inline void selfLoop(const int i, const T p = One<T>()) { insert(i,i,p); }


  /** @brief Compute the transposed matrix
   *  @return the transposed matrix
   *  @todo{FIXME : final states/rows are not set}
   */
  matrix<T> * Transpose() const;

  /** @brief Compute the composition of two matrices
   *  @param other is the matrix that will be composed with the current matrix this
   *  @return the composed matrix of the current matrix with the other matrix
   */
  matrix<T> * Compose(const matrix<T> &other) const;

  /** @brief Compute the probability/(min cost/count) to be at a final/(non final) row during the "nbSteps"th step
   *  @param nbSteps is the number of iterated products done on the matrix
   *  @param final is set when the computation is done on final states, otherwise, it is done on non final states (default, true, is on final states)
   *  @return the probability/(min cost) to be at a final/(non final) row during the "nbSteps"-th step
   *  @warning only work with "non window" matrices (must be "square" and "self-injecting") : otherwise final states are not set correctly
   *  @see Pr_transitive_final,Pr_one_step_from_one
   */
  const T Pr(const int nbSteps = gv_alignment_length, const bool final = true) const;

  /** @brief Compute the probability/(min cost/count) to be at a "transitive-sum of final values" during the "nbSteps"th step
   *  @param nbSteps is the number of iterated products done on the matrix
   *  @param max_final_value is the value that must not be reached by a "transitive-sum of finals states" (set to maximal \<int\> value by default)
   *  @param sub_final_value is the replacement value when the previous max_final_value is reached (set to 1 by default)
   *  @return the probability/(min cost) to be at a "transitive-sum of final values" during the "nbSteps"th step
   *  @warning only work "non window" matrices (must be "square" and "self-injecting") : otherwise final states are not set correctly
   *  @warning previous automata products must use the UNION_ADD operator only to keep "final values" greater than one
   *  @see Pr
   */
  std::vector<T> * Pr_transitive_final(const int nbSteps, const unsigned max_final_value = INT_INFINITY, const int sub_final_value = 1) const;

  /** @brief Compute a one step single walk from the initial to the final states marked by the matrix "m_final"
   *  @param m_final is the matrix used at the end to mark final states
   *    (matrices are otherwise represented as one final/row + a set of links, but nothing is known on the reaching states, so "m_final" is usefull here ...)
   *  @param final is set when the computation is done on final states, otherwise, it is done on non final states (default, true, is on final states)
   *  @return the probability/min cost/count to be at a final/(non final) row during the "nbSteps"-th step
   *  @warning only work with "full-window" computed matrices (not "self-injecting")
   *  @see Pr
   */
  const T Pr_one_step_from_one(const matrix<T> &m_final, const bool final = true) const;


  /**@name  IO methods
   * @brief IO streams to print/store/load matrix
   */
  // @{
  /// print matrix information
  template<typename U> friend ostream& operator<< (ostream& os, const matrix<U>& m);
  /** @brief print a matrix in "maple" format
   *  @param os is the outputstream
   *  @param separate_final indicates if final states ùmust be included in the matrix
   */

  void maple(ostream& os, bool separate_final = false) const;
  /** @brief print a row in "maple" recursive format
   *  @param os is the outputstream
   */
  void maple_recursive(ostream& os) const;

  /// load matrix information
  template<typename U> friend istream& operator>> (istream& is, matrix<U>& m);
  // @}

protected :
  /// vector of rows
  std::vector< row<T> > _rows;

  /// matrices_slicer is a friend class to ease access
  template<typename U> friend class matrices_slicer;
};


/// compare pairs of elements in a row only on their first index number (and not the second T element)
template<typename T> inline bool pairless(const std::pair<int, T> l, const std::pair<int, T> r) {
  return l.first < r.first;
}


/// output method for a row
template<typename T> inline ostream& operator<<(ostream& os, const row<T>& r) {
  os << "\t" << (r.size()) << "\t" << (r.final()) << endl;
  // for each row, display each cell
  if (r.sparse()) {
    for (typename std::vector< std::pair<int,T> >::const_iterator i = r._cells_sparse.begin(); i != r._cells_sparse.end(); i++) {
      os << "\t\t\t" << (i->first) << "\t" << (i->second) << endl;
    }
  } else {
    for (typename std::vector< T >::const_iterator i = r._cells_dense.begin(); i != r._cells_dense.end(); i++) {
      if ((*i) != Zero<T>()) {
        os << "\t\t\t" << (i - r._cells_dense.begin()) << "\t" << (*i) << endl;
      }
    }
  }
  return os;
}

/// output method as a "maple" row
template<typename T> inline void row<T>::maple(ostream& os, int max_size, std::vector<int> & ignore_indices_sorted, T & v_end_cell) const {
  os << "[";
  unsigned ignore_i = 0;

  // for each row, display each cell
  if (sparse()) {
    /* sparse implementation */
    int u_old = 0;
    bool first_element = false;
    for (typename std::vector< std::pair<int,T> >::const_iterator i = _cells_sparse.begin(); i != _cells_sparse.end(); i++) {
      while (u_old < (i->first)) {
        if (ignore_i < ignore_indices_sorted.size() && ignore_indices_sorted[ignore_i] == u_old) {
          ignore_i++;
        } else {
          if (first_element)
            os << ",";
          os << Zero<T>();
          first_element = true;
        }
        u_old++;
      }
      if (ignore_i < ignore_indices_sorted.size() && ignore_indices_sorted[ignore_i] == u_old) {
        v_end_cell = v_end_cell + i->second;
        ignore_i++;
      } else {
        if (first_element)
          os << ",";
        os << (i->second);
        first_element = true;
      }
      u_old++;
    }

    while (u_old < (max_size)) {
      if (ignore_i < ignore_indices_sorted.size() && ignore_indices_sorted[ignore_i] == u_old) {
        ignore_i++;
      } else {
        if (first_element)
          os << ",";
        os << Zero<T>();
        first_element = true;
      }
      u_old++;
    }

  } else {
    /* dense implementation */
    bool first_element = false;
    for (typename std::vector< T >::const_iterator i = _cells_dense.begin(); i != _cells_dense.end(); i++) {
      if (ignore_i < ignore_indices_sorted.size() && ignore_indices_sorted[ignore_i] == (i-_cells_dense.begin())) {
        v_end_cell = v_end_cell + (*i);
        ignore_i++;
      } else {
        if (first_element)
          os << ",";
        os << (*i);
        first_element = true;
      }
      os << (*i);
      first_element = true;
    }
  }
  os << "]  # final : " << final() << endl;
}

/// output method as a "maple" recursive definition
template<typename T> inline void row<T>::maple_recursive(ostream& os) const {
  // for each row, display each cell
  if (sparse()) {
    /* sparse implementation */
    bool first_element = false;
    for (typename std::vector< std::pair<int,T> >::const_iterator i = _cells_sparse.begin(); i != _cells_sparse.end(); i++) {
      if ((i->second) != Zero<T>()) {
        if (first_element)
          os << " + ";
        os << "q" << (i->first) << "(n-1) * " << (i->second);
        first_element = true;
      }
    }

  } else {
    /* dense implementation */
    bool first_element = false;
    for (typename std::vector< T >::const_iterator i = _cells_dense.begin(); i != _cells_dense.end(); i++) {
      if ((*i) != Zero<T>()) {
        if (first_element)
          os << " + ";
        os << "q" << (i - _cells_dense.begin()) << "(n-1) * " << (*i);
        first_element = true;
      }
    }
  }
}



/// input method for a row
template<typename T> inline istream& operator>>(istream& is, row<T>& r) {
  // previous data removed if any
  r._cells_dense.clear();
  r._cells_sparse.clear();
  r._sparse = true;

  // a) read row length
  int nbcells = 0;
  is >> nbcells;
#ifdef DEBUGREADING
  cerr << nbcells << endl;
#endif

  if (nbcells <= 0) {
    cerr << "> when reading nbcells" << endl;
    _ERROR("row operator>>","incorrect size "<< nbcells);
  }

  // b) read row final
  int final;
  is >> final;
  r._final = final;

  // c) read row data
  for (int j = 0; j<nbcells; j++){
    int  to = 0;
    T    pr = T();
    is >> to >> pr;
#ifdef DEBUGREADING
    cerr << "\t\t" << to << "\t" << pr << endl;
#endif
    if (to < 0) {
      cerr << "> when reading cell" << endl;
      _ERROR("row operator>>","incorrect column value "<< to);
    }
    r.insert(to,pr);
  }
  return is;
}


/// output method for the current matrix
template<typename T> inline ostream& operator<<(ostream& os, const matrix<T>& m) {
  os << m._rows.size() << endl;
  // display each row
  for (typename std::vector< row<T> >::const_iterator i = m._rows.begin(); i != m._rows.end(); i++) {
    const row<T> & r = *i;
    os << "\t" << r;
  }
  return os;
}

/// output method as a maple "matrix" (must be square here)
template<typename T> inline void matrix<T>::maple(ostream& os, bool separate_final) const {
  os << "[" << endl;
  std::vector<int> ignore_indices_sorted(0);
  if (separate_final) {
    for (typename std::vector< row<T> >::const_iterator i = _rows.begin(); i != _rows.end(); i++) {
      const row<T> & r = *i;
      if (r.final())
        ignore_indices_sorted.push_back(i-_rows.begin());
    }
  }
  std::vector<T>  v_end(_rows.size()-ignore_indices_sorted.size(),Zero<T>());

  // display each row
  bool first_element = false;
  unsigned ignore_i = 0;
  for (typename std::vector< row<T> >::const_iterator i = _rows.begin(); i != _rows.end(); i++) {
    const row<T> & r = *i;
    if (separate_final && ignore_i < ignore_indices_sorted.size() && (ignore_indices_sorted[ignore_i] == i - _rows.begin())) {
      ignore_i++;
    } else {
      if (first_element)
        os << "," << endl;
      r.maple(os,_rows.size(),ignore_indices_sorted,v_end[(i - _rows.begin()) - ignore_i]);
      first_element = true;
    }
  }
  os << "];" << endl;

  // display the v_end vector
  if (separate_final) {
    os << "[";
    for (unsigned u = 0; u < v_end.size(); u++) {
      if (u)
        os << ",";
      os << (v_end[u]);
    }
    os << "];" << endl;
  }
}

/// output method as a "maple" recursive definition
template<typename T> inline void matrix<T>::maple_recursive(ostream& os) const {
  matrix<T> * m_t = Transpose();
  os << "{" << endl;
  // display each row
  bool first_element = false;
  for (typename std::vector< row<T> >::const_iterator i = m_t->_rows.begin(); i != m_t->_rows.end(); i++) {
    const row<T> & r = *i;
    if (first_element)
      os << "," << endl;
    os << "q" << (i-m_t->_rows.begin()) << "(n) = ";
    r.maple_recursive(os);
    first_element = true;
  }

  // display initial conditions
  for (typename std::vector< row<T> >::const_iterator i = m_t->_rows.begin(); i != m_t->_rows.end(); i++) {
    os << "," << endl;
    os << "q" << (i-m_t->_rows.begin()) << "(0) = " << (((i-m_t->_rows.begin()) == 1) ? 1 : 0);
  }
  os << endl << "}," << endl;

  // display variables
  os << "{" << endl;
  for (typename std::vector< row<T> >::const_iterator i = m_t->_rows.begin(); i != m_t->_rows.end(); i++) {
    if (i != m_t->_rows.begin())
      os << "," << endl;
    os << "q" << (i-m_t->_rows.begin());
  }
  os << "}" << endl;
  delete m_t;
}


/// input method for the current matrix
template<typename T> inline istream& operator>>(istream& is, matrix<T>& m) {
  // clear previous matrix
  m._rows.clear();

  // read matrix size
  int nbrows = 0;
  is >> nbrows;
#ifdef DEBUGREADING
  cerr << nbrows << endl;
#endif
  if (nbrows <= 0) {
    cerr << "> when reading matrix size" << endl;
    _ERROR("matrix operator>>","incorrect size "<< nbrows);
  }

  // add rows first
  for (int i = 0;i<nbrows;i++){
    m.addNewRow();
  }

  // get each row content
  int introw;
  for (int from = 0; from<nbrows; from++){
    is >> introw >> (m._rows[from]);

#ifdef DEBUGREADING
    cerr << "\t"  << introw << "\t" << (m._rows[from]);
#endif

    if (introw < 0 || introw >= nbrows) {
      cerr << "> when reading matrix cell line " << introw << " (line nb " << from << ")" << endl;
      _ERROR("matrix operator>>","incorrect row number " <<  introw << " (not in [0.."<<(nbrows-1)<<"])");
    }
  }
  return is;
}




/// insert an element to a row (or "add" it using the T "+" operator)
template<typename T> inline void row<T>::insert(const int j, const T v) {
  if (_sparse) {
    /* sparse implementation */
    std::pair<int,T> p(j,v);
    typename std::vector< std::pair<int,T> >::iterator it = lower_bound(_cells_sparse.begin(),_cells_sparse.end(),p,pairless<T>);
    if (it != _cells_sparse.end())
      if (it->first == j)
        it->second = it->second + v;
      else
        _cells_sparse.insert(it,p);
    else
      _cells_sparse.push_back(p);
  } else {
    /* dense implementation */
    if (_cells_dense.size() <= (unsigned)j)
      _cells_dense.resize(j+1, Zero<T>());
    _cells_dense[j] =  _cells_dense[j] + v;
  }
}



/// transpose a matrix
template<typename T> inline matrix<T> * matrix<T>::Transpose() const {
  matrix * result = new matrix();

  // find the max column on the original matrix
  int max_column = 0;

  for (unsigned i = 0; i < _rows.size(); i++) {
    max_column = MAX(max_column,_rows[i].max_index());
  }

  // build the transpose matrix in dense format first
  for (int i = 0; i <= max_column; i++) {
    result->addNewRow(false/* final does not mean anything here so set to false */, false /* dense */);
    result->_rows[i]._cells_dense.resize(_rows.size(), Zero<T>());
  }

  // fill the transpose matrix in dense format first then put the sparse row when needed ...
  for (unsigned i = 0; i < _rows.size(); i++) {
    if (_rows[i]._sparse) {
      for (typename std::vector< std::pair<int,T> >::const_iterator j=_rows[i]._cells_sparse.begin(); j != _rows[i]._cells_sparse.end(); j++) {
        result->_rows[j->first]._cells_dense[i] = j->second;
      }
    } else {
      for (unsigned j = 0 ; j < _rows[i]._cells_dense.size(); j++) {
        result->_rows[j]._cells_dense[i] = _rows[i]._cells_dense[j];
      }
    }
  }

  for (unsigned i = 0; i < result->_rows.size(); i++) {
    if (result->_rows[i].density() < MATRIX_SPARSE_ROW_DENSITY) {
      result->_rows[i].setsparse(true);
    }
  }
  return result;
}


/// matrix composition (or product)
template<typename T> inline matrix<T> * matrix<T>::Compose(const matrix<T> &other) const {
  matrix * result = new matrix();
  matrix * other_transpose = other.Transpose();

  for (unsigned i = 0; i < _rows.size(); i++) {
    result->addNewRow(_rows[i].final(),false);
    result->_rows[i]._cells_dense.resize(other_transpose->_rows.size(), Zero<T>()); // FIXME dense mode first
    for (unsigned j = 0; j < other_transpose->_rows.size(); j++) {
      if (_rows[i]._sparse) {
        if (other_transpose->_rows[j]._sparse) {
          /* 1/4 */
          typename std::vector< std::pair<int,T> >::const_iterator k_this =                   _rows[i]._cells_sparse.begin();
          typename std::vector< std::pair<int,T> >::const_iterator k_other = other_transpose->_rows[j]._cells_sparse.begin();
          while (k_this != _rows[i]._cells_sparse.end() && k_other != other_transpose->_rows[j]._cells_sparse.end()) {
            if (k_this->first == k_other->first) {
              result->_rows[i]._cells_dense[j] = result->_rows[i]._cells_dense[j] + (k_this->second * k_other->second);
              k_this++;
              k_other++;
            } else if (k_this->first < k_other->first) {
              k_this++;
            } else {
              k_other++;
            }
          }
        } else {
          /* 2/4 */
          for (typename std::vector< std::pair<int,T> >::const_iterator k_this = _rows[i]._cells_sparse.begin(); k_this < _rows[i]._cells_sparse.end(); k_this++) {
            if (other_transpose->_rows[j]._cells_dense[k_this->first] != Zero<T>()) {
              result->_rows[i]._cells_dense[j] = result->_rows[i]._cells_dense[j] + (k_this->second * other_transpose->_rows[j]._cells_dense[k_this->first]);
            }
          }
        }
      } else {
        if (other_transpose->_rows[j]._sparse) {
          /* 3/4 */
          for (typename std::vector< std::pair<int,T> >::const_iterator k_other = other_transpose->_rows[j]._cells_sparse.begin(); k_other < other_transpose->_rows[j]._cells_sparse.end(); k_other++) {
            if (_rows[i]._cells_dense[k_other->first] != Zero<T>()) {
              result->_rows[i]._cells_dense[j] = result->_rows[i]._cells_dense[j] + (_rows[i]._cells_dense[k_other->first] * k_other->second);
            }
          }
        } else {
          /* 4/4 */
          for (unsigned k=0; k < _rows[i]._cells_dense.size(); k++) {
            result->_rows[i]._cells_dense[j] = result->_rows[i]._cells_dense[j] + (_rows[i]._cells_dense[k] * other_transpose->_rows[j]._cells_dense[k]);
          }
        }
      }
    }
    if (result->_rows[i].density() < MATRIX_SPARSE_ROW_DENSITY) {
      result->_rows[i].setsparse(true);
    }
  }

  delete other_transpose;
  return result;
}


/// matrix "Pr" classical algorithm : compute the "T"-lity to be at a final or non final state after nbSteps on the same matrix
template<typename T> inline const T matrix<T>::Pr(const int nbSteps, const bool final) const {
  int i_min = 1;
  int i_max = 1;

  // set all the other probabilities or counts to 0 (PLUS neutral element) or cost to infinity (MIN neutral element)
  std::vector<T> v0(size(),Zero<T>());
  std::vector<T> v1(size(),Zero<T>());

  // set initial state probability or count to 1 (PRODUCT neutral element) or initial cost to 0 (PLUS neutral element)
  v1[1] = One<T>();

  // compute probability/cost/count at step i provided probabilities/costs/counts at step i-1
  for (int k = 0 ; k < nbSteps ; k++) {
    int j_min = size()-1;
    int j_max = 0;
    if (k&1) {
      for (int i = i_min; i <= i_max; i++) {
        if (v0[i] != Zero<T>()) { // p is not the "+" neutral element (0 for proba or count, or infinity for costs)
          if (_rows[i]._sparse) {
            for (typename std::vector< std::pair<int,T> >::const_iterator j=_rows[i]._cells_sparse.begin(); j != _rows[i]._cells_sparse.end(); j++) {
              j_min = MIN(j_min,j->first);
              j_max = MAX(j_max,j->first);
              v1[j->first] = v1[j->first] + ((j->second) * (v0[i]));
            }
          } else {
            for (unsigned j = 0 ; j < _rows[i]._cells_dense.size(); j++) {
              if (_rows[i]._cells_dense[j] != Zero<T>()) {
                j_min = MIN(j_min,(int)j);
                j_max = MAX(j_max,(int)j);
                v1[j] = v1[j] + ((_rows[i]._cells_dense[j]) * (v0[i]));
              }
            }
          }
          v0[i] = Zero<T>();
        }
      }
    } else {
      for (int i = i_min; i <= i_max; i++) {
        if (v1[i] != Zero<T>()) { // p is not the "+" neutral element (0 for proba or count, or infinity for costs)
          if (_rows[i]._sparse) {
            for (typename std::vector< std::pair<int,T> >::const_iterator j = _rows[i]._cells_sparse.begin(); j != _rows[i]._cells_sparse.end(); j++) {
              j_min = MIN(j_min,j->first);
              j_max = MAX(j_max,j->first);
              v0[j->first] = v0[j->first] + ((j->second) * (v1[i]));
            }
          } else {
            for (unsigned j = 0 ; j < _rows[i]._cells_dense.size(); j++) {
              if (_rows[i]._cells_dense[j] != Zero<T>()) {
                j_min = MIN(j_min,(int)j);
                j_max = MAX(j_max,(int)j);
                v0[j] = v0[j] + ((_rows[i]._cells_dense[j]) * (v1[i]));
              }
            }
          }
          v1[i] = Zero<T>();
        }
      }
    }
    i_min = j_min;
    i_max = j_max;
  } // for


  T result = Zero<T>();
  // Sum "final" states probabilities or counts, or find the Min "non final" cost (NB : (+) ~ min )
  if (final) {
    if (nbSteps & 1) {
      for (int i=i_min;i<=i_max;i++)
        if (_rows[i].final())
          result = (result + (v0[i]));
    } else {
      for (int i=i_min;i<=i_max;i++)
        if (_rows[i].final())
          result = (result + (v1[i]));
    }
  } else {
    if (nbSteps & 1) {
      for (int i=i_min;i<=i_max;i++)
        if (!(_rows[i].final()))
          result = (result + (v0[i]));
    } else {
      for (int i=i_min;i<=i_max;i++)
        if (!(_rows[i].final()))
          result = (result + (v1[i]));
    }
  }

  v0.clear();
  v1.clear();
  return result;
}


/// matrix "Pr" classical algorithm : compute the "T"-lity to be at a transitive-sum of the final values after nbSteps on the same matrix
template<typename T> inline std::vector<T> * matrix<T>::Pr_transitive_final(const int nbSteps, const unsigned max_final_value, const int sub_final_value) const {
  int i_min = 1;
  int i_max = 1;
  unsigned u_max = 0;
  // set all the other probabilities or counts to 0 (PLUS neutral element) or cost to infinity (MIN neutral element)
  std::vector< std::vector<T> > v0(size(), std::vector<T>(1,Zero<T>()));
  std::vector< std::vector<T> > v1(size(), std::vector<T>(1,Zero<T>()));

  // set initial state probability or count to 1 (PRODUCT neutral element) or initial cost to 0 (PLUS neutral element)
  v1[1][0] = One<T>();

  // compute probability/cost/count at step i provided probabilities/costs/counts at step i-1
  for (int k = 0 ; k < nbSteps ; k++) {
    int j_min = size()-1;
    int j_max = 0;
    u_max = 0;
    if (k&1) {
      for (int i = i_min; i <= i_max; i++) {
        for (unsigned u = 0; u < v0[i].size(); u++) {
          if (v0[i][u] != Zero<T>()) { // p is not the "+" neutral element (0 for proba or count, or infinity for costs)
            if (_rows[i]._sparse) {
              for (typename std::vector< std::pair<int,T> >::const_iterator j=_rows[i]._cells_sparse.begin(); j != _rows[i]._cells_sparse.end(); j++) {
                j_min   =  MIN(j_min,j->first);
                j_max   =  MAX(j_max,j->first);
                unsigned u_delta =  u  +  _rows[j->first].final();
                while (u_delta > max_final_value)
                  u_delta -= sub_final_value;
                if (v1[j->first].size() <= u_delta)
                  v1[j->first].resize(u_delta+1,Zero<T>());
                u_max = MAX(u_max,u_delta);
                v1[j->first][u_delta] = v1[j->first][u_delta] + ((j->second) * (v0[i][u]));
              }
            } else {
              for (unsigned j = 0 ; j < _rows[i]._cells_dense.size(); j++) {
                if (_rows[i]._cells_dense[j] != Zero<T>()) {
                  j_min = MIN(j_min,(int)j);
                  j_max = MAX(j_max,(int)j);
                  unsigned u_delta =  u  +  _rows[j].final();
                  while (u_delta > max_final_value)
                    u_delta -= sub_final_value;
                  if (v1[j].size() <= u_delta)
                    v1[j].resize(u_delta+1,Zero<T>());
                  u_max = MAX(u_max,u_delta);
                  v1[j][u_delta] = v1[j][u_delta] + ((_rows[i]._cells_dense[j]) * (v0[i][u]));
                }
              }
            }
            v0[i][u] = Zero<T>();
          }
        }
      }
    } else {
      for (int i = i_min; i <= i_max; i++) {
        for (unsigned u = 0; u < v1[i].size(); u++) {
          if (v1[i][u] != Zero<T>()) { // p is not the "+" neutral element (0 for proba or count, or infinity for costs)
            if (_rows[i]._sparse) {
              for (typename std::vector< std::pair<int,T> >::const_iterator j = _rows[i]._cells_sparse.begin(); j != _rows[i]._cells_sparse.end(); j++) {
                j_min = MIN(j_min,j->first);
                j_max = MAX(j_max,j->first);
                unsigned u_delta =  u  +  _rows[j->first].final();
                while (u_delta > max_final_value)
                  u_delta -= sub_final_value;
                if (v0[j->first].size() <= u_delta)
                  v0[j->first].resize(u_delta+1,Zero<T>());
                u_max = MAX(u_max,u_delta);
                v0[j->first][u_delta] = v0[j->first][u_delta] + ((j->second) * (v1[i][u]));
              }
            } else {
              for (unsigned j = 0 ; j < _rows[i]._cells_dense.size(); j++) {
                if (_rows[i]._cells_dense[j] != Zero<T>()) {
                  j_min = MIN(j_min,(int)j);
                  j_max = MAX(j_max,(int)j);
                  unsigned u_delta =  u  +  _rows[j].final();
                  while (u_delta > max_final_value)
                    u_delta -= sub_final_value;
                  if (v0[j].size() <= u_delta)
                    v0[j].resize(u_delta+1,Zero<T>());
                  u_max = MAX(u_max,u_delta);
                  v0[j][u_delta] = v0[j][u_delta] + ((_rows[i]._cells_dense[j]) * (v1[i][u]));
                }
              }
            }
          }
          v1[i][u] = Zero<T>();
        }
      }
    }
    i_min = j_min;
    i_max = j_max;
  } // for

  std::vector<T> * result = new std::vector<T>(u_max+1,Zero<T>());
  // Sum "final" states probabilities or counts, or find the Min "non final" cost (NB : (+) ~ min )
  if (nbSteps & 1) {
    for (int i=i_min;i<=i_max;i++)
      for (unsigned u = 0; u < v0[i].size(); u++)
        (*result)[u] = ((*result)[u] + (v0[i][u]));
  } else {
    for (int i=i_min;i<=i_max;i++)
      for (unsigned u = 0; u < v1[i].size(); u++)
        (*result)[u] = ((*result)[u] + (v1[i][u]));
  }
  v0.clear();
  v1.clear();
  return result;
}


/// matrix "Pr_one_step_from_one" algorithm : compute the "T"-lity to be at a final or non final state after one step from 1 in the original matrix to the m_final one
template<typename T> inline const T matrix<T>::Pr_one_step_from_one(const matrix<T> &m_final, const bool final) const {
  T result = Zero<T>();
  // Sum "final" states probabilities or counts, or find the Min "non final" cost (NB : (+) ~ min )
  if (final) {
    if (_rows[1]._sparse) {
      for (typename std::vector< std::pair<int,T> >::const_iterator j = _rows[1]._cells_sparse.begin(); j != _rows[1]._cells_sparse.end(); j++)
        if (m_final._rows[j->first].final())
          result = (result + (j->second));
    } else {
      for (unsigned int j = 0 ; j < _rows[1]._cells_dense.size(); j++) {
        if (m_final._rows[j].final())
          result = (result + (_rows[1]._cells_dense[j]));
      }
    }
  } else {
    if (_rows[1]._sparse) {
      for (typename std::vector< std::pair<int,T> >::const_iterator j = _rows[1]._cells_sparse.begin(); j != _rows[1]._cells_sparse.end(); j++)
        if (!(m_final._rows[j->first].final()))
          result = (result + (j->second));
    } else {
      for (unsigned int j = 0 ; j < _rows[1]._cells_dense.size(); j++) {
        if (!(m_final._rows[j].final()))
          result = (result + (_rows[1]._cells_dense[j]));
      }
    }
  }
  return result;
}
// @}







/** @addtogroup matrix
 */
// @{


/** @brief Compute the slicer structure needed to quickly computer overlapped windows sensitivity :
 *       @li _i and _j are the left/right indices of the window
 *       @li _matrices_data keep the matrices (and update them to keep a "binary-trie" structure that fasten the possible products between _i and _j)
 *       @li _left/_right matrices block and index are also there to help in the computation of successive windows :
 *           _left is maintained as a stack of matrices whereas _right is a single matrix.
*
 */

template<typename T> class matrices_slicer {

public :

  /** @brief constructor for a slicer from the set of matrices obtained with a "slice" (or "step") product
   *  @param data is the vector of sliced matrices been used (and erased when finished)
   *  @see  matrices_step_pr_product
   *  @see  matrices_step_cost_product
   */
  matrices_slicer(std::vector< matrix<T> * > * data): _matrices_data(data), _i(0), _j(0), _middle(0) {
    _left_matrices_block = stack<matrix<T>*>();
    _left_matrices_index = stack<int>();
    _right_matrix_block = NULL;
    _right_matrix_index = 0;
#ifdef MATRIX_SLICER_STATS
    nb_prod_a           = 0;
    nb_prod_b           = 0;
    nb_prod_c           = 0;
    nb_prod_d           = 0;
    nb_middle_change    = 0;
#endif
#ifndef NO_FULL_MATRIX_BLOCK
    _full_matrix_block   = NULL;
    _full_matrix_block_i = -1;
    _full_matrix_block_j = -1;
#ifdef MATRIX_SLICER_STATS
    nb_prod_e = 0;
#endif
#endif

  };


  /** destructor for a slicer
   */
  ~matrices_slicer() {
    // empty stacks
    empty_stacks();
    // delete input data
    for (unsigned i=(unsigned)_i; i < _matrices_data->size();i++) {
      (*_matrices_data)[i]->clear();
      delete (*_matrices_data)[i];
    }
    _matrices_data->clear();
    delete _matrices_data;
    _matrices_data = NULL;

    _middle        = 0;
    _i             = 0;
    _j             = 0;
  }



  /** @brief empty stacks
   */
  inline void empty_stacks() {
    while(!_left_matrices_index.empty()) {
      _left_matrices_index.pop();
      delete _left_matrices_block.top();
      _left_matrices_block.pop();
    }
    _right_matrix_index = 0;
    if (_right_matrix_block) {
      delete _right_matrix_block;
      _right_matrix_block = NULL;
    }
  };



  /** @brief right_size
   *  "right" size of a block starting at "s" in the dataset and that must be < "j"
   *  @param s is the starting position of the block
   *  @param j is the right position that must be not crossed (e.g. put "_j+1")
   *  @return the length of the block starting at "s"
   *  @todo{FIXME: check the function due to s}
   */
  inline int right_size(const int s, const int j) const {
    int u = s;
    int t = 1; // 1,2,4,8
    while ((u % 2 == 0) && (s + (t << 1) <= j )) {
      t  = t << 1;
      u  = u >> 1;
    }
    return t;
  }

  /** @brief left_size
   *  "right" size of a block finishing at "s-1" in the dataset and that must be >= "i"
   *  @param i is the left position that can be crossed but not sub-overlapped
   *  @param s is the ending + 1 position of the block (e.g. put "_j+1")
   *  @return the length of the block ending at "s-1"
   *  @todo{FIXME: check the function due to s}
   */

  inline int left_size(const int i, const int s) const {
    int u = s;
    int t = 1; // 1,2,4,8
    while ((u % 2 == 0) && (i + (t << 1) <= s)) {
      t     = t << 1;
      u     = u >> 1;
    }
    return t;
  }

  /** @brief middle_pos_right
   *  @param i is the left pos
   *  @param j is the right pos
   *  @return the end " + 1" of the maximal block included withing i and < j
   *  @todo{FIXME: check the function due to s}
   */
  inline int middle_pos_right(const int i, const int j) const {
    int s = i;
    int u = 0;
    do {
      int un = right_size(s,j);
      if (un > u) {
        u = un;
        s += u;
      } else {
        return s;
      }
    } while (1);
    return s;
  };

  /** @brief middle_pos_left
   *  @param i is the left pos
   *  @param j is the right pos
   *  @return the beginning of the maximal block included withing i and < j
   *  @todo{FIXME: check the function due to s}
   */
  inline int middle_pos_left(const int i, const int j) const {
    int s = j;
    int u = 0;
    do {
      int un = left_size(i,s);
      if (un > u) {
        u = un;
        s -= u;
      } else {
        return s;
      }
    } while(1);
    return s;
  };

  /** @brief increment _j, add the rightmost matrix _matrices_data[_j] (already done) and update matrices blocks that now include _matrices_data[_j]
   */
  inline void add_right() {
    _j++;
    {
      int u     = _j+1;
      int t_old = 0;
      int t     = 1; // 1,3,7,15 ...
      while ((u % 2 == 0) && (_j - t >= _i)) {
        matrix<T> * tmp_j_t       = (*_matrices_data)[_j - t];
        (*_matrices_data)[_j - t] = tmp_j_t->Compose(*((*_matrices_data)[_j - t_old]));
#ifdef MATRIX_SLICER_STATS
        nb_prod_a++;
#endif
        delete      tmp_j_t;
        t_old = t;
        t    += t + 1;
        u     = u >> 1;
      }
    }
    // erase when middle changes
    int n_middle = middle_pos_left(_i,_j+1);
    if (n_middle != _middle) {
      // FIXME dont empty stack by reusing them
      empty_stacks();
      _middle = n_middle;
#ifdef MATRIX_SLICER_STATS
      nb_middle_change++;
#endif
    }
  }


  /** @brief delete the leftmost  matrix and increment _i
   */
  inline void del_left() {
    if (_i < _j) {
      delete (*_matrices_data)[_i];
      (*_matrices_data)[_i] = NULL;
      _i++;
      // erase when middle changes
      int n_middle = middle_pos_left(_i,_j+1);
      if (n_middle != _middle) {
        empty_stacks();
        _middle = n_middle;
#ifdef MATRIX_SLICER_STATS
        nb_middle_change++;
#endif
      }
    } else {
      _WARNING("[matrix.hh] matrices_slicer.del_left()","_i cannot be > _j : del_left() is ignored.");
    }
  }


  /** @brief  compute the current windows "probability/minimal cost" obtained from
   *          the matrix product from _i to _j
   *  @param  final (default true) select if final or non final states have to be reached.
   *  @return a "probability/minimal cost" depending on the T template
   */

  inline const T current_Pr(const bool final = true) {
    matrix<T> *     left_result  = NULL;
    matrix<T> *     right_result;
    matrix<T> * left_right_result = NULL;

#ifdef MATRIX_SLICER_BUILDPROD
    cerr << (*this) << endl;
    cerr << "* current_product(" << _i << "," << _j << ")  {middle = " << _middle << "}" << endl;
#endif


    // 0) Full matrix block
#ifndef NO_FULL_MATRIX_BLOCK

    if ((_full_matrix_block != NULL) && (_i == _full_matrix_block_i)) {
      // a) use the previous computation when _i has not changed ...
      int lj = right_size(_full_matrix_block_j,_j+1);
      while (_full_matrix_block_j + lj <= _j+1) { // FIXME  : <=  +1?
#ifdef MATRIX_SLICER_BUILDPROD
        cerr << "    MB(compose) -> = " << _full_matrix_block_j << " + " << lj << endl;
#endif
        matrix<T> * old_full_matrix_block = _full_matrix_block;
                       _full_matrix_block = _full_matrix_block->Compose(*((*_matrices_data)[_full_matrix_block_j]));
        delete      old_full_matrix_block;

#ifdef MATRIX_SLICER_STATS
        nb_prod_e++;
#endif

        _full_matrix_block_j += lj;
        lj = right_size(_full_matrix_block_j,_j+1);
      }
      // and jump to te computation
      goto full_matrix_reuse;
    } else {
      // b) delete since full_block will not be usefull when _i has changed ...
      _full_matrix_block_i = -1;
      _full_matrix_block_j = -1;
      if (_full_matrix_block)
        delete _full_matrix_block;
      _full_matrix_block = NULL;
    }
#endif




    // I) Right
    {
      int j_stack;
      matrix<T> * old_right_result = NULL;

      // a) get right border already computed from _middle if it exists ...
      if (_right_matrix_block != NULL) {
        right_result = _right_matrix_block;
        j_stack      = _right_matrix_index;

#ifdef MATRIX_SLICER_BUILDPROD
        cerr << " R(from index) -> = " << j_stack << endl;
#endif

      } else {
        right_result     = new matrix<T>(*((*_matrices_data)[_middle]));
        j_stack          = _middle + right_size(_middle,_j+1);

        old_right_result = right_result;

#ifdef MATRIX_SLICER_BUILDPROD
        cerr << "        _middle -> = " << _middle << ", _middle_length -> = " << (right_size(_middle,_j+1)) << endl;
        cerr << " R(from middle) -> = " << j_stack << endl;
#endif
      }

      // b) compute the right part and push it to the right stack ...
      int lj = right_size(j_stack,_j+1);
      while (j_stack + lj <= _j+1) { // FIXME  : <=  ?

#ifdef MATRIX_SLICER_BUILDPROD
        cerr << " R(compose) -> = " << j_stack << " + " << lj << endl;
#endif

        right_result = right_result->Compose(*((*_matrices_data)[j_stack]));
#ifdef MATRIX_SLICER_STATS
        nb_prod_b++;
#endif

        if (old_right_result)
          delete old_right_result;
        old_right_result = right_result;

        j_stack += lj;
        lj = right_size(j_stack,_j+1);
      }

      // c) update right_matrix_block when needed
      if (_right_matrix_block == right_result) {
        // NO CHANGE !!
#ifdef MATRIX_SLICER_BUILDPROD
        cerr << " R(finalize) -> = NO CHANGE " << j_stack << "=" << _right_matrix_index << endl;
#endif
      } else {
        // SOME CHANGE !!
#ifdef MATRIX_SLICER_BUILDPROD
        cerr << " R(finalize) -> = CHANGE from " << _right_matrix_index << " to " << j_stack << endl;
#endif
        if (_right_matrix_block != NULL)
          delete _right_matrix_block;
        _right_matrix_block = right_result;
        _right_matrix_index = j_stack;
      }
    }

    // II) Left

    // a) get left border already computed if exists ...
    {
      int i_stack;
      while ((!_left_matrices_index.empty()) && ((i_stack = _left_matrices_index.top()) < _i)) {
        _left_matrices_index.pop();
        matrix<T> * t = _left_matrices_block.top();
        _left_matrices_block.pop();
        delete t;
#ifdef MATRIX_SLICER_BUILDPROD
        cerr << " L(stack) -> popping " << i_stack<< endl;
#endif
      }

      if (!_left_matrices_index.empty()) {
        left_result = _left_matrices_block.top();
        i_stack     = _left_matrices_index.top();
#ifdef MATRIX_SLICER_BUILDPROD
        cerr << " L(stack) -> Search success at " << i_stack << endl;
#endif
      } else {
        i_stack = _middle;
#ifdef MATRIX_SLICER_BUILDPROD
        cerr << " L(stack) -> Search failure (empty stack after search)" << endl;
#endif
      }


      // b) compute (and store) the rest of the left elements ...
      int li = left_size(_i,i_stack);
      while (i_stack - li >= _i) {

#ifdef MATRIX_SLICER_BUILDPROD
        cerr << " L(compose) -> = " << i_stack << "-" << li << " = " << (i_stack - li) << endl;
#endif

        i_stack -= li;

        if (left_result) {
          left_result = (*_matrices_data)[i_stack]->Compose(*left_result);
#ifdef MATRIX_SLICER_STATS
          nb_prod_c++;
#endif

        } else {
          left_result   = new matrix<T>(*((*_matrices_data)[i_stack]));
        }
#ifdef MATRIX_SLICER_BUILDPROD
        cerr << " L(push) -> = " << i_stack << endl;
#endif

        _left_matrices_index.push(i_stack);
        _left_matrices_block.push(left_result);

        li = left_size(_i,i_stack);
      }
    }


    // III) Combine Left and Right
    if (left_result) {

#ifdef MATRIX_SLICER_BUILDPROD
      cerr << " M(join)" << endl;
#endif
      left_right_result = left_result->Compose(*right_result);
#ifdef MATRIX_SLICER_STATS
      nb_prod_d++;
#endif

    } else {
      left_right_result = right_result;
#ifdef MATRIX_SLICER_BUILDPROD
      cerr << " M(free join)" << endl;
#endif
    }


    // IV) Full_Matrix "Reuse"
#ifndef NO_FULL_MATRIX_BLOCK
    if (_full_matrix_block != NULL) {
    full_matrix_reuse:
      left_right_result = _full_matrix_block;
    } else {
      if (left_result)
        _full_matrix_block = left_right_result;
      else
        _full_matrix_block = new matrix<T> (*left_right_result); // in fact "left_right_result"="right_result" here, so must be copied to avoid "double free"
      _full_matrix_block_i = _i;
      _full_matrix_block_j = _j+1;
    }
#endif



    // V) Compute sensitivity or count, or lossless property
    T result = Zero<T>();
    matrix<T> * _m_final = (*_matrices_data)[_j+1];

    VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  - sliced matrix (" << (_j+1) << ") size : " << (_m_final->size())););

    // Sum "final" states probabilities or find the Min "non final" cost (NB : (+) ~ min )
    if (final) {
      if (left_right_result->_rows[1]._sparse) {
        for (typename std::vector< std::pair<int,T> >::const_iterator j = left_right_result->_rows[1]._cells_sparse.begin(); j != left_right_result->_rows[1]._cells_sparse.end(); j++)
          if (_m_final->_rows[j->first].final())
            result = (result + (j->second));
      } else {
        for (unsigned j = 0 ; j <  left_right_result->_rows[1]._cells_dense.size(); j++) {
          if (_m_final->_rows[j].final())
            result = (result + (left_right_result->_rows[1]._cells_dense[j]));
        }
      }
    } else {
      if (left_right_result->_rows[1]._sparse) {
        for (typename std::vector< std::pair<int,T> >::const_iterator j = left_right_result->_rows[1]._cells_sparse.begin(); j != left_right_result->_rows[1]._cells_sparse.end(); j++)
          if (!(_m_final->_rows[j->first].final()))
            result = (result + (j->second));
      } else {
        for (unsigned j = 0 ; j <  left_right_result->_rows[1]._cells_dense.size(); j++) {
          if (!(_m_final->_rows[j].final()))
            result = (result + (left_right_result->_rows[1]._cells_dense[j]));
        }
      }
    }
// delete if not used
#ifdef NO_FULL_MATRIX_BLOCK
   if (left_result)
     delete left_right_result;
#endif
   return result;
}



  /**@name  IO methods
   * @brief IO streams to print/store/load matrix
   */

  // @{
  /// print matrices_slicer information
  template<typename U> friend ostream& operator<< (ostream& os, const matrices_slicer<U>& ms);
  // @}

#ifdef MATRIX_SLICER_STATS
  /** @name statistics of the number of products done within the matrix_slicer
   */
  // @{
  unsigned              nb_prod_a,nb_prod_b,nb_prod_c,nb_prod_d,nb_middle_change;
#ifndef NO_FULL_MATRIX_BLOCK
  unsigned              nb_prod_e;
#endif
  // @}
#endif

  /** @brief data to be sliced (will be modified so must not be !!!)
   *
   */
  std::vector< matrix<T> * >  * _matrices_data;
  /// left window index (del only)
  int _i;
  /// right window index (add only)
  int _j;
  /// separator left/right = defined by (2^k) * r such that k is the largest integer (and r is a odd integer obviously) : i <= (2^k) * r <= j
  int _middle;
  /// left matrices (when preprocessed before...)
  // @{
  stack< matrix<T> * >  _left_matrices_block;
  stack< int >          _left_matrices_index;
  // @}
  /// right matrix (when preprocessed before...)
  // @{
  matrix<T> *           _right_matrix_block;
  int                   _right_matrix_index;
  // }@
#ifndef NO_FULL_MATRIX_BLOCK
  /// full block (when preprocessed before...)
  // @{
  matrix<T> *           _full_matrix_block;
  int                   _full_matrix_block_i;
  int                   _full_matrix_block_j;
  // @}
#endif
};




/// output method for the current matrix
template<typename T> inline ostream& operator<<(ostream& os, const matrices_slicer<T>& ms) {
  os << "i = " << (ms._i) << "\t" << "j = " << (ms._j) << "\t" << "m = " << (ms._middle) << "{" << (ms.middle_pos_left(ms._i,ms._j+1)) << "," << (ms.middle_pos_right(ms._i,ms._j+1)) << "}" << endl;

#ifdef MATRIX_SLICER_STATS
  os << (ms.nb_prod_a) << "\t" << (ms.nb_prod_b) << "\t" << (ms.nb_prod_c) << "\t" << (ms.nb_prod_d) << "\t=\t" <<
    (
     (ms.nb_prod_a)
     +
     (ms.nb_prod_b)
     +
     (ms.nb_prod_c)
     +
     (ms.nb_prod_d)
#ifndef NO_FULL_MATRIX_BLOCK
     +
     (ms.nb_prod_e)
#endif
     ) << "\t" << (ms.nb_middle_change) << "\t" << endl;
#endif

  if (ms._right_matrix_block)
    os << "  _right : " << (ms._right_matrix_index) << endl;
  else
    os << "  _right : unset" << endl;

  if (ms._left_matrices_index.size() > 0)
    os << "  _left  : "  << (ms._left_matrices_index.top()) << endl;
  else
    os << "  _left  : empty" << endl;

  os << "  _matrices_data[] (size) : ";
  for (int i = ms._i; i <= ms._j; i++) {
    os << " m[" << i << "] : " << (((*(ms._matrices_data))[i])->size()) << "{" << (((*(ms._matrices_data))[i])->fullsize()) << "} ; ";
  }
  os << endl;
  // FIXME to be continued
  return os;
}




// @}

#endif
