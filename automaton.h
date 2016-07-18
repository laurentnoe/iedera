/** @page automata Automata
 *  @brief Automata description and functions
 *  @tableofcontents
 *
 *  @section automata-description Description
 *  This part describes an automaton : each @ref automaton is mainly a represented by a set of @ref state , each itself being represented by a set of @ref transition.
 *
 *  An automaton can be deterministic or not, probabilistic or not.
 *
 *  By default the automaton constructor is almost empty (it creates only a final state 0 and the init state 1), but several methods are proposed to construct @ref seed-automaton, @ref probabilistic-automaton, @ref structural-automaton (@ref automaton-construction). Several methods are also proposed to manipulate theses automata (@ref automaton-manipulate), compute properties (@ref automaton-computed-properties), convert them into matrices (@ref automaton-matrix-conversion),
 *
 *  @section automaton-construction Construction
 *
 *  @subsection seed-automaton seed automaton
 *
 *  Three seed models are supported with their associated automata : subset seeds (built with the method @ref automaton::Automaton_SeedPrefixesMatching for loops, with the method @ref automaton::Automaton_SeedPrefixesMatching_CoverageDM for global coverage measure, or with the method @ref automaton::Automaton_SeedLinearMatching for linear measure without loops), vectorized seeds (built with a Aho-Corasick like @ref automaton::Automaton_SeedBuhler), or more generally vectorized subset seeds (built with the combined @ref automaton::Automaton_SeedScore).
 *
 * The "Cost" variants of the subset seed model (@ref automaton::Automaton_SeedPrefixesMatchingCost) and the vectorized subset seed model (@ref automaton::Automaton_SeedScoreCost) are proposed variants that prune the automaton when Lossless construction set is used (@ref gv_lossless_flag).
 *
 *  @subsection probabilistic-automaton probabilistic automaton
 *
 *  Two models are directly supported : Bernoulli (@ref automaton::Automaton_Bernoulli) and Markov (@ref automaton::Automaton_Markov). Others can be read via the automaton::operator>>().
 *
 *  @subsection structural-automaton structural automaton
 *
 *  Three structural models are proposed for various objectives : The homogeneous alignments (@ref automaton::Automaton_Homogeneous) are proposed to guaranty that any sub-part of the alignment cannot have a higher score than the full one; To set seeds on a (cyclic) set of positions, the cycle automaton is used (@ref automaton::Automaton_Cycle) ; To get the automaton of alignments that must be rejected by the lossless constraints, the lossless automaton is used (@ref automaton::Automaton_Lossless) ; To get the automaton that counts a set of matching symbols, the symbols-counting automaton (@ref automaton::Automaton_CountAlphabetSymbols) can be used.

 *  @section automaton-manipulate Manipulation
 *
 *  Several methods are proposed to manipulate or combine automata. The more general is the product (@ref automaton::product), that realizes the union, intersection, exclusion, with possible effect on merging final states. The second one ask to compute the @f$m@f$ occurrence of a final state to accept the alignment, and is thus useful for multiple hits and also global coverage measure (@ref automaton::mhit). The last one builds the minimization of the current automaton with the Hopcroft algorithm (@ref automaton::Hopcroft), respecting different labelled final states.
 *
 *  @section automaton-computed-properties computing properties
 *
 *  From the probabilities attached to each transition (@ref automaton::Prob), it is possible to compute the probability to be at a final state after some @f$n@f$ steps without any restriction (@ref automaton::PrFinal), or with the restriction to be in the lossless accepted language (@ref automaton::PrLossless). It's also possible to know (when this probability is too close to 1.0 so that it is rounded to 1.0) if the settings are truly lossless (@ref automaton::Lossless).
 *
 *  @section automaton-matrix-conversion converting into matrices
 *  Several methods are proposed to convert the @f$ Automaton \times ProbabilisticModel @f$, @f$ Automaton \times CostModel @f$, or @f$ Automaton \times CountModel @f$ into matrices for more convenient computations. The probability @ref automaton::matrix_pr_product , the cost @ref automaton::matrix_cost_product and the counting @ref automaton::matrix_count_product give the product of two automata into a resulting matrix @f$M@f$ (that store either probabilities, costs, or counts) so that @f$M^n@f$ usually compute the needed properties. There are also three "stepwise equivalent" methods @ref automaton::matrices_step_pr_product,  @ref automaton::matrices_step_cost_product, and @ref automaton::matrices_step_count_product : these three methods give the "breadth first" product as an ordered set of matrices  @f$M_1,M_2,M_3\ldots,M_l@f$, thus enabling any computation @f$M_i,M_{i+1}\ldots,M_{j}@f$ @f$\forall 0 \leq i < j \leq l@f$ @see matrix @see matrices_slicer .
 *
 *  @todo{FIXME : to be continued}
 */

#ifndef __AUTOMATON_H__
#define __AUTOMATON_H__

//STL
#include <vector>
#include <queue>
#include <stack>
#include <functional>
#include <map>
//STD
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
//STR
using namespace std;
#include "macro.h"
#include "seed.h"

#ifndef NOMATRIX
#include "matrix.hh"
#endif

#ifdef USEINFINT
#include "infint.hh"
#endif

class automaton;

/** Define the product families of operators
 */

typedef enum ProductSetFinalType {
  PRODUCT_UNION_FINAL_LOOP,
  PRODUCT_UNION_NO_FINAL_LOOP,
  PRODUCT_UNION_NO_FINAL_LOOP_ADD,
  PRODUCT_INTERSECTION_FINAL_LOOP,
  PRODUCT_INTERSECTION_NO_FINAL_LOOP,
  PRODUCT_BUTNOT_FINAL_LOOP,
  PRODUCT_BUTNOT_NO_FINAL_LOOP,
  PRODUCT_NOTBUT_FINAL_LOOP,
  PRODUCT_NOTBUT_NO_FINAL_LOOP,
  PRODUCT_ADDHOC_FINAL_LOOP,
  PRODUCT_ADDHOC_NO_FINAL_LOOP,
} ProductSetFinalType;


/** Define the probability family of operators
 */

typedef enum ProductProbabilityType {
  PRODUCT_NONE_IS_PROBABILIST,
  PRODUCT_THIS_IS_PROBABILIST,
  PRODUCT_OTHER_IS_PROBABILIST,
  PRODUCT_BOTH_ARE_PROBABILIST,
} ProductProbabilityType;


typedef int (*AddHoc_Final_Func)(int, int);


/// final
#define TRUE  1
/// non final
#define FALSE 0

/** @defgroup automaton automaton class
 *  @brief automaton (deterministic or non-deterministic, probabilistic or not)
 */
// @{

/**
 * @class transition
 *
 * @brief transition to a given state (on a given letter)
 */

class transition {
 public:
  /** @brief build a transition object
   *  @param state gives the state number reached by this transition
   *  @param prob  gives the probability of this transition
   */
 transition(int state = 0, double prob = 1e+0) : _state(state), _prob(prob) {};

 protected :
  /// state number to be reached
  int     _state;
  /// transition probability (1.0 if deterministic automaton)
  double  _prob;

  /// state is a friend class to ease access
  friend class state;
  /// automaton is a friend class to ease access
  friend class automaton;
  /// print automaton information is friend
  friend ostream& operator<<(ostream& os, const automaton& automaton);
  /// read  automaton information
  friend istream& operator>>(istream& is, automaton& automaton);
};

/**
 * @class state
 *
 * @brief state represented as a table of list of transitions ("next[gv_align_alphabet_size]")
 *
 */

class state {
 public:
  /// Build an empty state (empty transition lists)
 state(int final = 0): _final(final) { _next  =  vector< vector < transition > >(gv_align_alphabet_size, vector < transition > ());};


  /// Erase a state (clear transition lists first)
  ~state() { for(int a = 0; a < gv_align_alphabet_size; a++) _next[a].clear();  _next.clear(); };

  /** @brief Clear transition lists
   */
  void clear() { for(int a = 0 ; a < gv_align_alphabet_size ; a++) _next[a].clear();  _next.clear(); };

 protected:
  /// forward  transitions list on letter A
  vector < vector < transition > >  _next;
  /// final state
  int _final;

  /// automaton is a friend class to ease access
  friend class automaton;
  /// print automaton information is friend
  friend ostream& operator<<(ostream& os, const automaton& automaton);
  /// read  automaton information
  friend istream& operator>>(istream& is, automaton& automaton);

};




/**
 * @class automaton
 *
 * @brief automaton (deterministic or non-deterministic, probabilistic or not) roughly represented as a vector of states
 *
 *  this one can either be affected to:
 *      @li a seed automaton  (deterministic, initial state is 1, final states have the final tag)
 *      @li a probability model (Bernoulli, Markov, and old M3,M14,HMM models, can either be deterministic or non-deterministic)
 *      @li a target alignment model (deterministic)
 *
 */

class automaton {
 public:

  /// constructor for an empty automaton
  automaton() {
         _states = vector<state>(0);
    _init_states = vector<int>(0);
  };

  /// destructor for all possibles automata build by this class
  ~automaton() {
    for(int i = 0 ; i < (int)_states.size() ; i++){
      _states[i].clear();
    }
         _states.clear();
    _init_states.clear();
  };

  /** @name  Build automata
   *  @brief each one builds a special automaton
   */

  // @{

  /** @brief simple algorithm to build a linear automaton of the current seed "s", that only matches if the seed match from its very beginning position
   *         the very beginning of the sequence : only useful for covariance computation
   *  @param s is the seed descriptor
   *  @param matchingmatrix is a boolean matrix that gives for alignment letter "a", and seed letter "b", the matching with "matrix[a][b]"
   *  @return 0
   */
  int Automaton_SeedLinearMatching (const seed & s,
                                    const vector< vector <int> > & matchingmatrix);

  /** @brief Old algorithm in @f$ \mathcal{O}(w 2^{s-w})@f$ for subset seed matching (use on index of size @f$ \mathcal{O}(w 2^{s-w})@f$ )
   *  @param s is the seed descriptor
   *  @param matchingmatrix is a boolean matrix that gives for alignment letter "a", and seed letter "b", the matching with "matrix[a][b]"
   *  @param nomerge if final states dont have to be merged together into a single state (default behaviour is merging)
   *  @return 0
   */

  int Automaton_SeedPrefixesMatching_old (const seed & s,
                                          const vector< vector <int> > & matchingmatrix,
                                          const bool nomerge=false);


  /** @brief New linear algorithm for subset seed matching (CIAA 2007)
   *  @param s is the seed descriptor
   *  @param matchingmatrix is a boolean matrix that gives for alignment letter "a", and seed letter "b", the matching with "matrix[a][b]"
   *  @param nomerge if final states dont have to be merged together into a single state (default behaviour is merging)
   *  @return 0
   */

  int Automaton_SeedPrefixesMatching (const seed & s,
                                      const vector< vector <int> > & matchingmatrix,
                                      const bool nomerge=false);

  /** @brief Subset seed matching automaton for multiple subset seeds, but here with the Donald Martin extension
   *         that measures the (non-ovelapping) number of '#' elements that cover a '1' position : extended so
   *          weight_seed_alphabet
   *
   *  @param s is the set of seeds descriptos that are used to build the automaton
   *  @param weight_seed_alphabet is the weight given for each seed symbol : when several symbols are under a position,
   *         then the "winner" of the overlap is the one with the biggest weight
   *  @param matchingmatrix is a boolean matrix that gives for alignment letter "a", and seed letter "b", the matching with "matrix[a][b]"
   *  @return 0
   *  @see Automaton_SeedPrefixesMatching
   */

  int Automaton_SeedPrefixesMatching_CoverageDM (const vector<seed *> & s,
                                                 const vector<int> weight_seed_alphabet,
                                                 const vector< vector <int> > & matchingmatrix);

  /** @brief Modified version that takes lossless costs into account
   *  @param s is the seed descriptor
   *  @param matchingmatrix is a boolean  matrix that gives for alignment letter "a", and seed letter "b", the matching with "matrix[a][b]"
   *  @param nomerge if final states dont have to be merged together into a single state (default behaviour is merging)
   *  @param costs gives a vector of costs for each align alphabet letter
   *  @param cost_threshold gives a cost threshold that must not be reached by any alignment
   *  @return 0
   *  @see Automaton_SeedPrefixesMatching
   */

  int Automaton_SeedPrefixesMatchingCost (const seed& s,
                                          const vector< vector <int> > & matchingmatrix,
                                          const bool nomerge,
                                          const vector<int> & costs,
                                          const int cost_threshold);


  /** @brief Aho-Corasick automaton used by J.Buhler
   *  @param s is the seed descriptor
   *  @param matchingmatrix is a boolean matrix that gives for alignment letter "a", and seed letter "b", the matching with "matrix[a][b]"
   *  @param nomerge if final states dont have to be merged together into a single state (default behaviour is merging)
   *  @return 0
   */
  int Automaton_SeedBuhler (const seed & s,
                            const vector< vector <int> > & matchingmatrix,
                            const bool nomerge=false);

  /** @brief Aho-Corasick automaton with scoring method to prune prefixes
   *  @param s is the seed descriptor
   *  @param matchingmatrix is a boolean  matrix that gives for alignment letter "a", and seed letter "b", the matching with "matrix[a][b]"
   *  @param scoringmatrix is an integer matrix that gives for alignment letter "a", and seed letter "b", the score    with "matrix[a][b]"
   *  @param scoringthreehold is the minimal score that has to be reached to have a match for the seed
   *  @param nomerge if final states dont have to be merged together into a single state (default behaviour is merging)
   *  @return 0
   */

  int Automaton_SeedScore    (const seed & s,
                              const vector< vector <int> > & matchingmatrix,
                              const vector< vector<int> >  & scoringmatrix,
                              const int scoringthreehold,
                              const bool nomerge=false);


  /** @brief Modified version that takes lossless costs into account
   *  @param s is the seed descriptor
   *  @param matchingmatrix is a boolean matrix that gives for alignment letter "a", and seed letter "b", the matching with "matrix[a][b]"
   *  @param scoringmatrix is an integer matrix that gives for alignment letter "a", and seed letter "b", the score    with "matrix[a][b]"
   *  @param scoringthreehold is the minimal score that has to be reached to have a match for the seed
   *  @param nomerge if final states dont have to be merged together into a single state (default behaviour is merging)
   *  @param costs gives a vector of costs for each align alphabet letter
   *  @param cost_threshold gives a cost threshold that must not be reached by any alignment
   *  @return 0
   *  @see Automaton_SeedScore
   */

  int Automaton_SeedScoreCost(const seed & s,
                              const vector< vector <int> > & matchingmatrix,
                              const vector< vector<int> > & scoringmatrix,
                              const int scoringthreehold,
                              const bool nomerge,
                              const vector<int> & costs,
                              const int cost_threshold);


  /** @brief build a probabilistic bernoulli model
   *  @param p gives the probability of each letter of A (sum must be equal to 1.0)
   *  @return 0
   */

  int Automaton_Bernoulli(const vector<double> & p);


  /** @brief build a probabilistic markov model of order k
   *  @param p gives the probability of each word of A^k (sum must be equal to 1.0)
   *  @param k gives the order
   *  @return 0
   */

  int Automaton_Markov(const vector<double> & p,
                       const int k);

  /** @brief build an Homogeneous sequence automaton
   *  @details It represents an alignment such that no substring of the alignment has a score greater than the full alignment
   *  @param scores gives a table of scores for each of the A letters
   *  @param length gives the alignment length
   *  @return 0
   */

  int Automaton_Homogeneous(const vector<int> & scores,
                            const int length);


  /** @brief build an Counting symbols automaton (automaton with two states)
   *  @details It gives a 1-final state if "a >= min_a", or the init state (non-final) otherwise
   *  @param min_a minimal symbol that trigger a count
   *  @return 0
   */
  int Automaton_CountAlphabetSymbols(const int min_a = gv_align_alphabet_size - 1);

  /** @brief build a cyclic automaton
   *  @param cycle gives the cycle size
   *  @param final_list gives a table of positions that are final ( 0 <= final < cycle)
   *  @param final_nb   gives the size of the previous table
   *  @return 0
   */

  int Automaton_Cycle(const int cycle, const int * const final_list, const int final_nb);

  /** @brief build a lossless automaton
   *  @param costs gives a table of costs for each of the A letters
   *  @param max_cost is the maximal cost that can be reached so that the automaton is not in an accepting state
   *  @return 0
   */

  int Automaton_Lossless(const vector<int> & costs, const int max_cost);


  // @}





  /** @name  Manipulate automata
   *  @brief reduction/product/and misc probabilities
   */
  // @{


  /** @brief Hopcroft automaton minimization
   *  @return the minimized automaton, does not affect the current automaton
   *  @todo{FIXME : set a compatible algorithm for "final" integer states (init at : final=0/1/2 and cross ?? )}
   */

  automaton * Hopcroft() const;

  /** @brief Isomorphism check
   *  @param  other is the automaton to be compared to
   *  @return true if the two automata are isomorph
   */

  bool       isIsomorphTo(const automaton & other) const;



  /** @brief Generic Automata Product
   *  @param other is the second automaton used for the product
   *  @param productSetFinalType indicates if product final states are the crossproduct of both automaton final states or only one of these
   *       @li PRODUCT_UNION_xxx           : automata "union" of final states,
   *       @li PRODUCT_INTERSECTION_xxx    : automata "intersection" of final states
   *       @li PRODUCT_BUTNOT_xxx          : automata "this".final BUT NOT "other".final
   *       @li PRODUCT_NOTBUT_xxx          : automata NOT "this".final BUT "other".final
   *       @li PRODUCT_UNION_NO_LOOP_ADD   : automata "union" and "sum value" of final states, (integer value)
   *       @li PRODUCT_ADDHOC_xxx          : automata addhoc final function aff does this work
   *          with
   *       @li xxx : LOOP / NO_LOOP        : indicates if the final state is a single one absorbant (force it, otherwise keep it as in the "true" product)
   *  @param thisOrOtherIsProbabilist indicates that either one or both of the automata represents a probabilistic model (false by default)
   *       @li PRODUCT_NONE_IS_PROBABILIST  : no affectation for probabilities
   *       @li PRODUCT_THIS_IS_PROBABILIST  : only "this" probability is taken for each transition
   *       @li PRODUCT_OTHER_IS_PROBABILIST : only "other" probability is taken for each transition
   *       @li PRODUCT_BOTH_ARE_PROBABILIST : product of both "this" x "other" probability is taken for each transtion (sum is normalized to 1.0 per state)
   *  @param depth indicates the maximal depth that must be reached : extra states are non final selflooping states
   *         (by default, this value is greater than 2 Billions, but the given alignment length should be enought in most cases)
   *  @param aff function indicates the final value to be used, is is used when @param productSetFinalType = PRODUCT_ADDHOC_xxx
   *  @param shift is a positive or negative value that reset the first or second automaton on its initial state during a defined depth given by this parameter
   *  @return a new automaton that only gives reachable states of the product
   *  @see matrix_product
   */

  automaton * product(const automaton & other,
                      const ProductSetFinalType productSetFinalType,
                      const ProductProbabilityType thisOrOtherIsProbabilist = PRODUCT_NONE_IS_PROBABILIST,
                      const int depth = INT_INFINITY,
                      const AddHoc_Final_Func aff = NULL,
                      const int shift = 0) const;


  /** @brief Compute the multiple hit automaton
   *  @param m is the number of hits needed to hit an alignment
   *  @param depth indicates the maximal depth that must be reached
   *  @return the automaton that take at least @f$m@f$ hits of the
   *          original automaton to hit the alignment
   */

  automaton * mhit(const unsigned int m,
                   const int depth) const;



  /** @name  Matrix/Matrices Products
   *  @brief matrices probabilities/cost to speed up computation (no alphabet, step product functions)
   */

#ifndef NOMATRIX

  /** @addtogroup matrix */
  // @{
#include "automaton.hh"

  // @{
  /// probablistic specialized product
  matrix<double> * matrix_pr_product(const automaton & other,
                                     const ProductSetFinalType productSetFinalType,
                                     const int depth,
                                     const AddHoc_Final_Func aff = NULL) const {
    return matrix_product<double>(other,productSetFinalType,depth,aff);
  };

  /// cost specialized product
  matrix<cost<int> > * matrix_cost_product(const automaton & other,
                                           const ProductSetFinalType productSetFinalType,
                                           const int depth,
                                           const AddHoc_Final_Func aff = NULL) const {
    return matrix_product<cost<int> >(other,productSetFinalType,depth,aff);
  };

  /// count specialized product
  matrix<long long> * matrix_count_product(const automaton & other,
                                           const ProductSetFinalType productSetFinalType,
                                           const int depth,
                                           const AddHoc_Final_Func aff = NULL) const {
    return matrix_product<long long>(other,productSetFinalType,depth,aff);
  };

#ifdef USEINFINT
  /// count specialized product
  matrix<infint<long long int> > * matrix_count_infint_product(const automaton & other,
                                                  const ProductSetFinalType productSetFinalType,
                                                  const int depth,
                                                  const AddHoc_Final_Func aff = NULL) const {
    return matrix_product<infint<long long int> >(other,productSetFinalType,depth,aff);
  };
#endif


  // }@

  // @{
  /// probablistic specialized step product
  vector< matrix<double> * > * matrices_step_pr_product(const automaton & other,
                                                        const ProductSetFinalType productSetFinalType,
                                                        const int depth,
                                                        const AddHoc_Final_Func aff = NULL) const {
    return matrices_step_product<double>(other,productSetFinalType,depth,aff);
  };

  /// cost specialized step product
  vector< matrix<cost<int> > * > * matrices_step_cost_product(const automaton & other,
                                                              const ProductSetFinalType productSetFinalType,
                                                              const int depth,
                                                              const AddHoc_Final_Func aff = NULL) const {
    return matrices_step_product<cost<int> >(other,productSetFinalType,depth,aff);
  };

  /// count specialized step product
  vector< matrix<long long> * > * matrices_step_count_product(const automaton & other,
                                                              const ProductSetFinalType productSetFinalType,
                                                              const int depth,
                                                              const AddHoc_Final_Func aff = NULL) const {
    return matrices_step_product<long long>(other,productSetFinalType,depth,aff);
  };
  // @}
  // @}
#endif

  /** @brief Gives the transition probability
   *  @param a is the transition letter
   *  @param startingState is the starting state
   *  @param endingState is the ending state
   *  @return either 0 if there is no link, otherwise > 0
   */

  double    Prob(const int a, const int startingState, const int endingState) const;


  /** @brief Compute probability to be at a final state during the "nbSteps"th step
   *  @param nbSteps is the number of iterations done on the automaton to compute the probability
   *  @return the probabilty to be at a final state at the "nbSteps"-th step
   */

  double    PrFinal(const int nbSteps = gv_alignment_length) const;


  /** @brief Gives the probability to hit an alignment of the lossless set.
   *  @param nbSteps is the number of iterations done on the automaton to compute the probability
   *  @param costs is the vector of costs used on alignment alphabets
   *  @param cost_threshold is the maximal cost allowed
   *  @return the probabilty to hit such alignemnt (1.0 if lossless)
   *  @todo{FIXME : check if it does correct job for overlaps of 2 different seeds}
   *  @todo{FIXME : extends the algorithm for several hits of one/several seeds}
   *  @see Lossless
   */

  double PrLossless( const int nbSteps, const vector<int> & costs, const int cost_threshold) const;


  /** @brief Gives the lossless property of a seed (set of seeds)
   *  @param nbSteps is the number of iterations done on the automaton to compute the probability
   *  @param costs is the vector of costs used on alignment alphabets
   *  @param cost_threshold is the maximal cost allowed
   *  @param Nocc is the minimal number of seed occurences
   *  @return a boolean that indicates if the set of seeds is lossless (at least Nocc occurences everytime),
   *  @todo{FIXME : check if it does correct job for overlaps of 2 different seeds}
   *  @todo{FIXME : extends the algorithm for several hits of one/several seeds}
   *  @see PrLossless
   */

  bool  Lossless(const int nbSteps, const vector<int> & costs, const int cost_threshold, const int Nocc = 1) const;


  // @}




  /**@name  IO methods
   * @brief IO streams to print/store/load automaton
   */

  // @{
  /// print the graphviz/dot form of the automaton
  void dot(ostream& os) const;
  /// print the Mealy form of the automaton (gap - FR package), see http://www.gap-system.org/Manuals/pkg/fr-2.1.1/doc/chap0.html
  void gapFR(ostream& os) const;
  /// print the Moore form of the automaton (gap - Automata package), see http://www.gap-system.org/Manuals/pkg/automata/doc/chap0.html
  void gapAutomata(ostream& os) const;
  /// print automaton information
  friend ostream& operator<<(ostream& os, const automaton& automaton);
  /// read  automaton information
  friend istream& operator>>(istream& is, automaton& automaton);
  // @}

  /** @name  Miscellaneous
   *  @brief miscellaneous methods
   */

  // @{
  ///  return the size (number of states) of the current automaton
  inline int size() const { return _states.size(); }
  ///  make the state "stateNb" loop on itselft for any letter
  inline void selfLoop(const int stateNb);
  ///  generate an alignment sequence of given length from the current probabilistic automaton
  void GenSeq(const int len) const;
  void GenSeq(vector<int> & alignment) const;
  // @}


 protected:
  /** @name automaton manipulation routines
   */
  // @{

  /** @brief add a new state
   *  @param final indicates if this is a final state
   *  @return the index where the state has been added (= size() - 1)
   */

  inline int addNewState(const int final = 0);


  /** @brief add a new transition between two states
   *  @param a is the transition letter read
   *  @param startingState is the initial transition state number
   *  @param endingState is the final transition state number
   *  @param prob is the probability assigned to the current transition
   *  @return 0
   */

  inline void addNewTransition(const int a , const int startingState , const int endingState, const double prob = (1.00/gv_align_alphabet_size));


  /** @brief change the transition probability
   *  @param "a" is the transition letter
   *  @param startingState is the starting state
   *  @param endingState is the ending state
   *  @param prob is the new probability
   *  @return 0
   */

  inline int changeTransitionProb(const int a, const int startingState, const int endingState, const double prob);


  /** @brief check if there is a transition from the startingState labeled with "a"
   *  @param a is the transition letter
   *  @param startingState is the starting state
   *  @return true if there is at least one transition
   */

  inline bool hasTransition(const int a, const int startingState) const;


  // @}

  /// vector of states currently assigned
  vector<state> _states;
  /// vector of init states for each cycle (if empty then always state 1 is the init state)
  vector<int>   _init_states;
};
#endif


// @}
