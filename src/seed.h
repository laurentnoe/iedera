/** @page seeds Seeds
 *  @brief Seeds description and functions
 *  @tableofcontents
 *
 *  @section seed-description Description
 *  This part describes a seed element : each single seed element is represented by a @ref seed : each seed is thus
 *  represented first by a table of integer (@ref seed::_seedMotif_int) of given length (@ref seed::_span).
 *
 *  @see  seed::table()
 *  @see  seed::span()
 *
 *  Each value within this table must be a seed alphabet letter, and as thus, in the range
 *  [0.. @ref gv_seed_alphabet_size - 1].
 *
 *  It is possible to add some constraints on the positions where the seed will be applied (by default, it is applied
 *  a any possible position it fits).
 *  If any, a table of cyclic positions (seed::_seedCyclePos_int) of given length (@ref seed::_seedNbCyclePos) is
 *  allocated and filled. Each value within this table must be in  the range [0.. @ref seed::_seedMaxCyclePos - 1]
 *
 *  @see  seed::cycled()
 *  @see  seed::pos()
 *  @see  seed::nbpos()
 *  @see  seed::maxpos()
 *
 *  @section seed-enumeration Enumeration
 *  Enumerating seeds of given span, weight, signature is the main tool proposed in the @ref seed.
 *  Signature is a fixed set of number, for each non-null seed elements that are needed for the seed.
 *  Note that null seed elements can be added to increase the span without changing the signature validity.
 *
 *  @subsection full-enum Full enumeration
 *
 *  This is the more time consuming enumeration option, as it fully enumerate the seeds of given span
 *  (between @ref gv_minspan and @ref gv_maxspan), weight (between @ref gv_minweight and @ref gv_maxweight) and possible signature
 *  (@ref gv_signature in activated with @ref gv_signature_flag) constraints.
 *  It also enumerate the set of positions where the seed is applied (when cyclic positions seed:_seedCyclePos_int
 *  is set)
 *
 *  @see  seed::next()
 *
 *  @subsection random-enum Random enumeration
 *  This is the least costly, but also the least efficient way of enumeration seeds of given span
 *  (@ref gv_minspan and @ref gv_maxspan), weight (@ref gv_minweight and @ref gv_maxweight) and signature
 *  (@ref gv_signature in activated with @ref gv_signature_flag) constraints.
 *  since it randomly draw some seeds with these constraints.
 *
 *  @see  seed::random()
 *
 *  @subsection hill-climbing-enum Hill-climbing enumeration
 *  This is a classical local search optimization : when activated, it enumerates all the possible swaps of two
 *  elements, all the one "jocker symbol" add/delete events, all the "one position" changes for positioned seeds
 *  within the current seed until one gives a better result : in this case, this seed is set; otherwise if none is
 *  better, the original seed is reset. Several variables are provided to keep track of local modifications of
 *  the seed "jocker" last edit (@ref seed::_edit_operation, seed::_edit_from, seed::_edit_from_symbol), the seed
 *  last swap (@ref seed::_swap_from, seed::_swap_to) and for positioned seeds, the last edited position
 *  (@ref seed::_edit_cycle_pos_index, seed::_edit_cycle_pos_old_value).
 *
 *  @see  seed::next_hmove()
 *  @see  seed::set_hmove();
 *  @see  seed::reset_hmove();
 *
 *
 */

#ifndef __SEED_H__
#define __SEED_H__

//STL
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <algorithm>
//STD
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
using namespace std;
//STR
#include "macro.h"


/** @defgroup seed seed class
 *  @brief seed representation together with the set of cycle
 *         positions where it should be applied
 */
// @{


/**
 * @class seed
 *
 * @brief represents a subset seed (can also handle a list of positions where the seed is meant to apply)
 */

class seed {
 public:

  /** @name  seed(int nbcyclespos, int maxcyclepos)
   *  @brief build a seed of given span but full of jokers
   *  @param nbcyclepos is the number of positions in a cycle where the seed can be applied (0 = all positions),
   *  @param maxcyclepos is the size of the cycle.
   */

  seed(const int nbcyclepos = 0, const int maxcyclepos = 1);


  /** @name  seed(string  str_desc)
   *  @brief build a seed according to a descriptor
   *  @param str_desc is a string on ['0'-'9''A'-'Z''a'-'z']* or [(#[0-9]*)*]
   *         as "90A" or "#9#0#10" in equivalent form.
   *         if the Bsymbols option is selected, then the seed must respect
   *         the new symbols definition (eg for "-spaced", "##-##-###" ...)
   *  @param update_gv_span_weight checks if the seed given as str_desc respects the global seed constraints,
   *         otherwise changes these constraints (with a warning for the user)
   *  @see gv_min_span,gv_max_span,,gv_min_weight,gv_max_weight,
   */

  seed(string str_desc, const bool update_gv_span_weight = true);


 /** @name  ~seed()
   *  @brief delete a seed and associed structures (cyclepos and seed motif)
   */
  ~seed();


  /** @name Get functions
   *  @brief these methods are provided to get direct information from the current seed
   */

  //@{

  /** @brief gives the integer representation of the seed
   *  @return an integer table of lenght span, each element must
   *          be in the range 0 .. |B|-1
   */

  int  *   table() const {return _seedMotif_int;};


  /** @brief gives the span of the seed
   *  @return the span of the seed
   */

  int      span() const  {return _span;};


  /** @brief gives the integer positions where the seed should match
   *  @return an integer table of lenght nbpos(), each element must
   *          be in the range 0..(maxpos()-1)
   */

  int *    pos() const   {return _seedCyclePos_int;};


  /** @brief gives the number of positions where the seed should match (on a cycle)
   *  @return the number of positions
   */

  int    nbpos() const   {return _seedNbCyclePos;};


  /** @brief gives the size of the cycle
   *  @return the size of the cycle
   */

  int   maxpos()    {return _seedMaxCyclePos;};


  /** @brief  check if the seed has cycle positions (otherwise, is set at any valid position)
   *  @return true is there is cycle information associated with the current seed
   */

  bool  cycled()    {return _seedCyclePos_int != NULL;};

  //@}

  /** @name Check functions
   *  @brief these methods are provided to compute properties from the current seed
   */

  //@{

  /** @brief check if the seed is symetric
   *  @return true or false
   */
  bool symetric();

  /** @brief check if the seed is equal to another one
   *  @param other is the seed to be compared with this
   *  @return a boolean (only compare shapes, not positions)
   *  @todo{FIXME : do the same on position restricted seeds}
   */
  bool equal(seed * other);

  /** @brief gives the current seed shape (a string
   *         in [0-9A-Za-z]* or (#[0-9]*)* with a possible position restricted part ":1,3,5/7"
   *  @return a string object that represents this seed.
   */
  string str();

  /** @brief compute the weight of a given seed
   *  @return the weight as a double value
   */
  double weight();

  /** @brief compute seed selectivity according to the weight (bernoulli model estimation)
   *  @return the selectivity as a double value
   */
  double selectivityFromWeight();

  /** @brief check if the current seed is acceptable.
   *  This implies a correct weight  (according to gloval variables gv_minweight and gv_maxweight),
   *  no jocker on the first and last symbol of the seed, and a symetric shape when requested.
   *  @return true when the current seed is acceptable
   */
  bool acceptable();

  //@}



  //@{

  /** @brief check if the signature is realisable according to the span : stop the program otherwise
   */
  void checksignature();

  /** @brief set the very first seed signature
   */
  void setsignature();

  //@}



  //@{

  /** @brief set all the elements of position >= first_pos that are <= b in the seed in the decreasing order
   *         this function is used to enumerate signatures and inner used.
   *  @param first_pos gives the first position to be sorted
   *  @param b gives the maximal value that can be sorted (all others remain unchanged)
   *  @see next()
   */
  void reorder(int first_pos,int b);


  /** @brief this function does an enumerative generation of the acceptable seeds (and cycle positions when required),
   *         according to @ref gv_minspan and @ref gv_maxspan, @ref gv_minweight and @ref gv_maxweight
   *         with or without signature constraints
   *  @return 1 if next seed has been generated, 0 if all the seeds in the set have already been generated
   *  @see acceptable()
   *  @see setsignature(),checksignature()
   *  @todo{FIXME : avoid the possible infinity loop of this function}
   */
  int next();


  /** @brief this function generates a random acceptable seed (and random cycle positions when required)
   *         according to @ref gv_minspan and @ref gv_maxspan, @ref gv_minweight and @ref gv_maxweight
   *         with or without signature constraint
   *
   *  @return 1 all the time (but may stop the program if signature contradicts weight)
   *  @see acceptable()
   *  @see setsignature(),checksignature()
   */
  int random();

  //@}


  //@{

  /** @brief does a local move for the hill climbing function (one of the three posibiities, either :
   * - edit one of the cycle positions (if any),
   * - swap a pair of (different type) seeds elements,
   * - insert one jocker anywhere, or delete a jocker from the seed (if any).
   *  @return 1 if the next local move has been generated, 0 if all possible local moves have already been enumerated
   */
  int next_hmove();

  /** @brief set the last local move done previously for the hill climbing function, and keep its last pattern
   *  (all local will then be processed from this last pattern)
   */
  void set_hmove();

  /** @brief reset definitely the last local move done previously for the hill climbing function, and go back to the
   *  initial pattern. (all local will then be processed from this initial pattern, but it make nonsense doing that :
   *  better starting with another seed ...)
   */
  void reset_hmove();

  //@}

  /** @brief set cycle properties of the seed
   *  @param cycle_pos is a vector of positions where the seed is allowed to match
   *  @param cycle_size is the maximal position where a position can be set
   */
  void setCyclePos(vector<int> cycle_pos, int cycle_size);

  //@{

  /** @brief return the first position of a hit for the current (possibly positioned) subset seed, -1 otherwise
   *  @param alignment is the alignment to be checked
   *  @param matchingmatrix is a boolean matrix that gives for alignment letter "a", and seed letter "b", the matching with "matrix[a][b]"
   *  @return the first position of a hit for the current (possibly positioned) subset seed, -1 otherwise
   */
  int Hit(const vector<int> & alignment, const vector< vector <int> > & matchingmatrix);


  /** @brief return the number of hits for the current (possibly positioned) subset seed
   *  @param alignment is the alignment to be checked
   *  @param matchingmatrix is a boolean matrix that gives for alignment letter "a", and seed letter "b", the matching with "matrix[a][b]"
   *  @return the number of hits for the current (possibly positioned) subset seed
   */
  int mHits(const vector<int> & alignment, const vector< vector <int> > & matchingmatrix);

  //@}


  /** @brief output a seed object
   *  @param os is the outputstream
   *  @param s is the seed to output
   *  @return the outputsteam
   *  @see seed
   *  @todo{FIXME : implement also the input operator for a seed (set of seeds ?) with cycle reading and strong parsing}
   */

  friend ostream& operator<<(ostream& os, seed& s);


 protected:
  /// current seed motif represented as an int array {0,1,2, ...} of lenght "_span"
  int  * _seedMotif_int;
  /// fixed seed span (unless when an edit operation is ongoing)
  int _span;


  /// swap_and_edit_shift (random but fixed value)
  int _swap_and_edit_shift;
  /// swapped char cursor (from)
  int _swap_from;
  /// swapped char cursor (to)
  int _swap_to;

  /// edit the cycle pos index (-1 if no edition)
  int _edit_cycle_pos_index;
  /// edit the cycle pos old value (0 if no edition)
  int _edit_cycle_pos_old_value;

  /// edit operation : (0 if no edit yet, 1 for insertion, 2 for deletion, 3 if all edit operations have been done)
  int _edit_operation;
  /// edit position : where the symbol has been added removed
  int _edit_from;
  /// edit symbol : if deletion : which symbol has been removed ; if insertion : which symbol is tried to be added.
  int _edit_from_symbol;

  /// current seed positions to be set (NULL if applied on all positions)
  int * _seedCyclePos_int;
  /// size of the "_seedCyclePos_int" table
  int _seedNbCyclePos;
  /// max value in the  "_seedCyclePos_int" table
  int _seedMaxCyclePos;

  /**
   * @brief edit the cycle positions of the seeds
   * @return 1 if the next edit has been generated, 0 if all possible edits have already been enumerated
   */
  int next_cycleposedit();

  /**
   * @brief swap 2 different seed elements
   * @return 1 if the next swap has been generated, 0 if all possible swaps have already been enumerated
   */
  int next_swap();

  /**
   * @brief delete one jocker
   * @return 1 if the next deletion has been generated, 0 if all possible deletions have already been enumerated
   */
  int next_delete();

  /**
   * @brief insert one jocker
   * @return 1 if the next insertion has been generated, 0 if all possible insertions have already been enumerated
   */
  int next_insert();

  /**
   * @brief edit the seed by adding/removing any jocker like symbol
   * @return 1 if the next edit has been generated, 0 if all possible edits have already been enumerated
   * @see next_delete,next_insert
   */
  int next_edit();

  /**
   * @brief convert a char into a number
   *        '0' -> 0, '1' -> 1 ... '9' -> 9, 'A' -> 10 ... 'Z' -> 35, 'a' -> 36 ...
   * @param s is the string where conversion is currently processed
   * @param pos is the first char inside s where the conversion has to be done
   * @return the converted value (integer that represents a seed letter) and increments i to the next symbol
   */
  int convert(string & s, int * pos);

  /**
   *  @brief convert an integer into a string
   *        0 -> '0', 1 -> '1' ... 9 -> '9', 10 -> 'A' ... 35 -> 'Z', 36 -> 'a' ...
   *  @param i is the integer currently processed
   *  @return the converted value (string that represents a seed letter)
   */
  string rconvert(int i);

};

#endif

// @}
