#ifndef __MACRO_H__
#define __MACRO_H__
#include <iostream>
#include <string>
#include <vector>

using namespace std;


/** @addtogroup main */
// @{


/** @name compiling options
 *  @brief used to activate/unactivate code in iedera
 */
// @{

/* debug purposes */
/// some assertions that may be not checked to speedup the process
#define ASSERTB
#undef ASSERTB
/// some build debug
#define BUILD
#undef BUILD

/* map (or direct index) table */
/// use a map for the product and not an array (better for sparse results)
#define USEMAPPRODUCT


/// use a sampling on some alignments before computing the sensitivity (better when many seeds are designed together)
#define SAMPLING_STRATEGY_MF
#undef SAMPLING_STRATEGY_MF
/// keep product of seeds that are less subject to change in future (usefull when many seeds are designed together: faster for hillclimbing and full enumeration, but requires more memory)
#define KEEP_PRODUCT_MF
/// compute correlation coefficients for coverage or multi-hit automaton
#define COMPUTE_CORRELATION_LN
#undef COMPUTE_CORRELATION_LN
/// use the InfInt Header
#define USEINFINT
#undef USEINFINT

/// dont use the matrix slicer (slower than the slicer version)
#define NOSLICER
#undef NOSLICER
/// print matrix slices statitics (slow, for debugging purposes)
#define MATRIX_SLICER_STATS
#undef MATRIX_SLICER_STATS
/// compute lossless probabilities even for non hillclimbing
#define LOSSLESS_PROB
/// set the matrix sparsity level to be represended as a sparse row matrix
#define MATRIX_SPARSE_ROW_DENSITY 0.2
/// disable full matrix block computation (faster for i,j independent increases, but seems slower for i,j window increases, and keep the same algorithm proposed in the paper)
#define NO_FULL_MATRIX_BLOCK
/// use a queue for the product (otherwise it is a stack that disables some optimizations)
#define USEQUEUEPRODUCT

// @}


/** @name subalignments computation
 *  @brief parameters used for the sub-alignment computation and selection of the sensitivity/lossless value
 */
// @{
/// flag that enables the "subalignment" computation
extern bool                   gv_subalignment_flag;
/// length of the "subalignment"
extern int                    gv_subalignment_length;

/// index of the function used in @ref gv_subalignment_functions_names
extern int                    gv_subalignment_function_index;
/// number of functions proposed in the "subalignment" computation
#define SUBALIGNMENT_FUNCTIONS_NUMBER 4
/// type of the int functions proposed in the "subalignment" computation
typedef double (*subalignment_function_int)(vector<int> &);
/// vector of int functions proposed in the "subalignment" computation, see @ref gv_subalignment_functions_names for their names
extern subalignment_function_int    gv_subalignment_functions_int[SUBALIGNMENT_FUNCTIONS_NUMBER];
/// type of the double functions proposed in the "subalignment" computation
typedef double (*subalignment_function_double)(vector<double> &);
/// vector of double functions proposed in the "subalignment" computation, see @ref gv_subalignment_functions_names for their names
extern subalignment_function_double gv_subalignment_functions_double[SUBALIGNMENT_FUNCTIONS_NUMBER];
/// names of the functions proposed in the "subalignment" computation, see @ref gv_subalignment_functions_int and @ref gv_subalignment_functions_double for their respective implementation
extern char *                 gv_subalignment_functions_names[SUBALIGNMENT_FUNCTIONS_NUMBER];
// @}

/** @name verbosity
 *  @brief verbosity level and macros
 */
// @{
/// vebosity mode
extern int gv_verbose;

/// messages (code adapted from the Marta SToRM read mapper :-) )
/// colored messages (helpfull)
#define INFO__(message)     {cerr << "\033[32;1m" << message; cerr << "\033[0m" << endl;}
#define DEBUG__(message)    {cerr << "\033[35;1m" << message; cerr << "\033[0m" << endl;}
#define MESSAGE__(message)  {cerr << "\033[0m"    << message; cerr << "\033[0m" << endl;}
/// verbosity levels
#define VERBOSITY_NONE      0
#define VERBOSITY_LIGHT     1
#define VERBOSITY_MODERATE  2
#define VERBOSITY_HIGH      3
#define VERBOSITY_ANNOYING  4
#define VERBOSITY_DEBUGGING 5
#define VERBOSITY_MAX      VERBOSITY_DEBUGGING
#define VERBOSITY_DEFAULT  VERBOSITY_MODERATE
/// verbosity macro filter
#define VERB_FILTER(verbosity_level, f) if ((verbosity_level) <= gv_verbose) { f; }
// @}

/** @name seed enumeration
 *  @brief parameters used to enumerate seeds
 */
// @{

/// number of seed runs (0 when full enumeration of all possible seeds)
extern int                    gv_nbruns;
/// signature vector (number of seed letters per seed)
extern vector<int>            gv_signature;
/// signature flag that enable the signature selection @see gv_signature
extern bool                   gv_signature_flag;
/// signature shuffle flag that enable only shuffle from the -m pattern and is not compatible with signature @see gv_signature_flag
extern bool                   gv_signature_shuffle_from_m_pattern_flag;
/// min seed span
extern int                    gv_minspan;
/// max seed span
extern int                    gv_maxspan;

/// min seed weight
extern double                 gv_minweight;
/// max seed weight
extern double                 gv_maxweight;

/// flag that enables min and max seed weight @see gv_minweight,gv_maxweight
extern bool                   gv_weight_interval_flag;
/// flag that enables vectorized seeds @see gv_vectorizedsubsetseed_scoring_matrix,gv_vectorizedsubsetseed_scoring_threshold
extern bool                   gv_vectorized_flag;
/// set of symbols when "overdubbed"  @see gv_bsymbols_flag
extern char *                 gv_bsymbols_array;
/// flag that enables "overdubbed" symbols @see gv_bsymbols_array
extern bool                   gv_bsymbols_flag;

/// hmove choice is set randomly to help seed optimisation choose a way to generate candidate seeds
extern int                    gv_hmove_choice;
/// flag to select symetric seeds only
extern bool                   gv_symetric;
// @}

/** @name alphabet characteritics
 */
// @{
/// flag that enables the "overdubbed" symbols
extern bool                   gv_matching_symbol_flag;
/// length of the alignment on which seeds are optimized
extern int                    gv_alignment_length;
/// alignment alphabet size
extern int                    gv_align_alphabet_size;
/// seed alphabet size
extern int                    gv_seed_alphabet_size;
/// min seed letter selectivity for the seed alphabet
extern double                 gv_bsel_minprob;
/// max seed letter selectivity for the seed alphabet
extern double                 gv_bsel_maxprob;
/// weight as @f$\frac{ \log(\mathrm{seed\ letter\ selectivity})}{log (max (\mathrm{seed\ letter\ selectivity}))}@f$  for each seed  letter
extern vector<double>         gv_bsel_weight;
/// computer the weight (gv_bsel_weight) according to background probabilities (gv_bsel)
extern void computeWeight();
// @}


/** @name lossless seeds
 */
// @{
/// flag that activates lossless computation (default lossy computation)
extern bool                   gv_lossless_flag;
/// the lossless cost for A when lossless flag is activated @see gv_lossless_flag,gv_lossless_cost_threshold
extern vector<int>            gv_lossless_costs_vector;
/// lossless cost for an alignment (must be at most equal to this value) @see gv_lossless_costs_vector
extern int                    gv_lossless_cost_threshold;
// @}

/** @name helpfull functions and macros
 *  @brief classical macros and error/warning reporting
 */
/// ERROR macro
#define _ERROR(message1,message2) {                                               \
  cerr << "\033[31;1m*ERROR :\033[0m " << message1 << " : " << message2 << endl; \
  exit(-1);                                                                       \
}

/// WARNING macro
#define _WARNING(message1,message2) {                                               \
  cerr << "\033[33;1m*WARNING :\033[0m " << message1 << " : " << message2 << endl; \
}

/// MAX macro
#define MAX(a,b) ((a)>(b)?(a):(b))
/// MIN macro
#define MIN(a,b) ((a)<(b)?(a):(b))
/// INFINITY macro for signed int
#define INT_INFINITY   (0x7fffffff)

/** qpow functions
 */
// @{

/**
 *@brief quick power on integer base
 *@param x is the base of the power
 *@param n is the exponent of the power
 *@return \f$ x^n \f$
 */
int power(int x, int n);

/**
 *@brief quick power on double base
 *@param x is the base of the power
 *@param n is the exponent of the power
 *@return \f$ x^n \f$
 */
double power(double x, int n);

// @}

// @}
#endif

