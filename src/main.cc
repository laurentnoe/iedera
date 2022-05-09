/**
 * @mainpage iedera
 * @version 1.07
 * @section Introduction
 *   Iedera is a program to select and design subset seeds and vectorized
 *   subset seeds. This program heavily relies on automaton theory to compute the
 *   sensitivity of seed that are represented by the subset seed model, or the
 *   vectorized subset seed model.
 *
 *
 *
 * @section Download
 * Download the iedera source code and binaries on the following web page:\n
 * <A HREF="http://bioinfo.cristal.univ-lille.fr/yass/iedera.php">http://bioinfo.cristal.univ-lille.fr/yass/iedera.php</A>\n
 * A small tutorial to use iedera is also provided on this web page.\n
 *
 * @author <A HREF="http://www.cristal.univ-lille.fr/~noe" TARGET="_top">Laurent Noe</A>
 * @author <A HREF="http://www.mccme.ru/lifr/pers/roytberg.htm">Mikhail Roytberg</A>
 * @author <A HREF="http://www-igm.univ-mlv.fr/~koutcher/" TARGET="_top">Gregory Kucherov</A>
 *
 *
 */

/** @defgroup main main program
 *  @brief main program, with pareto seed/sens/sel attributes, and the set of global variables (gv_...)
 *         used to parametrise the computation
 */
// @{


//C stuff
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <limits>
//STL
#include <vector>
#include <list>
#include <utility>
#include <algorithm>
#include <iterator>
//STD
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
using namespace std;
//STR
#include "macro.h"
#include "seed.h"
#include "automaton.hh"
#include "matrix.hh"


#ifdef USEINFINT
#include "infint.hh"
/// BIGINT polynomial coefs computation defined as Infinite Integer (much slower ... but overflow secure)
#define BIGINT infint<long long>
#else
/// BIGINT polynomial coefs computation defined as long long (quite fast ... but overflows quickly ...)
#define BIGINT long long
#endif


/** @name program name and version number
 *  @brief defined if not provided by the compiler
 */
// @{

#ifndef PACKAGE
/// iedera program name
#define PACKAGE "iedera"
#endif

#ifndef VERSION
/// iedera program version
#define VERSION "1.07"
#endif
// @}

/** @name subset matching matrix (uninitialized)
 *  @brief must be initialized in the main function
 */
// @{
/// boolean matching matrix
std::vector< std::vector<int> >  gv_subsetseed_matching_matrix;
/// build an initial matrix gv_subsetseed_matching_matrix @see gv_subsetseed_matching_matrix
void build_default_subsetseed_matching_matrix();
// @}

/** @name background and foreground probabilities (uninitialized)
 *  @brief must be initialized in the main function
 */
// @{
/// foreground sensitivities
std::vector< double >       gv_bsens;
/// k-order model when Markov is activated (0 for Bernoulli) @see gv_bsens_k
int                         gv_bsens_k = 0;
/// probabilistic automaton associated with the sensitivity computation @see gv_bsens
automaton <double> *        gv_bsens_automaton = NULL;
/// background seletivities  @see gv_bsel_k
std::vector<double>         gv_bsel;
/// k-order model when Markov is activated (0 for Bernoulli) @see gv_bsel
int                         gv_bsel_k = 0;
/// build probabilities @see gv_bsens,gv_bsens_k,gv_bsel,gv_bsel_k
void build_default_probabilities();
// @}

/** @name vectorized seeds (uninitialized)
 *  @brief generate subset seeds with an additional scoring constraint
 */
// @{
/// vector scoring matrix when vectorized subset seeds are activated @see gv_vectorized_flag,gv_vectorizedsubsetseed_scoring_threshold
std::vector< std::vector<int> >  gv_vectorizedsubsetseed_scoring_matrix;
/// scoring threshold for a vectorized subset seed @see gv_vectorizedsubsetseed_scoring_matrix
int                              gv_vectorizedsubsetseed_scoring_threshold = 0;
/// flag that activates vectorized subset seeds (default subset seed)
bool                             gv_vectorized_flag = false;
/// build an initial matrix gv_vectorizedsubsetseed_scoring_matrix @see gv_vectorizedsubsetseed_scoring_matrix
void build_default_vectorizedsubsetseed_scoring_matrix();
// @}


/** @name seed enumeration
 *  @brief parameters used to enumerate seeds
 */
// @{
/// Hopcroft minimization applied on each deterministic automaton
bool     gv_minimize_flag      = true;
/// hillclimbing local optimization for each seed and set of positions
bool     gv_hillclimbing_flag  = false;
/// hillclimbing local optimization frequency
double   gv_hillclimbing_alpha = 5e-3;
/// flag set when a seed is given on command line
bool     gv_motif_flag         = 0;
/// set of seeds being computed or optimized
std::vector< seed*>  gv_seeds;
/// maximal difference of weight (given as a ratio) between seeds being optimized together
double   gv_jive               = 0.0;
// @}

/** @name cycle parameters
 *  @brief generates position restricted seeds of a cycle
 */
// @{
/// cycle flag @see gv_cycles
bool                             gv_cycles_flag = false;
/// cycle size (one for each seed) @see gv_cycles_flag
std::vector< int >               gv_cycles;
/// cycle pos (number of pos per cycle) @see gv_cycles_pos_list,gv_cycles_flag
std::vector< int >               gv_cycles_pos_nb;
/// cycle pos (pos for each cycle) @see gv_cycles_flag,gv_cycles_pos_nb
std::vector< std::vector<int>  > gv_cycles_pos_list;
// @}

/* @name multihit criterion
 * @brief consider a hit if at least $m$-hits (of any of the seeds) have been counted
 */
// @{
/// flag that activates the mhits @see gv_multihit_nb
bool gv_multihit_flag = false;
/// number of mhits (must be >= 1 if the flag is activated) @see gv_multihit_flag
int  gv_multihit_nb   = 0;
// @}

/* @name DM coverage criterion
 * @brief consider a hit if the DM (Donald Martin) coverage criterion is reached
 */
// @{
/// flag that activates the coverage criterion
bool gv_global_coverage_flag = false;
/// number of element to be covered (must be >= 1 if the flag is activated) @see gv_multihit_flag
int  gv_global_coverage_nb   = 0;
/// coverage cost for each seed alphabet element
std::vector<int> gv_global_coverage_cost;
// @}

/** @name excluded seeds
 *  @brief exclude seeds motif from alignments generated
 */
// @{
/// excluded seeds
std::vector<seed*>   gv_xseeds;
// @}

/** @name excluded seed cycle parameters
 *  @brief generates position restricted seeds of a cycle
 */
// @{
/// cycle flag @see gv_xseeds_cycles
bool                            gv_xseeds_cycles_flag = false;
/// cycle size (one for each seed) @see gv_xseeds_cycles_flag
std::vector<int>                gv_xseeds_cycles;
/// cycle pos (number of pos per cycle) @see gv_xseeds_cycles_pos_list,gv_xseeds_cycles_flag
std::vector<int>                gv_xseeds_cycles_pos_nb;
/// cycle pos (pos for each cycle) @see gv_xseeds_cycles_flag,gv_xseeds_cycles_pos_nb
std::vector< std::vector<int> > gv_xseeds_cycles_pos_list;
// @}

/* @name excluded seed multihit criterion
 * @brief consider a hit if at least $m$-hits (of any of the excluded seeds) have been counted
 */
// @{
/// flag that activates the mhits @see gv_xseeds_multihit_nb
bool gv_xseeds_multihit_flag = false;
/// number of mhits (must be >= 1 if the flag is activated) @see gv_xseeds_multihit_flag
int gv_xseeds_multihit_nb    = 0;
// @}



/** @name homogeneous alignments
 *  @brief generates homogeneous alignments according to scoring parameters
 */
// @{
/// flag that activates computation on homogeneous alignments
bool                gv_homogeneous_flag = false;
/// score for each letter of the alignment alphabet
std::vector<int>    gv_homogeneous_scores;
// @}



/** @name correlation computation
 *  @brief consider correlation computation on the set of alignments (and not sensitivity)
 */
// @{
/// flag that activates the correlation computation
bool gv_correlation_flag           = false;
/// function choosen for the correlation computation @see gv_correlation_flag
int  gv_correlation_function_index = 0;
/// function numbers for the correlation computation @see gv_correlation_function_index
#define CORRELATION_FUNCTIONS_ENUM_PEARSON  0
#define CORRELATION_FUNCTIONS_ENUM_SPEARMAN 1
#define CORRELATION_FUNCTIONS_ENUM_NUMBERS  2
/// function symbols for the correlation computation
char * gv_correlation_functions_names[CORRELATION_FUNCTIONS_ENUM_NUMBERS] = {(char *)"PEARSON",(char *)"SPEARMAN"};
// @}


/** @name simple covariance optimization
 *  @brief consider computing the covariance as done in "Morgenstern, Zhu et al. AMB 2015"
 */
// @{
bool gv_covariance_flag            = false;
// @}




/** @name polynomial output
 *  @brief consider the coefficients on the set of alignments
 */
// @{
/// flag that activate the output of the "Mak Benson 2009" coefficients (number of alignments with @f$ x @f$ mismatches that are not detected, or detected, or with coverage/multihit value ...). Note that it can be used with "Chung Park 2010" definition of integral hits, but you need to replace the  @f$ p^{l-x} \cdot (1-p)^{x} @f$ factor preceeding each coefficient for each mismatch rate of @f$ \frac{x}{l} @f$, by its integral  @f$ \int_{p_a}^{p_b} p^{l-x} \cdot (1-p)^{x} dp = \Big[ p^{l-x+1} \sum_{j=0}^{x} {x \choose j} \frac{(-p)^j}{l-x+j} \Big]_{p_a}^{p_b} @f$. Fast and precise computation can be obtained with a fast power of the coefficients of the initial polynomial @f$ p^n (1-p)^x @f$ to avoid some small p "vs" big binomial problem.

/// flag that activate the selection of dominant seeds as in "Mak Benson 2009", which is also true for "Chung Park 2010"
/// (and not the selection according to their sensitivity only: its a more general form that may select more seeds)
bool gv_polynomial_dominant_selection_flag = false;



/// flag that activate multinomial evaluation (note : this is independent from the previous dominance selection). Note that this is only a polynomial evaluation, and does not have any other effect than ... FIXME
bool gv_multipoly_file_flag = false;
automaton<polynomial<BIGINT > > * gv_multipoly_bsens_automaton = NULL;



// @}

/** @name data management
 *  @brief pareto input/output files
 */
// @{
/// pareto output filename @see outputPareto
char   * gv_output_filename = NULL;
/// pareto input filenames @see inputPareto
char   * gv_input_filenames[256] = {NULL};
/// pareto input filenames (number of files) @see gv_input_filenames
int      gv_nb_input_filenames =   0;
/// pareto selection every nbruns (print statitics also)
int      gv_pareto_select_runs = 1000;
// @}




/** @name correlation coefficients computations
 *  @brief binomial coefficients and associated weight used by correlation coefficients computations
 */
// @{

/// binomial coefficients, must be pre-computed for correlations algorithms
std::vector< std::vector<BIGINT> > binomials;

/// binomials computed with the Pascal Triangle (overflow are checked and stop the program)
void generate_binomials(int N) {
  binomials = std::vector< std::vector<BIGINT> > (N+1,  std::vector<BIGINT>(N+1,0) );
  int k, n;
  for (k = 1; k <= N; k++) binomials[0][k] = 0;
  for (n = 0; n <= N; n++) binomials[n][0] = 1;
  for (n = 1; n <= N; n++) {
    for (k = 1; k <= N; k++) {
      binomials[n][k] = binomials[n-1][k-1] + binomials[n-1][k];
#ifndef USEINFINT
      if (binomials[n][k] < 0) {
        _ERROR("this binary has been compiled with undefined USEINFINT (no infinite precision integer)"," generate_binomials DOES OVERFLOW ...");
      }
#endif
    }
  }
}

/// tools used to reduce size of binomials weight : greatest common divisor of all of them
BIGINT gcd(BIGINT a, BIGINT b) {
  if (a < b) {
    BIGINT t = b; b = a; a = t;
  }
  while (b > 0) {
    BIGINT t = b; b = a % b; a = t;
  }
  return a;
}

/// tools used to reduce size of binomials weight : lowest common multiple of all of them
BIGINT lcm(BIGINT a, BIGINT b) {
  return a/gcd(a,b) * b;
}


/// binomial weight (inverted values of binomials coeeficients), must be pre-computed for correlations algorithms
#ifdef USEINFINT
std::vector< std::vector<BIGINT> > binomial_weights;
#else
std::vector< std::vector<double> > binomial_weights;
#endif

/// binomial weight flag, must be set to false when already computed
bool gp_binomial_weights_not_computed_flag = true;

void generate_binomials_weights(int N) {
#ifdef USEINFINT
  binomial_weights = std::vector< std::vector<BIGINT> > (N+1,  std::vector<BIGINT>(N+1,0) );
#else
  binomial_weights = std::vector< std::vector<double> > (N+1,  std::vector<double>(N+1,0) );
#endif

  generate_binomials(N);

  /*
   * lcm { binomial(k,0) ; binomial(k,1) ; ... ; binomial(k,k) } = lcm {1;2;...;k;k+1} / (k+1)
   * "An identity involving the least common multiple of binomial coefficients and its application" Bakir FARHI
   */
  BIGINT lcmc = 1;
  for (long long k_plus_1 = 1; k_plus_1 <= N+1 ; k_plus_1++) {
    lcmc = lcm(lcmc,k_plus_1);
#ifndef USEINFINT
    if (lcmc < 0) {
      _ERROR("this binary has been compiled with undefined USEINFINT (no infinite precision integer)"," generate_binomials_weights DOES OVERFLOW ...");
    }
#endif
    BIGINT lcmc_div_k_plus_1 = lcmc / (BIGINT) k_plus_1;
    long long k = k_plus_1 - 1;
    for (long long i = 0; i <= k ; i++) {
      BIGINT lcmc_div_k_plus_1_div_binomial_k_i = lcmc_div_k_plus_1/(binomials[k][i]);
#ifdef USEINFINT
      binomial_weights[k][i] =         lcmc_div_k_plus_1_div_binomial_k_i;
#else
      binomial_weights[k][i] = (double)lcmc_div_k_plus_1_div_binomial_k_i;
#endif
    }
  }
}
// @}




/** @name command line functions
 *  @brief set of command line management functions
 */
// @{

/// display a 2d matrix (a vector of vectors of int)
void DISPLAYMATRIX( std::ostream & os, vector< vector<int> > & matrix );

/// describe the command line parameters
void USAGE() {
  cerr << "* USAGE ("<< PACKAGE << " v" << VERSION << ") :" << endl;
  cerr << "  " << PACKAGE << " [options]" << endl;
  cerr << "      -h                      : display this Help screen" << endl;
  cerr << "      -v <int>                : set verbose mode (default = " << gv_verbose << ")" << endl;
  cerr << "                                                                                      " << endl;
  cerr << "   1) Seed Model :" << endl;
  cerr << "     a) Data " << endl;
  cerr << "      -A <int>                : align alphabet size \'A\' (default = " << gv_align_alphabet_size << ")" << endl;
  cerr << "      -B <int>                :  seed alphabet size \'B\' (default = " <<  gv_seed_alphabet_size << ")" << endl;
  cerr << "      -M \"{{<bool>}}\" table   : subset seed Matching matrix (default = \""; DISPLAYMATRIX(cerr, gv_subsetseed_matching_matrix); cerr << "\")" << endl;
  cerr << "         * note    : parameter is a C-like matrix of 'A'x'B' cols x rows " << endl;
  cerr << "         * example : -M \"{{1,0,1},{0,1,1}}\"    (for A=3 and B=2) " << endl;
  cerr << "      -MF <filename>          : gives the matrix inside a file (when too large for command line)" << endl;
  cerr << "                                                                               " << endl;
  cerr << "     b) Subset Seed / Vectorized Subset Seed " << endl;
  cerr << "      -V                      : select the Vectorized subset seed model (default = disabled)" << endl;
  cerr << "      -S \"{{<int>}}\" table    : vectorized subset seed Scoring table (default = \""; DISPLAYMATRIX(cerr, gv_vectorizedsubsetseed_scoring_matrix); cerr << "\")" << endl;
  cerr << "      -SF <filename>          : gives the table inside a file (when too large for command line)" << endl;
  cerr << "      -T <int>                : vectorized subset seed scoring minimal Threshold (default = " << gv_vectorizedsubsetseed_scoring_threshold << ")" << endl;
  cerr << "                                                                               " << endl;
  cerr << "     c) Subset Seed / Lossless Subset Seed " << endl;
  cerr << "      -L <int>,<int>,...      : select the Lossless mode and set and the lossless alignment costs" << endl;
  cerr << "      -X <int>                : set the alignment maX cost to be considered in the alignment set" << endl;
  cerr << "                                                                               " << endl;
  cerr << "   2) Alignments:" << endl;
  cerr << "      -b <dbl>,<dbl>,...      : set the Bernoulli/Markov background distribution (size |'A'|^k)" << endl;
  cerr << "      -f <dbl>,<dbl>,...      : set the Bernoulli/Markov foreground distribution (size |'A'|^k)" << endl;
  cerr << "      -fF <filename>          : set the foreground model (as a probabilistic NFA)" << endl;
  cerr << "                                for more details, see http://bioinfo.cristal.univ-lille.fr/yass/iedera.php#fFFormat" << endl;
  cerr << "      -l <int>                : set the alignment length (default = " <<  gv_alignment_length << ") " << endl;
  cerr << "      -ll <int>               : select the sub-align window computation and set its length (default = disabled)" << endl;
  cerr << "      -llf <function>         : set the sub-align function used to merge the results (default = \"min\")" << endl;
  cerr << "         * note    : available functions are ";
  for (int i=0; i<SUBALIGNMENT_FUNCTIONS_NUMBER; i++) {
    cerr << " \"" << gv_subalignment_functions_names[i] << "\"";
  }
  cerr << endl;
  cerr << "      -u <int>,<int>,...      : select the homogeneous alignment model and set the scores (default = disabled)" << endl;
  cerr << "                                                                               " << endl;
  cerr << "   3) Seed Selection / Enumeration :" << endl;
  cerr << "      -n <int>                : number of consecutive seeds to design (default = " << gv_seeds.size() << ")" << endl;
  cerr << "      -s <int>,<int>          : minimal/maximal seed span (default min = " << gv_minspan << ",max = " << gv_maxspan << ")" << endl;
  cerr << "      -w <dbl>,<dbl>          : minimal/maximal seed weight (default min = " << gv_minweight << ", max = " << gv_maxweight << ")" << endl;
  cerr << "         * symbol # (if defined, otherwise the max of all symbols) is of weight 1.0 " << endl;
  cerr << "         * symbol _ is of weight 0.0 " << endl;
  cerr << "         * symbol e weight is given by 'log_#(e matches background distribution)'" << endl;
  cerr << "      -i <int>,<int>,...      : signature id (number of '0','1',...'B-1' elements inside a subset seed)" << endl;
  cerr << "         * note    : the first number can be 0 to be adjusted on the fly depending of seed span." << endl;
  cerr << "         * example : for transitive seeds, -i 0,2,8 enumerate seeds with 8 \'#\', 2 \'@\', and any number of \'-\'." << endl;
  cerr << "         * note     : -i shuffle set after -m <patterns> shuffles each seed independently without changing the weight/span" << endl;
  cerr << "      -j <dbl>                : difference of weight (given as a ratio) allowed between consecutive seeds." << endl;
  cerr << "         * note    : between 0.0 for miminal change and 1.0 for maximal change (default = " << gv_jive << ")" << endl;
  cerr << "      -r <int>                : random enumeration (default = " << gv_nbruns << " is complete enumeration)" << endl;
  cerr << "         * note    : must be used with -k hill climbing for better performances. " << endl;
  cerr << "      -k                      : hill climbing heuristic (for random enumeration)." << endl;
  cerr << "         * note    : it checks all permutations of two different elements inside one seed," << endl;
  cerr << "                     it checks the edit (add/removal) of one jocker '-' inside the seed." << endl;
  cerr << "      -a <dbl>                : hill climbing heuristic frequency: good seeds are more likely to be optimized" << endl;
  cerr << "      -x                      : only symetric seeds are selected" << endl;
  cerr << "      -y <int>                : consider sensitivity when y hits are required (multihit sensitivity)" << endl;
  cerr << "      -g <int>                : consider sensitivity when g coverage is required (coverage sensitivity)" << endl;
  cerr << "         * note    : when a seed hits, its coverage is computed as the sum of its elements symbols per position" << endl;
  cerr << "                     and for several possibly overlapping seed hits, the \"best\" element per position is kept." << endl;
  cerr << "         * example : \"##-#\"  gives 5 on a spaced seed alphabet, \"#@-#\"  gives 9 (4x2+1) on a transitive alphabet," << endl;
  cerr << "                      \"##-#\" because 5'#' with '#'=1.            \"#@-#\" because 4'#' wt '#'=2, one '@' wt '@'=1."<< endl;
  cerr << "         * note    : only available with subset seeds (not with vectorized subset seeds)." << endl;
  cerr << "                                                                               " << endl;
  cerr << "      * -y/g note  : -y/-g can be used with SPEARMAN/PEARSON correlation with the alignment %%id (only for -A 2)." << endl;
  cerr << "                     you can select the minimal number of matches (%%id) by setting -y <int> before -y <SPEARMAN/PEARSON>" << endl;
  cerr << "         * example : \"-spaced -l 64 -y 32 -y PEARSON\" for a %%id of 50%% and PEARSON correlation computation" << endl;
  cerr << "      -p                      : activate \"Mak Benson 2009\" dominant selection and output polynomial." << endl;
  cerr << "         * note    : this is useful for iid models only (Bernoulli,...) " << endl;
  cerr << "         * example : \"-spaced -l 8 -m \"##-#\" -p\" must output the additional values :"<< endl;
  cerr << "                          0,0=1;1,0=8;2,0=28;3,0=51;4,0=45;5,0=15;6,0=1;"<< endl;
  cerr << "                          3,1=5;4,1=25;5,1=41;6,1=27;7,1=8;8,1=1;" << endl;
  cerr << "                      for the (matching) polynomial (second line: terms contain \",1\") where q ~ 1-p:" << endl;
  cerr << "                             3 8-3    4 8-4    5 8-5    6 8-6    7 8-7   8 8-8" << endl;
  cerr << "                          5.p.q + 25.p.q + 41.p.q + 27.p.q  + 8.p.q + 1.p.q" << endl;
  cerr << "                                                                               " << endl;
  cerr << "                      or complementary (non-matching) polynomial (first line: terms contain \",0\"):" << endl;
  cerr << "                             0 8-0   1 8-1    2 8-2    3 8-3    4 8-4    5 8-5   6 8-6" << endl;
  cerr << "                          1.p.q + 8.p.q + 28.p.q + 51.p.q + 45.p.q + 15.p.q + 1.p.q" << endl;
  cerr << "                                                                                " << endl;
  cerr << "      -pF <filename>          : activate general polynomial evaluation and load the associated file automaton." << endl;
  cerr << "         * example : \"-spaced -l 5 -m \"##-#\" -pF _file_\" where the _file_ is:" << endl;
  cerr << "                                                                                " << endl;
  cerr << "                             |3   0 1    0 0" << endl;
  cerr << "                             |           1 0" << endl;
  cerr << "                             |    1 0    0 1     1 x" << endl;
  cerr << "                             |           1 1     2 1 - x" << endl;
  cerr << "                             |    2 0    0 1     1 1 - y" << endl;
  cerr << "                             |           1 1     2 y" << endl;
  cerr << "                                                                                " << endl;
  cerr << "                     will give the following result (with spaces between each ^*+- symbols):" << endl;
  cerr << "                                                                                " << endl;
  cerr << "                             [y - x*y - x*y^2 + 2*x*y^3 - x^2*y + 2*x^2*y^2 - 2*x^2*y^3 + x^3*y - x^3*y^2]" << endl;
  cerr << "                                                                                " << endl;
  cerr << "      -c <int>,<int>,...      : consider sensitivity when each seed is indexed on 1/c of its positions" << endl;
  cerr << "         * note    : the position is enumerated or choosen randomly depending on -r parameter" << endl;
  cerr << "         * example : \"##-#:1/5\" means that the seed \"##-#\" will be placed at 1st,6th,11th... positions" << endl;
  cerr << "      -q <int>,<int>,...      : consider sensitivity when each seed is indexed on q/c of its positions" << endl;
  cerr << "         * note    : parameter -c must be specified before, positions are enumerated or randomly drawn (-r)" << endl;
  cerr << "         * example : \"##-#:2,5/5\" means that the seed \"##-#\" will be placed at 2nd,5th,7th,10th... positions" << endl;
  cerr << "      -d                      : disable Hopcroft minimization (default is enabled) " << endl;
  cerr << "      -m \"<seeds>\"            : check the given seed selectivity/sensitivity " << endl;
  cerr << "         * note    : possibility to add cycle content as shown before (see -c -q \"* example :\" )" << endl;
  cerr << "      -mx \"<seeds>\"           : generate alignments that exclude the given seed pattern" << endl;
  cerr << "         * note    : possibility to add cycle content as shown before (see -c -q \"* example :\" )" << endl;
  cerr << "      -mxy <int>              : as -y, consider -mx when y hits are required (multihit criterion)" << endl;
  cerr << "                                                                               " << endl;
  cerr << "   4) Input / Ouput :" << endl;
  cerr << "      -BSymbols '<char>'      : seed alphabet symbols to display (default is '01..9AB..Zab..z')" << endl;
  cerr << "      -e <filename>           : input an initial \"pareto set\" file before optimizing" << endl;
  cerr << "      -o <filename>           : output the \"pareto set\" to this file after optimizing" << endl;
  cerr << "      -z <int>                : print statistics every \"z\" seeds (default = " << gv_pareto_select_runs << ", 0 to disable)" <<  endl;
  cerr << "                                                                               " << endl;
  cerr << "   5) Nucleic Spaced Seeds : (shortcuts : you may overload the weight -w and span -s after)" << endl;
  cerr << "      -transitive             : \"-A 3 -B 3 -f .15,.15,.70 -b 0.5,0.25,0.25 -BSymbols '-@#' -l 64 -s 1,8\"" << endl;
  cerr << "      -spaced                 : \"-A 2 -B 2 -f .30,.70     -b 0.75,0.25     -BSymbols '-#'  -l 64 -s 1,8\"" << endl;
  cerr << "      -iupac                  : \"-A 16 -B 27 -f <TAM30,gc50,kappa1> -b ... -BSymbols 'ACGTRrYySsWwKkMmBbDdHhVvn@N'\""<< endl;
  cerr << "                                \"-M {<iupac 16 x 27 specific matrix>} -l 64 -s 1,4\"" << endl;

  exit(-1);
}
/// display a 2d matrix (a vector of vectors of int)
void DISPLAYMATRIX(std::ostream & os, std::vector< std::vector<int> > & matrix ) {
  os << "{";
  for (int b = 0; b < gv_seed_alphabet_size; b++) {
    os << "{";
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      os << (matrix[a][b]);
      if (a < (gv_align_alphabet_size-1)) os << ",";
    }
    os << "}";
    if (b < (gv_seed_alphabet_size-1)) os << ",";
  }
  os << "}";
}

/// display a 1d vector (vector of int)
void DISPLAYTABLE(std::ostream & os, std::vector<double> & table ) {
  os << "{";
  for (int a = 0; a < (int)table.size(); a++) {
    os << (table[a]);
    if (a < (int)(table.size() - 1))
      os << ",";
  }
  os << "}";
}

/// check if the sum of probabilities is ~ 1.0
void CHECKPROB(std::vector<double> & table) {
  double result = 0.0;
  for (int j = 0; j < (int)table.size(); j++)
    result += (table)[j];
  if (result < 0.99 || result > 1.01)
    _ERROR("CHECKPROB"," sum of " << (table)[0] << " + ... + " << (table)[table.size()-1] << " must be 1.0 and not " << result);
}

/// parse one flag
void PARSEFLAG(int & i, char ** argv, int argc, bool & var, bool value) {
  var = value;
}

/// parse and check one integer
void PARSEINT(int & i,  char ** argv, int argc, int & var, int min, int max) {
  i++;
  if (i >= argc)
    _ERROR("PARSEINT","\"" << argv[i-1] << "\" found without argument");
  char * next = argv[i];
  int i_tmp = strtol(argv[i], &next, 10);
  if ((*argv[i] == '\0')  || ((*next) != '\0'))
    _ERROR("PARSEINT","\"" << argv[i-1] << "\" is not followed by an integer");
  if (i_tmp>max || i_tmp <min)
    _ERROR("PARSEINT","\"" << argv[i-1] << "\" is followed by an integer outside the valid range");
  var = i_tmp;
}

/// parse and check one integer or a set of options from an enum
void PARSEINTORENUM(int & i,  char ** argv, int argc, int & var, int min, int max, char * enum_string_values[], int enum_string_values_size, bool &enum_set, int &enum_set_value) {
  i++;
  if (i >= argc)
    _ERROR("PARSEINTORENUM","\"" << argv[i-1] << "\" found without argument");
  char * next = argv[i];
  /* check in allowed strings before */
  for (int j = 0; j < enum_string_values_size; j++) {
    if (strcmp(enum_string_values[j],next) == 0) {
      enum_set       = true;
      enum_set_value = j;
      return;
    }
  }
  /* check in integers in not a string */
  int i_tmp = strtol(argv[i], &next, 10);
  if ((*argv[i] == '\0')  || ((*next) != '\0')) {
    cerr << "[*] = (";
    for (int j = 0; j < enum_string_values_size; j++) {
      if (j)
        cerr << ",";
      cerr << (enum_string_values[j]);
    }
    cerr << ")" << endl;
    _ERROR("PARSEINTORENUM","\"" << argv[i-1] << "\" is not followed by an integer or a correct value [*]");
  }
  /* FIXME HERE : more detail */
  if (i_tmp > max || i_tmp < min)
    _ERROR("PARSEINTORENUM","\"" << argv[i-1] << "\" is followed by an integer outside the valid range");
  var = i_tmp;
}


/// parse and check a pair of integers
void PARSEDINT(int & i, char ** argv, int argc,
               int & var1, int min1, int max1,
               int & var2, int min2, int max2) {
  i++;
  if (i >= argc)
    _ERROR("PARSEDINT","\"" << argv[i-1] << "\" found without argument");
  int i_tmp = sscanf(argv[i],"%i,%i",&var1,&var2);
  if (i_tmp != 2 )
    _ERROR("PARSEDINT","\"" << argv[i-1] << "\" is not followed by an <int>,<int>");
  if (var1 < min1 || var1 > max1 || var2 < min2 || var2 > max2)
    _ERROR("PARSEDINT","\"" << argv[i-1] << "\" is followed by integers out of the range");
  if (var1 > var2)
    _ERROR("PARSEDINT","\"" << var1 << "\" must be <= \"" << var2 << "\"");
}

/// parse and check a list of integers
void PARSEINTS(int & i, char ** argv, int argc,
               std::vector<int> & r_table, int neededsize = 0, bool minmax=false, int vmin=0, int vmax=1000, bool complete=false) {
  i++;
  if (i >= argc)
    _ERROR("PARSEINTS","\"" << argv[i-1] << "\" found without argument");

  // read table
  r_table.clear();
  char * pch = strtok(argv[i],",;");
  while (pch){
    int value = 0;
    int i_tmp = sscanf(pch,"%d",&value);
    if (i_tmp != 1)
      _ERROR("PARSEINTS","\"" << pch << "\" is not a correct <int>");
    if (minmax) {
      if (value < vmin || value > vmax)
        _ERROR("PARSEINTS","\"" << value << "\" is by an integer out of the range");
    }
    r_table.push_back(value);
    pch = strtok(NULL,",;");
  }

  // check table size if needed
  if (neededsize != 0) {
    if (complete && (int) r_table.size() < neededsize) {
      while ((int)r_table.size() < neededsize)
        r_table.push_back(r_table.back());
    }
    if ((int)r_table.size() != neededsize)
      _ERROR("PARSEINTS", " there is not a correct number of <int> values (" << (r_table.size()) << ") given by " << argv[i-1] << " when compared with the needed size " << neededsize);
  }
}

/// parse and check a double
void PARSEDOUBLE(int & i, char ** argv, int argc, double & var , double min, double max) {
  i++;
  if (i >= argc)
    _ERROR("PARSEDOUBLE","\"" << argv[i-1] << "\" found without argument");
  char * next = argv[i];
  double d_tmp = strtod(argv[i],&next);
  if ((*argv[i] == '\0')  || ((*next) != '\0'))
    _ERROR("PARSEDOUBLE","\"" << argv[i-1] << "\" is not followed by a double");
  if (d_tmp>max || d_tmp <min)
    _ERROR("PARSEDOUBLE","\"" << argv[i-1] << "\" is followed by a double outside the valid range");
  var = d_tmp;
}

/// parse and check a pair of double
void PARSEDDOUBLE(int & i, char ** argv, int argc,
                  double & var1, double min1, double max1,
                  double & var2, double min2, double max2) {
  i++;
  if (i >= argc)
    _ERROR("PARSEDDOUBLE","\"" << argv[i-1] << "\" found without argument");
  int i_tmp = sscanf(argv[i],"%lf,%lf",&var1,&var2);
  if (i_tmp != 2 )
    _ERROR("PARSEDDOUBLE","\"" << argv[i-1] << "\" is not followed by an <dbl>,<dbl>");
  if (var1 < min1 || var1 > max1 || var2 < min2 || var2 > max2)
    _ERROR("PARSEDDOUBLE","\"" << argv[i-1] << "\" is followed by doubles out of the range");
  if (var1 > var2)
    _ERROR("PARSEDDOUBLE","\"" << var1 << "\" must be <= \"" << var2 << "\"");
}

/// parse and check a list of double
void PARSEDOUBLES(int & i, char ** argv, int argc,
                  std::vector<double> & r_table, int neededsize = 0, bool positive_values=false) {
  i++;
  if (i >= argc)
    _ERROR("PARSEDOUBLES","\"" << argv[i-1] << "\" found without argument");

  // read table
  r_table.clear();
  char * pch = strtok(argv[i],",;");
  while (pch){
    double value = 0.0;
    int i_tmp = sscanf(pch,"%lf",&value);
    if (i_tmp != 1)
      _ERROR("PARSEDOUBLES","\"" << pch << "\" is not a correct <dbl>");
    if (positive_values) {
      if (value < 0)
        _ERROR("PARSEDOUBLES","\"" << value << "\" is a negative <dbl>");
    }
    r_table.push_back(value);
    pch = strtok(NULL,",;");
  }

  // check table size if needed
  if (neededsize != 0) {
    if ((int)r_table.size() != neededsize)
      _ERROR("PARSEDOUBLES", " there is not a correct number of <dbl> values (" << (r_table.size()) << ") given by " << argv[i-1] << " when compared with the needed size " << neededsize);
  }
}

/// parse and check if an options commes from a list of strings (retunr the index of the string)
void PARSESTRING(int & i, char ** argv, int argc,
                 int & value,
                 char * values[],
                 int values_size
                 ) {
  i++;
  if (i >= argc)
    _ERROR("PARSESTRING","\"" << argv[i-1] << "\" found without argument");
  for (int j=0; j<values_size; j++) {
    if (strcmp(argv[i],values[j]) == 0) {
      value = j;
      return;
    }
  }
  _ERROR("PARSESTRING","\"" << argv[i-1] << "\" argument \"" << argv[i] << "\" is not accepted by this option");
}

/// check if a signature is compatible with a given weight range (this must be tested when the weight is changed, or the signature changed)
void CHECKSIGNATURE() {
  if (gv_weight_interval_flag && gv_signature_flag) {
    double signature_weight = 0.0;
    for (int b = 0; b < gv_seed_alphabet_size; b++)
      signature_weight += gv_signature[b] * gv_bsel_weight[b];
    if (signature_weight < gv_minweight)
      _ERROR("CHECKSIGNATURE", " the \"-i\" signature computed weight " << signature_weight << " is too small when compared with the current weight range " << gv_minweight << "," << gv_maxweight);
    if (signature_weight > gv_maxweight)
      _ERROR("CHECKSIGNATURE", " the \"-i\" signature computed weight " << signature_weight << " is too large when compared with the current weight range " << gv_minweight << "," << gv_maxweight);
  }
}

/// parse and check a signature (number of seed elements inside a seed)
void PARSESIGNATURE(int & i, char ** argv, int argc, std::vector<int> & r_table) {
  i++;
  if (i >= argc)
    _ERROR("PARSESIGNATURE","\"" << argv[i-1] << "\" found without argument");

  // read table
  char * pch = strtok(argv[i],",;");
  r_table.clear();
  while (pch){
    int value = 0;
    int i_tmp = sscanf(pch,"%d",&value);
    if (i_tmp != 1)
      _ERROR("PARSESIGNATURE","\"" << pch << "\" is not a correct integer");
    if (value < 0 || value >= 64)
      _ERROR("PARSESIGNATURE","\"" << pch << "\" is out of the range [0..N] (N=64)");
    r_table.push_back(value);
    pch = strtok(NULL,",;");
  }

  // check table size
  if ((int)r_table.size() != gv_seed_alphabet_size)
    _ERROR("PARSESIGNATURE", " there is not a correct number of <int> values given by " << argv[i-1] << " when compared with the current seed alphabet size B = " << gv_seed_alphabet_size);
  // check signature and weight interval
  CHECKSIGNATURE();
}

/// parse and check a homogeneous scoring system
void PARSEHOMOGENEOUS(int & i, char ** argv, int argc, std::vector<int> & r_table) {
  i++;
  if (i >= argc)
    _ERROR("PARSEHOMOGENEOUS","\"" << argv[i-1] << "\" found without argument");

  // read table
  char * pch = strtok(argv[i],",;");
  r_table.clear();
  while (pch){
    int value = 0;
    int i_tmp = sscanf(pch,"%d",&value);
    if (i_tmp != 1)
      _ERROR("PARSEHOMOGENEOUS","\"" << pch << "\" is not a correct integer");
    if (value < -1000 || value > 1000)
      _ERROR("PARSEHOMOGENEOUS","\"" << pch << "\" is out of the range [-1000..1000]");
    r_table.push_back(value);
    pch = strtok(NULL,",;");
  }

  // check table size
  if ((int)r_table.size() != gv_align_alphabet_size)
    _ERROR("PARSEHOMOGENEOUS", " there is not a correct number of <int> values given by " << argv[i-1] << " when compared with the current align alphabet size A = " << gv_align_alphabet_size);
}

/// parse and check a set of probabilities
void PARSEPROBS(int & i, char ** argv, int argc,
                std::vector<double> & r_table, int & k) {
  i++;
  if (i >= argc)
    _ERROR("PARSEPROBS", "\"" << argv[i-1] << "\" found without argument");

  // read table
  char * pch = strtok(argv[i], ",;");
  r_table.clear();
  while (pch){
    double value = 0.0;
    int i_tmp = sscanf(pch, "%lf", &value);
    if (i_tmp != 1)
      _ERROR("PARSEPROBS","\"" << pch << "\" is not a correct <dbl>");
    if (value < 0 || value > 1 )
      _ERROR("PARSEPROBS","\"" << pch << "\" is  out of the range [0,1]");
    r_table.push_back(value);
    pch = strtok(NULL,",;");
  }

  // check table size
  int a = gv_align_alphabet_size;
  k = 0;
  while(1) {
    if ((int)r_table.size() == a)
      break;
    else if ((int)r_table.size() < a)
      _ERROR("PARSEPROBS", " there is not a correct number of <dbl> values (" << (r_table.size()) << ") given by " << argv[i-1] << " when compared with the current alphabet size A = " << gv_align_alphabet_size);
    a *= gv_align_alphabet_size;
    k++;
  }
  CHECKPROB(r_table);
}

/// parse and check a set of probabilities as an automaton file
void PARSEPROBSAUTOMATONFILE(int & i, char ** argv, int argc, automaton<double> * &r_automaton) {
  i++;
  if (i >= argc)
    _ERROR("PARSEPROBSAUTOMATONFILE","\"" << argv[i-1] << "\" found without argument");

  // read file
  ifstream ifs_file;
  ifs_file.open(argv[i]);
  if (!ifs_file){
    _ERROR("PARSEPROBSAUTOMATONFILE","unreadable file \"" << argv[i] << "\" ");
  }

  // read the content and set the automaton
  if (r_automaton)
    delete r_automaton;
  r_automaton = new automaton<double>();

  ifs_file >> (*r_automaton);
  ifs_file.close();
}


/// parse and check a set of polynomial probabilities as an automaton file
void PARSEMULTIPOLYAUTOMATONFILE(int & i, char ** argv, int argc, automaton<polynomial<BIGINT > > * &r_automaton) {
  i++;
  if (i >= argc)
    _ERROR("PARSEPOLYPROBSAUTOMATONFILE","\"" << argv[i-1] << "\" found without argument");

  // read file
  ifstream ifs_file;
  ifs_file.open(argv[i]);
  if (!ifs_file){
    _ERROR("PARSEPOLYPROBSAUTOMATONFILE","unreadable file \"" << argv[i] << "\" ");
  }

  // read the content and set the automaton
  if (r_automaton)
    delete r_automaton;
  r_automaton = new automaton<polynomial<BIGINT > >();

  ifs_file >> (*r_automaton);
  ifs_file.close();
}



/// set of states on a matrix input @see PARSEMATRIX
typedef enum {
  dummy,
  readfirstopenbracket,
  readsecondopenbracket,
  readint,
  readintseparator,
  readsecondclosedbracket,
  readsecondbracketseparator,
  readfirstclosedbracket,
  readfirstbracketseparator
} automatareadmatrixstates;



/// parse and check a matrix input
void PARSEMATRIX(int & i, char ** argv, int argc, std::vector< std::vector<int> > & matrix, int min, int max, int gnbcolumns, int gnbrows) {

  i++;
  if (i >= argc)
    _ERROR("PARSEMATRIX","\"" << argv[i-1] << "\" found without argument");
  char * pch = argv[i];
  int nbcolumns = 0, nbrows = 0;

  automatareadmatrixstates state = dummy;
  // automata
  while(*pch) {

    switch (*pch){

    case '{':
      switch (state){
      case dummy:
        nbrows = 0;
        state = readfirstopenbracket;
        break;
      case readfirstopenbracket:
      case readsecondclosedbracket:
      case readsecondbracketseparator:
        nbcolumns = 0;
        nbrows++;
        if (nbrows > gnbrows)
          _ERROR("PARSEMATRIX","incorrect number of lines : " << nbrows << " != " << gnbrows << " in \"" << argv[i] << "\" string parameter ");

        state = readsecondopenbracket;
        break;
      default:
        _ERROR("PARSEMATRIX","incorrect \'" << *pch << "\' in \"" << argv[i] << "\" string parameter ");
      }
      pch++;
      break;

    case '}':
      switch (state){
      case readint:
      case readintseparator:
        state = readsecondclosedbracket;
        if (nbcolumns != gnbcolumns)
          _ERROR("PARSEMATRIX","incorrect number of columns :  " << nbcolumns << " != " << gnbcolumns << " in \"" << argv[i] << "\" string parameter ");
        break;
      case readsecondclosedbracket:
        state = readfirstclosedbracket;
        if (nbrows != gnbrows )
          _ERROR("PARSEMATRIX","incorrect number of lines : " << nbrows << " != " << gnbrows << " in \"" << argv[i] << "\" string parameter ");
        break;
      default:
        _ERROR("PARSEMATRIX","incorrect \'" << *pch << "\' in \"" << argv[i] << "\" string parameter ");
      }
      pch++;
      break;
    case ' ':
    case '\t':
    case '\n':
      pch++;
      break;

    case ';':
    case ',':
      switch (state){
      case readint:
        state = readintseparator;
        break;
      case readsecondclosedbracket:
        state = readsecondbracketseparator;
        break;
      case readfirstclosedbracket:
        state = readfirstbracketseparator;
        break;
      default:
        _ERROR("PARSEMATRIX","incorrect \'" << *pch << "\' in \"" << argv[i] << "\" string parameter ");
      }
      pch++;
      break;

    default: // read value
      if (state == readsecondopenbracket || state == readintseparator) {
        state = readint;
        int value = strtol (pch, &pch, 10);
        if (value>max || value<min)
          _ERROR("PARSEMATRIX","\"" << argv[i] << "\" contains at least one integer outside the valid range");
        if (nbcolumns >= gnbcolumns)
          _ERROR("PARSEMATRIX","incorrect number of columns : " << (nbcolumns+1) << ">" <<  gnbcolumns << " in \"" << argv[i] << "\" string parameter ");
        if (nbrows > gnbrows)
          _ERROR("PARSEMATRIX","incorrect number of lines : " << nbrows << ">" <<  gnbrows << " in \"" << argv[i] << "\" string parameter ");
        matrix[nbcolumns][nbrows-1] = value;
        nbcolumns++;

      } else {
        _ERROR("PARSEMATRIX","incorrect \'" << *pch << "\' in \"" << argv[i] << "\" string parameter ");
      }
    }
  } // while
  if (state == readfirstclosedbracket || state == readfirstbracketseparator)
    return;
  else
    _ERROR("PARSEMATRIX","\" unfinished string parameter \"" << argv[i] << "\" ");
}


/// parse and check a matrix input as a file
void PARSEMATRIXFILE(int & i, char ** argv, int argc, std::vector< std::vector<int> > & matrix, int min, int max, int gnbcolumns, int gnbrows) {

  i++;
  if (i >= argc)
    _ERROR("PARSEMATRIXFILE","\"" << argv[i-1] << "\" found without argument");

  // read file
  ifstream ifs_file;
  ifs_file.open(argv[i]);
  if (!ifs_file){
    _ERROR("PARSEMATRIXFILE","unreadable file \"" << argv[i] << "\" ");
  }

  // transform it simply into a string
  char c;
  string s;
  while (ifs_file.get(c)) {
    s += c;
  }
  const char * sc = s.c_str();


  char * pch = (char *) sc;
  int nbcolumns = 0, nbrows = 0;

  automatareadmatrixstates state = dummy;
  // automata
  while(*pch) {

    switch (*pch){

    case '{':
      switch (state){
      case dummy:
        nbrows = 0;
        state = readfirstopenbracket;
        break;
      case readfirstopenbracket:
      case readsecondclosedbracket:
      case readsecondbracketseparator:
        nbcolumns = 0;
        nbrows++;
        if (nbrows > gnbrows)
          _ERROR("PARSEMATRIXFILE","incorrect number of lines : " << nbrows << " != " << gnbrows << " in \"" << argv[i] << "\" string parameter ");

        state = readsecondopenbracket;
        break;
      default:
        _ERROR("PARSEMATRIXFILE","incorrect \'" << *pch << "\' in \"" << argv[i] << "\" string parameter ");
      }
      pch++;
      break;

    case '}':
      switch (state){
      case readint:
      case readintseparator:
        state = readsecondclosedbracket;
        if (nbcolumns != gnbcolumns)
          _ERROR("PARSEMATRIXFILE","incorrect number of columns :  " << nbcolumns << " != " << gnbcolumns << " in \"" << argv[i] << "\" string parameter ");
        break;
      case readsecondclosedbracket:
        state = readfirstclosedbracket;
        if (nbrows != gnbrows )
          _ERROR("PARSEMATRIXFILE","incorrect number of lines : " << nbrows << " != " << gnbrows << " in \"" << argv[i] << "\" string parameter ");
        break;
      default:
        _ERROR("PARSEMATRIXFILE","incorrect \'" << *pch << "\' in \"" << argv[i] << "\" string parameter ");
      }
      pch++;
      break;

    case ' ':
    case '\t':
    case '\n':
      pch++;
      break;

    case ';':
    case ',':
      switch (state){
      case readint:
        state = readintseparator;
        break;
      case readsecondclosedbracket:
        state = readsecondbracketseparator;
        break;
      case readfirstclosedbracket:
        state = readfirstbracketseparator;
        break;
      default:
        _ERROR("PARSEMATRIX","incorrect \'" << *pch << "\' in \"" << argv[i] << "\" string parameter ");
      }
      pch++;
      break;

    default: // read value
      if (state == readsecondopenbracket || state == readintseparator) {
        state = readint;
        int value = strtol (pch, &pch, 10);
        if (value>max || value<min)
          _ERROR("PARSEMATRIXFILE","\"" << argv[i] << "\" contains at least one integer outside the valid range");
        if (nbcolumns >= gnbcolumns)
          _ERROR("PARSEMATRIXFILE","incorrect number of columns : " << (nbcolumns+1) << ">" <<  gnbcolumns << " in \"" << argv[i] << "\" string parameter ");
        if (nbrows > gnbrows)
          _ERROR("PARSEMATRIXFILE","incorrect number of lines : " << nbrows << ">" <<  gnbrows << " in \"" << argv[i] << "\" string parameter ");
        matrix[nbcolumns][nbrows-1] = value;
        nbcolumns++;

      } else {
        _ERROR("PARSEMATRIXFILE","incorrect \'" << *pch << "\' in \"" << argv[i] << "\" string parameter ");
      }
    }
  } // while
  if (state == readfirstclosedbracket || state == readfirstbracketseparator)
    return;
  else
    _ERROR("PARSEMATRIXFILE","\" unfinished string parameter \"" << argv[i] << "\" ");
  ifs_file.close();
}


/// check a matrix input as a subset seed matching set
void CHECKMATCHINGMATRIX(std::vector< std::vector<int> > & matrix ) {
  gv_matching_symbol_flag = true;

  // check the '1' symbol existance
  for (int b=0; b<gv_seed_alphabet_size; b++) {
    if ( matrix[gv_align_alphabet_size-1][b] == 0 ) {
      gv_matching_symbol_flag = false;
    }
  }
  // check the '#' symbol existance
  for (int a = 0; a < gv_align_alphabet_size - 1; a++) {
    if (matrix[a][gv_seed_alphabet_size-1] == 1) {
      gv_matching_symbol_flag = false;
    }
  }
  if ( matrix[gv_align_alphabet_size-1][gv_seed_alphabet_size-1] == 0 ) {
    gv_matching_symbol_flag = false;
  }
}

/**
 * @brief parse new seed symbols when "overdubbed"
 * @param i is the command line argument number being processed
 * @param argv is the command line arguments being processed
 * @param argc is the command line maximal number of arguments being processed
 * @param size is the size of the seed alphabet
 * @param r_str is the new str where seed alphabet symbols will be stored
 * @see gv_bsymbols_flag, gv_bsymbols_array
 */
void PARSESYMBOLS(int & i, char ** argv, int argc, int size, char * &r_str) {
  i++;
  if (i >= argc)
    _ERROR("PARSESYMBOLS","\"" << argv[i-1] << "\" found without argument");
  if (strlen(argv[i]) != (unsigned) size)
    _ERROR("PARSESYMBOLS","\"" << argv[i] << "\" is not of required length " << size);
  if (r_str)
    delete[] r_str;
  r_str = new char[size];
  for (char * d = r_str, *s = argv[i]; *s; d++, s++){
    if (*s == ':') {
      _ERROR("PARSESYMBOLS","\'" << *s << "\' is a reserved symbol that must not be used in a seed");
    }
    *d = *s;
  }
}


/// set of states on a seed input @see PARSESEEDS,PARSEXSEEDS
typedef enum {
  readseed,
  readphase,
  readphasecycle
} automatareadseedstates;


/// parse a seed or a set of seeds (possibly with a set of restricted positions in a cycle)
void PARSESEEDS(int & i, char ** argv, int argc) {
  i++;
  if (i >= argc)
    _ERROR("PARSESEEDS","" << argv[i-1] << "\" found without argument");
  gv_motif_flag = true;

  for (unsigned v = 0; v < gv_seeds.size(); v++)
    delete gv_seeds[v];
  gv_seeds.clear();
  gv_cycles_flag = false;
  gv_cycles.clear();
  gv_cycles_pos_nb.clear();
  gv_cycles_pos_list.clear();

  char * pch = argv[i];
  int    num = 0;
  char * pch_seed = NULL;
  char * pch_last_seed = NULL;
  automatareadseedstates state = readseed;

  while (pch != NULL) {
    switch (state) {
      // a) reading a seed
    case readseed:
      if (!pch_seed)
        pch_last_seed = pch_seed = pch;
      switch (*pch) {
      case ':':
        if (gv_global_coverage_flag) {
          _ERROR("\"-g\" and cycles within \"-m\" are not compatible together","<to be implemented>");
        }
        state = readphase;
      case ',':
      case ';':
        *pch = '\0';
        {
          string s(pch_seed);
          gv_seeds.push_back(new seed(s));
          pch_seed = NULL;
        }
        break;
      case '\0':
        {
          string s(pch_seed);
          gv_seeds.push_back(new seed(s));
          pch_seed = NULL;
        }
        goto eowhile;

      } // switch(*pch)
      break;

      // b) reading phases
    case readphase:
      // cycle allowed

      // >>
      gv_cycles_flag = true;
      while (gv_cycles.size() < gv_seeds.size())
        gv_cycles.push_back(1);
      while (gv_cycles_pos_nb.size() < gv_seeds.size())
        gv_cycles_pos_nb.push_back(0);
      while (gv_cycles_pos_list.size() < gv_seeds.size())
        gv_cycles_pos_list.push_back(std::vector<int>(0));
      // <<

      switch(*pch) {
      case '/':
        state = readphasecycle;
      case ',':
      case ';':
        gv_cycles_pos_nb[gv_seeds.size()-1]++;
        gv_cycles_pos_list[gv_seeds.size()-1].push_back(num);
        num = 0;
        break;
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        num = num * 10 + (*pch - '0');
        break;
      default:
        _ERROR("PARSESEEDS","[between ':' and '/'] incorrect \'" << *pch << "\'  after the seed \"" << pch_last_seed << "\" string parameter ");
      } // switch(*pch)
      break;

      // c) cycle
    case readphasecycle:
      switch(*pch) {
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        num = num * 10 + (*pch - '0');
        break;
      case ',':
      case ';':
        if (num <= 0)
          _ERROR("PARSESEEDS","[after '/'] incorrect cycle size \'" << num << "\' after the seed \"" << pch_last_seed << "\" string parameter ");
        gv_cycles[gv_seeds.size()-1] = num;
        state = readseed;
        num = 0;
        break;
      case '\0':
        if (num <= 0)
          _ERROR("PARSESEEDS","[after '/'] incorrect cycle size \'" << num << "\'  after the seed \"" << pch_last_seed << "\" string parameter ");
        gv_cycles[gv_seeds.size()-1] = num;
        num = 0;
        goto eowhile;
      default:
        _ERROR("PARSESEEDS","[after '/'] incorrect \'" << *pch << "\'  after the seed \"" << pch_last_seed << "\" string parameter ");
      } // switch(*pch)
      break;
    } // switch(state)
    pch++;
  }
 eowhile:
  // >>
  if (gv_cycles_flag) {
    while (gv_cycles.size() < gv_seeds.size())
      gv_cycles.push_back(1);
    while (gv_cycles_pos_nb.size() < gv_seeds.size())
      gv_cycles_pos_nb.push_back(0);
    while (gv_cycles_pos_list.size() < gv_seeds.size())
      gv_cycles_pos_list.push_back(std::vector<int>(0));


    for (unsigned int i = 0; i < gv_seeds.size(); i++) {

      if (gv_cycles_pos_nb[i] == 0) {
        gv_cycles_pos_nb[i] = 1;
        gv_cycles_pos_list[i].push_back(1);
      }

      for (int j = 0; j < gv_cycles_pos_nb[i]; j++) {
        if (gv_cycles_pos_list[i][j] > gv_cycles[i]) {
          _ERROR("PARSESEEDS"," incorrect cycle position \"" << (gv_cycles_pos_list[i][j]) << "\" > cycle size \"" << (gv_cycles[i]) << "\"");
        }
        if (gv_cycles_pos_list[i][j] < 1) {
          _ERROR("PARSESEEDS"," incorrect cycle position \"" << (gv_cycles_pos_list[i][j]) << "\" < 1 \"");
        }
        gv_cycles_pos_list[i][j] = (gv_cycles_pos_list[i][j] + gv_seeds[i]->span() - 1) % gv_cycles[i];
      }
      sort(gv_cycles_pos_list[i].begin(), gv_cycles_pos_list[i].end());
      gv_seeds[i]->setCyclePos(gv_cycles_pos_list[i], gv_cycles[i]);
    }
  }
  // <<
}


/// parse a seed or a set of seeds to exclude
void PARSEXSEEDS(int & i, char ** argv, int argc) {
  i++;
  if (i >= argc)
    _ERROR("PARSESXEEDS","" << argv[i-1] << "\" found without argument");

  gv_xseeds.clear();
  gv_xseeds_cycles_flag = true;
  gv_xseeds_cycles.clear();
  gv_xseeds_cycles_pos_nb.clear();
  gv_xseeds_cycles_pos_list.clear();

  char * pch = argv[i];
  int    num = 0;
  char * pch_seed = NULL;
  char * pch_last_seed = NULL;
  automatareadseedstates state = readseed;

  while (pch != NULL) {
    switch (state) {
      // a) reading a seed
    case readseed:
      if (!pch_seed)
        pch_last_seed = pch_seed = pch;
      switch (*pch) {
      case ':':
        state = readphase;
      case ',':
      case ';':
        *pch = '\0';
        {
          string s (pch_seed);
          gv_xseeds.push_back(new seed(s,false));
          pch_seed = NULL;
        }
        break;
      case '\0':
        {
          string s (pch_seed);
          gv_xseeds.push_back(new seed(s,false));
          pch_seed = NULL;
        }
        goto eowhile;

      } // switch(*pch)
      break;

      // b) reading phases
    case readphase:
      // cycle allowed

      // >>
      gv_xseeds_cycles_flag = true;
      while (gv_xseeds_cycles.size() < gv_xseeds.size())
        gv_xseeds_cycles.push_back(1);
      while (gv_xseeds_cycles_pos_nb.size() < gv_xseeds.size())
        gv_xseeds_cycles_pos_nb.push_back(0);
      while (gv_xseeds_cycles_pos_list.size()  < gv_xseeds.size())
        gv_xseeds_cycles_pos_list.push_back(std::vector<int>(0));
      // <<

      switch(*pch) {
      case '/':
        state = readphasecycle;
      case ',':
      case ';':
        gv_xseeds_cycles_pos_nb[gv_xseeds.size()-1]++;
        gv_xseeds_cycles_pos_list[gv_xseeds.size()-1].push_back(num);
        num = 0;
        break;
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        num = num * 10 + (*pch - '0');
        break;
      default:
        _ERROR("PARSEXSEEDS","[between ':' and '/'] incorrect \'" << *pch << "\'  after the seed \"" << pch_last_seed << "\" string parameter ");
      } // switch(*pch)
      break;

      // c) cycle
    case readphasecycle:
      switch(*pch) {
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        num = num * 10 + (*pch - '0');
        break;
      case ',':
      case ';':
        if (num <= 0)
          _ERROR("PARSEXSEEDS","[after '/'] incorrect cycle size \'" << num << "\' after the seed \"" << pch_last_seed << "\" string parameter ");
        gv_xseeds_cycles[gv_xseeds.size()-1] = num;
        state = readseed;
        num = 0;
        break;
      case '\0':
        if (num <= 0)
          _ERROR("PARSEXSEEDS","[after '/'] incorrect cycle size \'" << num << "\'  after the seed \"" << pch_last_seed << "\" string parameter ");
        gv_xseeds_cycles[gv_xseeds.size()-1] = num;
        num = 0;
        goto eowhile;
      default:
        _ERROR("PARSEXSEEDS","[after '/'] incorrect \'" << *pch << "\'  after the seed \"" << pch_last_seed << "\" string parameter ");
      } // switch(*pch)
      break;
    } // switch(state)
    pch++;
  }
 eowhile:
  // >>
  if (gv_xseeds_cycles_flag) {
    while (gv_xseeds_cycles.size() < gv_xseeds.size())
      gv_xseeds_cycles.push_back(1);
    while (gv_xseeds_cycles_pos_nb.size() < gv_xseeds.size())
      gv_xseeds_cycles_pos_nb.push_back(0);
    while (gv_xseeds_cycles_pos_list.size() < gv_xseeds.size())
      gv_xseeds_cycles_pos_list.push_back(std::vector<int>(0));


    for (unsigned i = 0; i < gv_xseeds.size(); i++) {

      if (gv_xseeds_cycles_pos_nb[i] == 0) {
        gv_xseeds_cycles_pos_nb[i] = 1;
        gv_xseeds_cycles_pos_list[i].push_back(1);
      }

      for (int j = 0; j < gv_xseeds_cycles_pos_nb[i]; j++) {
        if (gv_xseeds_cycles_pos_list[i][j] > gv_xseeds_cycles[i]) {
          _ERROR("PARSESEEDS"," incorrect cycle position \"" << (gv_xseeds_cycles_pos_list[i][j]) << "\" > cycle size \"" << (gv_xseeds_cycles[i]) << "\"");
        }
        if (gv_xseeds_cycles_pos_list[i][j] < 1) {
          _ERROR("PARSESEEDS"," incorrect cycle position \"" << (gv_xseeds_cycles_pos_list[i][j]) << "\" < 1 \"");
        }
        gv_xseeds_cycles_pos_list[i][j] = (gv_xseeds_cycles_pos_list[i][j] + gv_xseeds[i]->span() - 1) % gv_xseeds_cycles[i];
      }
      sort(gv_xseeds_cycles_pos_list[i].begin(), gv_xseeds_cycles_pos_list[i].end());
      gv_xseeds[i]->setCyclePos(gv_xseeds_cycles_pos_list[i], gv_xseeds_cycles[i]);
    }
  }
  // <<
}


/// parse a pareto input file
void PARSEINPUT(int & i, char ** argv, int argc) {
  i++;
  if (i >= argc)
    _ERROR("PARSEINPUT","\"" << argv[i-1] << "\" found without argument");

  char * pch = strtok(argv[i],",;");
  while (pch) {
    char * file = strdup(pch);
    gv_input_filenames[gv_nb_input_filenames++] = file;
    pch = strtok(NULL,",;");
  }
}

/// parse the pareto output file
void PARSEOUTPUT(int & i, char ** argv, int argc) {
  i++;
  if (i >= argc)
    _ERROR("PARSEINPUT","\"" << argv[i-1] << "\" found without argument");
  gv_output_filename = strdup(argv[i]);
}

/// main function that scan command line arguments
void SCANARG(int argc , char ** argv) {
  if (argc <= 1) {
    USAGE();
  }

  for (int i = 1; i < argc; i++){
    // 1) Seed model
    if (!strcmp(argv[i],"-A")||!strcmp(argv[i],"--A")) {
      PARSEINT(i, argv, argc, gv_align_alphabet_size, 2, 1024);
      build_default_subsetseed_matching_matrix();
      build_default_vectorizedsubsetseed_scoring_matrix();
      build_default_probabilities();
      if (gv_homogeneous_flag) {
        gv_homogeneous_flag = false;
        _WARNING("\"-u\" OPTION DISABLED","align alphabet size was changed \"after\" setting the \"-u\" option");
      }
      if (gv_lossless_flag) {
        gv_lossless_flag    = false;
        _WARNING("\"-L\" OPTION DISABLED","align alphabet size was changed \"after\" setting the \"-L\" option");
      }
      if (gv_align_alphabet_size > 2 && gv_correlation_flag) {
        _ERROR("\"Alignment alphabet of size greater than 2 is not compatible with correlation computation","<not implemented yet>");
      }
    } else if (!strcmp(argv[i],"-B")||!strcmp(argv[i],"--B")) {
      PARSEINT(i, argv, argc, gv_seed_alphabet_size, 1, 1024);
      build_default_subsetseed_matching_matrix();
      build_default_vectorizedsubsetseed_scoring_matrix();
      build_default_probabilities();
      computeWeight();
      gv_global_coverage_cost = std::vector<int>(gv_seed_alphabet_size);
      for (int i = 0; i < gv_seed_alphabet_size; i++)
        gv_global_coverage_cost[i] = i;
      if (gv_signature_flag || gv_signature_shuffle_from_m_pattern_flag) {
        gv_signature_flag = false;
        gv_signature.clear();
        gv_signature_shuffle_from_m_pattern_flag = false;
        _WARNING("\"-i\" OPTION DISABLED","seed alphabet size was changed \"after\" setting the \"-i\" option");
      }
      if (gv_bsymbols_flag) {
        gv_bsymbols_flag  = false;
        _WARNING("\"-BSymbols\" OPTION DISABLED","seed alphabet size was changed \"after\" setting the \"-BSymbols\" option");
      }
      if (gv_xseeds.size()) {
        gv_xseeds.clear();
        gv_xseeds_cycles_flag = false;
        _WARNING("\"-mx\" OPTION DISABLED","seed alphabet size was changed \"after\" setting the \"-mx\" option");
      }
      if (gv_motif_flag) {
        gv_motif_flag = false;
        for (unsigned v = 0; v < gv_seeds.size(); v++) {
          delete gv_seeds[v];
          gv_seeds[v] = NULL;
        }
        _WARNING("\"-m\" OPTION DISABLED","seed alphabet size was changed \"after\" setting the \"-m\" option");
      }
      if (gv_lossless_flag) {
        gv_lossless_flag = false;
        _WARNING("\"-L\" OPTION DISABLED","seed alphabet size was changed \"after\" setting the \"-L\" option");
      }
      if (gv_homogeneous_flag) {
        gv_homogeneous_flag = false;
        gv_homogeneous_scores.clear();
        _WARNING("\"-u\" OPTION DISABLED","seed alphabet size was changed \"after\" setting the \"-u\" option");
      }
      // 1.1) subset seeds
    } else if (!strcmp(argv[i],"-M")||!strcmp(argv[i],"--Matching")) {
      PARSEMATRIX(i, argv, argc, gv_subsetseed_matching_matrix, 0, 1, gv_align_alphabet_size, gv_seed_alphabet_size);
      CHECKMATCHINGMATRIX(gv_subsetseed_matching_matrix);
    } else if (!strcmp(argv[i],"-MF")||!strcmp(argv[i],"--MatchingFile")) {
      PARSEMATRIXFILE(i, argv, argc, gv_subsetseed_matching_matrix, 0, 1, gv_align_alphabet_size, gv_seed_alphabet_size);
      CHECKMATCHINGMATRIX(gv_subsetseed_matching_matrix);

      // 1.2) vectorized subset seeds
    } else if (!strcmp(argv[i],"-S")||!strcmp(argv[i],"--Scoring")) {
      PARSEMATRIX(i, argv, argc, gv_vectorizedsubsetseed_scoring_matrix,-10000, 10000, gv_align_alphabet_size, gv_seed_alphabet_size);
    } else if (!strcmp(argv[i],"-SF")||!strcmp(argv[i],"--ScoringFile")) {
      PARSEMATRIXFILE(i, argv, argc, gv_vectorizedsubsetseed_scoring_matrix,-10000, 10000, gv_align_alphabet_size, gv_seed_alphabet_size);
    } else if (!strcmp(argv[i],"-T")||!strcmp(argv[i],"--Threshold")) {
      PARSEINT(i, argv, argc, gv_vectorizedsubsetseed_scoring_threshold,-10000,+10000);
    } else if (!strcmp(argv[i],"-V")||!strcmp(argv[i],"--Vectorized")) {
      if (gv_global_coverage_flag) {
        _ERROR("\"-V\" and \"-g\" are not compatible together","<to be defined before implemented>");
      }
      PARSEFLAG(i, argv, argc, gv_vectorized_flag, true);

      // 1.3) lossless seeds
    } else if (!strcmp(argv[i],"-L")||!strcmp(argv[i],"--LosslessCosts")) {
      if (gv_correlation_flag) {
        _ERROR("lossless mode and \"-g /or/ -y <CORRELATION>\" are not compatible together","<not implemented yet, but do you need it ?>");
      }
      if (gv_homogeneous_flag) {
        _ERROR("lossless mode and \"-u\" homogeneous are not compatible together","<not implemented yet, but do you need it ?>");
      }
      if (gv_subalignment_flag) {
        _ERROR("lossless mode and \"-ll\" are not compatible together","<not implemented yet, and no simple way of doing this ?>");
      }
      PARSEINTS(i, argv, argc, gv_lossless_costs_vector, gv_align_alphabet_size, true, 0, 1000);
      gv_lossless_flag = true;

    } else if (!strcmp(argv[i],"-X")||!strcmp(argv[i],"--LosslessCostThreshold")) {
      PARSEINT(i, argv, argc, gv_lossless_cost_threshold, 0, 1000);

      // 2) Alignment model
    } else if (!strcmp(argv[i],"-b")||!strcmp(argv[i],"--background")) {
      PARSEPROBS(i, argv, argc, gv_bsel, gv_bsel_k);
      computeWeight();
    } else if (!strcmp(argv[i],"-f")||!strcmp(argv[i],"--foreground")) {
      PARSEPROBS(i, argv, argc, gv_bsens, gv_bsens_k);
    } else if (!strcmp(argv[i],"-fF")||!strcmp(argv[i],"--foregroundFile")) {
      PARSEPROBSAUTOMATONFILE(i, argv, argc, gv_bsens_automaton);
    } else if (!strcmp(argv[i],"-l")||!strcmp(argv[i],"--length")) {
      PARSEINT(i, argv, argc, gv_alignment_length, 1, 1000000);
      if (gv_subalignment_flag && gv_subalignment_length >= gv_alignment_length) {
        _ERROR(" \"-l\" alignment length value is set UNDER the sub-alignment length \"-ll\"" ,"you must provide alignment length GREATER THAN sub-alignment length");
      }
    } else if (!strcmp(argv[i],"-ll")||!strcmp(argv[i],"--sublength")) {
      PARSEINT(i, argv, argc, gv_subalignment_length, 1, gv_alignment_length);
      if (gv_correlation_flag) {
        _ERROR("\"-ll\" and \"-g /or/ -y <CORRELATION>\"  are not compatible together","<many ways to make sense with this combination ?>");
      }
      if (gv_homogeneous_flag) {
        _ERROR("\"-ll\" and \"-u\" homogeneous are not compatible together","<not implemented yet, and no simple way of doing this ?>");
      }
      if (gv_lossless_flag) {
        _ERROR("\"-ll\" and lossless mode are not compatible together","<not implemented yet, and no simple way of doing this ?>");
      }
      gv_subalignment_flag = true;
    } else if (!strcmp(argv[i],"-llf")||!strcmp(argv[i],"--sublengthfunction")) {
      PARSESTRING(i, argv, argc, gv_subalignment_function_index, gv_subalignment_functions_names,
                  sizeof(gv_subalignment_functions_names) / sizeof(*gv_subalignment_functions_names));
    } else if (!strcmp(argv[i],"-u")||!strcmp(argv[i],"--homogeneous")) {
      if (gv_correlation_flag) {
        _ERROR("\"-u\" homogeneous and \"-g /or/ -y <CORRELATION>\" are not compatible together","<not implemented yet, but do you need it ?>");
      }
      if (gv_subalignment_flag) {
        _ERROR("\"-u\" homogeneous and \"-ll\" are not compatible together","<not implemented yet, and no simple way of doing this ?>");
      }
      if (gv_lossless_flag) {
        _ERROR("\"-u\" homogeneous and lossless mode are not compatible together","<not implemented yet, but do you need it ?>");
      }
      PARSEHOMOGENEOUS(i, argv, argc, gv_homogeneous_scores);
      gv_homogeneous_flag = true;

      // 3) Seed enumeration
    } else if (!strcmp(argv[i],"-s")||!strcmp(argv[i],"--span")) {
      PARSEDINT(i, argv, argc, gv_minspan, 1, 1000, gv_maxspan, 1, 1000);
    } else if (!strcmp(argv[i],"-w")||!strcmp(argv[i],"--weight")) {
      gv_weight_interval_flag = true;
      PARSEDDOUBLE(i, argv, argc, gv_minweight, 0, 1e32, gv_maxweight, 0, 1e32);
      CHECKSIGNATURE();
    } else if (!strcmp(argv[i],"-n")||!strcmp(argv[i],"--number")) {
      int gv_nbseeds;
      PARSEINT(i, argv, argc, gv_nbseeds, 1, 1000);
      if (gv_motif_flag) {
        _ERROR(" \"-m\" seed pattern set BEFORE changing the number \"-n\" of seeds is not correct","you must provide the number of seeds inside \"-m\"");
      }
      if (gv_cycles_flag) {
        _ERROR(" \"-c\" cycle flag set BEFORE changing the number \"-n\" of seeds is not correct","you must provide the \"-c\"  and \"-q\" parameters AFTER the \"-n\"");
      }
      gv_seeds = std::vector<seed*>(gv_nbseeds);
    } else if (!strcmp(argv[i],"-r")||!strcmp(argv[i],"--random")) {
      PARSEINT(i, argv, argc, gv_nbruns, 1, 1000000000);
    } else if (!strcmp(argv[i],"-k")||!strcmp(argv[i],"--hillclimbing")) {
      if (gv_nbruns == 0) {
        _ERROR(" \"-k\" hillclimbing flag set but \"-r\" has not been set before (still in full enumeration mode)","always provide the \"-k\" parameter AFTER changing the default \"-r\" behavior");
      }
      PARSEFLAG(i, argv, argc, gv_hillclimbing_flag, true);
    } else if (!strcmp(argv[i],"-a")||!strcmp(argv[i],"--hillclimbing-alpha")) {
      PARSEDOUBLE(i, argv, argc, gv_hillclimbing_alpha, 0, 1);
    } else if (!strcmp(argv[i],"-d")||!strcmp(argv[i],"--disable")) {
      PARSEFLAG(i, argv, argc, gv_minimize_flag, false);
    } else if (!strcmp(argv[i],"-m")||!strcmp(argv[i],"--motifs")) {
      if (gv_cycles_flag) {
        _ERROR("\"-m\" pattern and \"-c\" are not compatible together","you must provide the cycle positions on the shape (e.g  \" -m 100101001:1, 3, 5/6\") ");
      }
      if (gv_signature_flag) {
        _ERROR(" \"-i "<< gv_signature[0] <<",...\" and \"-m\" pattern options are not compatible together", "you must provide only one of them, or use \"-i shuffle\" after setting \"-m\"");
      }
      unsigned gv_nbseeds = gv_seeds.size();
      PARSESEEDS(i, argv, argc);
      if (gv_nbseeds > 1 && gv_nbseeds != gv_seeds.size()) {
        _WARNING("\"-n " << gv_nbseeds << "\" OPTION DISABLED","\"-m\" option was set \"after\" setting the \"-n\" option");
      }
    } else if (!strcmp(argv[i],"-mx")||!strcmp(argv[i],"--motifs-excluded")) {
      if (gv_correlation_flag) {
        _ERROR("\"-g /or/ -y <CORRELATION>\" and \"-mx\" excluded seeds are not compatible together","<to be implemented>");
      }
      PARSEXSEEDS(i, argv, argc);
    } else if (!strcmp(argv[i],"-mxy")||!strcmp(argv[i],"--motifs-excluded-multihit")) {
      if (gv_xseeds.size() == 0) {
        _ERROR(" \"-mxy\" value set but \"-mx\" has not been set before","always provide the \"-mxy\" parameter AFTER setting \"-mx\"");
      }
      PARSEINT(i, argv, argc, gv_xseeds_multihit_nb, 1, 64);
      gv_xseeds_multihit_flag = true;
    } else if (!strcmp(argv[i],"-i")||!strcmp(argv[i],"--idsign")) {
      if (argc > i+1) {
        if (!strcmp(argv[i+1],"shuffle")) {
          // FIXME WARNING IF OTHER was TRUE
          if (!gv_motif_flag) {
            _ERROR(" \"-i shuffle\" set BUT \"-m <pattern>\" pattern not set BEFORE","always provide the \"-i shuffle\" AFTER setting \"-m\"");
          }
          if (gv_signature_flag) {
            _WARNING("\"-i <int>,<int>,...\" OPTION DISABLED","\"-i shuffle\" option was set \"after\" setting the \"-i <int>,<int>,...\" option");
            gv_signature_flag = false;
            gv_signature.clear();
          }
          gv_signature_shuffle_from_m_pattern_flag = true;
          i++;
        } else {
          if (gv_signature_shuffle_from_m_pattern_flag) {
            _WARNING("\"-i shuffle\" OPTION DISABLED","\"-i <int>,<int>,...\" option was set \"after\" setting the \"-i shuffle\" option");
            gv_signature_shuffle_from_m_pattern_flag = false;
          }
          if (gv_motif_flag) {
            _ERROR("\"-m\" pattern and \"-i <int>,<int>,...\" options are not compatible together", "you must provide only one of them, or use \"-i shuffle\" after setting -m");
          }
          gv_signature_flag = true;
          PARSESIGNATURE(i, argv, argc, gv_signature);
        }
      } else {
        _ERROR(" \"-i\" found without argument","always provide \"-m <pattern> -i shuffle\", or \"-i <int>,<int>,...\" with the number of seed elements required");
      }
    } else if (!strcmp(argv[i],"-x")||!strcmp(argv[i],"--symetric")) {
      PARSEFLAG(i, argv, argc, gv_symetric, true);
    } else if (!strcmp(argv[i],"-y")||!strcmp(argv[i],"--multihit")) {
      PARSEINTORENUM(i, argv, argc, gv_multihit_nb, 1, 1024, gv_correlation_functions_names, CORRELATION_FUNCTIONS_ENUM_NUMBERS, gv_correlation_flag, gv_correlation_function_index);
      if (gv_global_coverage_flag) {
        _ERROR("\"-y\" and \"-g\" are not compatible together","it does not make sense to use both.");
      }
      if (gv_lossless_flag && gv_correlation_flag) {
        _ERROR("\"-y <CORRELATION>\" and lossless mode are not compatible together","<not implemented yet, but do you need it ?>");
      }
      if (gv_subalignment_flag && gv_correlation_flag) {
        _ERROR("\"-y <CORRELATION>\" and \"-ll\" mode are not compatible together","<many ways to make sense with this combination ?>");
      }
      if (gv_homogeneous_flag && gv_correlation_flag) {
        _ERROR("\"-y <CORRELATION>\" and \"-u\" homogeneous are not compatible together","<not implemented yet, but do you need it ?>");
      }
      if (gv_xseeds.size() && gv_correlation_flag) {
        _ERROR("\"-y <CORRELATION>\" and \"-mx\" excluded seeds are not compatible together","<to be implemented>");
      }
      if (gv_align_alphabet_size > 2 && gv_correlation_flag) {
        _ERROR("\"Alignment alphabet of size greater than 2 is not compatible with correlation computation","<not implemented yet>");
      }
#ifndef USEINFINT
      if (gv_correlation_flag) {
        _WARNING("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","these specific functions : Pearson/Spearman correlation count on <long long> and computation on <double> may overflow ...\n you can compile this program with USEINFINT defined (-DUSEINFINT or see inside \"macro.hh\") but it will be much slower");
      }
#endif
      gv_multihit_flag = true;
    } else if (!strcmp(argv[i],"-g")||!strcmp(argv[i],"--global-coverage")) {
      PARSEINTORENUM(i, argv, argc, gv_global_coverage_nb, 1, 1024, gv_correlation_functions_names, CORRELATION_FUNCTIONS_ENUM_NUMBERS, gv_correlation_flag, gv_correlation_function_index);
      if (gv_align_alphabet_size > 2 && gv_correlation_flag) {
        _ERROR("\"Alignment alphabet of size greater than 2 is not compatible with correlation computation","<not implemented yet>");
      }
      if (gv_cycles_flag) {
        if (gv_motif_flag) {
          _ERROR("\"-g\" and cycles within \"-m\" are not compatible together","<to be implemented>");
        } else {
          _ERROR("\"-g\" and \"-c\" are not compatible together","<to be implemented>");
        }
      }
      if (gv_vectorized_flag) {
        _ERROR("\"-g\" and \"-V\" are not compatible together","<to be defined before implemented>");
      }
      if (gv_multihit_flag) {
        _ERROR("\"-g\" and \"-y\" are not compatible together","it does not make sense to use both.");
      }
      if (gv_lossless_flag && gv_correlation_flag) {
        _ERROR("\"-g <CORRELATION>\" and lossless mode are not compatible together.","<not implemented yet, but do you need it ?>");
      }
      if (gv_subalignment_flag && gv_correlation_flag) {
        _ERROR("\"-g <CORRELATION>\" and \"-ll\" mode are not compatible together","<many ways to make sense to this combination ?>");
      }
      if (gv_homogeneous_flag && gv_correlation_flag) {
        _ERROR("\"-g <CORRELATION>\" and \"-u\" homogeneous are not compatible together","<not implemented yet, but do you need it ?>");
      }
      if (gv_xseeds.size() && gv_correlation_flag) {
        _ERROR("\"-g <CORRELATION>\" and \"-mx\" excluded seeds are not compatible together","<to be implemented>");
      }
#ifndef USEINFINT
      if (gv_correlation_flag) {
        _WARNING("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","these specific functions : Pearson/Spearman correlation count on <long long> and computation on <double> may overflow ...\n you can compile this program with USEINFINT defined (-DUSEINFINT or see inside \"macro.hh\") but it will be much slower");
      }
#endif
      gv_global_coverage_flag = true;
      gv_global_coverage_cost = std::vector<int>(gv_seed_alphabet_size);
      for (int i = 0; i < gv_seed_alphabet_size; i++)
        gv_global_coverage_cost[i] = i;
    } else if (!strcmp(argv[i],"-p")||!strcmp(argv[i],"--polynomial-dominance")) {
      ///@todo{FIXME : check several parameters incompatible with dominant selection and output}
      gv_polynomial_dominant_selection_flag = true;
#ifndef USEINFINT
      if (gv_polynomial_dominant_selection_flag) {
        _WARNING("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","these specific functions : Polynomial coefficients count on <long long> may overflow ...\n you can compile this program with USEINFINT defined (-DUSEINFINT or see inside \"macro.hh\") but it will be much slower");
      }
#endif
    } else if (!strcmp(argv[i],"-pF")||!strcmp(argv[i],"--multipolynomial-file")) {
      ///@todo{FIXME : check several parameters incompatible with dominant selection and output}
#ifndef USEINFINT
      _WARNING("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","these specific functions :  Multivariate polynomial evaluation on <long long> may overflow ...\n you can compile this program with USEINFINT defined (-DUSEINFINT or see inside \"macro.hh\") but it will be much slower");
#endif
      gv_multipoly_file_flag = true;
      PARSEMULTIPOLYAUTOMATONFILE(i, argv, argc,  gv_multipoly_bsens_automaton);
    } else if (!strcmp(argv[i],"-c")||!strcmp(argv[i],"--cycles")) {
      if (gv_motif_flag) {
        _ERROR("\"-c\" pattern and \"-m\" are not compatible together", "you must provide the cycle positions on the shape (e.g  \" -m 100101001:1, 3, 5/6\") ");
      }
      if (gv_global_coverage_flag) {
        _ERROR("\"-c\" and \"-g\" are not compatible together","<to be implemented>");
      }
      gv_cycles_flag = true;
      PARSEINTS(i, argv, argc, gv_cycles, gv_seeds.size(), true, 1, 10000, true);
      // set the init cycles positions
      gv_cycles_pos_nb    = std::vector<int> (gv_seeds.size(), 1);
      gv_cycles_pos_list  = std::vector< std::vector<int> >(gv_seeds.size(), std::vector<int>(1, 0));
    } else if (!strcmp(argv[i],"-q")||!strcmp(argv[i],"--cyclespos")) {
      if (gv_cycles_flag) {
        if (gv_motif_flag) {
          _ERROR("\"-m\" pattern and \"-c\" are not compatible together", "you must provide the cycle positions on the shape (e.g  \" -m 100101001:1,3,5/6\") ");
        }
        // set the init cycles positions
        PARSEINTS(i, argv, argc, gv_cycles_pos_nb, gv_seeds.size(), true, 1, 10000, true);
        gv_cycles_pos_list = std::vector< std::vector<int> >(gv_seeds.size(), std::vector<int>(0));
        for (unsigned i = 0; i < gv_seeds.size(); i++) {
          gv_cycles_pos_nb[i] = MIN(gv_cycles_pos_nb[i], gv_cycles[i]);
          for (int j = 0; j < gv_cycles_pos_nb[i]; j++)
            gv_cycles_pos_list[i].push_back(j);
        }
      } else {
        _ERROR(" -q ","\"" << argv[i] << "\" must be preceded by -c ... ");
      }
    } else if (!strcmp(argv[i],"-j")||!strcmp(argv[i],"--jive")) {
      PARSEDOUBLE(i, argv, argc, gv_jive, 0, 1);
      // 4) misc
    } else if (!strcmp(argv[i],"-e")||!strcmp(argv[i],"--enter")) {
      PARSEINPUT(i, argv, argc);
    } else if (!strcmp(argv[i],"-o")||!strcmp(argv[i],"--output")) {
      PARSEOUTPUT(i, argv, argc);
    } else if (!strcmp(argv[i],"-z")||!strcmp(argv[i],"--statistics")) {
      PARSEINT(i, argv, argc, gv_pareto_select_runs, 0, 1000000000);
    } else if (!strcmp(argv[i],"-h")||!strcmp(argv[i],"--help")) {
      USAGE();
    } else if (!strcmp(argv[i],"-v")||!strcmp(argv[i],"--verbose")) {
      PARSEINT(i, argv, argc, gv_verbose, VERBOSITY_NONE, VERBOSITY_MAX);
    } else if (!strcmp(argv[i],"-BSymbols")) {
      PARSESYMBOLS(i, argv, argc, gv_seed_alphabet_size, gv_bsymbols_array);
      gv_bsymbols_flag = true;
    } else if (!strcmp(argv[i],"-transitive")) {
      // transitives params
      gv_align_alphabet_size = 3;
      gv_seed_alphabet_size  = 3;
      if (gv_signature_flag || gv_signature_shuffle_from_m_pattern_flag) {
        gv_signature.clear();
        gv_signature_flag = false;
        gv_signature_shuffle_from_m_pattern_flag = false;
        _WARNING("\"-i\" OPTION DISABLED","\"-transitive\" option was set \"after\" setting the \"-i\" option");
      }
      if (gv_xseeds.size()) {
        for (unsigned v = 0; v < gv_xseeds.size(); v++)
          delete gv_xseeds[v];
        gv_xseeds.clear();
        gv_xseeds_cycles.clear();
        gv_xseeds_cycles_pos_nb.clear();
        gv_xseeds_cycles_pos_list.clear();
        _WARNING("\"-mx\" OPTION DISABLED","\"-transitive\" option was set \"after\" setting the \"-mx\" option");
      }
      gv_xseeds_cycles_flag = false;
      gv_xseeds_multihit_flag = false;
      if (gv_motif_flag) {
        gv_motif_flag = false;
        for (unsigned v = 0; v < gv_seeds.size(); v++) {
          delete gv_seeds[v];
          gv_seeds[v] = NULL;
        }
        _WARNING("\"-m\" OPTION DISABLED","\"-transitive\" option was set \"after\" setting the \"-m\" option");
      }
      if (gv_lossless_flag) {
        gv_lossless_flag = false;
        _WARNING("\"-L\" OPTION DISABLED","\"-transitive\" option was set \"after\" setting the \"-L\" option");
      }
      if (gv_homogeneous_flag) {
        gv_homogeneous_flag = false;
        gv_homogeneous_scores.clear();
        _WARNING("\"-u\" OPTION DISABLED","\"-transitive\" option was set \"after\" setting the \"-u\" option");
      }
      if (gv_align_alphabet_size > 2 && gv_correlation_flag) {
        _ERROR("\"Alignment alphabet of size greater than 2 is not compatible with correlation computation","<not implemented yet>");
      }
      gv_bsens.clear();
      gv_bsens = std::vector<double>(3); gv_bsens[0] = 0.15; gv_bsens[1] = 0.15; gv_bsens[2] = 0.7;
      gv_bsens_k = 0;
      gv_bsel.clear();
      gv_bsel  = std::vector<double>(3); gv_bsel[0]  = 0.50; gv_bsel[1]  = 0.25; gv_bsel[2]  = 0.25;
      gv_bsel_k = 0;
      gv_bsel_weight = std::vector<double>(3); gv_bsel_weight[0] = 0.0; gv_bsel_weight[1] = 0.5; gv_bsel_weight[2] = 1.0;
      gv_bsel_minprob = 0.25;
      gv_bsel_maxprob = 1.0;
      gv_minspan = 1;
      gv_maxspan = 8;
      gv_minweight = 1; gv_maxweight = 8;
      gv_weight_interval_flag = false;
      gv_vectorized_flag = false;
      if (gv_bsymbols_flag)
        free(gv_bsymbols_array);
      gv_bsymbols_array = strdup(string("-@#").c_str());
      gv_bsymbols_flag = true;
      build_default_subsetseed_matching_matrix();
      build_default_vectorizedsubsetseed_scoring_matrix();
      gv_alignment_length = 64;
    } else if (!strcmp(argv[i],"-spaced")) {
      // spaced params
      gv_align_alphabet_size = 2;
      gv_seed_alphabet_size  = 2;
      if (gv_signature_flag || gv_signature_shuffle_from_m_pattern_flag) {
        gv_signature.clear();
        gv_signature_flag = false;
        gv_signature_shuffle_from_m_pattern_flag = false;
        _WARNING("\"-i\" OPTION DISABLED","\"-spaced\" option was set \"after\" setting the \"-i\" option");
      }
      if (gv_xseeds.size()) {
        for (unsigned v = 0; v < gv_xseeds.size(); v++)
          delete gv_xseeds[v];
        gv_xseeds.clear();
        gv_xseeds_cycles.clear();
        gv_xseeds_cycles_pos_nb.clear();
        gv_xseeds_cycles_pos_list.clear();
        _WARNING("\"-mx\" OPTION DISABLED","\"-spaced\" option was set \"after\" setting the \"-mx\" option");
      }
      gv_xseeds_cycles_flag = false;
      gv_xseeds_multihit_flag = false;
      if (gv_motif_flag) {
        gv_motif_flag = false;
        for (unsigned v = 0; v < gv_seeds.size(); v++) {
          delete gv_seeds[v];
          gv_seeds[v] = NULL;
        }
        _WARNING("\"-m\" OPTION DISABLED","\"-spaced\" option was set \"after\" setting the \"-m\" option");
      }
      if (gv_lossless_flag) {
        gv_lossless_flag = false;
        _WARNING("\"-L\" OPTION DISABLED","\"-spaced\" option was set \"after\" setting the \"-L\" option");
      }
      if (gv_homogeneous_flag) {
        gv_homogeneous_flag = false;
        gv_homogeneous_scores.clear();
        _WARNING("\"-u\" OPTION DISABLED","\"-spaced\" option was set \"after\" setting the \"-u\" option");
      }
      gv_bsens.clear();
      gv_bsens = std::vector<double>(2); gv_bsens[0] = 0.30; gv_bsens[1] = 0.70;
      gv_bsens_k = 0;
      gv_bsel.clear();
      gv_bsel  = std::vector<double>(2); gv_bsel[0]  = 0.75; gv_bsel[1]  = 0.25;
      gv_bsel_k = 0;
      gv_bsel_weight = std::vector<double>(2); gv_bsel_weight[0] = 0.0; gv_bsel_weight[1] = 1.0;
      gv_bsel_minprob = 0.25;
      gv_bsel_maxprob = 1.0;
      gv_minspan = 1;
      gv_maxspan = 8;
      gv_minweight = 1; gv_maxweight = 8;
      gv_weight_interval_flag = false;
      gv_vectorized_flag = false;
      if (gv_bsymbols_flag)
        free(gv_bsymbols_array);
      gv_bsymbols_array = strdup(string("-#").c_str());
      gv_bsymbols_flag  = true;
      build_default_subsetseed_matching_matrix();
      build_default_vectorizedsubsetseed_scoring_matrix();
      gv_alignment_length = 64;
    } else if (!strcmp(argv[i],"-iupac")) {
      // iupac 16x16 params
      gv_align_alphabet_size = 16;
      gv_seed_alphabet_size  = 27;
      if (gv_signature_flag || gv_signature_shuffle_from_m_pattern_flag) {
        gv_signature.clear();
        gv_signature_flag = false;
        gv_signature_shuffle_from_m_pattern_flag = false;
        _WARNING("\"-i\" OPTION DISABLED","\"-iupac\" option was set \"after\" setting the \"-i\" option");
      }
      if (gv_xseeds.size()) {
        for (unsigned v = 0; v < gv_xseeds.size(); v++)
          delete gv_xseeds[v];
        gv_xseeds.clear();
        gv_xseeds_cycles.clear();
        gv_xseeds_cycles_pos_nb.clear();
        gv_xseeds_cycles_pos_list.clear();
        _WARNING("\"-mx\" OPTION DISABLED","\"-iupac\" option was set \"after\" setting the \"-mx\" option");
      }
      gv_xseeds_cycles_flag = false;
      gv_xseeds_multihit_flag = false;
      if (gv_motif_flag) {
        gv_motif_flag = false;
        for (unsigned v = 0; v < gv_seeds.size(); v++) {
          delete gv_seeds[v];
          gv_seeds[v] = NULL;
        }
        _WARNING("\"-m\" OPTION DISABLED","\"-iupac\" option was set \"after\" setting the \"-m\" option");
      }
      if (gv_lossless_flag) {
        gv_lossless_flag = false;
        _WARNING("\"-L\" OPTION DISABLED","\"-iupac\" option was set \"after\" setting the \"-L\" option");
      }
      if (gv_homogeneous_flag) {
        gv_homogeneous_flag = false;
        gv_homogeneous_scores.clear();
        _WARNING("\"-u\" OPTION DISABLED","\"-iupac\" option was set \"after\" setting the \"-u\" option");
      }
      gv_bsens.clear();
      gv_bsens = std::vector<double>(16);
      // Doing an add-hoc alignment matrix out of TAM matrices with pam=30,gc=50%,kappa=1.0
      // * matrix filling (ACGT x ACGT order)
      // * TAM matrix set with  PAM level at 30, gc percent at 50%, kappa=1.0 (no transition bias)
      double gv_bsens_TAM_pam30_gc50_kappa1[16] = {0.188185008631682360,0.020604997122772542,0.020604997122772542,0.020604997122772542,
                                                   0.020604997122772542,0.188185008631682360,0.020604997122772542,0.020604997122772542,
                                                   0.020604997122772542,0.020604997122772542,0.188185008631682360,0.020604997122772542,
                                                   0.020604997122772542,0.020604997122772542,0.020604997122772542,0.188185008631682360};
      for (int a = 0; a < gv_align_alphabet_size; a++)
        gv_bsens[a] = gv_bsens_TAM_pam30_gc50_kappa1[a];
      gv_bsens_k = 0;
      gv_bsel.clear();
      gv_bsel  = std::vector<double>(16);
      for (int a = 0; a < gv_align_alphabet_size; a++)
        gv_bsel[a]  = 0.0625; // when gc=50%
      gv_bsel_k = 0;
      double gv_bsel_weight_gc50_tmp[27] = {1,1,1,1, // 'A','C','G','T' are only "x = x" matches in alignments, occuring with the best weight (set by default to 1)
                                            0.75,0.5,0.75,0.5,0.75,0.5,0.75,0.5,0.75,0.5,0.75,0.5,  // IUPAC symbols with 2 letters (e.g. S = {"G = G" or "C = C" matches},   and s =  {"G = G" or "C = C" OR "G = C" or "G = C" })
                                            0.60375937482, 0.20751874963,  0.60375937482, 0.20751874963,  0.60375937482, 0.20751874963,  0.60375937482, 0.20751874963, // IUPAC symbols with 3 letters
                                            0.0, 0.25, 0.5,// 'n' is a jocker, '@' is the transition acception symbol, 'N' is the match symbol

      };
      gv_bsel_weight = std::vector<double>(27); for (int b = 0; b < gv_seed_alphabet_size; b++ ) gv_bsel_weight[b] = gv_bsel_weight_gc50_tmp[b];
      gv_bsel_minprob = 0.0625;
      gv_bsel_maxprob = 1.0;
      gv_minspan = 1;
      gv_maxspan = 4;
      gv_minweight = 1; gv_maxweight = 4;
      gv_weight_interval_flag = false;
      gv_vectorized_flag = false;
      if (gv_bsymbols_flag)
        free(gv_bsymbols_array);
      gv_bsymbols_array = strdup(string("ACGTRrYySsWwKkMmBbDdHhVvn@N").c_str());
      gv_bsymbols_flag  = true;
      // THERE IS NOT MATCHING SYMBOL HERE :
      gv_matching_symbol_flag = false;
      // >>
      // Doing an add-hoc matching matrix, where "n" is the "#" symbol
      // . matrix allocation
      for (unsigned u = 0; u < gv_subsetseed_matching_matrix.size(); u++)
        gv_subsetseed_matching_matrix[u].clear();
      gv_subsetseed_matching_matrix.clear();
      for (int a = 0; a < gv_align_alphabet_size; a++) {
        gv_subsetseed_matching_matrix.push_back(vector<int>(gv_seed_alphabet_size, 0));
      }
      // . matrix filling (ACGT x ACGT order)
      int subsetseed_matching_matrix_tmp_reversed[27][16] = {
                                                    {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 'A' (1)
                                                    {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0}, // 'C' (1)
                                                    {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}, // 'G' (1)
                                                    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, // 'T' (1)
                                                    {1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}, //  iupac strict 'R' (2)
                                                    {1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0}, //  iupac lazy   'r' (4)
                                                    {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1}, //  iupac strict 'Y' (2)
                                                    {0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1}, //  iupac lazy   'y' (4)
                                                    {0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0}, //  iupac strict 'S' (2)
                                                    {0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0}, //  iupac lazy   's' (4)
                                                    {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, //  iupac strict 'W' (2)
                                                    {1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1}, //  iupac lazy   'w' (4)
                                                    {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1}, //  iupac strict 'K' (2)
                                                    {0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1}, //  iupac lazy   'k' (4)
                                                    {1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0}, //  iupac strict 'M' (2)
                                                    {1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0}, //  iupac lazy   'm' (4)
                                                    {0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1}, //  iupac strict 'B' (3)
                                                    {0,0,0,0,0,1,1,1,0,1,1,1,0,1,1,1}, //  iupac lazy   'b' (3)
                                                    {1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1}, //  iupac strict 'D' (3)
                                                    {1,0,1,1,0,0,0,0,1,0,1,1,1,0,1,1}, //  iupac lazy   'd' (3)
                                                    {1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1}, //  iupac strict 'H' (3)
                                                    {1,1,0,1,1,1,0,1,0,0,0,0,1,1,0,1}, //  iupac lazy   'h' (3)
                                                    {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0}, //  iupac strict 'V' (3)
                                                    {1,1,1,0,1,1,1,0,1,1,1,0,0,0,0,0}, //  iupac lazy   'v' (3)
                                                    {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}, // 'n' (or '-')
                                                    {1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1}, //     (   '@') transition accepting symbol
                                                    {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1}, // 'N' (or '#')

      };

      for (int a = 0; a < gv_align_alphabet_size; a++ )
        for (int b = 0; b < gv_seed_alphabet_size; b++ )
          gv_subsetseed_matching_matrix[a][b] = subsetseed_matching_matrix_tmp_reversed[b][a];
      // <<



      build_default_vectorizedsubsetseed_scoring_matrix(); // FIXME ?? NOT LEAVE IT
      gv_alignment_length = 64;
    } else _ERROR("PARSE","\"" << argv[i] << "\" is not a valid argument");
  }
}
// @}

void computeWeight() {
  gv_bsel_weight  = std::vector<double>(gv_seed_alphabet_size);
  gv_bsel_minprob =  1e32;
  gv_bsel_maxprob = -1e32;

  // compute background probs
  for (int b = 0; b < gv_seed_alphabet_size; b++) {
    double pr = 0.0;
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      if (gv_subsetseed_matching_matrix[a][b]) {
        for (int i = 0; i < power(gv_align_alphabet_size, gv_bsel_k); i++)
          pr += gv_bsel[a+i*gv_align_alphabet_size];
      }
    }
    gv_bsel_weight[b] = pr;
    gv_bsel_minprob   = MIN(pr, gv_bsel_minprob);
    gv_bsel_maxprob   = MAX(pr, gv_bsel_maxprob);

  }

  // take the log_#(b)
  for (int b = 0; b < gv_seed_alphabet_size; b++) {
    gv_bsel_weight[b] = log(gv_bsel_weight[b])/log(gv_bsel_minprob);
  }
}


/// number of classes for the (from "0" to max percent of identity = gv_alignment_length)
#define CLASSES (gv_alignment_length+1)
/// minimal percent of identity from which to compute the correlation coefficient
#define MIN_PVAL (gv_multihit_flag?gv_multihit_nb:(gv_global_coverage_flag?gv_global_coverage_nb:0))

/** @brief compute the PEARSON/SPEARMAN correlation coefficient from the 2D array
 *   (i%CLASSES) = percent %id;  (i/CLASSES) = value (number of hits or coverage)
 *
 * SPEARMAN was proposed by "Brinda Sykulski Kucherov 2015"
 *   whereas
 * PEARSON was used in "Noe Martin 2014"
 * @param polynom is the set of coefficients as a vector of pair\< pair\<first_var_count, second_var_count\>, BIGINT_count\>
 * @param correlation_enum is the index of the correlation function being used @see gv_correlation_function_index
 * @return the correlation computed accordingly
 */
double compute_correlation(vector< pair<pair<int,int>,BIGINT> > * polynom, int correlation_enum) {

  int min_p = MIN_PVAL;
  if (min_p >= gv_alignment_length) {
    _WARNING("\"-y/-g <int>\" OPTION set to alignment length - 1","(was too large before)");
    min_p = gv_alignment_length-1;
  }

  if (gp_binomial_weights_not_computed_flag) {
    generate_binomials_weights(CLASSES);
    gp_binomial_weights_not_computed_flag = false;
  }

  switch (correlation_enum) {

  case CORRELATION_FUNCTIONS_ENUM_PEARSON :
    {
      vector<BIGINT> y  = vector<BIGINT>(CLASSES,0);
      vector<BIGINT> p  = vector<BIGINT>(CLASSES,0);
      vector<BIGINT> y2 = vector<BIGINT>(CLASSES,0);
      vector<BIGINT> p2 = vector<BIGINT>(CLASSES,0);
      vector<BIGINT> yp = vector<BIGINT>(CLASSES,0);
      vector<BIGINT> n  = vector<BIGINT>(CLASSES,0);

      vector<BIGINT>  y_sum = vector<BIGINT>(CLASSES,0);
      vector<BIGINT>  p_sum = vector<BIGINT>(CLASSES,0);
      vector<BIGINT> y2_sum = vector<BIGINT>(CLASSES,0);
      vector<BIGINT> p2_sum = vector<BIGINT>(CLASSES,0);
      vector<BIGINT> yp_sum = vector<BIGINT>(CLASSES,0);
      vector<BIGINT>  n_sum = vector<BIGINT>(CLASSES,0);

      // fill initial tables
      for (unsigned i = 0; i < polynom->size(); i++) {
        int p_count = (*polynom)[i].first.first;
        int y_count = (*polynom)[i].first.second;
#ifdef USEINFINT
        BIGINT number         = (*polynom)[i].second;
        BIGINT y_count_number = number * y_count;
        BIGINT p_count_number = number * p_count;
#else
#warning "undefined USEINFINT (no infinite precision integer) : very specific functions as Pearson/Spearman correlation (or Polynomial coefficient) count on <long long> and computation on <double> may overflow ... you can compile this program with USEINFINT defined (-DUSEINFINT) but it will be much slower"
        double number         = (*polynom)[i].second;
        double y_count_number = number * y_count;
        double p_count_number = number * p_count;
#endif
         y[p_count] =  y[p_count] + y_count_number;
        y2[p_count] = y2[p_count] + y_count_number * y_count;
         p[p_count] =  p[p_count] + p_count_number;
        p2[p_count] = p2[p_count] + p_count_number * p_count;
        yp[p_count] = yp[p_count] + y_count_number * p_count;
         n[p_count] =  n[p_count] + number;
      }

      // "p"/gv_alignment_length partial sums to full sums >>
      {
        int u;
        // normalize sets by their binomial weight
        for (u = 0; u <= gv_alignment_length; u++) {
#ifdef USEINFINT
           y_sum[u] = binomial_weights[gv_alignment_length][u] *  y[u];
          y2_sum[u] = binomial_weights[gv_alignment_length][u] * y2[u];
           p_sum[u] = binomial_weights[gv_alignment_length][u] *  p[u];
          p2_sum[u] = binomial_weights[gv_alignment_length][u] * p2[u];
          yp_sum[u] = binomial_weights[gv_alignment_length][u] * yp[u];
           n_sum[u] = binomial_weights[gv_alignment_length][u] *  n[u];
#else
           y_sum[u] = (double) binomial_weights[gv_alignment_length][u] * (double)  y[u];
          y2_sum[u] = (double) binomial_weights[gv_alignment_length][u] * (double) y2[u];
           p_sum[u] = (double) binomial_weights[gv_alignment_length][u] * (double)  p[u];
          p2_sum[u] = (double) binomial_weights[gv_alignment_length][u] * (double) p2[u];
          yp_sum[u] = (double) binomial_weights[gv_alignment_length][u] * (double) yp[u];
           n_sum[u] = (double) binomial_weights[gv_alignment_length][u] * (double)  n[u];
#endif
        }

        // overlap counts
        for (u = gv_alignment_length; u > 0; u--) {
           y_sum[u-1] =  y_sum[u-1] +  y_sum[u];
          y2_sum[u-1] = y2_sum[u-1] + y2_sum[u];
           p_sum[u-1] =  p_sum[u-1] +  p_sum[u];
          p2_sum[u-1] = p2_sum[u-1] + p2_sum[u];
          yp_sum[u-1] = yp_sum[u-1] + yp_sum[u];
           n_sum[u-1] =  n_sum[u-1] +  n_sum[u];
        }
      }
      // << "p"/gv_alignment_length

      // computing and printing correlations
      {
#ifdef USEINFINT
        BIGINT covyp  =       yp_sum[min_p] * n_sum[min_p] - y_sum[min_p] * p_sum[min_p];
        BIGINT sigy2  =       y2_sum[min_p] * n_sum[min_p] - y_sum[min_p] * y_sum[min_p];
        BIGINT sigp2  =       p2_sum[min_p] * n_sum[min_p] - p_sum[min_p] * p_sum[min_p];
        BIGINT sqrt_sigy2_sigp2 = (sigy2 * sigp2).intSqrt();
        while (covyp > BIGINT(~((long long)1<<(sizeof(long long)*8-1))) || sqrt_sigy2_sigp2 > BIGINT(~((long long)1<<(sizeof(long long)*8-1)))) {
          covyp /= 2;
          sqrt_sigy2_sigp2 /= 2;
        }
        double result = (double)(covyp.toLongLong())/(double)(sqrt_sigy2_sigp2.toLongLong());
        return result;
#else
        double covyp  =       yp_sum[min_p] * n_sum[min_p] - y_sum[min_p] * p_sum[min_p];
        double sigy2  = sqrt( y2_sum[min_p] * n_sum[min_p] - y_sum[min_p] * y_sum[min_p]);
        double sigp2  = sqrt( p2_sum[min_p] * n_sum[min_p] - p_sum[min_p] * p_sum[min_p]);
        return covyp / (sigy2*sigp2);
#endif
      }
    }



  case CORRELATION_FUNCTIONS_ENUM_SPEARMAN :
    {
      // sum the number to sort in two arrays

#ifdef USEINFINT
      vector<BIGINT> y  = vector<BIGINT>(CLASSES,0);
      vector<BIGINT> p  = vector<BIGINT>(CLASSES,0);
      BIGINT full_sum = 0;
#else
      vector<double> y = vector<double>(CLASSES,0);
      vector<double> p = vector<double>(CLASSES,0);
      double full_sum = 0.0;
#endif
      for (unsigned i = 0; i < polynom->size(); i++) {
        int p_count = (*polynom)[i].first.first;
        int y_count = (*polynom)[i].first.second;
        BIGINT number = (*polynom)[i].second;
        if (p_count >= min_p) {
#ifdef USEINFINT
          p[p_count] += binomial_weights[gv_alignment_length][p_count] * number;
          y[y_count] += binomial_weights[gv_alignment_length][p_count] * number;
          full_sum   += binomial_weights[gv_alignment_length][p_count] * number;
#else
          p[p_count] += (double) binomial_weights[gv_alignment_length][p_count] * number;
          y[y_count] += (double) binomial_weights[gv_alignment_length][p_count] * number;
          full_sum   += (double) binomial_weights[gv_alignment_length][p_count] * number;
#endif
        }
      }

      // compute ex-aequo correction
#ifdef USEINFINT
      BIGINT twelve_mc_p = 0;
      BIGINT twelve_mc_y = 0;
#else
      double mc_p = 0;
      double mc_y = 0;
#endif

#ifdef USEINFINT
      for (int u = min_p; u <= gv_alignment_length; u++)
        twelve_mc_p += ((p[u]*p[u] - 1) * p[u]);
      for (int u = 0; u <= gv_alignment_length; u++)
        twelve_mc_y += ((y[u]*y[u] - 1) * y[u]);
#else
      for (int u = min_p; u <= gv_alignment_length; u++)
        mc_p += ((p[u]*p[u] - 1) * p[u]) / 12.0;
      for (int u = 0; u <= gv_alignment_length; u++)
        mc_y += ((y[u]*y[u] - 1) * y[u]) / 12.0;
#endif

      // measure rank difference
#ifdef USEINFINT
      BIGINT d_square = 0;
#else
      double d_square = 0;
#endif

      for (unsigned i = 0; i < polynom->size(); i++) {
        int p_count = (*polynom)[i].first.first;
        int y_count = (*polynom)[i].first.second;
        BIGINT number = (*polynom)[i].second;
        if (p_count >= min_p) {
#ifdef USEINFINT
          BIGINT rank_p = 0;
#else
          double rank_p = 0;
#endif
          for (int i = min_p; i < p_count; i++)
            rank_p += p[i];
          rank_p += p[p_count]/2;
#ifdef USEINFINT
          BIGINT rank_y = 0;
#else
          double rank_y = 0;
#endif
          for (int j = 0; j < y_count; j++)
            rank_y += y[j];
          rank_y += y[y_count]/2;
#ifdef USEINFINT
          d_square +=         binomial_weights[gv_alignment_length][p_count]  * number * (rank_p - rank_y) * (rank_p - rank_y);
#else
          d_square += (double)binomial_weights[gv_alignment_length][p_count]  * number * (rank_p - rank_y) * (rank_p - rank_y);
#endif
        }
      }

      // return coeff [ftp://ftp.biostat.wisc.edu/pub/chappell/800/hw/spearman.pdf]
#ifdef USEINFINT
      BIGINT n3_minus_n =  ((full_sum*full_sum - 1) * full_sum);
      BIGINT num = (n3_minus_n - (d_square*6 + twelve_mc_p/2 + twelve_mc_y/2));
      BIGINT den = ((n3_minus_n - twelve_mc_p)*(n3_minus_n - twelve_mc_y)).intSqrt();
      while (num > BIGINT(~((long long)1<<(sizeof(long long)*8-1))) || den > BIGINT(~((long long)1<<(sizeof(long long)*8-1)))) {
        num /= 2;
        den /= 2;
      }
      double result = (double)(num.toLongLong()) / (double)(den.toLongLong());
      return result;
#else
      double n3_minus_n_over_6 =  ((full_sum*full_sum - 1) * full_sum) / 6.0;
      return (n3_minus_n_over_6 - (d_square + (mc_p + mc_y))) / sqrt((n3_minus_n_over_6 - 2.0 * mc_p)*(n3_minus_n_over_6 - 2.0 * mc_y));
#endif
    }
  default:
    return 0.0;
  }
}


/** @name seedproperties functions
 *  @brief class and set of function that keep seed properties and the pareto set
 */
// @{

/** @class seedproperties
 *  @brief simple "inner class" to keep seed properties as
 *   sensitivity, selectivity ... and other properties
 */

class seedproperties {

public:
  /** @brief build a "seedproperties" object (keep current seed properties when needed)
   */
  seedproperties(double sel, double sens, double dist, std::string str, bool lossless = false,  vector< pair<pair<int,int>,BIGINT> > * polynom = NULL, polynomial<BIGINT > * multipoly = NULL) {

    this->sel   = sel;
    this->sens  = sens;
    this->dist  = dist;
    this->str   = string(str);
    this->lossless = lossless;
    if (polynom)
      this->polynom = vector< pair<pair<int,int>,BIGINT> >(polynom->begin(),polynom->end());
    else
      this->polynom = vector< pair<pair<int,int>,BIGINT> >(0);
    if (multipoly)
      this->multipoly = polynomial<BIGINT >(*multipoly);
    else
      this->multipoly = polynomial<BIGINT >();

  }

  /** @brief clone a "seedproperties" object
   */
  seedproperties(const seedproperties & other){
    this->sel   = other.sel;
    this->sens  = other.sens;
    this->dist  = other.dist;
    this->str   = string(other.str);
    this->lossless = other.lossless;
    this->polynom   = vector< pair<pair<int,int>,BIGINT> >(other.polynom.begin(),other.polynom.end());
    this->multipoly = polynomial<BIGINT >(other.multipoly);
  }


  /// seed selectivity
  double sel;
  /// seed sensitivity
  double sens;
  /// (sel, sens) euclidian distance to (1, 1)
  double dist;
  /// string representing the seed
  string str;
  /// is this seed lossless
  bool   lossless;
  /// keep polynomial factors for multihit / coverage hit  /vs/
  vector< pair<pair<int,int>,BIGINT> > polynom;
  /// keep multivariable polynomial for output only
  polynomial<BIGINT > multipoly;

  /** @brief delete a "seedproperties" object (this is just a polynom "check and erase")
   */
  ~seedproperties() {
    this->polynom.clear();
  }

  /** @brief dominance operator (done only on "p_count = 0" elements)
   */
  bool dominant(const seedproperties & other) const {
    unsigned i_this = 0, i_other = 0;
    while (i_this < this->polynom.size() && i_other < other.polynom.size()) {
      int    this_p_count = this->polynom[i_this].first.first;
      int    this_y_count = this->polynom[i_this].first.second;
      BIGINT  this_number = this->polynom[i_this].second;
      int   other_p_count = other.polynom[i_other].first.first;
      int   other_y_count = other.polynom[i_other].first.second;
      BIGINT other_number = other.polynom[i_other].second;

      if (this_y_count > 0)
        i_this++;
      else
        if (other_y_count > 0)
          i_other++;
        else
          if (this_p_count > other_p_count)
            i_other++;
          else if (this_p_count < other_p_count || this_number > other_number) {// more are "at zero" for "this" than for "other" : "this" is less sensitive than "other", so not dominant ...
            return false;
          } else {
            i_this++; i_other++;
          }
    }
    if (i_other == other.polynom.size()) {
      int   this_y_count = this->polynom[i_this].first.second;
      BIGINT this_number = this->polynom[i_this].second;

      while (i_this < this->polynom.size()) {
        if (this_y_count == 0 && this_number > 0) {// more are "at zero" for "this" than for "other" : "this" is less sensitive than "other", so not dominant ...
          return false;
        }
        i_this++;
      }
    }
    return true;
  }


  /** @brief comparison operator
   */
  friend bool operator<(const seedproperties & e1, const seedproperties & e2);

  /** @brief comparison operator
   */
  friend bool operator==(const seedproperties & e1, const seedproperties & e2);

  /** @brief output operator
   */
  friend std::ostream& operator<<(std::ostream& os, seedproperties& e);

};


/**  @brief comparison operator for a seedproperties object
 *   @param e1 is the first seed property to compare
 *   @param e2 is the second seed property to compare
 *   this function is used to sort on selectivity first,
 *   then on sensitivity and lossless property for an equal selectivity
 *   @return the comparison value (e1 < e2)
 */
bool operator<(const seedproperties & e1, const seedproperties & e2) {
  return (e1.sel < e2.sel) || ((e1.sel == e2.sel) && (!(e1.lossless) && e2.lossless)) || ((e1.sel == e2.sel) && (e1.lossless == e2.lossless) && (e1.sens < e2.sens));
}


/**  @brief comparison operator for a seedproperties object
 *   @param e1 is the first seed property to compare
 *   @param e2 is the second seed property to compare
 *   this function is used to avoid duplicates
 *   @return the comparison value (e1 == e2)
 */
bool operator==(const seedproperties & e1, const seedproperties & e2) {
  return /*(e1.sens == e2.sens) && (e1.sel == e2.sel) && (e1.lossless == e2.lossless) &&*/ ((e1.str.compare(e2.str)) == 0) /*&& (e1.polynom.size() == e2.polynom.size())*/;//FIMXE not true comparison
}

/**
 * @brief output a seedproperties object
 * @param os is the outputstream
 * @param e is the seedproperties to output
 * @return the outputsteam
 * @see seedproperties
 */
std::ostream& operator<<(std::ostream& os, seedproperties& e){
  os.setf(ios_base::fixed, ios_base::floatfield);
  os.precision(13);
  os <<  e.str << "\t" << e.sel << "\t";
  os.precision(6);
  if (gv_lossless_flag && e.lossless)
    os << "[lossless]";
  else
    os << e.sens;
  os <<  "\t" << e.dist;
  if (gv_polynomial_dominant_selection_flag) {
    os << "\t";
    for (vector< pair<pair<int,int>,BIGINT> >::iterator i = e.polynom.begin(); i != e.polynom.end(); i++) {
      os << (i->first.first) << "," << (i->first.second) << "=" << (i->second) << ";";
    }
  }
  if (gv_multipoly_file_flag) {
    os << "\t[" << e.multipoly << "]";
  }
  return os;
}


/**
 * @brief compute the euclidian distance of (sens, sel) to (1.0, 1.0)
 */

double dist(double sel, double sens) {
  return sqrt((1.0 - sel)*(1.0 - sel) + (1.0 - sens)*(1.0 - sens));
}


/**
 * @brief display a list of seeds properties
 */

int display (list<seedproperties> & l) {
  cerr << "#" << "\t" << "1th column : seed motif"
       << "\t" << "2nd column : selectivity"
       << "\t" << "3rd column : " << (gv_correlation_flag?"correlation":"sensitivity")
       << "\t" << "4th column : distance to (1, 1)";
  if (gv_polynomial_dominant_selection_flag) {
    cerr << "\t" << "5th column : nbmatches,property:count;... ";
  }

  l.sort();
  for (list<seedproperties>::iterator i = l.begin(); i != l.end(); i++) {
    if (i->str.length() > 0)
      cout << (*i) << endl;
  }
  return 0;
}



/**
 * @brief Select seedproperties that are on the "pareto set"
 * @brief usefull to keep best seeds only
 *
 * @param l is a (possibly unsorted) list of seedproperties
 * @return the area under the "pareto set"
 */

double selectPareto(list<seedproperties> & l) {

  // (1) push border seedproperties and sort
  l.sort();

  // (2) compute border hull
  list<seedproperties>::iterator i = l.begin();
  list<seedproperties>::iterator j = l.begin();
  i++;

  while(!gv_polynomial_dominant_selection_flag && i != l.end()) {

    // increased sensitivity
    if ((i->lossless && !(j->lossless)) || ((i->lossless == j->lossless) && (i->sens > j->sens))) { // note that i->sel >= j->sel
      if (j != l.begin()) {
        list<seedproperties>::iterator k = j;
        j--;
        l.erase(k);
      } else {
        l.erase(j);
        j = i;
        i++;
      }
    } else {
      // increased selectivity
      if (i->sel > j->sel) {
        j = i;
        i++;
      } else {
        // equal selectivity but decreased sensitivity
        list<seedproperties>::iterator k = i;
        i++;
        l.erase(k);
      }
    }
  } // while



  // (3) compute covered area (set of rectangles)
  double area = 0.0;
  j = l.begin();
  for ( i = l.begin(); i !=  l.end(); i++ ){
    area += (i->sel - j->sel) * (i->sens);
    j = i;
  }
  return area;
}


/**
 * @brief Compute the "pareto set" and "pareto set" area (set of rectangles)
 * @brief for the selectivity/sensitivity of a set of seeds
 * @brief [this definition is more suitable]
 *
 * @param l is a (possibly unsorted) list of seedproperties
 * @return the area under the "pareto set"
 */

double list_and_areaPareto(list<seedproperties> & l) {
  // (1) select pareto set and compute area
  double area = selectPareto(l);
  VERB_FILTER(VERBOSITY_MODERATE, INFO__("pareto area : " << area););
  // (2) display seeds
  for ( list<seedproperties>::iterator i = l.begin(); i !=  l.end(); i++ ){
    if (i->str.length() > 0)
      cout << (*i) << endl;
  }
  return area;
}


/**
 * @brief this method check if a given sel/sens is out of the pareto curve and thus does improve the overall selectivity/sensitivity
 * @param l is the list of seedproperties where e must be inserted
 * @param e is the seedproperties element that has to be inserted
 * @return a positive value if the e improve the pareto and has thus be inserted
 */

double insertPareto(list<seedproperties> & l, seedproperties & e) {
  // seach the position where to insert
  list<seedproperties>::iterator i = l.begin();
  double last_sens = 1.0;
  while(i != l.end() && i->sel + 1e-13 < e.sel) {
    last_sens = i->sens;
    i++;
  }
  // insert (front/back/between)
  if ( i ==  l.end()) {
    // "e" is the first seed, or has higher selectivity
    l.push_back(e);
    return (e.sens - last_sens);
  } else {
    // "e" has lower selectivity
    if (i->sel - 1e-13 > e.sel) {
      l.insert(i, e);
      return (e.sens - last_sens);
    } else {
      // "e" has approximately the same selectivity : (i->sel - 1e-13 <= e.sel && i->sel + 1e+13)
      bool lossless = e.lossless && (!(i->lossless));
      double dist   = e.sens - i->sens;
      if (!gv_polynomial_dominant_selection_flag && (lossless || (dist > 0))) {
        // better lossy seed, or "first" new lossless seed : replace i by e
        l.insert(i, e);
        l.erase(i);
        return dist;
      } else {
        if (gv_polynomial_dominant_selection_flag) {
          bool inserted = false;
          while (i != l.end() && i->sel - 1e-13 <= e.sel && e.sel <= i->sel + 1e-13) {
            if ((*i) == e)
              return 0;
            dist = MIN(dist, e.sens - i->sens);
            if (i->dominant(e)) {
              return dist;
            } else {
              if (e.dominant(*i)) {
                VERB_FILTER(VERBOSITY_MODERATE, MESSAGE__("\t[-] erased : " << (*i) << endl););
                i = l.erase(i);
                if (!inserted) {
                  dist = e.sens - i->sens;
                  VERB_FILTER(VERBOSITY_MODERATE, MESSAGE__("[+] inserted (after erasing seed(s)): " << (e) << endl););
                  i = l.insert(i, e);
                  i++;
                  inserted = true;
                }
              } else {
                i++;
              }
            }
          }
          if (!inserted) {
            VERB_FILTER(VERBOSITY_MODERATE, MESSAGE__("[+] inserted (no erased seed): " << (e) << endl););
            l.insert(i, e);
          }
          return dist;
        }
      }
    }
  }
  return 0;
}



/**
 * @brief read a list of initial seeds and merge them in the pareto set
 * @param l gives the list to complete the pareto part
 * @param filename gives the input file name
 * @return error code (always 0, even if the file is unreadable)
 */

int inputPareto(list<seedproperties> & l, char * filename) {
  ifstream in;
  in.open(filename, ifstream::in);
  if (in.good()) {
    string line;
    while (getline(in, line)) {
      istringstream row(line);
      string shape;
      string sens_or_lossless;
      double sel = 0.0;
      double sens = 1.0;
      double dist = 0.0;
      vector< pair<pair<int,int>,BIGINT> > polynom(0);

      // read begining
      row >> shape >> sel >> sens_or_lossless >> dist;

      // read lossless flag or sens
      bool lossless = false;
      if (!(sens_or_lossless.compare("[lossless]"))) {
        lossless =  true;
      } else {
        stringstream is(stringstream::in | stringstream::out);
        is << sens_or_lossless;
        is >> sens;
      }

      // read polynom if any
      if (gv_polynomial_dominant_selection_flag) {
        int ygm =  0,  p = 0;
        BIGINT c = 0;
        char x,y,z;
        while (row >> ygm >> x >> p >> y >> c >> z) {
          if (x != ',') {
            _ERROR("inputPareto"," missing , symbol");
          }
          if (y != '=') {
            _ERROR("inputPareto"," missing = symbol");
          }
          if (z != ';') {
            _ERROR("inputPareto"," missing ; symbol");
          }
          polynom.push_back(pair<pair<int,int>,BIGINT>(pair<int,int>(ygm,p),c));
        }
      }

      // build and insert seed property
      seedproperties sp(sel, sens, dist, shape, lossless, &polynom);
      insertPareto(l,sp);
    }
  } else {
    _WARNING("inputPareto", " unable to read  \"" << filename << "\"");
  }
  in.close();
  selectPareto(l);
  return 0;
}


/**
 * @brief write the pareto set to a file
 * @param l gives the list to write
 * @param filename gives the output file name
 * @return error code (always 0, but quit if the file in unwritable)
 */

int outputPareto(list<seedproperties> & l, char * filename) {
  // select pareto
  selectPareto(l);
  // write pareto
  ofstream out;
  out.open(filename, ifstream::out);
  if (out.good()) {
    for (list<seedproperties>::iterator i = l.begin(); i != l.end(); i++ ){
      if (i->str.length() > 0)
        out << (*i) << endl;
    }
  } else {
    _ERROR("outputPareto"," unable to write into \"" << filename << "\"");
  }
  out.close();
  return 0;
}

// @}




/**
 * @brief this method is needed when |A| or |B| size has changed
 */

void build_default_subsetseed_matching_matrix() {
  for (unsigned u = 0; u < gv_subsetseed_matching_matrix.size(); u++)
    gv_subsetseed_matching_matrix[u].clear();
  gv_subsetseed_matching_matrix.clear();

  for (int a = 0; a < gv_align_alphabet_size; a++) {
    gv_subsetseed_matching_matrix.push_back(vector<int>(gv_seed_alphabet_size, 0));
  }

  // fill it in
  for (int b = 0; b < gv_seed_alphabet_size; b++ ) {
    // sharp
    gv_subsetseed_matching_matrix[gv_align_alphabet_size-1][b] = 1;
    for (int a = b; a < gv_align_alphabet_size; a++ )
      gv_subsetseed_matching_matrix[a][b] = 1;
  }

  // set the flag to signal the '#' element
  CHECKMATCHINGMATRIX(gv_subsetseed_matching_matrix);
}




/**
 * @brief this method is needed when |A| or |B| size has changed
 */

void build_default_vectorizedsubsetseed_scoring_matrix() {
  for (unsigned u = 0; u < gv_vectorizedsubsetseed_scoring_matrix.size(); u++)
    gv_vectorizedsubsetseed_scoring_matrix[u].clear();
  gv_vectorizedsubsetseed_scoring_matrix.clear();

  for (int a = 0; a < gv_align_alphabet_size; a++) {
    gv_vectorizedsubsetseed_scoring_matrix.push_back(vector<int>(gv_seed_alphabet_size,-1));
  }

  // fill it in
  for (int a = 0; a < gv_align_alphabet_size; a++) {
    for (int b = 0; b < gv_seed_alphabet_size; b++ ) {
      if (gv_subsetseed_matching_matrix[a][b])
        gv_vectorizedsubsetseed_scoring_matrix[a][b] = 1;
    }
  }
}



/**
 * @brief this method is needed when |A| size has changed
 */

void build_default_probabilities()  {
  gv_bsens.clear();
  gv_bsens    = std::vector<double>(gv_align_alphabet_size, 0);
  gv_bsens_k  = 0;
  gv_bsel.clear();
  gv_bsel     = std::vector<double>(gv_align_alphabet_size, 0);
  gv_bsel_k   = 0;

  double remaining = 1.0;

  // bsel  :   0.5, 0.25,  0.125, 0.125,
  // bsens : 0.125, 0.125, 0.25,  0.5
  for (int i = 0; i < gv_align_alphabet_size; i++) {
    if ( i < (gv_align_alphabet_size - 1) ) {
      remaining *= .5;
    }

    gv_bsel[i] = remaining;
    gv_bsens[gv_align_alphabet_size-i-1] = remaining;
  }
}


/** @name system functions and macros
 *  @brief systems functions and macros that are windows/unix compatible,
 *  and a terminal signal handler if something goes wrong
 */
// @{

#if defined(WIN32) || defined(WIN64)
// windows part
/// SIGNAL terminal flag (disable terminal handler)
#undef SIGNAL
#include <windows.h>
/// get process id (used to init pseudo-random function)
#define GETPID() GetCurrentProcessId()
#else

// unix part
/// SIGNAL terminal flag (enable terminal handler)
#define SIGNAL
#include <unistd.h>
#include <pthread.h>
/// get process id (used to init pseudo-random function)
#define GETPID() getpid()
#endif

#ifdef SIGNAL
#include <signal.h>
#endif

#if defined(WIN32) || defined(WIN64)
/// windows getmstime (milliseconds) function (used to init pseudo-random function)
typedef unsigned __int64 QWORD;
QWORD GETMSTIME() {
  FILETIME filetime;
  GetSystemTimeAsFileTime(&filetime);
  return filetime.dwLowDateTime + filetime.dwHighDateTime * 1000000;
}
#else
/// unix getmstime (milliseconds) function (used to init pseudo-random function)
unsigned long long GETMSTIME() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_usec + tv.tv_sec * 1000000;
}
#endif



/**
 * @brief l (list of seedproperties) defined global to be saved by the "terminate handler"
 * @see termin_handler
 */

list<seedproperties> l;


#ifdef SIGNAL

/**
 * @brief this handler writes output file on premature ending
 */

void termin_handler( int signal ) {
  VERB_FILTER(VERBOSITY_MODERATE, INFO__( "Signal handler"););
  // security input
  if (gv_nb_input_filenames > 0) {
    for (int i = 0; i < gv_nb_input_filenames; i++)
      inputPareto(l, gv_input_filenames[i]);
  }

  // final output
  VERB_FILTER(VERBOSITY_LIGHT, MESSAGE__("# seeds\t sel\t " << (gv_correlation_flag?"corr":"sens") << "\t dist" << (gv_polynomial_dominant_selection_flag?"\tpol":"")););
  if (!gv_output_filename)
    list_and_areaPareto(l);
  else
    outputPareto(l, gv_output_filename);
  exit(0);
}

#endif

// @}




/// choose between a subset seed and a vectorized subset seed (when @ref gv_vectorized_flag is activated) , with lossless simplification (when @ref gv_lossless_flag is activated)
#define SEEDAUTOMATON(automaton, seed, nomerge)                                         \
  ( gv_vectorized_flag ?                                                                \
    ( gv_lossless_flag ?                                                                \
      ((automaton)->Automaton_SeedScoreCost(                                            \
                                            *(seed),                                    \
                                            gv_subsetseed_matching_matrix,              \
                                            gv_vectorizedsubsetseed_scoring_matrix,     \
                                            gv_vectorizedsubsetseed_scoring_threshold,  \
                                            nomerge,                                    \
                                            gv_lossless_costs_vector,                   \
                                            gv_lossless_cost_threshold                  \
                                                                        ))              \
      :                                                                                 \
      ((automaton)->Automaton_SeedScore(                                                \
                                        *(seed),                                        \
                                        gv_subsetseed_matching_matrix,                  \
                                        gv_vectorizedsubsetseed_scoring_matrix,         \
                                        gv_vectorizedsubsetseed_scoring_threshold,      \
                                        nomerge                                         \
                                                                        ))              \
      )                                                                                 \
    :                                                                                   \
    ( gv_lossless_flag ?                                                                \
      ((automaton)->Automaton_SeedPrefixesMatchingCost(                                 \
                                                       *(seed),                         \
                                                       gv_subsetseed_matching_matrix,   \
                                                       nomerge,                         \
                                                       gv_lossless_costs_vector,        \
                                                       gv_lossless_cost_threshold       \
                                                                        ))              \
      :                                                                                 \
      ((automaton)->Automaton_SeedPrefixesMatching(                                     \
                                                   *(seed),                             \
                                                   gv_subsetseed_matching_matrix,       \
                                                   nomerge                              \
                                                                        ))              \
      )                                                                                 \
    )



/**
 *@brief iedera main program
 *@param argc is the number of arguments (including program name)
 *@param argv is the value of all the arguments (including program name)
 *@return error code
 */

int main(int argc, char * argv[]) {
  int nbruns = 0;

  build_default_subsetseed_matching_matrix();
  build_default_vectorizedsubsetseed_scoring_matrix();
  build_default_probabilities();
  computeWeight();
  gv_seeds  = std::vector<seed*>(1);
  gv_xseeds = std::vector<seed*>(0);

  SCANARG(argc, argv);

  l = list<seedproperties>();
  vector< pair<pair<int,int>,BIGINT> > v1_one = vector< pair<pair<int,int>,BIGINT> >(1,pair<pair<int,int>,BIGINT>(pair<int,int>(0,1),1));
  vector< pair<pair<int,int>,BIGINT> > v0_one = vector< pair<pair<int,int>,BIGINT> >(1,pair<pair<int,int>,BIGINT>(pair<int,int>(0,0),1));
  l.push_front(seedproperties(0.0, 1.0, 1.0, string(""), gv_lossless_flag, &v1_one));
  l.push_back( seedproperties(1.0, 0.0, 1.0, string(""), false,            &v0_one));

  // randomize
  srand(GETMSTIME()*GETPID());

  // initial pareto set
  if (gv_nb_input_filenames > 0) {
    for (int i = 0; i < gv_nb_input_filenames; i++)
      inputPareto(l, gv_input_filenames[i]);
  }

#ifdef SIGNAL
  signal(SIGTERM, &termin_handler );
  signal(SIGSTOP, &termin_handler );
#endif


  //
  // (O) DISPLAY PARAMETERS
  //

  VERB_FILTER(VERBOSITY_MODERATE, INFO__("* alignment alphabet size |A| : " << gv_align_alphabet_size << " , " << "seed alphabet size |B| : "  << gv_seed_alphabet_size ););
  VERB_FILTER(VERBOSITY_MODERATE, INFO__("* subset seed matching matrix : "; DISPLAYMATRIX(cerr, gv_subsetseed_matching_matrix);););
  if (gv_vectorized_flag) {
    VERB_FILTER(VERBOSITY_MODERATE, INFO__("* subset seed scoring matrix : "; DISPLAYMATRIX(cerr, gv_vectorizedsubsetseed_scoring_matrix);););
    VERB_FILTER(VERBOSITY_MODERATE, INFO__("* subset seed scoring threshold : " << gv_vectorizedsubsetseed_scoring_threshold););
  }
  VERB_FILTER(VERBOSITY_MODERATE, INFO__("* seed alphabet weights "; DISPLAYTABLE(cerr, gv_bsel_weight)););

  // build the foreground/background models
  automaton<double> a_sel;
  if (gv_bsel_k == 0) {
    a_sel.Automaton_Bernoulli(gv_bsel);
    VERB_FILTER(VERBOSITY_MODERATE, INFO__("* background Bernoulli model "; DISPLAYTABLE(cerr, gv_bsel)););
  } else {
    a_sel.Automaton_Markov(gv_bsel, gv_bsel_k);
    VERB_FILTER(VERBOSITY_MODERATE, INFO__("* background Markov model "; DISPLAYTABLE(cerr, gv_bsel)););
    VERB_FILTER(VERBOSITY_MODERATE, INFO__(" (order : " << (gv_bsel_k) << ")" ););
  }

  automaton<double> a_sens;
  if (!gv_bsens_automaton) {
    if (gv_bsens_k == 0) {
      a_sens.Automaton_Bernoulli(gv_bsens);
      VERB_FILTER(VERBOSITY_MODERATE, INFO__("* foreground Bernoulli model "; DISPLAYTABLE(cerr, gv_bsens);););
    } else {
      a_sens.Automaton_Markov(gv_bsens, gv_bsens_k);
      VERB_FILTER(VERBOSITY_MODERATE, INFO__("* foreground Markov model "; DISPLAYTABLE(cerr, gv_bsens)););
      VERB_FILTER(VERBOSITY_MODERATE, INFO__( " (order : " << (gv_bsens_k) << ")" ););
    }
  } else {
    VERB_FILTER(VERBOSITY_MODERATE, INFO__("* foreground Ad hoc model " ););
    a_sens = *gv_bsens_automaton;
  }

  // build the lossless automaton
  automaton<void> a_lossless;
  double    lossless_set_sens = 1.0;
  if (gv_lossless_flag) {
    a_lossless.Automaton_Lossless(gv_lossless_costs_vector, gv_lossless_cost_threshold);
    matrix<double> *    m_lossless_pr_sens = a_lossless.matrix_pr_product(a_sens, PRODUCT_UNION_NO_FINAL_LOOP, gv_alignment_length);
    lossless_set_sens = m_lossless_pr_sens->Pr(gv_alignment_length, false);
    delete              m_lossless_pr_sens;
    VERB_FILTER(VERBOSITY_MODERATE, INFO__(" - pr_div_lossless_set : " << lossless_set_sens ););
  }


  // build the cost automaton
  automaton<cost<int> > a_cost;
  if (gv_lossless_flag) {
    vector<cost<int> > lossless_costs_vector(gv_seed_alphabet_size,cost<int>(0));
    for (int i = 0; i < gv_seed_alphabet_size; i++)
      lossless_costs_vector[i] = cost<int>(gv_lossless_costs_vector[i]);
    a_cost.Automaton_Bernoulli(lossless_costs_vector);
    //[FIXME]
  }

  //
  // build the homogeneous automaton if requested and
  // compute the foreground probability of homogeneous paths
  //

  double pr_div_homogeneous = 1.00;
  automaton<void> * a_homogeneous = NULL;
  if (gv_homogeneous_flag) {
    a_homogeneous = new automaton<void>();
    VERB_FILTER(VERBOSITY_MODERATE, INFO__("* Homogeneous automaton : {";
                                           for (int a = 0; a < gv_align_alphabet_size; a++) {
                                             if  (a>0)  cerr << ",";
                                             cerr << (gv_homogeneous_scores[a]);
                                           }
                                           cerr << "}";
                                           ););

    a_homogeneous->Automaton_Homogeneous(gv_homogeneous_scores, gv_alignment_length);
    VERB_FILTER(VERBOSITY_ANNOYING, INFO__("   - size : " << (a_homogeneous->size())););
    if (gv_minimize_flag) {
      automaton<void> * na = a_homogeneous->Hopcroft();
      delete a_homogeneous;
      a_homogeneous = na;
      VERB_FILTER(VERBOSITY_ANNOYING, INFO__("   - reduced size : " << (a_homogeneous->size())););
    }
    automaton<double> * a_product_homogeneous_sens = a_homogeneous->product(a_sens, PRODUCT_BUTNOT_FINAL_LOOP, gv_alignment_length);

    if (gv_lossless_flag) {
      pr_div_homogeneous = a_product_homogeneous_sens->PrLossless(gv_alignment_length, gv_lossless_costs_vector, gv_lossless_cost_threshold);
    } else {
      pr_div_homogeneous = a_product_homogeneous_sens->Pr(gv_alignment_length);
    }

    delete a_product_homogeneous_sens;
    VERB_FILTER(VERBOSITY_MODERATE, INFO__(" - pr_div_homogenous : " << pr_div_homogeneous););
  }


  //
  // build the excluded seed automaton if requested
  // compute the foreground probability of excluded paths
  //

  double pr_div_excluded = 0.00;
  automaton<void> * a_excluded = NULL;
  if (gv_xseeds.size()) {

    VERB_FILTER(VERBOSITY_MODERATE, INFO__("* Excluded automaton : ";
                                           for (unsigned i = 0; i < gv_xseeds.size(); i++) {
                                             if  (i>0)  cerr << ",";
                                             cerr << (gv_xseeds[i])->str();
                                           }
                                           ););

    // a) build the excluded seed automaton
    for (unsigned i = 0; i < gv_xseeds.size(); i++){

      automaton<void> * a_se = new automaton<void>();
      VERB_FILTER(VERBOSITY_MODERATE, INFO__(" = seed : " << (*gv_xseeds[i])););
      SEEDAUTOMATON(a_se, gv_xseeds[i],  gv_xseeds[i]->cycled() || gv_xseeds_multihit_flag );
      VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  - size : " << (a_se->size())););
      if (gv_minimize_flag) {
        automaton<void> * na = a_se->Hopcroft();
        delete a_se;
        a_se = na;
        VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  - reduced size : " << (a_se->size())););
      }



      // excluded seed cycle
      if (gv_xseeds[i]->cycled()) {
        automaton<void> * na = a_se;
        automaton<void> * a_cycle = new automaton<void>();
        a_cycle->Automaton_Cycle(gv_xseeds[i]->maxpos(), gv_xseeds[i]->pos(), gv_xseeds[i]->nbpos());
        a_se = (a_se)->product(*a_cycle, gv_xseeds_multihit_flag?PRODUCT_INTERSECTION_NO_FINAL_LOOP:PRODUCT_INTERSECTION_FINAL_LOOP, gv_alignment_length);
        delete a_cycle;
        delete na;
        VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  - cycled size : " << (a_se->size())););
        if (gv_minimize_flag) {
          automaton<void> * na = a_se->Hopcroft();
          delete a_se;
          a_se = na;
          VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  - reduced cycled size : " << (a_se->size())););
        }
      }

      // excluded automaton
      if (i > 0) {
        automaton<void> * a_excluded_temp = a_excluded;
        a_excluded = a_excluded->product(*a_se,  gv_xseeds_multihit_flag?PRODUCT_UNION_NO_FINAL_LOOP_ADD:PRODUCT_UNION_FINAL_LOOP, gv_alignment_length);
        VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" = excluded product size : " << (a_excluded->size())););
        if (gv_minimize_flag) {
          automaton<void> * na = a_excluded->Hopcroft();
          delete a_excluded;
          a_excluded = na;
          VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" - reduced excluded product size : " << (a_excluded->size())););
        }
        delete a_excluded_temp;
        delete a_se;
      } else {
        a_excluded = a_se;
        VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" = excluded size : " << (a_excluded->size())););
      }
    }// for each excluded seed


    // excluded automaton multihit
    if (gv_xseeds_multihit_flag) {
      automaton<void> * a_excluded_temp = a_excluded;
      a_excluded = a_excluded->mhit(gv_xseeds_multihit_nb, gv_alignment_length);
      delete a_excluded_temp;
      VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" = mhits excluded size : " << (a_excluded->size())););
      if (gv_minimize_flag) {
        automaton<void> * na = a_excluded->Hopcroft();
        delete a_excluded;
        a_excluded = na;
        VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" - reduced mhits excluded size : " << (a_excluded->size())););
      }
    }




    // b) compute the foreground probability
    if (gv_homogeneous_flag) {
      automaton<void> * na = a_excluded->product(*a_homogeneous, PRODUCT_INTERSECTION_FINAL_LOOP, gv_alignment_length);
      delete a_excluded;
      a_excluded = na;
      VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" = mhits excluded x homogeneous size : " << (a_excluded->size())););
      if (gv_minimize_flag) {
        automaton<void> * na = a_excluded->Hopcroft();
        delete a_excluded;
        a_excluded = na;
        VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" - reduced mhits excluded x homogeneous size : " << (a_excluded->size())););
      }
    }


    automaton<double> * a_xpr_sens = a_excluded->product(a_sens, PRODUCT_UNION_NO_FINAL_LOOP, gv_alignment_length);
    if (gv_lossless_flag) {
      pr_div_excluded = a_xpr_sens->PrLossless(gv_alignment_length, gv_lossless_costs_vector, gv_lossless_cost_threshold);
    } else {
      pr_div_excluded = a_xpr_sens->Pr(gv_alignment_length);
    }
    delete      a_xpr_sens;
    if (pr_div_excluded == 1.0) {
      _ERROR("EXCLUDED SEEDS"," the seed(s) provided by the \"-mx\" option reaches a full sensitivity of 1.0");
    }
    VERB_FILTER(VERBOSITY_MODERATE, INFO__(" - pr_div_excluded : " << pr_div_excluded););
  }

  // check if some options are not compatible with the current sampling (cannot reimplement everything)
#ifdef SAMPLING_STRATEGY_MF
  if (gv_subalignment_flag)         {_ERROR("SAMPLING_STRATEGY_MF : subalignment functions are not supported by sampling","<to be implemented>");};
  if (gv_global_coverage_flag)      {_ERROR("SAMPLING_STRATEGY_MF : coverage sampling functions are not implemented","<to be implemented>");};
  if (gv_lossless_flag)             {_ERROR("SAMPLING_STRATEGY_MF : lossless functions are supported by sampling","<to be implemented>");};
  if (gv_xseeds.size())             {_ERROR("SAMPLING_STRATEGY_MF : excluded seeds are not supported by sampling","<to be implemented>");};
  if (gv_homogeneous_flag)          {_ERROR("SAMPLING_STRATEGY_MF : homogeneous aligments are not supported by sampling","<to be implemented>");};
  if (gv_vectorized_flag)           {_ERROR("SAMPLING_STRATEGY_MF : vectorized subset seeds are not supported by sampling","<to be implemented>");};
#endif

  VERB_FILTER(VERBOSITY_MODERATE, INFO__("** Main : "););

  //
  // (A) compute for all sets of seeds ...
  //

  bool   hillclimbing_flag      = false;
  std::vector<double>      sel  = std::vector<double>    (gv_seeds.size(),0.0);
  std::vector<automaton<void> *> a_s  = std::vector<automaton<void>*>(gv_seeds.size(),NULL);
  double hillclimbing_threshold = 1e-6;

  int    seed_to_hillclimbing      =                                   rand()%gv_seeds.size();
  int    last_seed_to_hillclimbing = (seed_to_hillclimbing+gv_seeds.size()-1)%gv_seeds.size();

  double sensitivity_threshold = 0.00;
  double sens_hillclimbing    = selectPareto(l);
#ifdef SAMPLING_STRATEGY_MF
  double max_sens             = selectPareto(l);
#endif
  double sens = 0.0;
  double selp = 0.0;

#ifdef KEEP_PRODUCT_MF
  std::vector<automaton<void>*> a_s_product(gv_seeds.size(), NULL);
  std::vector<int>        a_s_product_seed(gv_seeds.size(), -1);
#endif

  if (!gv_motif_flag) {
    // init the seeds and set the cycles
    for (unsigned i = 0; i < gv_seeds.size(); i++) {
      gv_seeds[i] = new seed();
      if (gv_cycles_flag)
        gv_seeds[i]->setCyclePos(gv_cycles_pos_list[i], gv_cycles[i]);
      if (gv_nbruns)
        (gv_seeds[i])->random();
    }
    // if not acceptable : generate a new set of seeds
    for (unsigned i = 0; i < gv_seeds.size(); i++)
      if (!(gv_seeds[i]->acceptable()))
        goto new_seeds;

    // discard twice same seed in the same set
    for (unsigned i = 0; i < gv_seeds.size(); i++)
      if (!gv_seeds[i]->cycled())
        for (unsigned j = 0; j < i; j++)
          if (!gv_seeds[j]->cycled())
            if (gv_seeds[i]->equal(gv_seeds[j]))
              goto new_seeds;
  }

  if (gv_motif_flag && gv_nbruns == 0 && gv_signature_shuffle_from_m_pattern_flag) {
    for (unsigned i = 0; i < gv_seeds.size(); i++) {
      gv_seeds[i]->reorder(0,/*b*/gv_seed_alphabet_size-1);
      while(!gv_seeds[i]->acceptable())
        if (!gv_seeds[i]->next())
          goto end_loop;
    }
  }


  while (1) {

    // for all seeds
#ifdef SAMPLING_STRATEGY_MF

    // Martin part >>
    {//(0)>>
      /* this part estimates the sensitivity by sampling before going to full computation if any improvement can be detected :
       * - for full enumeration, it can miss the best seed / set of seed
       * - for random enumeration, it can also miss some good seeds
       * - for random enumeration with hillclimbing, the hillclimbing can stop early due to the same reason
       */

#define TESTSET 1000
#define RUNSET  1000

#define AFFINE(p,q,t) ((p)*(1.0-(t)) + ((q)-(p))*(t))
#define COEF_9999_9(t) (AFFINE(3.29,1.64,t))
      int c = 0;
      for (int i = 0; i < TESTSET; i++)  {
        /* (0.1) do the sampling */
        for (int j = 0; j < RUNSET; j++)  {
          std::vector<int> alignment(gv_alignment_length,0);
          a_sens.GenSeq(alignment);
          if (gv_multihit_flag) {
            int t = 0;
            for (unsigned u = 0; u < gv_seeds.size(); u++) {
              t += gv_seeds[u]->mHits(alignment,gv_subsetseed_matching_matrix);
            }
            if (t >= gv_multihit_nb) {
              c++;
              break;
            }
          } else {
            for (unsigned u = 0; u < gv_seeds.size(); u++) {
              if (gv_seeds[u]->Hit(alignment,gv_subsetseed_matching_matrix) >= 0) {
                c++;
                break;
              }
            }
          }
        } // RUNSET
        /* (0.2) : [1/3] are we in hillclimbing mode activated and ongoing ? */
        if (hillclimbing_flag) {
          if (c < sens_hillclimbing * (i+1) * RUNSET - COEF_9999_9((double)i/RUNSET) * sqrt(sens_hillclimbing * (1 - sens_hillclimbing) * (i+1) * RUNSET)) {
            /* hillclimbing seed estimated bad ... so generate the next one */
            VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("{hx " << ((double)c/ ((i+1) * RUNSET)) << " after " << ((i+1) * RUNSET) << " runs)"););
            /* fake "sens" and "sensitivity_threshold" to update then go to 8 */
            sens = ((double)c/ ((i+1) * RUNSET));
            sensitivity_threshold = sens - max_sens;
            goto next_progress_seed;
          } else {
            /* hillclimbing seed estimated promissing ... so estimate its true probablity */
            if (c  > sens_hillclimbing * (i+1) * RUNSET + COEF_9999_9((double)i/RUNSET) * sqrt(sens_hillclimbing * (1 - sens_hillclimbing) * (i+1) * RUNSET)) {
              VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("{h? " << ((double)c/ ((i+1) * RUNSET)) << " after " << ((i+1) * RUNSET) << " runs)"););
              break;
            }
          }
        } else {
          /* (0.2) : [2/3] or [3/3] : we are in random generation mode, either in hillclimbing (but not activated), or in full one */
          if (gv_hillclimbing_flag) {
            /* (0.2) : [2/3] we are in random generation mode in hillclimbing (but not activated) */
            if (c < (max_sens - hillclimbing_threshold) * (i+1) * RUNSET - COEF_9999_9((double)i/RUNSET) * sqrt((max_sens - hillclimbing_threshold) * (1 - (max_sens - hillclimbing_threshold)) * (i+1) * RUNSET)) {
              /* random seed estimated bad ... so generate a random one */
              VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("{rx " << ((double)c/ ((i+1) * RUNSET)) << " after " << ((i+1) * RUNSET) << " runs)"););
              /* fake "sens" and "sensitivity_threshold" to update then go to 8 */
              sens = ((double)c/ ((i+1) * RUNSET));
              sensitivity_threshold =  sens - max_sens;
              goto next_progress_seed;
            } else {
              /* hillclimbing seed estimated promissing ... so estimate its true probablity */
              if (c  > (max_sens - hillclimbing_threshold) * (i+1) * RUNSET + COEF_9999_9((double)i/RUNSET) * sqrt((max_sens - hillclimbing_threshold) * (1 - (max_sens - hillclimbing_threshold)) * (i+1) * RUNSET)) {
                VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("{r? " << ((double)c/ ((i+1) * RUNSET)) << " after " << ((i+1) * RUNSET) << " runs)"););
                break;
              }
            }
          } else {
            /* (0.2) : [3/3] we are in random generation mode full one */
            if (c < max_sens * (i+1) * RUNSET - COEF_9999_9((double)i/RUNSET) * sqrt(max_sens * (1 - max_sens) * (i+1) * RUNSET)) {
              /* random seed estimated bad ... so generate a random one */
              VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("{rx " << ((double)c/ ((i+1) * RUNSET)) << " after " << ((i+1) * RUNSET) << " runs)"););
              /* fake "sens" to update then go to 9 */
              sens = ((double)c/ ((i+1) * RUNSET));
              goto next_progress_seed;
            } else {
              /* hillclimbing seed estimated promissing ... so estimate its true probablity */
              if (c  > max_sens * (i+1) * RUNSET + COEF_9999_9((double)i/RUNSET) * sqrt(max_sens * (1 - max_sens) * (i+1) * RUNSET)) {
                VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("{r? " << ((double)c/ ((i+1) * RUNSET)) << " after " << ((i+1) * RUNSET) << " runs)"););
                break;
              }
            }
          }
        }
      } //TESTSET
      VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("going to true sensitivity "););
    } //<<(0)
    //<< Martin part



#endif



    //(1-6)>>
    {
      // (1) build seed automata
      for (unsigned i = 0; i < gv_seeds.size(); i++) {
        VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" = seed : " << (*gv_seeds[i])););
        if (gv_global_coverage_flag || gv_covariance_flag) {//FIXMECOV
          sel[i] = gv_seeds[i]->selectivityFromWeight();
        } else {
          if (!a_s[i]) {
#ifdef KEEP_PRODUCT_MF
            for (unsigned u = 0; u < gv_seeds.size(); u++) {
              if (a_s_product_seed[u] == (int) i) {
                if (u < 8 * gv_seeds.size() / 10) {
                  u = 0;
                  VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  * cleaning all the already computed products ( < 80%% reuse criterion)" ););
                } else {
                  VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  * cleaning already computed products in range [ " << (u+1) << " .. " << (gv_seeds.size()) << " ]" ););
                }
                for (unsigned v = u; v < gv_seeds.size(); v++) {
                  if (v>0)
                    delete a_s_product[v];
                  a_s_product[v] = NULL;
                  a_s_product_seed[v] = -1;
                }
                break;
              }
            }
#endif
            a_s[i] = new automaton<void>();
            SEEDAUTOMATON(a_s[i], gv_seeds[i], gv_seeds[i]->cycled() || gv_multihit_flag);

            VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  - automaton size : " << (a_s[i]->size())););

            if (gv_minimize_flag) {
              automaton<void> * na = a_s[i]->Hopcroft();
              delete a_s[i];
              a_s[i] = na;
              VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  - automaton reduced size : " << (a_s[i]->size())););
            }

            // compute partial selectivity
            if (gv_lossless_flag) {
              sel[i] = gv_seeds[i]->selectivityFromWeight();
            } else {
              automaton<double> * a_pr_s_sel = (a_s[i])->product(a_sel, PRODUCT_UNION_FINAL_LOOP, gv_alignment_length);
              sel[i]                 = a_pr_s_sel->Pr(gv_seeds[i]->span());
              delete a_pr_s_sel;
            }
            if (gv_seeds[i]->cycled()) {
              sel[i] *= (double)gv_seeds[i]->nbpos()/gv_seeds[i]->maxpos();
            }

            // cycle
            if (gv_seeds[i]->cycled()) {
              automaton<void> * na = a_s[i];
              automaton<void> * a_cycle = new automaton<void>();
              a_cycle->Automaton_Cycle(gv_seeds[i]->maxpos(), gv_seeds[i]->pos(), gv_seeds[i]->nbpos());
              a_s[i] = (a_s[i])->product(*a_cycle, gv_multihit_flag?PRODUCT_INTERSECTION_NO_FINAL_LOOP:PRODUCT_INTERSECTION_FINAL_LOOP, gv_alignment_length);
              VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  - automaton cycled size : " << (a_s[i]->size())););

              delete a_cycle;
              delete na;
              if (gv_minimize_flag) {
                automaton<void> * na = a_s[i]->Hopcroft();
                delete a_s[i];
                a_s[i] = na;
                VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  - automaton reduced cycled size : " << (a_s[i]->size())););
              }
            }
          } // "if (!a_s[i])"
        } // "if (gv_global_coverage_size)"
      } // for each seed

      selp = 0.0;
      for (unsigned i = 0; i < gv_seeds.size(); i++) {
        selp += sel[i];
      } // for each seed


      // (2) compute sensitivity (and lossless property when needed ...)
      int lossless = 0;
      automaton<void> * a_spr = NULL;

      // FIXMECOV >>
      if (gv_covariance_flag) {
        for (unsigned i = 0; i < gv_seeds.size(); i++) {
          if (!a_s[i]) {
            a_s[i] = new automaton<void>();
            a_s[i]->Automaton_SeedLinearMatching(*(gv_seeds[i]),gv_subsetseed_matching_matrix);
          }
        }
        goto gv_covariance_flag_0;
      }
      // FIXMECOV<<


      if (gv_global_coverage_flag) {

        // (2.1) global coverage constraint
        a_spr = new automaton<void>();
        a_spr->Automaton_SeedPrefixesMatching_CoverageDM(gv_seeds,
                                                         gv_global_coverage_cost,
                                                         gv_subsetseed_matching_matrix);
        // Donald part >>
        /*
         * { std::ofstream out("_moore_automaton.dot"); a_spr->dot(out); out.close();}
         * { std::ofstream out("_moore_automaton_-_mealy_like_conversion.gapfr"); a_spr->gapFR(out); out.close();}
         */
        // <<
        VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" = automaton global coverage size : " << (a_spr->size())););
        if (gv_minimize_flag) {
          automaton<void> * na = a_spr->Hopcroft();
          delete a_spr;
          a_spr = na;
          VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" - automaton reduced global coverage size : " << (a_spr->size())););
        }
        // Donald part >>
        /*
         * { std::ofstream out("_moore_automaton_reduced.dot"); a_spr->dot(out); out.close();}
         * { std::ofstream out("_moore_automaton_reduced_-_mealy_like_conversion.gapfr"); a_spr->gapFR(out); out.close();}
         */
        // <<

      } else { // "if (gv_global_coverage_flag)"

        // (2.2) classical hit/multihit constraint

#ifdef KEEP_PRODUCT_MF
        // get the last product index
        int last_product_index = -1;
        for (unsigned u = 0; u < gv_seeds.size(); u++) {
          if (a_s_product_seed[u] >= 0) {
            last_product_index = u;
            a_spr = a_s_product[u];
            VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" * multiseed [";
                                                   for (int s = 0; s <= last_product_index; s++) {
                                                     if (s)  cerr << " x "; cerr << (a_s_product_seed[s]+1);
                                                   }
                                                   cerr << "] preloaded product size : " << (a_spr->size());
                                                   );
                        );
          }
        }

        // do the product for non-found seeds either from nothing or from the last product index
        int seed_index_start;
        if (gv_hillclimbing_flag)
          seed_index_start = (seed_to_hillclimbing + gv_seeds.size() - 1) % gv_seeds.size();
        else
          seed_index_start = gv_seeds.size()-1;
        for (unsigned i = seed_index_start , j = 0; j < gv_seeds.size(); i += gv_seeds.size() - 1 , i %= gv_seeds.size(), j++) {
          // search the seed
          for (int v = 0; v <= last_product_index; v++) {
            if (a_s_product_seed[v] == (int) i)
              goto seed_found;
          }
          // seed not found : do the product and store it
          if (last_product_index >= 0) {
            a_spr = a_spr->product(*(a_s[i]), gv_multihit_flag?PRODUCT_UNION_NO_FINAL_LOOP_ADD:PRODUCT_UNION_FINAL_LOOP, gv_alignment_length);
            VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" + multiseed [";
                                                   for (int s = 0; s <= last_product_index; s++) {
                                                     if (s)  cerr << " x "; cerr << (a_s_product_seed[s]+1);
                                                   }
                                                   cerr << "] x [" << (i+1) << "] product size : " << (a_spr->size());
                                                   );
                        );

            if (gv_minimize_flag) {
              automaton<void> * na = a_spr->Hopcroft();
              delete a_spr;
              a_spr = na;
              VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" + reduced multiseed [";
                                                     for (int s = 0; s <= last_product_index; s++) {
                                                       if (s)  cerr << " x "; cerr << (a_s_product_seed[s]+1);
                                                     }
                                                     cerr << "] x [" << (i+1) << "] product size : " << (a_spr->size());
                                                     );
                          );

            }
          } else {
            a_spr = a_s[i];
          }
          last_product_index++;
          a_s_product[last_product_index] = a_spr;
          a_s_product_seed[last_product_index] = i;
        seed_found:;

        }
#else
        for (unsigned i = 0; i < gv_seeds.size(); i++) {
          if (i != 0) {
            automaton<void> * a_spr_temp = a_spr;
            a_spr = a_spr->product(*(a_s[i]), gv_multihit_flag?PRODUCT_UNION_NO_FINAL_LOOP_ADD:PRODUCT_UNION_FINAL_LOOP, gv_alignment_length);

            VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" = multiseed 1.." << (i) << " x " << (i+1) << " product size : " << (a_spr->size())););

            if (gv_minimize_flag) {
              automaton<void> * na = a_spr->Hopcroft();
              delete a_spr;
              a_spr = na;
              VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" - reduced multiseed 1.." << (i) << " x " << (i+1) << " product size : " << (a_spr->size())););
            }

            if (i > 1)
              delete a_spr_temp;
          } else {
            a_spr = (a_s[0]);
          }
        } // for each seed
#endif
      } // "else" part of "if (gv_global_coverage_flag)"

      //FIXMECOV>>
    gv_covariance_flag_0:
      //<<

      // multihit (xor) global_coverage
      automaton<void> * a_spr_mhits_or_gcov_res = a_spr;
      if (gv_multihit_flag || gv_global_coverage_flag) {

        a_spr_mhits_or_gcov_res = a_spr->mhit(gv_multihit_flag?gv_multihit_nb:gv_global_coverage_nb, gv_alignment_length);

        VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" = mhits/gcov size : " << (a_spr_mhits_or_gcov_res->size())););

        if (gv_minimize_flag) {
          automaton<void> * na = a_spr_mhits_or_gcov_res->Hopcroft();
          delete a_spr_mhits_or_gcov_res;
          a_spr_mhits_or_gcov_res = na;
          VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" - reduced mhits/gcov size : " << (a_spr_mhits_or_gcov_res->size())););
        }
      }

      // homogeneous
      automaton<void> * a_spr_h_res = a_spr_mhits_or_gcov_res;
      if (gv_homogeneous_flag) {
        a_spr_h_res = a_spr_mhits_or_gcov_res->product(*a_homogeneous, PRODUCT_INTERSECTION_FINAL_LOOP, gv_alignment_length);

        VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" = homogeneous product size : " << (a_spr_h_res->size())););

        if (gv_minimize_flag) {
          automaton<void> * na = a_spr_h_res->Hopcroft();
          delete a_spr_h_res;
          a_spr_h_res = na;
          VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" - reduced homogeneous product size : " << (a_spr_h_res->size())););
        }
      }

      // excluded seeds
      automaton<void> * a_spr_mx_h_res = a_spr_h_res;
      if (gv_xseeds.size()) {
        a_spr_mx_h_res = a_spr_h_res->product(*a_excluded, PRODUCT_BUTNOT_NO_FINAL_LOOP, gv_alignment_length);

        VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" = mx product size : " << (a_spr_mx_h_res->size())););

        if (gv_minimize_flag) {
          automaton<void> * na = a_spr_mx_h_res->Hopcroft();
          delete a_spr_mx_h_res;
          a_spr_mx_h_res = na;
          VERB_FILTER(VERBOSITY_ANNOYING, INFO__(" - reduced mx product size : " << (a_spr_mx_h_res->size())););
        }
      }

      std::vector< pair<pair<int,int>,BIGINT> > *  polynom   = NULL;
      polynomial<BIGINT > * multipoly = NULL;
      //FIXMECOV>>
      if (gv_covariance_flag) goto gv_covariance_flag_1;
      //<<

      // precompute the polynom when needed (for correlation computation, or for mere polynomial output)
      if (gv_correlation_flag || gv_polynomial_dominant_selection_flag) {

        // use percentage of identity Automaton to count "CLASSES"
        automaton<BIGINT> a_sens_count_ones = automaton<BIGINT>();
        a_sens_count_ones.Automaton_CountAlphabetSymbols();

        // use a "COUNT" Semi-Ring mixed with the "CLASS" of Percentage
        struct MyAFF { static int funaff (int a, int b) { return (gv_multihit_flag||gv_global_coverage_flag) ? (a*CLASSES + b) : ((a>0)?(CLASSES+b):(b)); } };

        matrix<BIGINT>   * m_ct_sens_dist = a_spr->matrix_product(a_sens_count_ones, PRODUCT_ADDHOC_NO_FINAL_LOOP, gv_alignment_length, (AddHoc_Final_Func) MyAFF::funaff);
        VERB_FILTER(VERBOSITY_MODERATE, INFO__("- count matrix product size : " << (m_ct_sens_dist->size())););
        std::vector<BIGINT> * v_ct_sens_dist = m_ct_sens_dist->Pr_transitive_final(gv_alignment_length, (gv_multihit_flag||gv_global_coverage_flag)?(INT_INFINITY):(2*CLASSES-1),(gv_multihit_flag||gv_global_coverage_flag)?(1):(CLASSES));
        // transfor as a polynom (64 bits numbers !! warning on overflow)
        polynom = new vector< pair<pair<int,int>,BIGINT> >(0);
        for (unsigned i = 0; i < v_ct_sens_dist->size(); i++) {
          BIGINT number = (*v_ct_sens_dist)[i];
          if (number > 0) {
            int p_count = (i % CLASSES), y_count = (i / CLASSES);
            polynom->push_back(pair<pair<int,int>,BIGINT>(pair<int,int>(p_count,y_count),(BIGINT)number));
            //cerr << "[" << i << "] -> (" << p_count << "," << y_count << ") : " << number << endl;
          }
#ifndef USEINFINT
          else {
            if (number < 0) {
              _ERROR("this binary has been compiled with undefined USEINFINT (no infinite precision integer)","(gv_correlation_flag || gv_polynomial_dominant_selection_flag) p_count/y_count DOES OVERFLOW ...");
            }
          }
#endif
        }
        delete m_ct_sens_dist;
        delete v_ct_sens_dist;
      }



      // (3) Computation (correlation // sensitivity (-l/-ll) // losslessness)

      // (3.x) correlation mode computation (no probability or lossless : all alignments are considered
      if (gv_correlation_flag) {
        // (3.x.3) compute correlation value between percent of identity and "y" single/multi-hit/coverage value for each "CLASS"
        sens = compute_correlation(polynom,gv_correlation_function_index);

      } else {

        // (3.a) with matrix
        if (gv_lossless_flag) {
          // (3.a.1) lossless with subalignment matrix
          if (gv_subalignment_flag) {

            std::vector< matrix<cost<int> > * > * vm = a_spr_mx_h_res->matrices_step_cost_product(a_cost, PRODUCT_UNION_NO_FINAL_LOOP, gv_alignment_length);
            std::vector<int> results(gv_alignment_length-gv_subalignment_length+1);
#ifndef NOSLICER
            // Sliced version
            // >>
            matrices_slicer< cost<int> > vm_slicer(vm);

            // set the sub-alignment length
            for (int i = 1; i < gv_subalignment_length; i++) vm_slicer.add_right();

            // move along the sub-alignment windows and compute Pr
            for (int i = gv_subalignment_length; i <= gv_alignment_length; i++) {
              cost<int> c = vm_slicer.current_Pr(false);
              results[i - gv_subalignment_length] = (c > cost<int>(gv_lossless_cost_threshold)) ? 1 : 0;
              vm_slicer.add_right();
              vm_slicer.del_left();
            }
            // <<
#else
            // Non sliced version
            // >>
            for (int i = 0; i <= gv_alignment_length - gv_subalignment_length; i++) {
              matrix< cost<int> > * m = (*vm)[i];
              matrix< cost<int> > * m_old = NULL;
              m->resize(2);
              for (int k = 1; k < gv_subalignment_length; k++) {
                m = m->Compose(*((*vm)[i+k]));
                if (m_old) delete m_old;
                m_old = m;
              }

              VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  - non-sliced cost matrix (" << i << ") size : " << (m->size())););

              cost<int> c = m->Pr_one_step_from_one(*((*vm)[i+gv_subalignment_length]),false);
              results[i] = (c > cost<int>(gv_lossless_cost_threshold)) ? 1 : 0;
              if (m_old) delete m_old;
              delete (*vm)[i];
            }
            for (int i = gv_alignment_length - gv_subalignment_length + 1; i<(int)vm->size(); i++)
              delete (*vm)[i];
            delete vm;
            // <<
#endif
            /// @todo{FIXME: complete this with probability computation}
            lossless = (int)(gv_subalignment_functions_int[gv_subalignment_function_index](results));
            sens     =       gv_subalignment_functions_int[gv_subalignment_function_index](results);

          } else {
            // (3.a.2) lossless without subalignment with matrix
            matrix<cost<int> > * m_cost_sens = a_spr_mx_h_res->matrix_cost_product(a_cost, PRODUCT_UNION_NO_FINAL_LOOP, gv_alignment_length);

            VERB_FILTER(VERBOSITY_ANNOYING, INFO__("= lossless cost matrix product size : " << (m_cost_sens->size())););

            // @note{NOTE : check if "non final states" (i.e. the ones that are "not rejected" or "non final") have all more than "k" costs}
            lossless = m_cost_sens->Pr(gv_alignment_length, false) > cost<int>(gv_lossless_cost_threshold);
            delete   m_cost_sens;
            // probability is needed when hill climbing is activated
#ifndef LOSSLESS_PROB
            if (gv_hillclimbing_flag) {
#endif
              automaton<void> * a_spr_mx_h_res_loss = a_spr_mx_h_res->product(a_lossless, PRODUCT_UNION_NO_FINAL_LOOP, gv_alignment_length);
              // @note{NOTE : both "reject" and "final" states are final so "final" should not be use to mesure probs}
              if (gv_minimize_flag) {
                automaton<void> * na = a_spr_mx_h_res_loss->Hopcroft();
                delete a_spr_mx_h_res_loss;
                a_spr_mx_h_res_loss = na;
              }
              matrix<double>     * m_pr_sens = a_spr_mx_h_res_loss->matrix_pr_product(a_sens, PRODUCT_UNION_NO_FINAL_LOOP, gv_alignment_length);

              VERB_FILTER(VERBOSITY_ANNOYING, INFO__("= lossless prob matrix product size : " << (m_pr_sens->size())););

              sens = (lossless_set_sens - (m_pr_sens->Pr(gv_alignment_length, false/* measure on non final (see upper) */)))/lossless_set_sens;
              delete               m_pr_sens;
              delete       a_spr_mx_h_res_loss;
#ifndef LOSSLESS_PROB
            }
#endif
          }
        } else {

          // (3.a.3) lossy with subalignment matrix
          if (gv_subalignment_flag) {

            std::vector< matrix<double> * > * vm = a_spr_mx_h_res->matrices_step_pr_product(a_sens, PRODUCT_UNION_NO_FINAL_LOOP, gv_alignment_length);
            std::vector<double> results(gv_alignment_length-gv_subalignment_length+1);
#ifndef NOSLICER
            // Sliced version
            // >>
            matrices_slicer<double> vm_slicer(vm);

            // set the sub-alignment length
            for (int i = 1; i < gv_subalignment_length; i++) vm_slicer.add_right();

            // move along the sub-alignment windows and compute Pr
            for (int i = gv_subalignment_length; i <= gv_alignment_length; i++) {
              results[i - gv_subalignment_length] = vm_slicer.current_Pr(true);
              vm_slicer.add_right();
              vm_slicer.del_left();
            }
            // <<
#else
            // Non sliced version
            // >>
            for (int i = 0; i <= gv_alignment_length - gv_subalignment_length; i++) {
              matrix<double> * m = (*vm)[i];
              matrix<double> * m_old = NULL;
              m->resize(2);
              for (int k = 1; k < gv_subalignment_length; k++) {
                m = m->Compose(*((*vm)[i+k]));
                if (m_old) delete m_old;
                m_old = m;
              }

              VERB_FILTER(VERBOSITY_ANNOYING, INFO__("  - non-sliced prob matrix (" << i << ") size : " << (m->size())););

              results[i] = m->Pr_one_step_from_one(*((*vm)[i+gv_subalignment_length]),true);
              if (m_old) delete m_old;
              delete (*vm)[i];
            }
            for (int i = gv_alignment_length - gv_subalignment_length + 1; i<(int)vm->size(); i++)
              delete (*vm)[i];
            delete vm;
            // <<
#endif
            sens = gv_subalignment_functions_double[gv_subalignment_function_index](results);

          } else {
            // (3.a.4) lossy without subalignment matrix
            matrix<double> * m_pr_sens = a_spr_mx_h_res->matrix_pr_product(a_sens, PRODUCT_UNION_NO_FINAL_LOOP, gv_alignment_length);
            VERB_FILTER(VERBOSITY_ANNOYING, INFO__("= prob matrix product size : " << (m_pr_sens->size())););

            // Multinomial evaluation can be enabled in that case
            if (gv_multipoly_file_flag) {
              /*
              // test 1 on polynomials
              automaton<polynomial<BIGINT > > * pr = a_spr_mx_h_res->product(*gv_multipoly_bsens_automaton, PRODUCT_UNION_NO_FINAL_LOOP, PRODUCT_OTHER_IS_PROBABILIST, gv_alignment_length);
              polynomial<BIGINT > pol1  = pr->Pr(gv_alignment_length,true);
              cout << endl << "(a) [" << pol1 << "]" << endl;
              polynomial<BIGINT > inv_pol1  = pr->Pr(gv_alignment_length,false);
              cout << endl << "(a) {" << inv_pol1 << "}" << endl;
              cout << endl << "(a) <" << (pol1 + inv_pol1) << ">" << endl;
              delete pr;

              // test 2 on matrices
              matrix<polynomial<BIGINT > > * m_pr = a_spr_mx_h_res->matrix_product(*gv_multipoly_bsens_automaton, PRODUCT_UNION_NO_FINAL_LOOP, gv_alignment_length);
              polynomial<BIGINT > pol2 = m_pr->Pr(gv_alignment_length,true);
              cout << endl << "(b) [" << pol2 << "]" << endl;
              polynomial<BIGINT > inv_pol2  = m_pr->Pr(gv_alignment_length,false);
              cout << endl << "(b) {" << inv_pol2 << "}" << endl;
              cout << endl << "(b) <" << (pol2 + inv_pol2) << ">" << endl;
              delete m_pr;
              */
              // implemented on matrices
              matrix<polynomial<BIGINT > > * m_pr = a_spr_mx_h_res->matrix_product(*gv_multipoly_bsens_automaton, PRODUCT_UNION_NO_FINAL_LOOP, gv_alignment_length);
              polynomial<BIGINT > mpol = m_pr->Pr(gv_alignment_length,true);
              multipoly = new polynomial<BIGINT >(mpol);
              delete m_pr;
            }

            sens                       = m_pr_sens->Pr(gv_alignment_length, true);
            delete m_pr_sens;
          }
        } // "else" part of "if (gv_subalignment_flag)"
      } // "else" part of "if (gv_correlation_flag)"


      //FIXMECOV>>
    gv_covariance_flag_1:
      //<<

      if (lossless) {
        sens = 1.0;
      }

      if (
#ifndef KEEP_PRODUCT_MF
          gv_seeds.size() > 1 ||
#endif
          gv_global_coverage_flag) {
        delete a_spr;
      }

      if (gv_multihit_flag || gv_global_coverage_flag) {
        delete a_spr_mhits_or_gcov_res;
      }
      if (gv_homogeneous_flag) {
        delete a_spr_h_res;
      }
      if (gv_xseeds.size()) {
        delete a_spr_mx_h_res;
      }
      sens /= pr_div_homogeneous - pr_div_excluded;

      //FIXMECOV>>
      if (gv_covariance_flag) {
        double covariance = 0.0;
        for (unsigned i = 0; i < gv_seeds.size(); i++) {
          int span_i = gv_seeds[i]->span();
          automaton<double> * a_x_v_i = a_s[i]->product(a_sens, PRODUCT_UNION_NO_FINAL_LOOP, span_i);
          double x_v_i        = a_x_v_i->Pr(span_i);
          delete a_x_v_i;
          for (unsigned j = i; j < gv_seeds.size(); j++) {
            int span_j = gv_seeds[j]->span();
            automaton<double> * a_x_v_j = a_s[j]->product(a_sens, PRODUCT_UNION_NO_FINAL_LOOP, span_j);
            double x_v_j        = a_x_v_j->Pr(span_j);
            delete a_x_v_j;
            for (int shift = -span_j+1; shift <= span_i-1; shift++) {
              int len = MAX(span_i - MIN(shift,0) , span_j + MAX(shift,0));
              automaton<void> * p = a_s[i]->product(*(a_s[j]),PRODUCT_INTERSECTION_FINAL_LOOP, len, NULL, shift);
              if (gv_minimize_flag) {
                automaton<void> * na = p->Hopcroft();
                delete p;
                p = na;
              }
              automaton<double> * a_x_cov = p->product(a_sens, PRODUCT_UNION_NO_FINAL_LOOP, len);
              double x_cov        = a_x_cov->Pr(len);
              covariance += x_cov - x_v_i*x_v_j;
              //cout << "i:" << i << ",j:"<< j << ",shift:" << shift << ":" << (x_cov) << endl;
              delete a_x_cov;
              delete p;
            }
          }
        }
        //cout << "covariance : " << covariance << endl;
        sens = 1.0 / (1.0 + covariance);
      }
      //FIXMECOV<<
      //cout << "sens : " << sens << endl;

      // compute selectivity
      selp = 1.0 - (double)MIN(1.0, selp);

#ifdef SAMPLING_STRATEGY_MF
      VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__("" << (sens)););
#endif

      // (4) build the memorizing structure
      double distance = dist(selp, sens);
      ostringstream outs;
      for (unsigned i = 0; i < gv_seeds.size(); i++) {
        if (i > 0)
          outs << ",";
        outs << (gv_seeds[i])->str();
      }
      seedproperties e = seedproperties(selp, sens, distance, outs.str(), lossless, polynom, multipoly);

      if (gv_correlation_flag || gv_polynomial_dominant_selection_flag) {
        polynom->clear();
        delete polynom;
      }

      if (gv_multipoly_file_flag) {
        delete multipoly;
      }

      //
      // (5) insertion inside pareto set
      //
      //cerr << e.polynom << endl;
      sensitivity_threshold = insertPareto(l, e);

#ifdef SAMPLING_STRATEGY_MF
      max_sens = MAX(max_sens,sens);
#endif
      if (gv_motif_flag && gv_nbruns == 0 && !(gv_signature_flag || gv_signature_shuffle_from_m_pattern_flag)) {
        goto end_loop;
      }

      //
      // (6)  Sleep if requested
      //
      // not needed anymore ...

    } //<<(1-6)

#ifdef SAMPLING_STRATEGY_MF
  next_progress_seed:
#endif

    { //(7)>>

      //
      // (7) Display progress
      //
      if (gv_pareto_select_runs && !(nbruns % gv_pareto_select_runs)) {
        double area    = selectPareto(l);
        time_t    t    = time(NULL);
        struct tm * tm = localtime(&t);
        ostringstream outs;
        for (unsigned i = 0; i < gv_seeds.size(); i++) {
          if (i > 0)
            outs << ",";
          outs << (gv_seeds[i])->str();
        }
        VERB_FILTER(VERBOSITY_LIGHT, MESSAGE__(
                                               "#runs:" << nbruns << "," << "\t" << "set size:" << l.size() << "," << "\t" << "area:" << area
                                               << "\t"  << (tm->tm_year + 1900) << "-" << setw(2) << setfill('0') << (tm->tm_mon+1) << "-"  << setw(2) << setfill('0') << (tm->tm_mday) << ","
                                               << "\t"  << setw(2) << setfill('0') << (tm->tm_hour)  << ":"  << setw(2) << setfill('0') << (tm->tm_min) << ":" <<  setw(2) << setfill('0') << (tm->tm_sec)
                                               << "\t"  << "last motif checked:" << outs.str()
                                               << "\t"  << (gv_correlation_flag?"last corr:":"last sens:") << sens;
                                               if (gv_hillclimbing_flag) {
                                                 cerr << "," << "\t" << "hillclimbing_threshold:" << hillclimbing_threshold;
                                               }
                                               ););
      }
    } //<<(7)

    { //(8)>>

      nbruns++;
      if (gv_nbruns && nbruns >= gv_nbruns)
        goto end_loop;

      // (8) hill climbing method
      if (gv_hillclimbing_flag) {

        // (8.1) set hill_climbing_flag
        if (!hillclimbing_flag) {
          double threshold = sensitivity_threshold + hillclimbing_threshold;
          if (threshold  > 0) {
            VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__(
                                                      "\t [+] : " << sens << " (" << threshold << " = " <<  sensitivity_threshold << " + " << hillclimbing_threshold << ")" << endl;
                                                      ););
            hillclimbing_flag = true;
            hillclimbing_threshold = MAX(1e-6, hillclimbing_threshold - threshold); // more stringent
          } else {
            VERB_FILTER(VERBOSITY_ANNOYING, MESSAGE__(
                                                      "\t [-] : " << sens << " (" << threshold << " = " <<  sensitivity_threshold << " + " << hillclimbing_threshold << ")" << endl;
                                                      ););
            hillclimbing_threshold = MIN(0.999999, hillclimbing_threshold - (threshold*power(gv_hillclimbing_alpha,log(log(gv_seeds.size()+1)+1)))); // little less stringent
          }
        }

        // (8.2) launch hill_climbing procedure
        if (hillclimbing_flag) {

          // hillclimb phenomenom detected
          if (sens > sens_hillclimbing) {

            // (8.0) printing global improvements if any
            if (sensitivity_threshold > 0) {
              VERB_FILTER(VERBOSITY_MODERATE, MESSAGE__(
                                                        "- global hill-climbing improvement : ";
                                                        for (unsigned i = 0; i < gv_seeds.size(); i++) {
                                                          cerr << *(gv_seeds[i]) << ",";
                                                        }
                                                        cerr << "\tsel:" << selp << "\t" << (gv_correlation_flag?"corr:":"sens:") << sens << "\t [" << sensitivity_threshold << "]" << endl;
                                                        ););
            } else {
              VERB_FILTER(VERBOSITY_MODERATE, MESSAGE__(
                                                        "- local hill-climbing improvement : ";
                                                        for (unsigned i = 0; i < gv_seeds.size(); i++) {
                                                          cerr << *(gv_seeds[i]) << ",";
                                                        }
                                                        cerr << "\tsel:" << selp << "\t" << (gv_correlation_flag?"corr:":"sens:") << sens << "\t [" << sensitivity_threshold << "]" << endl;
                                                        ););
              sens_hillclimbing = sens;
              for (unsigned i = 0; i < gv_seeds.size(); i++) {
                (gv_seeds[i])->set_hmove();
              }
              gv_hmove_choice = rand();

              last_seed_to_hillclimbing = seed_to_hillclimbing;
              seed_to_hillclimbing ++;
              seed_to_hillclimbing %= gv_seeds.size();
            }
          }

          {
          next_hmove:

            // next seed in hillclimbing order
            while ((gv_seeds[seed_to_hillclimbing])->next_hmove() == 0 ) {

              // delete automaton associated with the modified seed
              if (a_s[seed_to_hillclimbing]) {
                delete a_s[seed_to_hillclimbing];
                a_s[seed_to_hillclimbing] = NULL;
              }
              if (seed_to_hillclimbing == last_seed_to_hillclimbing) {
                goto eof_hmove;
              }
              seed_to_hillclimbing ++;
              seed_to_hillclimbing %= gv_seeds.size();
            }


            // delete automaton associated with the modified seed
            if (a_s[seed_to_hillclimbing]) {
              delete a_s[seed_to_hillclimbing];
              a_s[seed_to_hillclimbing] = NULL;
            }

            // (8.3) : discard twice seeds in the same set
            for (unsigned i = 0; i < gv_seeds.size(); i++)
              if (!gv_seeds[i]->cycled())
                for (unsigned j = 0; j < i; j++)
                  if (!gv_seeds[j]->cycled())
                    if (gv_seeds[i]->equal(gv_seeds[j]))
                      goto next_hmove;
            continue;

          eof_hmove:
            // end of hill climbing
            hillclimbing_flag = false;
            sens = sens_hillclimbing;
            sens_hillclimbing = 0.0;

            for (unsigned i = 0; i < gv_seeds.size(); i++) {
              (gv_seeds[seed_to_hillclimbing])->reset_hmove();
            }
            gv_hmove_choice = rand();

            last_seed_to_hillclimbing = seed_to_hillclimbing;
            seed_to_hillclimbing ++;
            seed_to_hillclimbing %= gv_seeds.size();

            VERB_FILTER(VERBOSITY_MODERATE, MESSAGE__(
                                                      "- final local hill-climbing search for : ";
                                                      for (unsigned i = 0; i < gv_seeds.size(); i++) {
                                                        cerr << *(gv_seeds[i]) << ",";
                                                      }
                                                      cerr << "\tsel:" << selp << "\t" <<  (gv_correlation_flag?"corr:":"sens:") << sens << "\t [" << sensitivity_threshold << "]" << endl;
                                                      ););
          }
        } // if (hill_climbing_flag)
      } // if (gv_hill_climbing_flag)
    } //<<(8)


    { //>>(9)

      // (9) : enumerate the next seed(s)
    new_seeds:
      if (gv_nbruns > 0) {
        double first_seed_weight = 0;

        // (9.1) : random enumeration
        for (unsigned i = 0; i < gv_seeds.size(); i++ ) {

        next_rand:
          (gv_seeds[i])->random();

          // delete automaton associated with the modified seed
          if (a_s[i]) {
            delete a_s[i];
            a_s[i] = NULL;
          }


          // jive selection activated of not  ?
          if (gv_jive > 0){
            if (i > 0) {
              double weight = (gv_seeds[i])->weight();
              if ( ( weight / first_seed_weight < gv_jive ) || ( first_seed_weight / weight < gv_jive ) )
                goto next_rand;
            } else {
              first_seed_weight = (gv_seeds[i])->weight();
              if (first_seed_weight <= 0)
                goto next_rand;
            }
          }
        }

      } else { // gv_nbruns == 0

        // (9.2) : complete enumeration
        unsigned j = 0;

        while ( (gv_seeds[j])->next() == 0 ) {

          // delete automaton associated with the modified seed
          if (a_s[j]) {
            delete a_s[j];
            a_s[j] = NULL;
          }

          j++;

          if (j == gv_seeds.size())
            goto end_loop;
        } // while (seeds->next())


        // delete automaton associated with the modified seed (a_s[j]->next == 1)
        if (j < gv_seeds.size() && a_s[j]) {
          delete a_s[j];
          a_s[j] = NULL;
        }

        for (unsigned i = 0; i < j; i++) { // reset the span of previous j-th seeds to min
          if (gv_signature_shuffle_from_m_pattern_flag || gv_signature_flag) {
            gv_seeds[i]->reorder(0,/*b*/ gv_seed_alphabet_size-1);
            while(!gv_seeds[i]->acceptable())
              if (!gv_seeds[i]->next())
                goto end_loop;
          } else {
            delete gv_seeds[i];
            gv_seeds[i] = new seed();
            if (gv_cycles_flag)
              gv_seeds[i]->setCyclePos(gv_cycles_pos_list[i], gv_cycles[i]);
          }
        }
      } // if (gv_nbruns > 0)

      // (9.3) : discard twice same seed in the same set
      for (unsigned i = 0; i < gv_seeds.size(); i++)
        if (!gv_seeds[i]->cycled())
          for (unsigned j = 0; j < i; j++)
            if (!gv_seeds[j]->cycled())
              if (gv_seeds[i]->equal(gv_seeds[j]))
                goto new_seeds;

    } //<<(9)

  }// while (1)
 end_loop:;

  // deleting several global allocated structures
  if (gv_homogeneous_flag)
    delete a_homogeneous;

  if (gv_xseeds.size()) {
    delete a_excluded;
    for (unsigned i = 0; i < gv_xseeds.size(); i++)
      delete gv_xseeds[i];
  }

  for (unsigned i = 0; i < gv_seeds.size(); i++) {
    delete gv_seeds[i];
    if (a_s[i])
      delete a_s[i];
  }

  if (gv_bsens_automaton)
    delete gv_bsens_automaton;

  if (gv_multipoly_file_flag)
    delete gv_multipoly_bsens_automaton;

  if (gv_bsymbols_array)
    free(gv_bsymbols_array);

#ifdef KEEP_PRODUCT_MF
  for (unsigned v = 0; v < gv_seeds.size(); v++)
    if (v > 0 && a_s_product[v])
      delete a_s_product[v];
#endif


  //
  // (10) Display pareto set and area
  //

  // update pareto set according to other processes (helpful on parallel computation)
  if (gv_nb_input_filenames > 0) {
    for (int i = 0; i < gv_nb_input_filenames; i++)
      inputPareto(l, gv_input_filenames[i]);
  }

  // final output
  VERB_FILTER(VERBOSITY_LIGHT, MESSAGE__("# seeds\t sel\t " << (gv_correlation_flag?"corr":"sens") << "\t dist" << (gv_polynomial_dominant_selection_flag?"\tpol":"")););
  if (!gv_output_filename)
    list_and_areaPareto(l);
  else
    outputPareto(l, gv_output_filename);

  // cleaning filename memory
  for (int i = 0; i < gv_nb_input_filenames; i++)
    free(gv_input_filenames[i]);
  if (gv_output_filename)
    free(gv_output_filename);
  return 0;
}

// @}
