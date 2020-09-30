#include "macro.h"

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>




int                    gv_hmove_choice = 0;
bool                   gv_symetric = false;
int                    gv_nbruns = 0;
bool                   gv_matching_symbol_flag = true;
int                    gv_alignment_length = 64;
int                    gv_verbose = 1;
bool                   gv_subalignment_flag = false;
int                    gv_subalignment_length = 0;
int                    gv_subalignment_function_index = 0;
char *                 gv_subalignment_functions_names[SUBALIGNMENT_FUNCTIONS_NUMBER] = {(char *)"min",(char *)"max",(char *)"avg", (char*)"med"};

template<typename T> double min_vector(vector<T> & v) {return (double)(*(std::min_element(v.begin(),v.end())));}
template<typename T> double max_vector(vector<T> & v) {return (double)(*(std::max_element(v.begin(),v.end())));}
template<typename T> double avg_vector(vector<T> & v) {return (double)(std::accumulate(v.begin(), v.end(), 0.0)) / v.size();}
template<typename T> double med_vector(vector<T> & v) {size_t n = v.size() / 2; nth_element(v.begin(), v.begin()+n, v.end()); return (double)v[n];}

subalignment_function_int gv_subalignment_functions_int[SUBALIGNMENT_FUNCTIONS_NUMBER] = {
  (subalignment_function_int)  &min_vector<int>,
  (subalignment_function_int)  &max_vector<int>,
  (subalignment_function_int)  &avg_vector<int>,
  (subalignment_function_int)  &med_vector<int>,
};

subalignment_function_double gv_subalignment_functions_double[SUBALIGNMENT_FUNCTIONS_NUMBER] = {
  (subalignment_function_double)  &min_vector<double>,
  (subalignment_function_double)  &max_vector<double>,
  (subalignment_function_double)  &avg_vector<double>,
  (subalignment_function_double)  &med_vector<double>
};

int                    gv_align_alphabet_size = 2;
int                    gv_seed_alphabet_size = 2;
vector< double >       gv_bsel_weight;
vector< int >          gv_signature;
bool                   gv_signature_flag = false;
bool                   gv_signature_shuffle_from_m_pattern_flag = false;
double                 gv_bsel_minprob,gv_bsel_maxprob;
int                    gv_minspan = 1, gv_maxspan = 16;
double                 gv_minweight = -1e32, gv_maxweight = 1e32;
bool                   gv_weight_interval_flag = false;
bool                   gv_lossless_flag = false;
vector<int>            gv_lossless_costs_vector;
int                    gv_lossless_cost_threshold = 0;
char *                 gv_bsymbols_array = NULL;
bool                   gv_bsymbols_flag  = false;


int power(int x, int n) {
  int res = 1;
  int u   = x;
  for(int i = 1; i <= n; i<<=1){
    if (i&n)
      res *= u;
    u *= u;
  }
  return res;
}


double power(double x, int n) {
  double res = 1.0;
  double u   = x;
  for(int i = 1; i <= n; i<<=1){
    if (i&n)
      res *= u;
    u *= u;
  }
  return res;
}


