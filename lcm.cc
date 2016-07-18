#include <iostream>
#include <iomanip>
using namespace std;

#define USEINFINT

#ifdef USEINFINT
#include "infint.hh"
#define BIGINT InfInt<long long int>
#else
#define BIGINT __uint128_t
#endif


/*
 * Please change this value if you need to increase the size of the table (this is quadratic in size according to N ...)
 */

#define N 64
BIGINT BINOMIAL[N+1][N+1];


void generate_binomials() {
  int k, n;
  for (k = 1; k <= N; k++)
    BINOMIAL[0][k] = 0;
  for (n = 0; n <= N; n++)
    BINOMIAL[n][0] = 1;
  for (n = 1; n <= N; n++)
    for (k = 1; k <= N; k++)
      BINOMIAL[n][k] = BINOMIAL[n-1][k-1] + BINOMIAL[n-1][k];
}

BIGINT GCD(BIGINT a, BIGINT b) {
  if (a < b) {
    BIGINT t = b;
    b = a;
    a = t;
  }
  while (b > 0) {
    BIGINT t = b;
    b = a % b;
    a = t;
  }
  return a;
}

BIGINT LCM(BIGINT a, BIGINT b) {
  return a/GCD(a,b) * b;
}

ostream& DOUBLE(ostream & os, const BIGINT&b) {
  BIGINT bcopy = BIGINT(b);
  int    exp   = 0;
  while (bcopy >= 100000000) {
    bcopy = bcopy / 10;
    exp++;
  }
#ifdef USEINFINT
  os << (bcopy) << "e" << exp;
#else
  os << ((long long int)(bcopy)) << "e" << exp;
#endif
  return os;
}




#ifdef USEINFINT
#define INT128_LONG_LONG_CODE(a) "(((__uint128_t)" << (((BIGINT(a)/BIGINT("18446744073709551616"))%BIGINT("18446744073709551616")).toUnsignedLongLong()) <<  "ULL)<<64) | (__uint128_t)" << ((BIGINT(a)%BIGINT("18446744073709551616")).toUnsignedLongLong()) << "ULL"
#else
#define INT128_LONG_LONG_CODE(a) "(((__uint128_t)" << (long long unsigned)((a)>>64) << "ULL)<<64) | (__uint128_t)" << (long long unsigned)((a)&((((__uint128_t)1)<<64)-1)) << "ULL"
#endif


int main() {

  generate_binomials();
  BIGINT lcm = 1;
  long long int k_plus_1;
  cout << "#ifndef __BINOMIAL_WEIGHT__" << endl;
  cout << "#define __BINOMIAL_WEIGHT__" << endl;
  cout << "#define N_binomial_weight " << N << endl;

  /*
   * Generate this code if USEITINF activated >>
   */

#ifdef USEINFINT
  cout << "#ifdef USEINFINT" << endl;
  cout << "InfInt<long long int> binomial_weight[N_binomial_weight+1][N_binomial_weight+1] = {" << endl;
  lcm = 1;
  for (k_plus_1 = 1; k_plus_1 <= N+1 ; k_plus_1++) {
    lcm = LCM(lcm,k_plus_1);

    /*
     * lcm { binomial(k,0) ; binomial(k,1) ; ... ; binomial(k,k) } = lcm {1;2;...;k;k+1} / (k+1)
     * "An identity involving the least common multiple of binomial coefficients and its application" Bakir FARHI
     */

    BIGINT lcm_div_k_plus_1 = lcm / (k_plus_1);
    cout << "\t{" << endl;
    long long int i;
    long long int k = k_plus_1 - 1;
    for (i = 0; i <= k ; i++) {
      BIGINT lcm_div_k_plus_1_div_binomial_k_i = lcm_div_k_plus_1/(BINOMIAL[k][i]);
      cout << "\t\t /* llcm / Binomial{(k=" << k << "),(i=" << i << ")} = */ InfInt<long long int>(\"" << ((lcm_div_k_plus_1_div_binomial_k_i)) <<"\") ," << endl;
    }
    cout << "\t}," << endl;
  }
  cout << "};" << endl;
  cout << "#else" << endl;
#endif

  /*
   * << Generate this code if USEITINF activated
   */
  cout << "#ifdef __SIZEOF_INT128__" << endl;
  cout << "__uint128_t binomial_weight[N_binomial_weight+1][N_binomial_weight+1] = {" << endl;
  lcm = 1;
  for (k_plus_1 = 1; k_plus_1 <= N+1 ; k_plus_1++) {
    lcm = LCM(lcm,k_plus_1);

    /*
     * lcm { binomial(k,0) ; binomial(k,1) ; ... ; binomial(k,k) } = lcm {1;2;...;k;k+1} / (k+1)
     * "An identity involving the least common multiple of binomial coefficients and its application" Bakir FARHI
     */

    BIGINT lcm_div_k_plus_1 = lcm / (k_plus_1);
    cout << "\t{" << endl;
    long long int i;
    long long int k = k_plus_1 - 1;
    for (i = 0; i <= k ; i++) {
      BIGINT lcm_div_k_plus_1_div_binomial_k_i = lcm_div_k_plus_1/(BINOMIAL[k][i]);
      cout << "\t\t /* llcm / Binomial{(k=" << k << "),(i=" << i << ")} = */" <<  INT128_LONG_LONG_CODE(lcm_div_k_plus_1_div_binomial_k_i)<< "," << endl;
    }
    cout << "\t}," << endl;
  }
  cout << "};" << endl;

  cout << "#else" << endl;
  cout << "double binomial_weight[N_binomial_weight+1][N_binomial_weight+1] = {" << endl;
  lcm = 1;
  for (k_plus_1 = 1; k_plus_1 <= N+1 ; k_plus_1++) {
    lcm = LCM(lcm,k_plus_1);

    /*
     * lcm { binomial(k,0) ; binomial(k,1) ; ... ; binomial(k,k) } = lcm {1;2;...;k;k+1} / (k+1)
     * "An identity involving the least common multiple of binomial coefficients and its application" Bakir FARHI
     */

    BIGINT lcm_div_k_plus_1 = lcm / (k_plus_1);
    cout << "\t{" << endl;
    long long int i;
    long long int k = k_plus_1 - 1;
    for (i = 0; i <= k ; i++) {
      BIGINT lcm_div_k_plus_1_div_binomial_k_i = lcm_div_k_plus_1/(BINOMIAL[k][i]);
      cout.precision(15);
      cout << "\t\t /* llcm / Binomial{(k=" << k << "),(i=" << i << ")} = */ (double)" << DOUBLE(cout,lcm_div_k_plus_1_div_binomial_k_i) << "," << endl;
    }
    cout << "\t}," << endl;
  }
  cout << "};" << endl;
  cout << "#endif" << endl;
#ifdef USEINFINT
  cout << "#endif" << endl;
#endif
  cout << "#endif" << endl;
  return 0;
}


