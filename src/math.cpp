/*!
 * \file    math.cpp
 * \ingroup math
 *
 * \brief   The definition of basic mathematical functions.
 */

#include "math.h"

#include "common.h"

#include <gmp.h>

#include <stdint.h>

uint32_t abs_i(const int32_t x) {
  return (uint32_t)((x < 0) ? -x : x);
}

double abs_d(const double x) {
  return (x < 0) ? -x : x;
}

long double abs_ld(const long double x) {
  return (x < 0) ? -x : x;
}

int32_t sgn_i(const int32_t x) {
  return (x < 0) ? -1 : 1;
}

int32_t sgn_d(const double x) {
  return (x < 0) ? -1 : 1;
}

int32_t sgn_ld(const long double x) {
  return (x < 0) ? -1 : 1;
}

int32_t max_i(const int32_t a, const int32_t b) {
  return (a > b) ? a : b;
}

uint32_t max_ui(const uint32_t a, const uint32_t b) {
  return (a > b) ? a : b;
}

void mod_reduce(mpz_t x, const mpz_t n) {
  mpz_mod(x, x, n); /* x = x mod n */

  mpz_t tmp;
  mpz_init(tmp);
  mpz_mul_ui(tmp, x, 2); /* tmp = 2 * x */

  if (mpz_cmp(tmp, n) >= 0) { /* if 2 * x >= n  <=> x >= n/2 */
    mpz_sub(x, x, n); /* subtract n from x */
  }

  /* Clear memory. */
  mpz_clear(tmp);
}

uint32_t kappa(const mpz_t x)
{
  mpz_t acc;
  mpz_init(acc);
  mpz_abs(acc, x);

  mpz_t result;
  mpz_init(result);

  uint32_t kappa_value = 0;

  while (mpz_cmp_ui(acc, 0) > 0) {
    mpz_mod_ui(result, acc, 2);
    if (mpz_cmp_ui(result, 1) == 0) {
      break;
    }

    mpz_div_ui(acc, acc, 2);

    kappa_value++;
  }

  /* Clear memory. */
  mpz_clear(acc);
  mpz_clear(result);

  /* Return kappa. */
  return kappa_value;
}

bool is_pow2(const uint32_t x) {
  for (uint32_t pow = 1; 0 != pow; pow <<= 1) {
    if (x == pow) {
      return TRUE;
    }
  }

  return FALSE;
}
