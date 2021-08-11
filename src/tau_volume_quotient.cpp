/*!
 * \file    tau_volume_quotient.cpp
 * \ingroup estimating_volume_quotients
 *
 * \brief   The definition of functions for converting tau estimates into volume
 *          quotient estimates.
 */

#include "tau_volume_quotient.h"

#include "common.h"

#include <gmp.h>
#include <mpfr.h>

#include <stddef.h>
#include <stdint.h>

bool tau_volume_quotient(
  const uint32_t m,
  const uint32_t l,
  const uint32_t n,
  const mpz_t d_or_r,
  const long double tau,
  mpfr_t v)
{
  /* Temporary variables. */
  mpfr_t tmp;
  mpfr_init2(tmp, PRECISION);

  mpfr_t radius;
  mpfr_init2(radius, PRECISION);

  mpfr_t gamma;
  mpfr_init2(gamma, PRECISION);

  /*
   * Given
   * 
   *    tau = log_( (1 / n) \sum_{i = 1}^{n} alpha^2_i ) / 2 - m 
   * 
   * we first compute the radius
   * 
   *    R = sqrt( \sum_{i = 1}^{n} alpha^2_i + ( d^2 or r^2 ) ),
   * 
   * and then the volume quotient
   * 
   *    v = pi^(D / 2) R^D / Gamma(D / 2 + 1) / 2^((l + m) n),
   * 
   * where D = n + 1 is the dimension of the lattice.
   * 
   * Below, R is computed from tau in in the following steps:
   * 
   * 1. 2^(2 (tau + m)) = (1 / n) \sum_{i = 1}^{n} alpha^2_i
   * 
   * 2. n 2^(2 (tau + m)) = \sum_{i = 1}^{n} alpha^2_i
   * 
   * 3. n 2^(2 (tau + m)) + d^2 = \sum_{i = 1}^{n} alpha^2_i + d^2
   * 
   * 4. sqrt(n 2^(2 (tau + m))) + d^2 = sqrt( \sum_{i = 1}^{n} alpha^2_i + d^2 )
   * 
   * The steps required to compute v given R are rather self-explanatory.
   */
  
  mpfr_set_ld(radius, 2.0f * (tau + m), MPFR_RNDN); /* radius = 2 (tau + m) */
  mpfr_exp2(radius, radius, MPFR_RNDN);
    /* radius = 2^(2 (tau + m)) = (1 / n) \sum_{i = 1}^{n} alpha^2_i */
  mpfr_mul_ui(radius, radius, n, MPFR_RNDN); 
    /* radius = \sum_{i = 1}^{n} alpha^2_i */
  
  if (NULL == d_or_r) {
    mpfr_set_ui_2exp(tmp, 1, (mpfr_exp_t)m, MPFR_RNDN); /* tmp = 2^m */
  } else {
    mpfr_set_z(tmp, d_or_r, MPFR_RNDN); /* tmp = d or r */
  }

  mpfr_sqr(tmp, tmp, MPFR_RNDN); /* tmp = d^2 or r^2 */
  mpfr_add(radius, radius, tmp, MPFR_RNDN); 
    /* radius = \sum_{i = 1}^{n} alpha^2_i + (d^2 or r^2) */
  mpfr_sqrt(radius, radius, MPFR_RNDN);
    /* radius = sqrt(\sum_{i = 1}^{n} alpha^2_i + (d^2 or r^2)) */

  mpfr_set_ld(gamma, (long double)(n + 1) / 2 + 1, MPFR_RNDN);
    /* gamma = (n + 1) / 2 + 1 = D / 2 + 1 */
  mpfr_gamma(gamma, gamma, MPFR_RNDN);
    /* gamma = Gamma(D / 2 + 1) */

  /* Compute the volume of a D-dimensional ball of radius R. */
  mpfr_const_pi(v, MPFR_RNDN); /* v = pi */
  mpfr_set_ld(tmp, (long double)(n + 1) / (long double)2, MPFR_RNDN);
    /* tmp = (n + 1) / 2 = D / 2 */
  mpfr_pow(v, v, tmp, MPFR_RNDN); /* v = pi^(D / 2) */

  mpfr_set_ui(tmp, n + 1, MPFR_RNDN); /* tmp = n + 1 = D */
  mpfr_pow(tmp, radius, tmp, MPFR_RNDN); /* tmp = radius^D */
  mpfr_mul(v, v, tmp, MPFR_RNDN); /* v = pi^(D / 2) radius^D */
  mpfr_div(v, v, gamma, MPFR_RNDN);
    /* v = pi^(D / 2) radius^D / Gamma(D / 2 + 1) */

  /* Divide by the determinant. */
  mpfr_set_ui_2exp(tmp, 1, (mpfr_exp_t)((l + m) * n), MPFR_RNDN);
    /* tmp = 2^((l + m) n) */
  mpfr_div(v, v, tmp, MPFR_RNDN);
    /* v = pi^(D / 2) radius^D / Gamma(D / 2 + 1) / 2^((l + m) n) */

  /* Clear memory. */
  mpfr_clear(tmp);
  mpfr_clear(radius);
  mpfr_clear(gamma);

  /* Return. */
  return mpfr_cmp_ui(v, 2) < 0;
}
