/*!
 * \file    sample.cpp
 * \ingroup sample_distribution
 *
 * \brief   The definition of functions for sampling probability distributions.
 */

#include "sample.h"

#include "parameters.h"
#include "diagonal_parameters.h"
#include "random.h"
#include "errors.h"
#include "math.h"
#include "common.h"

#include <gmp.h>
#include <mpfr.h>

#include <math.h>

#include <stdint.h>

void sample_approximate_alpha_from_region(
  mpfr_t alpha,
  const double min_log_alpha,
  const double max_log_alpha,
  Random_State * const random_state)
{
  /* Sanity checks. */
  if (sgn_d(min_log_alpha) != sgn_d(max_log_alpha)) {
    critical("sample_approximate_alpha_from_region(): "
      "Incompatible signs for min_log_alpha and max_log_alpha.");
  }

  if (abs_d(min_log_alpha) >= abs_d(max_log_alpha)) {
    critical("sample_approximate_alpha_from_region(): "
      "Incompatible absolute values for min_log_alpha and max_log_alpha.");
  }

  /* Construct the minimum and maximum absolute alpha. */
  mpfr_t alpha_min;
  mpfr_init2(alpha_min, PRECISION);

  mpfr_t alpha_max;
  mpfr_init2(alpha_max, PRECISION);

  mpfr_t alpha_diff;
  mpfr_init2(alpha_diff, PRECISION);

  /* Sample alpha. */
  mpfr_set_d(alpha_min, abs_d(min_log_alpha), MPFR_RNDN);
  mpfr_exp2(alpha_min, alpha_min, MPFR_RNDN);
  mpfr_round(alpha_min, alpha_min);

  mpfr_set_d(alpha_max, abs_d(max_log_alpha), MPFR_RNDN);
  mpfr_exp2(alpha_max, alpha_max, MPFR_RNDN);
  mpfr_round(alpha_max, alpha_max);

  mpfr_sub(alpha_diff, alpha_max, alpha_min, MPFR_RNDN);

  const long double fraction = random_generate_pivot_exclusive(random_state);

  mpfr_mul_d(alpha, alpha_diff, fraction, MPFR_RNDN);
  mpfr_add(alpha, alpha, alpha_min, MPFR_RNDN);

  /* Select the sign of alpha. */
  if (sgn_d(min_log_alpha) == -1) {
    mpfr_neg(alpha, alpha, MPFR_RNDN);
  }

  /* Clear. */
  mpfr_clear(alpha_min);
  mpfr_clear(alpha_max);
  mpfr_clear(alpha_diff);
}

void sample_alpha_from_region(
  mpz_t alpha,
  const double min_log_alpha,
  const double max_log_alpha,
  const uint32_t kappa,
  Random_State * const random_state)
{
  /* Sanity checks. */
  if (sgn_d(min_log_alpha) != sgn_d(max_log_alpha)) {
    critical("sample_alpha_from_region(): "
      "Incompatible signs for min_log_alpha and max_log_alpha.");
  }

  if (abs_d(min_log_alpha) >= abs_d(max_log_alpha)) {
    critical("sample_alpha_from_region(): "
      "Incompatible absolute values for min_log_alpha and max_log_alpha.");
  }

  /* Setup precision. */
  const uint32_t m = ceil(abs_d(max_log_alpha));

  const uint32_t precision = 3 * m;

  /* Construct min_alpha. */
  mpfr_t min_alpha;
  mpfr_init2(min_alpha, precision);

  mpfr_set_d(min_alpha, abs_d(min_log_alpha), MPFR_RNDN);
  mpfr_exp2(min_alpha, min_alpha, MPFR_RNDN);
  mpfr_round(min_alpha, min_alpha);

  /* Construct max_alpha. */
  mpfr_t max_alpha;
  mpfr_init2(max_alpha, precision);

  mpfr_set_d(max_alpha, abs_d(max_log_alpha), MPFR_RNDN);
  mpfr_exp2(max_alpha, max_alpha, MPFR_RNDN);
  mpfr_round(max_alpha, max_alpha);

  /* Randomize alpha given min_alpha and max_alpha. */
  mpz_t min_alpha_z;
  mpz_init(min_alpha_z);
  mpfr_get_z(min_alpha_z, min_alpha, MPFR_RNDN);

  mpz_t max_alpha_z;
  mpz_init(max_alpha_z);
  mpfr_get_z(max_alpha_z, max_alpha, MPFR_RNDN);

  mpz_t alpha_modulus_z;
  mpz_init(alpha_modulus_z);
  mpz_sub(alpha_modulus_z, max_alpha_z, min_alpha_z);

  random_generate_mpz(alpha, alpha_modulus_z, random_state);
  mpz_add(alpha, min_alpha_z, alpha);

  if (kappa > 0) {
    /* Force the least significant 2^kappa bits to zero. */
    mpz_t x;
    mpz_init(x);
    mpz_set_ui(x, 0);

    mpz_setbit(x, kappa); /* x = 2^kappa */
    mpz_mod(x, alpha, x); /* x = alpha mod 2^kappa */
    mpz_sub(alpha, alpha, x); /* alpha = alpha - (alpha mod 2^kappa) */

    mpz_clear(x);
  }

  /* Adjust the sign. */
  if (sgn_d(min_log_alpha) == -1) {
    mpz_neg(alpha, alpha);
  }

  /* Clear memory. */
  mpfr_clear(min_alpha);
  mpfr_clear(max_alpha);

  mpz_clear(min_alpha_z);
  mpz_clear(max_alpha_z);
  mpz_clear(alpha_modulus_z);
}

void sample_j_from_alpha_r(
  mpz_t j,
  const mpz_t alpha_r,
  const Parameters * const parameters,
  Random_State * const random_state)
{
  /* Declare variables. */
  mpz_t t_r;
  mpz_init(t_r);

  mpz_t pow2;
  mpz_init(pow2);

  /* Compute kappa_r. */
  uint32_t kappa_r = kappa(parameters->r);

  if (kappa_r > 0) {
    mpz_set_ui(pow2, 0);
    mpz_setbit(pow2, kappa_r);
    random_generate_mpz(t_r, pow2, random_state);
  }

  /* Compute the integer j. */
  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, kappa_r);
  mpz_div(j, parameters->r, pow2);

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->l);
  mpz_invert(j, j, pow2);

  mpz_mul(j, j, alpha_r);

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, kappa_r);
  mpz_div(j, j, pow2);

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->l - kappa_r);
  mpz_mul(pow2, pow2, t_r);
  mpz_add(j, j, pow2);

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->l);
  mpz_mod(j, j, pow2);

  /* Clear memory. */
  mpz_clear(pow2);
}

void sample_j_k_from_alpha_d(
  mpz_t j,
  mpz_t k,
  const mpz_t alpha_d,
  const Parameters * const parameters,
  Random_State * const random_state)
{
  /* Declare variables. */
  mpz_t t_d;
  mpz_init(t_d);

  mpz_t pow2;
  mpz_init(pow2);

  mpz_t tmp;
  mpz_init(tmp);

  /* Compute kappa_d. */
  uint32_t kappa_d = kappa(parameters->d);

  if (kappa_d > 0) {
    mpz_set_ui(pow2, 0);
    mpz_setbit(pow2, kappa_d);
    random_generate_mpz(t_d, pow2, random_state);
  }

  /* Randomize the integer k. */
  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->l);
  random_generate_mpz(k, pow2, random_state);

  /* Compute the integer j. */
  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, kappa_d);
  mpz_div(j, parameters->d, pow2);

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->l);
  mpz_invert(j, j, pow2);

  mpz_set_ui(tmp, 0);
  mpz_setbit(tmp, parameters->m);
  mpz_mul(tmp, tmp, k);
  mpz_sub(tmp, alpha_d, tmp);

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, kappa_d);
  mpz_div(tmp, tmp, pow2);

  mpz_mul(j, j, tmp);

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->l - kappa_d);
  mpz_mul(pow2, pow2, t_d);
  mpz_add(j, j, pow2);

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->l);
  mpz_mod(j, j, pow2);

  /* Clear memory. */
  mpz_clear(pow2);
  mpz_clear(tmp);
}

void sample_j_k_from_alpha_d_r(
  mpz_t j,
  mpz_t k,
  const mpz_t alpha_d,
  const mpz_t alpha_r,
  const Parameters * const parameters,
  Random_State * const random_state)
{
  /* Declare variables. */
  mpz_t t_r;
  mpz_init(t_r);

  mpz_t pow2;
  mpz_init(pow2);

  /* Compute kappa_d, kappa_r and kappa_t_r. */
  uint32_t kappa_d = kappa(parameters->d);
  uint32_t kappa_r = kappa(parameters->r);

  uint32_t kappa_t_r =
    max_i(0, (int32_t)kappa_r - (int32_t)kappa_d - (int32_t)(parameters->l));

  /* Sample t_r at random on the interval 0 <= t_r < 2^kappa_r respecting the
   * above requirement that t_r must be a multiple of 2^(kappa_t_r). */
  if (kappa_r > 0) {
    mpz_set_ui(pow2, 0);
    mpz_setbit(pow2, kappa_r - kappa_t_r); /* pow2 = 2^(kappa_r - kappa_t_r) */
    random_generate_mpz(t_r, pow2, random_state);

    /* Multiply by 2^(kappa_t_r) to obtain 0 <= t_r < 2^kappa_r respecting the
     * above requirement that t_r must be a multiple of 2^(kappa_t_r). */
    mpz_set_ui(pow2, 0);
    mpz_setbit(pow2, kappa_t_r); /* pow2 = 2^kappa_t_r */
    mpz_mul(t_r, t_r, pow2);
  }

  /* Compute the integer j. */
  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, kappa_r); /* pow2 = 2^kappa_r */
  mpz_div(j, parameters->r, pow2); /* j = r / 2^kappa_r */

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->l); /* pow2 = 2^(m+l) */
  mpz_invert(j, j, pow2); /* j = (r / 2^kappa_r)^-1 mod 2^(m + l) */

  mpz_mul(j, j, alpha_r); /* j = alpha_r ((r / 2^kappa_r)^-1 mod 2^(m + l)) */

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, kappa_r); /* pow2 = 2^kappa_r */
  mpz_div(j, j, pow2);
    /* j = j_0 = (alpha_r / 2^kappa_r) ((r / 2^kappa_r)^-1 mod 2^(m + l)) */

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->l - kappa_r);
    /* pow2 = 2^(m + l - kappa_r) */
  mpz_mul(pow2, pow2, t_r); /* pow2 = t_r * 2^(m + l - kappa_r) */
  mpz_add(j, j, pow2); /* j = j_0 + t_r * 2^(m + l - kappa_r) */

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->l); /* pow2 = 2^(m + l) */
  mpz_mod(j, j, pow2); /* j = (j_0 + t_r * 2^(m + l - kappa_r)) mod 2^(m + l) */

  /* Compute the integer k given j to form the pair (j, k). */
  mpz_mul(k, parameters->d, j); /* k = dj */
  mpz_sub(k, alpha_d, k); /* k = alpha_d - dj */

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m); /* pow2 = 2^m */
  mpz_div(k, k, pow2); /* k = (alpha_d - dj) / 2^m */

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->l); /* pow2 = 2^l */
  mpz_mod(k, k, pow2); /* k = ((alpha_d - dj) / 2^m) mod 2^l */

  /* Clear memory. */
  mpz_clear(t_r);
  mpz_clear(pow2);
}

void sample_j_from_diagonal_alpha_r(
  mpz_t j,
  const mpz_t alpha_r,
  const Diagonal_Parameters * const parameters,
  Random_State * const random_state)
{
  /* Declare variables. */
  mpz_t t_r;
  mpz_init(t_r);

  mpz_t pow2;
  mpz_init(pow2);

  /* Compute kappa_r. */
  uint32_t kappa_r = kappa(parameters->r);

  /* Sample t_r at random on the interval 0 <= t_r < 2^kappa_r. */
  if (kappa_r > 0) {
    mpz_set_ui(pow2, 0);
    mpz_setbit(pow2, kappa_r); /* pow2 = 2^kappa_r */
    random_generate_mpz(t_r, pow2, random_state);
  }

  /* Compute the integer j. */
  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, kappa_r); /* pow2 = 2^kappa_r */
  mpz_div(j, parameters->r, pow2); /* j = r / 2^kappa_r */

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->sigma);
    /* pow2 = 2^(m + sigma) */
  mpz_invert(j, j, pow2); /* j = (r / 2^kappa_r)^-1 mod 2^(m + sigma) */

  mpz_mul(j, j, alpha_r);
    /* j = alpha_r ((r / 2^kappa_r)^-1 mod 2^(m + sigma)) */

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, kappa_r); /* pow2 = 2^kappa_r */
  mpz_div(j, j, pow2);
    /* j = j_0 = (alpha_r / 2^kappa_r) ((r / 2^kappa_r)^-1 mod 2^(m + sigma)) */

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->sigma - kappa_r);
    /* pow2 = 2^(m + sigma - kappa_r) */
  mpz_mul(pow2, pow2, t_r); /* pow2 = t_r * 2^(m + sigma - kappa_r) */
  mpz_add(j, j, pow2); /* j = j_0 + t_r * 2^(m + sigma - kappa_r) */

  mpz_set_ui(pow2, 0);
  mpz_setbit(pow2, parameters->m + parameters->sigma);
    /* pow2 = 2^(m + sigma) */
  mpz_mod(j, j, pow2);
    /* j = (j_0 + t_r * 2^(m + sigma - kappa_r)) mod 2^(m + sigma) */

  /* Clear memory. */
  mpz_clear(t_r);
  mpz_clear(pow2);
}
