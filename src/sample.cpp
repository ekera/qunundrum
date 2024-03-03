/*!
 * \file    sample.cpp
 * \ingroup sample_distribution
 *
 * \brief   The definition of functions for sampling probability distributions.
 */

#include "sample.h"

#include "common.h"
#include "diagonal_parameters.h"
#include "diagonal_probability.h"
#include "errors.h"
#include "math.h"
#include "parameters.h"
#include "random.h"

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

bool sample_k_from_diagonal_j_eta_pivot(
  const Diagonal_Parameters * const parameters,
  long double pivot,
  const mpz_t j,
  const int32_t eta,
  const uint32_t delta_bound,
  mpz_t k,
  mpfr_t alpha_phi)
{
  /* Sanity checks. */
  if ((pivot < 0) || (pivot > 1)) {
    critical("sample_k_from_diagonal_j_eta_pivot(): "
      "The pivot is out of bounds.");
  }

  /* Setup precision. */
  const uint32_t precision =
    3 * (max_ui(parameters->m + parameters->sigma, PRECISION));

  /* Declare variables. */
  mpz_t alpha_r;
  mpz_init(alpha_r);

  mpz_t k0;
  mpz_init(k0);

  mpz_t tmp_z;
  mpz_init(tmp_z);

  mpfr_t tmp;
  mpfr_init2(tmp, precision);

  /* Declare constants. */
  mpz_t pow2_m_sigma;
  mpz_init(pow2_m_sigma);

  mpz_set_ui(pow2_m_sigma, 0); /* pow2_m_sigma = 0 */
  mpz_setbit(pow2_m_sigma,
    parameters->m +
    parameters->sigma); /* pow2_m_sigma = 2^(m + sigma) */

  mpz_t pow2_m_sigma_l;
  mpz_init(pow2_m_sigma_l);

  mpz_set_ui(pow2_m_sigma_l, 0); /* pow2_m_sigma_l = 0 */
  mpz_setbit(pow2_m_sigma_l,
    parameters->m +
    parameters->sigma -
    parameters->l); /* pow2_m_sigma_l = 2^(m + sigma - l) */

  mpz_t pow2_l;
  mpz_init(pow2_l);

  mpz_set_ui(pow2_l, 0); /* pow2_l = 0 */
  mpz_setbit(pow2_l, parameters->l); /* pow2_l = 2^l */

  /* Compute alpha_r. */
  mpz_mul(alpha_r, parameters->r, j); /* alpha_r = rj */
  mod_reduce(alpha_r, pow2_m_sigma); /* alpha_r = {rj}_{2^(m + sigma)} */

  /* Compute k0. */
  mpz_mul_si(tmp_z, pow2_m_sigma, eta); /* tmp_z = 2^(m + sigma) eta */
  mpz_sub(tmp_z, alpha_r, tmp_z); /* tmp_z = alpha_r - 2^(m + sigma) eta */

  mpfr_set_z(tmp, tmp_z, MPFR_RNDN); /* tmp = alpha_r - 2^(m + sigma) eta */
  mpfr_mul_z(tmp, tmp, parameters->d, MPFR_RNDN);
    /* tmp = d * (alpha_r - 2^(m + sigma) eta) */
  mpfr_div_z(tmp, tmp, parameters->r, MPFR_RNDN);
    /* tmp = (d / r) * (alpha_r - 2^(m + sigma) eta) */

  mpz_mul(tmp_z, parameters->d, j); /* tmp_z = d * j */
  mpfr_sub_z(tmp, tmp, tmp_z, MPFR_RNDN);
    /* tmp = (d / r) * (alpha_r - 2^(m + sigma) eta) - d * j */

  mpfr_div_z(tmp, tmp, pow2_m_sigma_l, MPFR_RNDN);
    /* tmp = ((d / r) * (alpha_r - 2^(m + sigma) eta) - d * j) /
               2^(m + sigma - l) */

  mpfr_round(tmp, tmp);
  mpfr_get_z(k0, tmp, MPFR_RNDN); /* k0 = round(tmp) */
  mpz_mod(k0, k0, pow2_l); /* k0 = k0 mod 2^l */

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("sample_k_from_diagonal_j_eta_pivot(): "
    "Debug: Computed k0: %Zd\n", k0);
  #endif

  /* Compute phi and evaluate h(phi). Pre-compute some values. */
  mpfr_t phi;
  mpfr_init2(phi, precision);

  mpfr_t scale;
  mpfr_init2(scale, precision);

  mpfr_const_pi(scale, MPFR_RNDN); /* scale = pi */
  mpfr_mul_ui(scale, scale, 2, MPFR_RNDN); /* scale = 2 pi */
  mpfr_div_z(scale, scale, pow2_m_sigma, MPFR_RNDN);
    /* scale = 2 pi / 2^(m + sigma) */

  mpfr_t term_j_eta;
  mpfr_init2(term_j_eta, precision);

  mpz_mul_si(tmp_z, pow2_m_sigma, eta); /* tmp_z = 2^(m + sigma) eta */
  mpz_sub(tmp_z, alpha_r, tmp_z); /* tmp_z = alpha_r - 2^(m + sigma) eta */

  mpfr_set_z(term_j_eta, tmp_z, MPFR_RNDN);
    /* tmp = alpha_r - 2^(m + sigma) eta */
  mpfr_mul_z(term_j_eta, term_j_eta, parameters->d, MPFR_RNDN);
    /* term_j_eta = d * (alpha_r - 2^(m + sigma) eta) */
  mpfr_div_z(term_j_eta, term_j_eta, parameters->r, MPFR_RNDN);
    /* term_j_eta = (d / r) * (alpha_r - 2^(m + sigma) eta) */

  mpz_mul(tmp_z, parameters->d, j); /* tmp_z = d j */
  mpfr_sub_z(term_j_eta, term_j_eta, tmp_z, MPFR_RNDN);
    /* term_j_eta = (d / r) * (alpha_r - 2^(m + sigma) eta) - d j */
  mpfr_neg(term_j_eta, term_j_eta, MPFR_RNDN);
    /* term_j_eta = d j - (d / r) * (alpha_r - 2^(m + sigma) eta) */

  #ifdef DEBUG_TRACE_SAMPLING
  mpfr_printf("sample_k_from_diagonal_j_eta_pivot(): "
    "Debug: Computed term_j_eta: %Rg\n", term_j_eta);
  #endif

  for (uint32_t delta_abs = 0; delta_abs <= delta_bound; delta_abs++) {
    for (int32_t delta_sgn = 1; delta_sgn >= -1; delta_sgn -= 2) {

      if ((0 == delta_abs) && (-1 == delta_sgn)) {
        continue;
      }

      /* Compute k. */
      if (1 == delta_sgn) {
        mpz_add_ui(k, k0, delta_abs); /* k = k0 + delta = k0 + delta_abs */
      } else {
        mpz_sub_ui(k, k0, delta_abs); /* k = k0 + delta = k0 - delta_abs */
      }

      mpz_mod(k, k, pow2_l); /* k = k mod 2^l */

      /* Compute phi. */
      mpz_mul(tmp_z, pow2_m_sigma_l, k); /* tmp_z = 2^(m + sigma - l) k */
      mpfr_add_z(phi, term_j_eta, tmp_z, MPFR_RNDN);
        /* phi = dj + 2^(m + sigma - l) k -
                   (d / r) (alpha_r - 2^(m + sigma) eta) */

      mpfr_set_z(tmp, pow2_m_sigma, MPFR_RNDN); /* tmp = 2^(m + sigma) */
      mpfr_fmod(phi, phi, tmp, MPFR_RNDN);
        /* phi = (dj + 2^(m + sigma - l) k -
                   (d / r) (alpha_r - 2^(m + sigma) eta)) mod 2^(m + sigma) */
      mpfr_div_ui(tmp, tmp, 2, MPFR_RNDN); /* tmp = 2^(m + sigma) / 2 */
      if (mpfr_cmp(phi, tmp) >= 0) {
        mpfr_sub_z(phi, phi, pow2_m_sigma, MPFR_RNDN);
      }

      /* phi = {dj + 2^(m + sigma - l) k -
                 (d / r) (alpha_r - 2^(m + sigma) eta)}_{2^(m + sigma)} */

      if (NULL != alpha_phi) {
        mpfr_set(alpha_phi, phi, MPFR_RNDN);
      }

      mpfr_mul(phi, phi, scale, MPFR_RNDN);
        /* phi = (2 pi / 2^(m + sigma)) {dj + 2^(m + sigma - l) k -
                   (d / r) (alpha_r - 2^(m + sigma) eta)}_{2^(m + sigma)} */

      #ifdef DEBUG_TRACE_SAMPLING
      mpfr_printf("sample_k_from_diagonal_j_eta_pivot(): "
        "Debug: For Delta = %d, computed phi: %Rg\n",
          ((int32_t)delta_abs) * delta_sgn, phi);
      #endif

      diagonal_probability_approx_h(tmp, phi, parameters);

      pivot -= mpfr_get_ld(tmp, MPFR_RNDN);

      #ifdef DEBUG_TRACE_SAMPLING
      printf("sample_k_from_diagonal_j_eta_pivot(): "
        "Debug: Decremented pivot to: %Lf\n", pivot);
      #endif

      if (pivot <= 0) {
        break;
      }
    }

    if (pivot <= 0) {
      break;
    }
  }

  /* Clear memory. */
  mpz_clear(alpha_r);

  mpz_clear(k0);
  mpz_clear(tmp_z);

  mpz_clear(pow2_m_sigma);
  mpz_clear(pow2_m_sigma_l);
  mpz_clear(pow2_l);

  mpfr_clear(phi);

  mpfr_clear(scale);
  mpfr_clear(term_j_eta);

  mpfr_clear(tmp);

  if (pivot > 0) {
    #ifdef DEBUG_TRACE_SAMPLING
    printf("sample_k_from_diagonal_j_eta_pivot(): "
      "Debug: Sampled k given j and eta: Out of bounds.\n");
    #endif

    mpz_set_ui(k, 0);

    if (NULL != alpha_phi) {
      mpfr_set_ui(alpha_phi, 0, MPFR_RNDN);
    }

    /* Signal sampling error in h. */
    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("sample_k_from_diagonal_j_eta_pivot(): "
    "Debug: Sampled k: %Zd\n", k);
  #endif

  /* Signal success. */
  return TRUE;
}

bool sample_k_from_diagonal_j_eta(
  const Diagonal_Parameters * const parameters,
  Random_State * const random_state,
  const mpz_t j,
  const int32_t eta,
  const uint32_t delta_bound,
  mpz_t k,
  mpfr_t alpha_phi)
{
  /* Select a pivot uniformly at random. */
  long double pivot = random_generate_pivot_inclusive(random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("sample_k_from_diagonal_j_eta(): Debug: Sampled pivot: %Lf\n", pivot);
  #endif

  /* Call and signal result. */
  return sample_k_from_diagonal_j_eta_pivot(
            parameters,
            pivot,
            j,
            eta,
            delta_bound,
            k,
            alpha_phi);
}