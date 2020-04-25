/*!
 * \file    tau_estimate.cpp
 * \ingroup estimating_tau
 *
 * \brief   The definition of functions for computing tau estimates that may be
 *          converted into volume quotient estimates.
 */

#include "tau_estimate.h"

#include "distribution.h"
#include "linear_distribution.h"
#include "random.h"
#include "errors.h"
#include "common.h"

#include <mpfr.h>
#include <gmp.h>

#include <float.h>

#include <stdint.h>
#include <stdlib.h>

bool tau_estimate(
  const Distribution * const distribution,
  Random_State * const random_state,
  const uint32_t n,
  long double &tau_d,
  long double &tau_r)
{
  mpfr_t alpha_d;
  mpfr_init2(alpha_d, PRECISION);

  mpfr_t alpha_r;
  mpfr_init2(alpha_r, PRECISION);

  mpfr_t sum_d;
  mpfr_init2(sum_d, PRECISION);
  mpfr_set_ui(sum_d, 0, MPFR_RNDN);

  mpfr_t sum_r;
  mpfr_init2(sum_r, PRECISION);
  mpfr_set_ui(sum_r, 0, MPFR_RNDN);

  bool result = FALSE;

  for (uint32_t i = 0; i < n; i++) {
    result = distribution_sample_approximate_alpha_d_r(
                distribution,
                random_state,
                alpha_d,
                alpha_r);
    if (FALSE == result) {
      /* Break and return maximal values for tau_d and tau_r below. */
      break;
    }

    mpfr_sqr(alpha_d, alpha_d, MPFR_RNDN);
    mpfr_add(sum_d, sum_d, alpha_d, MPFR_RNDN);

    mpfr_sqr(alpha_r, alpha_r, MPFR_RNDN);
    mpfr_add(sum_r, sum_r, alpha_r, MPFR_RNDN);
  }

  if (FALSE == result) {
    /* Return maximal values for tau_d and tau_r. */
    tau_d = DBL_MAX;
    tau_r = DBL_MAX;
  } else {
    mpfr_div_ui(sum_d, sum_d, n, MPFR_RNDN);
    mpfr_log2(sum_d, sum_d, MPFR_RNDN);
    tau_d = mpfr_get_ld(sum_d, MPFR_RNDN) / 2 - distribution->parameters.m;

    mpfr_div_ui(sum_r, sum_r, n, MPFR_RNDN);
    mpfr_log2(sum_r, sum_r, MPFR_RNDN);
    tau_r = mpfr_get_ld(sum_r, MPFR_RNDN) / 2 - distribution->parameters.m;
  }

  /* Clear memory. */
  mpfr_clear(alpha_d);
  mpfr_clear(alpha_r);

  mpfr_clear(sum_d);
  mpfr_clear(sum_r);

  /* Signal success or failure. */
  return result;
}

bool tau_estimate_linear(
  const Linear_Distribution * const distribution,
  Random_State * const random_state,
  const uint32_t n,
  long double &tau)
{
  mpfr_t alpha;
  mpfr_init2(alpha, PRECISION);

  mpfr_t sum;
  mpfr_init2(sum, PRECISION);
  mpfr_set_ui(sum, 0, MPFR_RNDN);

  bool result = FALSE;

  for (uint32_t i = 0; i < n; i++) {
    result = linear_distribution_sample_approximate_alpha(
      distribution,
      random_state,
      alpha);
    if (FALSE == result) {
      /* Break and return the maximal value for tau below. */
      break;
    }

    mpfr_sqr(alpha, alpha, MPFR_RNDN);
    mpfr_add(sum, sum, alpha, MPFR_RNDN);
  }

  if (FALSE == result) {
    /* Return the maximal value for tau. */
    tau = DBL_MAX;
  } else {
    mpfr_div_ui(sum, sum, n, MPFR_RNDN);
    mpfr_log2(sum, sum, MPFR_RNDN);
    tau = mpfr_get_ld(sum, MPFR_RNDN) / 2.0f - distribution->parameters.m;
  }

  /* Clear memory. */
  mpfr_clear(alpha);
  mpfr_clear(sum);

  /* Signal success or failure. */
  return result;
}
