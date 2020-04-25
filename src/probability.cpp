/*!
 * \file    probability.cpp
 * \ingroup two_dimensional_distribution_probability
 *
 * \brief   The definition of functions for computing the probability of
 *          observing a pair (j, k) with angle pair (theta_d, theta_r).
 */

#include "probability.h"

#include "parameters.h"
#include "math.h"

#include <mpfr.h>

#include <stdlib.h>

bool probability_approx_adjust_sigma(
  mpfr_t best_norm,
  mpfr_t best_error,
  uint32_t &best_sigma,
  const mpfr_t theta_d,
  const mpfr_t theta_r,
  const Parameters * const parameters)
{
  mpfr_t norm;
  mpfr_init2(norm, PRECISION);

  mpfr_t error;
  mpfr_init2(error, PRECISION);

  uint32_t sigma = best_sigma; /* Initially start with the best sigma. */

  /* Compute the norm and error at the starting point for sigma. */
  bool bounded_error = probability_approx(
                        norm,
                        error,
                        sigma,
                        theta_d,
                        theta_r,
                        parameters);
  mpfr_set(best_error, error, MPFR_RNDN);
  mpfr_set(best_norm, norm, MPFR_RNDN);
  
  bool best_bounded_error = bounded_error;

  /* Attempt to decrease sigma and see if we obtain better values. */
  for (uint32_t sigma = best_sigma - 1; sigma >= 1; sigma--) {
    bool bounded_error = probability_approx(
                          norm,
                          error,
                          sigma,
                          theta_d,
                          theta_r,
                          parameters);

    if (mpfr_cmp(error, best_error) >= 0) {
      if (best_sigma != sigma) {
        mpfr_clear(norm);
        mpfr_clear(error);

        return best_bounded_error;
      } else {
        break;
      }
    }

    mpfr_set(best_error, error, MPFR_RNDN);
    mpfr_set(best_norm, norm, MPFR_RNDN);
    best_bounded_error = bounded_error;
    best_sigma = sigma;
  }

  /* Attempt to increase sigma and see if we obtain better values. */
  for (uint32_t sigma = best_sigma + 1; sigma < (parameters->l - 1); sigma++) {
    bool bounded_error = probability_approx(
                          norm,
                          error,
                          sigma,
                          theta_d,
                          theta_r,
                          parameters);
    if (mpfr_cmp(error, best_error) >= 0) {
      break;
    }

    mpfr_set(best_error, error, MPFR_RNDN);
    mpfr_set(best_norm, norm, MPFR_RNDN);
    best_bounded_error = bounded_error;
    best_sigma = sigma;
  }

  /* Clear memory. */
  mpfr_clear(norm);
  mpfr_clear(error);

  return best_bounded_error;
}

bool probability_approx_optimal_sigma(
  mpfr_t norm,
  mpfr_t error,
  uint32_t &sigma,
  const mpfr_t theta_d,
  const mpfr_t theta_r,
  const Parameters * const parameters)
{
  mpfr_t best_norm;
  mpfr_init2(best_norm, PRECISION);

  mpfr_t best_error;
  mpfr_init2(best_error, PRECISION);

  uint32_t best_sigma = 0;

  bool best_bounded_error = FALSE;

  mpfr_set_ui(best_error, 1, MPFR_RNDN);

  for (sigma = 1; sigma < (parameters->l - 1); sigma++) {
    bool bounded_error = probability_approx(
                          norm,
                          error,
                          sigma,
                          theta_d,
                          theta_r,
                          parameters);

    if (mpfr_cmp(error, best_error) < 0) {
      mpfr_set(best_error, error, MPFR_RNDN);
      mpfr_set(best_norm, norm, MPFR_RNDN);
      best_sigma = sigma;
      best_bounded_error = bounded_error;
    }
  }

  mpfr_set(norm, best_norm, MPFR_RNDN);
  mpfr_set(error, best_error, MPFR_RNDN);
  sigma = best_sigma;

  /* Clear memory. */
  mpfr_clear(best_norm);
  mpfr_clear(best_error);

  return best_bounded_error;
}

bool probability_approx(
  mpfr_t norm,
  mpfr_t error,
  const uint32_t sigma,
  const mpfr_t theta_d,
  const mpfr_t theta_r,
  const Parameters * const parameters)
{
  mpfr_t tmp;
  mpfr_init2(tmp, PRECISION);

  mpfr_t tmp2;
  mpfr_init2(tmp2, PRECISION);

  /* Compute the probability estimate. */
  mpfr_set_ui_2exp(tmp, 1, (mpfr_exp_t)(sigma), MPFR_RNDN); /* tmp = 2^sigma */

  mpfr_neg(tmp2, tmp, MPFR_RNDN); /* tmp2 = -2^sigma */
  mpfr_mul_z(tmp2, tmp2, parameters->d, MPFR_RNDN); /* tmp2 = -2^sigma d */
  mpfr_div_z(tmp2, tmp2, parameters->r, MPFR_RNDN); /* tmp2 = -2^sigma d / r */
  mpfr_ceil(tmp2, tmp2); /* tmp2 = ceil(-2^sigma d / r) */
  mpfr_mul(tmp2, theta_r, tmp2, MPFR_RNDN);
    /* tmp2 = theta_r * ceil(-2^sigma d / r) */

  mpfr_mul(tmp, theta_d, tmp, MPFR_RNDN); /* tmp = theta_d 2^sigma */
  mpfr_add(tmp, tmp, tmp2, MPFR_RNDN); 
    /* tmp = theta_d 2^sigma + theta_r * ceil(-2^sigma d / r) */

  mpfr_set_ui_2exp(tmp2, 1, (mpfr_exp_t)(parameters->l - sigma), MPFR_RNDN);
    /* tmp2 = 2^(l - sigma) */
  
  if (0 == mpfr_cmp_ui(tmp, 0)) {
    mpfr_sqr(norm, tmp2, MPFR_RNDN); /* norm = 2^(2 (l - sigma)) */
  } else {
    /* We use that 2 * sin(phi / 2)^2 = 1 - cos(phi) to reduce the precision 
     * requirements when evaluating the the quotient
     * 
     * (1 - cos(phi L)) / (1 - cos(phi)) = (sin(phi L / 2) / sin(phi / 2))^2.
     * 
     * It is tricky to evaluate cos(phi) for very small phi, as cos(phi) is then
     * extremely close to 1, requiring extremely high precision in c when the 
     * value of cos(phi) is represented as c * 2^-e in the intermediary step 
     * prior to subtracting one. It is advantageous therefore to perform the 
     * computation using the expression in sin(phi / 2). */

    mpfr_div_ui(tmp, tmp, 2, MPFR_RNDN);
      /* tmp = (theta_d 2^sigma + theta_r * ceil(-2^sigma d / r)) / 2 */
    mpfr_mul(tmp2, tmp2, tmp, MPFR_RNDN);
      /* tmp2 = ((theta_d 2^sigma + 
       *            theta_r * ceil(-2^sigma d / r)) / 2) 2^(l-sigma) */
    mpfr_sin(tmp2, tmp2, MPFR_RNDN);
      /* tmp2 = sin(((theta_d 2^sigma + 
       *                theta_r * ceil(-2^sigma d / r)) / 2) 2^(l-sigma)) */
    mpfr_sin(tmp, tmp, MPFR_RNDN);
      /* tmp = sin((theta_d 2^sigma + theta_r * ceil(-2^sigma d / r)) / 2) */
    mpfr_div(tmp, tmp2, tmp, MPFR_RNDN);
      /* tmp = sin(((theta_d 2^sigma + 
       *                theta_r * ceil(-2^sigma d / r)) / 2) 2^(l-sigma)) /
       *          sin((theta_d 2^sigma + theta_r * ceil(-2^sigma d / r)) / 2) */

    mpfr_sqr(norm, tmp, MPFR_RNDN);
      /* norm = (sin(((theta_d 2^sigma + 
       *                 theta_r * ceil(-2^sigma d / r)) / 2) 2^(l-sigma)) /
       *       sin((theta_d 2^sigma + theta_r * ceil(-2^sigma d / r)) / 2))^2 */
  }

  mpfr_set_ui_2exp(tmp2, 1,
    (mpfr_exp_t)(parameters->m + parameters->l), MPFR_RNDN);
      /* tmp2 = 2^(m + l) */
  mpfr_div_z(tmp2, tmp2, parameters->r, MPFR_RNDN); /* tmp2 = 2^(m + l) / r */
  mpfr_ceil(tmp2, tmp2); /* tmp2 = ceil(2^(m + l) / r) */

  if (0 == mpfr_cmp_ui(theta_r, 0)) {
    mpfr_sqr(tmp2, tmp2, MPFR_RNDN); /* tmp2 = ceil(2^(m + l) / r)^2 */

    mpfr_mul(norm, norm, tmp2, MPFR_RNDN);
  } else {
    /* We use that 2 * sin(phi / 2)^2 = 1 - cos(phi) to reduce the precision 
     * requirements when evaluating the probability. */

    mpfr_div_ui(tmp, theta_r, 2, MPFR_RNDN); /* tmp = theta_r / 2 */
    mpfr_mul(tmp2, tmp, tmp2, MPFR_RNDN);
      /* tmp2 = ceil(2^(m + l) / r) * theta_r / 2 */
    mpfr_sin(tmp2, tmp2, MPFR_RNDN);
      /* tmp2 = sin(ceil(2^(m + l) / r) * theta_r / 2) */
    mpfr_sin(tmp, tmp, MPFR_RNDN); /* tmp = sin(theta_r / 2) */
    mpfr_div(tmp, tmp2, tmp, MPFR_RNDN);
      /* tmp = sin(ceil(2^(m + l) / r) * theta_r / 2) / sin(theta_r / 2) */
    mpfr_sqr(tmp, tmp, MPFR_RNDN);
      /* tmp = (sin(ceil(2^(m + l) / r) * theta_r / 2) / sin(theta_r / 2))^2 */

    mpfr_mul(norm, norm, tmp, MPFR_RNDN);
  }

  mpfr_set_ui_2exp(tmp2, 1, (mpfr_exp_t)(2 * sigma) - 
    (mpfr_exp_t)(2 * (parameters->m + 2 * parameters->l)), MPFR_RNDN);
      /* tmp2 = 2^(2 sigma) / 2^(2 (m + 2l)) */
  mpfr_mul_z(tmp2, tmp2, parameters->r, MPFR_RNDN);
    /* tmp2 = 2^(2 sigma) r / 2^(2 (m + 2l)) */
  mpfr_mul(norm, norm, tmp2, MPFR_RNDN);

  /* Compute the error bound given the probability estimate. */
  mpfr_abs(tmp, theta_d, MPFR_RNDN); /* tmp = abs(theta_d) */
  mpfr_abs(tmp2, theta_r, MPFR_RNDN); /* tmp2 = abs(theta_r) */
  mpfr_add(tmp, tmp, tmp2, MPFR_RNDN); /* tmp = abs(theta_d) + abs(theta_r) */

  mpfr_set_ui_2exp(tmp2, 1, (mpfr_exp_t)(sigma), MPFR_RNDN);
    /* tmp2 = 2^sigma */
  mpfr_div_ui(tmp2, tmp2, 2, MPFR_RNDN); /* tmp2 = 2^sigma / 2 */
  mpfr_mul(tmp, tmp2, tmp, MPFR_RNDN);
    /* tmp = (2^sigma / 2) (abs(theta_d) + abs(theta_r)) */
  mpfr_add_ui(tmp2, tmp, 2, MPFR_RNDN);
    /* tmp2 = 2 + (2^sigma / 2) (abs(theta_d) + abs(theta_r)) */
  mpfr_mul(tmp, tmp, tmp2, MPFR_RNDN);
    /* tmp2 = (2^sigma / 2) (abs(theta_d) + abs(theta_r)) *
     *          (2 + (2^sigma / 2) (abs(theta_d) + abs(theta_r)) */
  mpfr_mul(tmp, tmp, norm, MPFR_RNDN);
    /* tmp2 = (2^sigma / 2) (abs(theta_d) + abs(theta_r)) *
     *          (2 + (2^sigma / 2) (abs(theta_d) + abs(theta_r)) * norm */
  
  mpfr_set_ui_2exp(tmp2, 1,
    (mpfr_exp_t)(4) - (mpfr_exp_t)(parameters->m + sigma), MPFR_RNDN);
      /* tmp2 = 2^4 / 2^(m + sigma) */
  mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
  mpfr_set_ui_2exp(tmp2, 1, 
    (mpfr_exp_t)(3) - (mpfr_exp_t)(parameters->m + parameters->l), MPFR_RNDN);
      /* tmp2 = 2^3 / 2^(m + l) */
  mpfr_add(error, tmp, tmp2, MPFR_RNDN);
  
  /* Check the stated bound on the relative error. */
  mpfr_div(tmp, error, norm, MPFR_RNDN);
  const bool bounded_error = (mpfr_cmp_d(tmp, (double)0.01f) <= 0);

  /* Clear memory. */
  mpfr_clear(tmp);
  mpfr_clear(tmp2);
  
  return bounded_error ? TRUE : FALSE;
}

void probability_approx_quick(
  mpfr_t norm,
  const mpfr_t theta_d,
  const mpfr_t theta_r,
  const Parameters * const parameters)
{
  mpfr_t tmp;
  mpfr_init2(tmp, PRECISION);

  mpfr_t tmp2;
  mpfr_init2(tmp2, PRECISION);

  mpfr_mul_z(tmp, theta_r, parameters->d, MPFR_RNDN); /* tmp = theta_r d */
  mpfr_div_z(tmp, tmp, parameters->r, MPFR_RNDN); /* tmp = theta_r d / r */
  mpfr_sub(tmp, theta_d, tmp, MPFR_RNDN); /* tmp = theta_d - theta_r d / r */

  if (0 == mpfr_cmp_ui(tmp, 0)) {
    mpfr_set_ui_2exp(norm, 1, (mpfr_exp_t)(2 * parameters->l), MPFR_RNDN);
      /* norm = 2^(2l) */
  } else {
    mpfr_div_ui(tmp, tmp, 2, MPFR_RNDN);
      /* tmp = (theta_d - theta_r d / r) / 2 */

    mpfr_set_ui_2exp(tmp2, 1, parameters->l, MPFR_RNDN); /* tmp2 = 2^l */
    mpfr_mul(tmp2, tmp, tmp2, MPFR_RNDN);
      /* tmp2 = 2^l (theta_d - theta_r d / r) / 2 */
    
    /* We use that 2 * sin(phi / 2)^2 = 1 - cos(phi) to reduce the precision 
     * requirements when evaluating the probability. */

    mpfr_sin(tmp, tmp, MPFR_RNDN);
      /* tmp = sin((theta_d - theta_r d / r) / 2) */
    mpfr_sin(tmp2, tmp2, MPFR_RNDN);
      /* tmp2 = sin(2^l (theta_d - theta_r d / r) / 2) */
    mpfr_div(tmp, tmp2, tmp, MPFR_RNDN);
      /* tmp = sin(2^l (theta_d - theta_r d / r) / 2) / 
       *          sin((theta_d - theta_r d / r) / 2) */
    mpfr_sqr(norm, tmp, MPFR_RNDN);
      /* norm = (sin(2^l (theta_d - theta_r d / r) / 2) / 
       *            sin((theta_d - theta_r d / r) / 2))^2 */
  }

  mpfr_set_ui_2exp(tmp, 1,
    (mpfr_exp_t)(parameters->m + parameters->l), MPFR_RNDN);
      /* tmp = 2^(m + l) */
  mpfr_div_z(tmp, tmp, parameters->r, MPFR_RNDN); /* tmp = 2^(m + l) / r */
  mpfr_ceil(tmp, tmp); /* tmp = (ceil(2^(m + l) / r) */

  if (0 == mpfr_cmp_ui(theta_r, 0)) {
    mpfr_sqr(tmp2, tmp, MPFR_RNDN); /* tmp2 = (ceil(2^(m + l) / r))^2 */

    mpfr_mul(norm, norm, tmp2, MPFR_RNDN);
  } else {
    mpfr_div_ui(tmp2, theta_r, 2, MPFR_RNDN); /* tmp2 = theta_r / 2 */
    mpfr_mul(tmp, tmp, tmp2, MPFR_RNDN);
      /* tmp = ceil(2^(m + l) / r) theta_r / 2 */

    /* We use that 2 * sin(phi / 2)^2 = 1 - cos(phi) to reduce the precision 
     * requirements when evaluating the probability. */

    mpfr_sin(tmp, tmp, MPFR_RNDN);
      /* tmp = sin(ceil(2^(m + l) / r) theta_r / 2) */
    mpfr_sin(tmp2, tmp2, MPFR_RNDN);
      /* tmp2 = sin(theta_r / 2) */
    
    mpfr_div(tmp, tmp, tmp2, MPFR_RNDN);
      /* tmp = sin(ceil(2^(m + l) / r) theta_r / 2) / sin(theta_r / 2) */
    mpfr_sqr(tmp, tmp, MPFR_RNDN);
      /* tmp = (sin(ceil(2^(m + l) / r) theta_r / 2) / sin(theta_r / 2))^2 */

    mpfr_mul(norm, norm, tmp, MPFR_RNDN);
  }

  mpfr_mul_z(norm, norm, parameters->r, MPFR_RNDN);
  mpfr_set_ui_2exp(tmp, 1,
    (mpfr_exp_t)(2 * (parameters->m + 2 * parameters->l)), MPFR_RNDN);
      /* tmp = 2^(2(m + 2l)) */
  mpfr_div(norm, norm, tmp, MPFR_RNDN);
  
  /* Clear memory. */
  mpfr_clear(tmp);
  mpfr_clear(tmp2);
}