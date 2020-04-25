/*!
 * \file    diagonal_probability.cpp
 * \ingroup diagonal_distribution_probability
 *
 * \brief   The definition of functions for computing the probability of
 *          observing a pair (j, k) with angle pair (theta_d, theta_r).
 */

#include "diagonal_probability.h"

#include "parameters.h"
#include "errors.h"
#include "math.h"

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>

void diagonal_probability_approx(
  mpfr_t norm,
  const mpfr_t theta_d,
  const mpfr_t theta_r,
  const Parameters * const parameters)
{
  /* Since the probability is high when theta_r d/r is very close to theta_d, 
   * this function uses higher precision than the default precision of 
   * #PRECISION internally to ensure numeric stability. */

  const uint32_t precision =
    2 * (parameters->l + max_ui(parameters->m, PRECISION));
  
  mpfr_t tmp;
  mpfr_init2(tmp, precision);

  mpfr_t tmp2;
  mpfr_init2(tmp2, precision);

  mpfr_t tmp3;
  mpfr_init2(tmp3, precision);

  if (0 == mpfr_cmp_ui(theta_r, 0)) {
    critical("diagonal_probability(): "
      "Support for theta_r = 0 is not implemented.");
  }

  mpfr_set_ui_2exp(tmp, 1, 
    (mpfr_exp_t)(parameters->l + parameters->m), MPFR_RNDN);
      /* tmp = 2^(l + m) */

  mpfr_div_z(tmp2, tmp, parameters->r, MPFR_RNDN); /* tmp2 = 2^(l + m) / r */
  mpfr_mul(tmp2, tmp2, theta_r, MPFR_RNDN); /* tmp2 = (2^(l + m) / r) theta_r */
  mpfr_div_ui(tmp2, tmp2, 2, MPFR_RNDN);
    /* tmp2 = (2^(l + m) / r) theta_r / 2 */

  /* We use that 2 * sin(phi / 2)^2 = 1 - cos(phi) to reduce the precision 
   * requirements when evaluating the probability. */

  mpfr_sin(tmp2, tmp2, MPFR_RNDN); /* tmp = sin((2^(l + m) / r) theta_r / 2) */
  mpfr_sqr(tmp2, tmp2, MPFR_RNDN);
    /* tmp = sin((2^(l + m) / r) theta_r / 2)^2 */
  mpfr_mul_ui(tmp2, tmp2, 4, MPFR_RNDN);
    /* tmp = 4 sin((2^(l + m) / r) theta_r / 2)^2 =
     *       2 (1 - cos((2^(l + m) / r) theta_r)) */

  mpfr_sqr(tmp3, theta_r, MPFR_RNDN); /* tmp3 = theta_r^2 */
  mpfr_div(tmp3, tmp2, tmp3, MPFR_RNDN);
    /* tmp3 = 2 (1 - cos((2^(l + m) / r) theta_r)) / theta_r^2 */

  mpfr_mul_z(tmp2, theta_r, parameters->d, MPFR_RNDN); /* tmp2 = theta_r d */
  mpfr_div_z(tmp2, tmp2, parameters->r, MPFR_RNDN); /* tmp2 = theta_r d / r */
  mpfr_sub(tmp2, theta_d, tmp2, MPFR_RNDN); /* tmp2 = theta_d - theta_r d / r */

  if (0 == mpfr_cmp_ui(tmp2, 0)) {
    mpfr_set_ui_2exp(tmp, 1, 
      (mpfr_exp_t)(2 * (parameters->l + parameters->m)), MPFR_RNDN);
        /* tmp = 2^(2 (l + m)) */
  } else {
    mpfr_div_ui(tmp2, tmp2, 2, MPFR_RNDN);
      /* tmp2 = (theta_d - theta_r d / r) / 2 */

    mpfr_mul(tmp, tmp, tmp2, MPFR_RNDN);
      /* tmp = 2^(l + m) (theta_d - theta_r d / r) / 2 */
    
    /* We use that 2 * sin(phi / 2)^2 = 1 - cos(phi) to reduce the precision 
     * requirements when evaluating the probability. */

    mpfr_sin(tmp, tmp, MPFR_RNDN);
      /* tmp = sin(2^(l + m) (theta_d - theta_r d / r) / 2) */
    mpfr_sqr(tmp, tmp, MPFR_RNDN);
      /* tmp = sin(2^(l + m) (theta_d - theta_r d / r) / 2)^2 
      *     = 2 (1 - cos(2^(l + m) (theta_d - theta_r d / r) / 2)) */

    mpfr_sin(tmp2, tmp2, MPFR_RNDN);
      /* tmp2 = sin((theta_d - theta_r d / r) / 2) */
    mpfr_sqr(tmp2, tmp2, MPFR_RNDN);
      /* tmp2 = sin((theta_d - theta_r d / r) / 2)^2 
      *      = 2 (1 - cos((theta_d - theta_r d / r) / 2)) */
    
    mpfr_div(tmp, tmp, tmp2, MPFR_RNDN);
      /* tmp = (1 - cos(2^(l + m) (theta_d - theta_r d / r) / 2)) / 
      *          (1 - cos((theta_d - theta_r d / r) / 2)) */
  }

  mpfr_mul(tmp, tmp3, tmp, MPFR_RNDN);
  mpfr_mul_z(tmp, tmp, parameters->r, MPFR_RNDN);
  mpfr_set_ui_2exp(tmp2, 1, 
    (mpfr_exp_t)(4 * (parameters->l + parameters->m)), MPFR_RNDN);
  mpfr_div(norm, tmp, tmp2, MPFR_RNDN);

  /* Clear memory. */
  mpfr_clear(tmp);
  mpfr_clear(tmp2);
  mpfr_clear(tmp3);
}