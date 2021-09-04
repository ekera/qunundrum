/*!
 * \file    diagonal_probability.cpp
 * \ingroup diagonal_distribution_probability
 *
 * \brief   The definition of functions for computing the probability of
 *          observing a pair (j, k) with angle pair (theta_d, theta_r).
 */

#include "diagonal_probability.h"

#include "diagonal_parameters.h"
#include "math.h"

#include <mpfr.h>

#include <stdint.h>

void diagonal_probability_approx_f(
  mpfr_t norm,
  const mpfr_t theta_r,
  const Diagonal_Parameters * const parameters)
{
  if (0 == mpfr_cmp_ui(theta_r, 0)) {
    mpfr_set_ui(norm, 1, MPFR_RNDN);
    mpfr_div_z(norm, norm, parameters->r, MPFR_RNDN);

    return;
  }

  /* This function uses higher precision than the default precision of
   * #PRECISION internally to ensure numeric stability.
   *
   * It may be possible to reduce the precision, as the diagonal probability is
   * now computed in two steps, and as this the first step requires less
   * precision than the second step. */

  const uint32_t precision =
    2 * (max_ui(parameters->m + parameters->sigma, PRECISION));

  mpfr_t tmp;
  mpfr_init2(tmp, precision);

  mpfr_t tmp2;
  mpfr_init2(tmp2, precision);


  mpfr_set_ui_2exp(tmp, 1,
    (mpfr_exp_t)(parameters->m + parameters->sigma), MPFR_RNDN);
      /* tmp = 2^(m + sigma) */

  mpfr_div_z(tmp, tmp, parameters->r, MPFR_RNDN);
    /* tmp = 2^(m + sigma) / r */
  mpfr_mul(tmp, tmp, theta_r, MPFR_RNDN);
    /* tmp = (2^(m + sigma) / r) theta_r */
  mpfr_div_ui(tmp, tmp, 2, MPFR_RNDN);
    /* tmp = (2^(m + sigma) / r) theta_r / 2 */

  /* We use that 2 * sin(phi / 2)^2 = 1 - cos(phi) to reduce the precision
   * requirements when evaluating the probability. */

  mpfr_sin(tmp, tmp, MPFR_RNDN);
    /* tmp = sin((2^(m + sigma) / r) theta_r / 2) */
  mpfr_sqr(tmp, tmp, MPFR_RNDN);
    /* tmp = sin((2^(m + sigma) / r) theta_r / 2)^2 */
  mpfr_mul_ui(tmp, tmp, 4, MPFR_RNDN);
    /* tmp = 4 sin((2^(m + sigma) / r) theta_r / 2)^2 =
     *       2 (1 - cos((2^(m + sigma) / r) theta_r)) */

  mpfr_sqr(tmp2, theta_r, MPFR_RNDN); /* tmp2 = theta_r^2 */
  mpfr_div(tmp, tmp, tmp2, MPFR_RNDN);
    /* tmp = 2 (1 - cos((2^(m + sigma) / r) theta_r)) / theta_r^2 */

  mpfr_mul_z(tmp, tmp, parameters->r, MPFR_RNDN);
    /* tmp = r * 2 (1 - cos((2^(m + sigma) / r) theta_r)) / theta_r^2 */

  mpfr_set_ui_2exp(tmp2, 1,
    (mpfr_exp_t)(2 * (parameters->m + parameters->sigma)), MPFR_RNDN);
      /* tmp2 = 2^(2 (m + sigma)) */
  mpfr_div(norm, tmp, tmp2, MPFR_RNDN);
    /* norm = (r / 2^(2 (m + sigma))) *
                (2 (1 - cos((2^(m + sigma) / r) theta_r)) / theta_r^2) */

  /* Clear memory. */
  mpfr_clear(tmp);
  mpfr_clear(tmp2);
}

void diagonal_probability_approx_h(
  mpfr_t norm,
  const mpfr_t phi,
  const Diagonal_Parameters * const parameters)
{
  if (0 == mpfr_cmp_ui(phi, 0)) {
    mpfr_set_ui(norm, 1, MPFR_RNDN);

    return;
  }

  /* Since the probability is high when theta_r d/r is very close to theta_d,
   * this function uses higher precision than the default precision of
   * #PRECISION internally to ensure numeric stability. */

  const uint32_t precision =
    2 * (max_ui(parameters->m + parameters->sigma, PRECISION));

  mpfr_t tmp;
  mpfr_init2(tmp, precision);

  mpfr_t tmp2;
  mpfr_init2(tmp2, precision);


  mpfr_set_ui_2exp(tmp, 1, (mpfr_exp_t)(parameters->l), MPFR_RNDN);
    /* tmp = 2^l */

  mpfr_mul(tmp, phi, tmp, MPFR_RNDN); /* tmp = phi 2^l */
  mpfr_div_ui(tmp, tmp, 2, MPFR_RNDN); /* tmp = phi 2^l / 2 */
  mpfr_div_ui(tmp2, phi, 2, MPFR_RNDN); /* tmp2 = phi / 2 */


  /* We use that 2 * sin(omega / 2)^2 = 1 - cos(omega) to reduce the precision
   * requirements when evaluating the probability. */

  mpfr_sin(tmp, tmp, MPFR_RNDN);
    /* tmp = sin(phi 2^l / 2) */
  mpfr_sqr(tmp, tmp, MPFR_RNDN);
    /* tmp = sin(phi 2^l / 2)^2 */
  mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN);
    /* tmp = 2 sin(phi 2^l / 2)^2 = 1 - cos(phi 2^l) */

  mpfr_sin(tmp2, tmp2, MPFR_RNDN);
    /* tmp2 = sin(phi / 2) */
  mpfr_sqr(tmp2, tmp2, MPFR_RNDN);
    /* tmp2 = sin(phi / 2)^2 */
  mpfr_mul_ui(tmp2, tmp2, 2, MPFR_RNDN);
    /* tmp2 = 2 sin(phi / 2)^2 = 1 - cos(phi) */

  mpfr_div(tmp, tmp, tmp2, MPFR_RNDN);
    /* tmp = (1 - cos(phi 2^l)) / (1 - cos(phi 2^l)) */

  mpfr_set_ui_2exp(tmp2, 1, (mpfr_exp_t)(2 * parameters->l), MPFR_RNDN);
      /* tmp2 = 2^2l */
  mpfr_div(norm, tmp, tmp2, MPFR_RNDN);
    /* norm = tmp / 2^2l =
     *          (1 / 2^2l) * (1 - cos(phi 2^l)) / (1 - cos(phi 2^l)) */

  /* Clear memory. */
  mpfr_clear(tmp);
  mpfr_clear(tmp2);
}
