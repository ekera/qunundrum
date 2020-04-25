/*!
 * \file    linear_probability.cpp
 * \ingroup linear_distribution_probability
 *
 * \brief   The definition of functions for computing the probability of
 *          observing a pair (j, k) with angle theta_d, in the short discrete
 *          logarithm case, or an integer j with angle theta_r in the order
 *          case.
 */

#include "linear_probability.h"

#include "parameters.h"
#include "math.h"

#include <mpfr.h>
#include <gmp.h>

#include <stdint.h>

void linear_probability_d(
  mpfr_t norm,
  const mpfr_t theta_d,
  const Parameters * const parameters)
{
  const uint32_t precision = 3 * max_ui(parameters->m, PRECISION);

  mpfr_t tmp;
  mpfr_init2(tmp, precision);
  
  mpfr_t tmp2;
  mpfr_init2(tmp2, precision);

  mpfr_t tmp3;
  mpfr_init2(tmp3, precision);

  if (0 == mpfr_cmp_ui(theta_d, 0)) {
    mpfr_set_ui_2exp(tmp, 1, (mpfr_exp_t)(parameters->l), MPFR_RNDN);
      /* tmp = 2^l */
    mpfr_sub_ui(tmp, tmp, 1, MPFR_RNDN); /* tmp = 2^l - 1 */
    mpfr_mul_z(tmp, tmp, parameters->d, MPFR_RNDN); /* tmp = d (2^l - 1) */

    mpfr_set_ui_2exp(tmp2, 1, (mpfr_exp_t)(parameters->l + 1), MPFR_RNDN);
      /* tmp2 = 2^(l + 1) */
    mpfr_sub_ui(tmp2, tmp2, 1, MPFR_RNDN); /* tmp2 = 2^(l + 1) - 1 */
    mpfr_mul(tmp3, tmp, tmp2, MPFR_RNDN); 
      /* tmp3 = d (2^l - 1) (2^(l + 1) - 1) */
    mpfr_set_ui_2exp(tmp2, 1, (mpfr_exp_t)(parameters->l), MPFR_RNDN);
      /* tmp2 = 2^l */
    mpfr_mul(tmp3, tmp3, tmp2, MPFR_RNDN);
      /* tmp3 = 2^l d (2^l - 1) (2^(l + 1) - 1) */
    mpfr_div_ui(tmp3, tmp3, 3, MPFR_RNDN);
      /* tmp3 = (2^l d / 3) (2^l - 1) (2^(l + 1) - 1) */

    mpfr_set_ui_2exp(tmp2, 1,
      (mpfr_exp_t)(parameters->l + parameters->m), MPFR_RNDN);
        /* tmp2 = 2^(l + m) */
    mpfr_sub(tmp, tmp2, tmp, MPFR_RNDN);
      /* tmp = 2^(l + m) - d (2^l - 1) */
    mpfr_set_ui_2exp(tmp2, 1, (mpfr_exp_t)(2 * parameters->l), MPFR_RNDN);
      /* tmp2 = 2^(2l) */
    mpfr_mul(tmp2, tmp, tmp2, MPFR_RNDN);
      /* tmp2 = (2^(l + m) - d (2^l - 1)) 2^(2l) */
    mpfr_add(tmp3, tmp3, tmp2, MPFR_RNDN);
      /* tmp3 = (2^(l + m) - d (2^l - 1)) 2^(2l) + 
       *           (2^l d / 3) (2^l - 1) (2^(l + 1) - 1) */

    mpfr_set_ui_2exp(tmp, 1, 
      (mpfr_exp_t)(2 * (parameters->m + 2 * parameters->l)), MPFR_RNDN);
        /* tmp = 2^(2(m + 2l)) */
    mpfr_div(norm, tmp3, tmp, MPFR_RNDN);
      /* norm = (1 / 2^(2(m + 2l))) (2^(l + m) - d (2^l - 1)) 2^(2l) + 
       *                               (2^l d / 3) (2^l - 1) (2^(l + 1) - 1) */
    
    /* Clean up. */
    mpfr_clear(tmp);
    mpfr_clear(tmp2);
    mpfr_clear(tmp3);
  }

  mpfr_t one_minus_cos_theta; /* := 1 - cos(theta_d) */
  mpfr_init2(one_minus_cos_theta, precision);

  /* We use that 2 * sin(phi / 2)^2 = 1 - cos(phi) to reduce the precision 
   * requirements when evaluating the probability. */

  mpfr_div_ui(tmp, theta_d, 2, MPFR_RNDN); /* tmp = theta_d / 2 */
  mpfr_sin(tmp2, tmp, MPFR_RNDN); /* tmp2 = sin(theta_d / 2) */
  mpfr_sqr(tmp2, tmp2, MPFR_RNDN); /* tmp2 = sin(theta_d / 2)^2 */
  mpfr_mul_ui(one_minus_cos_theta, tmp2, 2, MPFR_RNDN);
    /* one_minus_cos_theta = 2 sin(theta_d / 2)^2 = 1 - cos(theta_d) */

  mpfr_t one_minus_cos_2l_theta; /* := 1 - cos(2^l theta_d) */
  mpfr_init2(one_minus_cos_2l_theta, precision);

  mpfr_set_ui_2exp(tmp3, 1, (mpfr_exp_t)(parameters->l), MPFR_RNDN);
    /* tmp3 = 2^l */
  mpfr_mul(tmp2, tmp3, tmp, MPFR_RNDN); /* tmp2 = 2^l (theta_d / 2) */
  mpfr_sin(tmp2, tmp2, MPFR_RNDN); /* tmp2 = sin(2^l (theta_d / 2)) */
  mpfr_sqr(tmp2, tmp2, MPFR_RNDN); /* tmp2 = sin(2^l (theta_d / 2))^2 */
  mpfr_mul_ui(one_minus_cos_2l_theta, tmp2, 2, MPFR_RNDN);
    /* one_minus_cos_2l_theta = 
     *    2 sin(2^l (theta_d / 2))^2 = 1 - cos(2^l theta_d) */

  mpfr_sub_ui(tmp3, tmp3, 1, MPFR_RNDN); /* tmp3 = 2^l - 1 */
  mpfr_mul(tmp2, tmp3, tmp, MPFR_RNDN); /* tmp2 = (2^l - 1) (theta_d / 2) */
  mpfr_sin(tmp2, tmp2, MPFR_RNDN); /* tmp2 = sin((2^l - 1) (theta_d / 2)) */
  mpfr_sqr(tmp2, tmp2, MPFR_RNDN); /* tmp2 = sin((2^l - 1) (theta_d / 2))^2 */
  mpfr_mul_ui(tmp2, tmp2, 2, MPFR_RNDN);
    /* tmp2 = 2 sin((2^l - 1) (theta_d / 2))^2 = 1 - cos((2^l - 1) theta_d) */

  mpfr_sub(tmp2, one_minus_cos_2l_theta, tmp2, MPFR_RNDN);
    /* tmp2 = cos((2^l - 1) theta_d) - cos(2^l theta_d) */
  mpfr_div(tmp2, tmp2, one_minus_cos_theta, MPFR_RNDN);
    /* tmp2 = (cos((2^l - 1) theta_d) - cos(2^l theta_d)) / 
     *           (1 - cos(theta_d)) */
  mpfr_sub_ui(tmp2, tmp2, 1, MPFR_RNDN);
    /* tmp2 = (cos((2^l - 1) theta_d) - cos(2^l theta_d)) / 
     *           (1 - cos(theta_d)) - 1 */
  mpfr_div_ui(tmp2, tmp2, 2, MPFR_RNDN);
    /* tmp2 = (1/2) ((cos((2^l - 1) theta_d) - cos(2^l theta_d)) / 
                      (1 - cos(theta_d)) - 1) */
  mpfr_sub(tmp2, tmp3, tmp2, MPFR_RNDN);
    /* tmp2 = (2^l - 1) - (1/2) ((cos((2^l - 1) theta_d) - cos(2^l theta_d)) / 
                                  (1 - cos(theta_d)) - 1) */
  mpfr_mul_ui(tmp2, tmp2, 2, MPFR_RNDN);
    /* tmp2 = 2 ( (2^l - 1) - (1/2) ((cos((2^l - 1) theta_d) - 
     *           cos(2^l theta_d)) / (1 - cos(theta_d)) - 1) ) */
  mpfr_mul_z(tmp2, tmp2, parameters->d, MPFR_RNDN);
    /* tmp2 = 2d ( (2^l - 1) - (1/2) ((cos((2^l - 1) theta_d) - 
     *           cos(2^l theta_d)) / (1 - cos(theta_d)) - 1) ) */

  mpfr_mul_z(tmp3, tmp3, parameters->d, MPFR_RNDN); /* tmp3 = (2^l - 1) d */
  mpfr_set_ui_2exp(tmp, 1, 
    (mpfr_exp_t)(parameters->l + parameters->m), MPFR_RNDN);
      /* tmp = 2^(l + m) */
  mpfr_sub(tmp3, tmp, tmp3, MPFR_RNDN); /* tmp2 = 2^(l + m) - (2^l - 1) d */
  mpfr_mul(tmp3, tmp3, one_minus_cos_2l_theta, MPFR_RNDN);
    /* tmp2 = (2^(l + m) - (2^l - 1) d) (1 - cos(2^l theta_d)) */
  
  mpfr_add(tmp, tmp3, tmp2, MPFR_RNDN);
    /* tmp = (2^(l + m) - (2^l - 1) d) (1 - cos(2^l theta_d)) + 
     *    2d ( (2^l - 1) - (1/2) ((cos((2^l - 1) theta_d) - cos(2^l theta_d)) / 
     *       (1 - cos(theta_d)) - 1) ) */

  /* Normalize. */
  mpfr_div(tmp, tmp, one_minus_cos_theta, MPFR_RNDN);

  mpfr_set_ui_2exp(tmp2, 1, 
    (mpfr_exp_t)(2 * (parameters->m + 2 * parameters->l)), MPFR_RNDN);
      /* tmp2 = 2^(2(m + 2l)) */
  mpfr_div(norm, tmp, tmp2, MPFR_RNDN);
    /* norm = (2^(l + m) - (2^l - 1) d) (1 - cos(2^l theta_d)) + 
     *    2d ( (2^l - 1) - (1/2) ((cos((2^l - 1) theta_d) - cos(2^l theta_d)) / 
     *       (1 - cos(theta_d)) - 1) ) / 2^(2(m + 2l)) */

  /* Clean up. */
  mpfr_clear(one_minus_cos_theta);
  mpfr_clear(one_minus_cos_2l_theta);

  mpfr_clear(tmp);
  mpfr_clear(tmp2);
  mpfr_clear(tmp3);
}

void linear_probability_r(
  mpfr_t norm,
  const mpfr_t theta_r,
  const Parameters * const parameters)
{
  mpz_t beta;
  mpz_init(beta);

  mpfr_t N;
  mpfr_init2(N, PRECISION);

  mpfr_t tmp;
  mpfr_init2(tmp, PRECISION);

  mpfr_t tmp2;
  mpfr_init2(tmp2, PRECISION);

  mpfr_t tmp3;
  mpfr_init2(tmp3, PRECISION);

  mpz_set_ui(beta, 0); /* beta = 0 */
  mpz_setbit(beta, parameters->l + parameters->m); /* beta = 2^(l + m) */
  mpz_mod(beta, beta, parameters->r); /* beta = 2^(l + m) mod r */

  mpfr_set_ui_2exp(N, 1, 
    (mpfr_exp_t)(parameters->l + parameters->m), MPFR_RNDN); /* N = 2^(l + m) */
  mpfr_div_z(N, N, parameters->r, MPFR_RNDN); /* N = 2^(l + m) / r */
  mpfr_floor(N, N); /* N = floor(2^(l + m) / r) */

  if (0 == mpfr_cmp_ui(theta_r, 0)) {
    mpfr_add_ui(tmp, N, 1, MPFR_RNDN); /* tmp = N + 1 */
    mpfr_sqr(tmp, tmp, MPFR_RNDN); /* tmp = (N + 1)^2 */
    mpfr_mul_z(tmp, tmp, beta, MPFR_RNDN); /* tmp = beta * (N + 1)^2 */

    mpz_sub(beta, parameters->r, beta); /* beta' = r - beta */
    
    mpfr_sqr(tmp2, N, MPFR_RNDN); /* tmp2 = N^2 */
    mpfr_mul_z(tmp2, tmp2, beta, MPFR_RNDN); /* tmp2 = (r - beta) N^2 */

    mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
      /* tmp = beta * (N + 1)^2 + (r - beta) N^2 */
  } else {
    mpfr_div_ui(tmp3, theta_r, 2, MPFR_RNDN); /* tmp3 = theta_r / 2 */

    /* We use that 2 * sin(phi / 2)^2 = 1 - cos(phi) to reduce the precision 
     * requirements when evaluating the probability. */

    mpfr_add_ui(tmp, N, 1, MPFR_RNDN); /* tmp = N + 1 */
    mpfr_mul(tmp, tmp, tmp3, MPFR_RNDN); /* tmp = (N + 1) * theta_r / 2 */
    mpfr_sin(tmp, tmp, MPFR_RNDN); /* tmp = sin((N + 1) * theta_r / 2) */

    mpfr_sin(tmp2, tmp3, MPFR_RNDN); /* tmp2 = sin(theta_r / 2) */

    mpfr_mul(tmp3, tmp3, N, MPFR_RNDN); /* tmp3 = N * theta_r / 2 */
    mpfr_sin(tmp3, tmp3, MPFR_RNDN); /* tmp3 = sin(N * theta_r / 2) */

    mpfr_div(tmp, tmp, tmp2, MPFR_RNDN);
      /* tmp = sin((N + 1) * theta_r / 2) / sin(theta_r / 2) */
    mpfr_sqr(tmp, tmp, MPFR_RNDN);
      /* tmp = (sin((N + 1) * theta_r / 2) / sin(theta_r / 2))^2 */

    mpfr_div(tmp3, tmp3, tmp2, MPFR_RNDN);
      /* tmp3 = sin(N * theta_r / 2) / sin(theta_r / 2) */
    mpfr_sqr(tmp3, tmp3, MPFR_RNDN);
      /* tmp3 = (sin(N * theta_r / 2) / sin(theta_r / 2))^2 */

    mpfr_mul_z(tmp, tmp, beta, MPFR_RNDN);
      /* tmp = beta * (sin((N + 1) * theta_r / 2) / sin(theta_r / 2))^2 */
    mpz_sub(beta, parameters->r, beta); /* beta' = r - beta */
    mpfr_mul_z(tmp3, tmp3, beta, MPFR_RNDN);
      /* tmp3 = (r - beta) * (sin(N * theta_r / 2) / sin(theta_r / 2))^2 */

    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
      /* tmp = beta * (sin((N + 1) * theta_r / 2) / sin(theta_r / 2))^2 + 
       *          (r - beta) * (sin(N * theta_r / 2) / sin(theta_r / 2))^2 */
  }

  mpfr_set_ui_2exp(tmp2, 1, 
    (mpfr_exp_t)(2 * (parameters->l + parameters->m)), MPFR_RNDN);
      /* tmp = 2^(2 (l + m)) */
  mpfr_div(norm, tmp, tmp2, MPFR_RNDN);
    /* norm = (beta * (sin((N + 1) * theta_r / 2) / sin(theta_r / 2))^2 + 
     *            (r - beta) * (sin(N * theta_r / 2) / sin(theta_r / 2))^2) /
     *               2^(2 (l + m)) */

  /* Clean up. */
  mpz_clear(beta);

  mpfr_clear(N);

  mpfr_clear(tmp);
  mpfr_clear(tmp2);
  mpfr_clear(tmp3);
}