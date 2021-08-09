/*!
 * \file    diagonal_distribution_slice_compute.cpp
 * \ingroup diagonal_distribution_slice
 *
 * \brief   The definition of functions for computing slices in diagonal
 *          probability distributions.
 */

#include "diagonal_distribution_slice.h"

#include "diagonal_probability.h"
#include "parameters.h"
#include "errors.h"
#include "math.h"
#include "common.h"

#include <mpfr.h>
#include <gmp.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Note that the implementation of this function is fairly explicit. Performance 
 * and memory is sacrificed for explicitness and ease of verification.
 * 
 * For instance, we first compute norms and store them in a vector to explicitly
 * iterate over sub-vectors when performing Simpson's method. It would be enough
 * to keep the current two endpoints and the midpoint in temporary registers, 
 * but we have kept this design for now to keep everything analogous with the 
 * two-dimensional implementation in distribution_slice_compute(). */
void diagonal_distribution_slice_compute(
  Diagonal_Distribution_Slice * const slice,
  const Parameters * const parameters,
  const int32_t min_log_alpha_r,
  const int32_t offset_alpha_d)
{
  const uint32_t dimension = slice->dimension;

  const double step = (double) 1 / (double)(dimension);

  /* Allocate memory for norms computed before interpolating. */
  mpfr_t * norm_vector =
    (mpfr_t *)malloc((2 * dimension + 1) * sizeof(mpfr_t));
  if (NULL == norm_vector) {
    critical("diagonal_distribution_slice_compute(): "
      "Failed to allocate memory.");
  }

  /* Constants. */
  const uint32_t precision =
    2 * (parameters->l + max_ui(parameters->m, PRECISION));

  mpfr_t scale_factor_theta;
  mpfr_init2(scale_factor_theta, precision);
  
  {
    /* Setup the scale factor when mapping alpha to theta. */
    mpfr_t tmp;
    mpfr_init2(tmp, precision);

    mpfr_const_pi(scale_factor_theta, MPFR_RNDN);
    mpfr_mul_ui(scale_factor_theta, scale_factor_theta, 2, MPFR_RNDN);

    mpfr_set_ui_2exp(tmp, 1, 
      (mpfr_exp_t)(parameters->l + parameters->m), MPFR_RNDN);
    mpfr_div(scale_factor_theta, scale_factor_theta, tmp, MPFR_RNDN);

    /* Clear memory. */
    mpfr_clear(tmp);
  }

  mpfr_t pow_2step;
  mpfr_init2(pow_2step, precision);
  mpfr_set_d(pow_2step, step, MPFR_RNDN);
  mpfr_exp2(pow_2step, pow_2step, MPFR_RNDN);

  /* Temporary variables. */
  mpfr_t avg_norm;
  mpfr_init2(avg_norm, PRECISION);

  mpfr_t alpha_d;
  mpfr_init2(alpha_d, precision);

  mpfr_t alpha_r;
  mpfr_init2(alpha_r, precision);

  mpfr_t theta_d;
  mpfr_init2(theta_d, precision);

  mpfr_t theta_r;
  mpfr_init2(theta_r, precision);

  mpfr_t min_alpha_r;
  mpfr_init2(min_alpha_r, precision);

  mpfr_t max_alpha_r;
  mpfr_init2(max_alpha_r, precision);

  /* Evaluate the norm in the dimension + 1 main points given by
   * 
   *   alpha_r = sgn(min_log_alpha_r) * 2^(abs(min_log_alpha_r) + i * step),
   * 
   * where 0 <= i <= dimension and step = 1 / dimension. Also evaluate the norm
   * in the average points inbetween the main points, given by
   * 
   *   alpha_r = sgn(min_log_alpha_r) * 
   *               (2^(abs(min_log_alpha_r) + i * step) + \
   *                2^(abs(min_log_alpha_r) + (i + 1) * step)) / 2
   * 
   * and store the results in an interleaved fashion in the norm vector. For X 
   * the main points, and A the averages inbetween, store as:
   * 
   *   X A X A X : X.
   * 
   * Note: To interleave we let i run through 0 <= i <= 2 * dimension below,
   * and use for the main points that
   * 
   *  alpha_r = sgn(min_log_alpha_r) 2^(abs(min_log_alpha_r) + (i / 2) * step).
   */

  /* Bootstrap the exponentiation in alpha_r. */
  mpfr_set_ui_2exp(max_alpha_r, 1, abs_i(min_log_alpha_r), MPFR_RNDN);

  for (uint32_t i = 0; i <= (2 * dimension); i++) {
    /* Compute the absolute value of alpha_r in a main or average point. */
    if (0 == (i % 2)) {
      /* Compute a main point. */
      mpfr_set(min_alpha_r, max_alpha_r, MPFR_RNDN);
      mpfr_round(alpha_r, min_alpha_r);
    } else {
      /* Compute an average point. */
      mpfr_mul(max_alpha_r, min_alpha_r, pow_2step, MPFR_RNDN);

      mpfr_add(alpha_r, min_alpha_r, max_alpha_r, MPFR_RNDN);
      mpfr_div_ui(alpha_r, alpha_r, 2, MPFR_RNDN);
      mpfr_round(alpha_r, alpha_r);
    }

    /* Adjust the sign. */
    if (sgn_d(min_log_alpha_r) == -1) {
      mpfr_neg(alpha_r, alpha_r, MPFR_RNDN);
    }

    /* Compute alpha_d. */
    mpfr_set(alpha_d, alpha_r, MPFR_RNDN);
    mpfr_mul_z(alpha_d, alpha_d, parameters->d, MPFR_RNDN);
    mpfr_div_z(alpha_d, alpha_d, parameters->r, MPFR_RNDN);
    mpfr_round(alpha_d, alpha_d);
    mpfr_add_si(alpha_d, alpha_d, offset_alpha_d, MPFR_RNDN);

    /* Compute theta_d and theta_r. */
    mpfr_mul(theta_d, scale_factor_theta, alpha_d, MPFR_RNDN);
    mpfr_mul(theta_r, scale_factor_theta, alpha_r, MPFR_RNDN);

    /* Compute and store the probability. */
    mpfr_init2(norm_vector[i], PRECISION);
    diagonal_probability_approx(norm_vector[i], theta_d, theta_r, parameters);
  }

  /* We can now drop the precision. The probability function is sensitive to 
   * the precision in alpha_r and alpha_d, not the summation that follows. */
  mpfr_set_prec(alpha_r, PRECISION);
  mpfr_set_prec(min_alpha_r, PRECISION);
  mpfr_set_prec(max_alpha_r, PRECISION);
  
  /* Apply Simpson's method to each 3 element sub-vector, with main points X at 
   * the edges, to sum up the probabilities in the norm vector.
   * 
   * This is easy as we have stored the points in an interleaved fashion
   * 
   *   X A X A X : X.
   * 
   * For each sub-vector X A X, we apply weights 1 4 1, sum up the results, and 
   * divide by 2 * 1 + 4 = 6. */

  slice->total_probability = 0;
  slice->total_error = 0; /* Not used but explicitly set to zero. */

  /* Bootstrap the exponentiation in alpha_r. */
  mpfr_set_ui_2exp(max_alpha_r, 1, abs_i(min_log_alpha_r), MPFR_RNDN);

  for (uint32_t i = 0; i < (2 * dimension); i += 2) {
    mpfr_set(min_alpha_r, max_alpha_r, MPFR_RNDN);
    mpfr_mul(max_alpha_r, min_alpha_r, pow_2step, MPFR_RNDN);
    mpfr_sub(alpha_r, max_alpha_r, min_alpha_r, MPFR_RNDN);

    /* Compute the norm using Simpson's rule. */
    mpfr_set_ui(avg_norm, 4, MPFR_RNDN);
    mpfr_mul(avg_norm, avg_norm, norm_vector[i + 1], MPFR_RNDN);
    mpfr_add(avg_norm, avg_norm, norm_vector[i    ], MPFR_RNDN);
    mpfr_add(avg_norm, avg_norm, norm_vector[i + 2], MPFR_RNDN);
    mpfr_div_ui(avg_norm, avg_norm, 6, MPFR_RNDN);

    /* Account for the number of alpha_r-values in this region. */
    mpfr_mul(avg_norm, avg_norm, alpha_r, MPFR_RNDN);
    
    const long double norm = mpfr_get_ld(avg_norm, MPFR_RNDN);
    slice->norm_vector[i / 2] = norm;
    slice->total_probability += norm;
  }

  slice->min_log_alpha_r = min_log_alpha_r;
  slice->offset_alpha_d = offset_alpha_d;

  /* Set the method flag. */
  slice->flags &= ~(SLICE_FLAGS_MASK_METHOD);
  slice->flags |= SLICE_FLAGS_METHOD_SIMPSON;

  /* Clear memory. */
  for (uint32_t i = 0; i <= (2 * dimension); i++) {
    mpfr_clear(norm_vector[i]);
  }
  free(norm_vector);

  mpfr_clear(scale_factor_theta);
  mpfr_clear(pow_2step);

  mpfr_clear(avg_norm);
  mpfr_clear(alpha_d);
  mpfr_clear(alpha_r);
  mpfr_clear(theta_d);
  mpfr_clear(theta_r);
  mpfr_clear(min_alpha_r);
  mpfr_clear(max_alpha_r);
}
