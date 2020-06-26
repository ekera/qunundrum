/*!
 * \file    linear_distribution_slice_compute.cpp
 * \ingroup linear_distribution_slice
 *
 * \brief   The definition of functions for computing slices in linear
 *          probability distributions.
 */

#include "linear_distribution_slice.h"

#include "linear_probability.h"
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
void linear_distribution_slice_compute(
  Linear_Distribution_Slice * const slice,
  const Parameters * const parameters,
  const Linear_Distribution_Slice_Compute_Target target,
  const int32_t min_log_alpha)
{
  /* Basic constants. */
  const uint32_t dimension = slice->dimension;

  const double step = (double) 1 / (double)(dimension);

  /* Allocate memory for norms computed before interpolating. */
  mpfr_t * norm_vector =
    (mpfr_t *)malloc((2 * dimension + 1) * sizeof(mpfr_t));
  if (NULL == norm_vector) {
    critical("linear_distribution_slice_compute_d(): "
      "Failed to allocate memory.");
  }

  /* Constants. */
  mpfr_t scale_factor_theta;
  mpfr_init2(scale_factor_theta, PRECISION);
  
  {
    /* Setup the scale factor when mapping alpha to theta. */
    mpfr_t tmp;
    mpfr_init2(tmp, PRECISION);

    mpfr_const_pi(scale_factor_theta, MPFR_RNDN);
    mpfr_mul_ui(scale_factor_theta, scale_factor_theta, 2, MPFR_RNDN);

    mpfr_set_ui_2exp(tmp, 1, 
      (mpfr_exp_t)(parameters->l + parameters->m), MPFR_RNDN);
    mpfr_div(scale_factor_theta, scale_factor_theta, tmp, MPFR_RNDN);

    /* Clean up temporary memory. */
    mpfr_clear(tmp);
  }

  mpfr_t pow_2l;
  mpfr_init2(pow_2l, PRECISION);
  mpfr_set_ui_2exp(pow_2l, 1, (mpfr_exp_t)(parameters->l), MPFR_RNDN);

  mpfr_t pow_2step;
  mpfr_init2(pow_2step, PRECISION);
  mpfr_set_d(pow_2step, step, MPFR_RNDN);
  mpfr_exp2(pow_2step, pow_2step, MPFR_RNDN);

  /* Temporary variables. */
  mpfr_t avg_norm;
  mpfr_init2(avg_norm, PRECISION);

  mpfr_t alpha;
  mpfr_init2(alpha, PRECISION);

  mpfr_t theta;
  mpfr_init2(theta, PRECISION);

  mpfr_t min_alpha;
  mpfr_init2(min_alpha, PRECISION);

  mpfr_t max_alpha;
  mpfr_init2(max_alpha, PRECISION);

  /* Evaluate the norm in the dimension + 1 main the points given by
   * 
   *      alpha = 2^(min_log_alpha + i * step,
   * 
   * where 0 <= i <= dimension and step = 1 / dimension. Also evaluate in the 
   * average points inbetween the main points, given by
   * 
   * (2^(min_log_alpha + i * step) + 2^(min_log_alpha + (i + 1) * step)) / 2
   * 
   * and store the results in an interleaved fashion in the norm vector. For X 
   * the main points, and A the averages inbetween, store as:
   * 
   *      X A X A X : X
   * 
   * Note: To interleave we let i run through 0 <= i <= 2 * dimension below,
   * and use for the main points that
   * 
   *      alpha = 2^(min_log_alpha_d + (i / 2) * step). */

  /* Bootstrap the exponentiation in alpha. */
  mpfr_set_ui_2exp(max_alpha, 1, abs_i(min_log_alpha), MPFR_RNDN);

  for (uint32_t i = 0; i <= (2 * dimension); i++) {
    /* Compute the absolute value of alpha in a main or average point. */
    if (0 == (i % 2)) {
      /* Compute a main point. */
      mpfr_set(min_alpha, max_alpha, MPFR_RNDN);
      mpfr_round(alpha, min_alpha);
    } else {
      /* Compute an average point. */
      mpfr_mul(max_alpha, min_alpha, pow_2step, MPFR_RNDN);

      mpfr_add(alpha, min_alpha, max_alpha, MPFR_RNDN);
      mpfr_div_ui(alpha, alpha, 2, MPFR_RNDN);
      mpfr_round(alpha, alpha);
    }

    /* Adjust the sign. */
    if (sgn_d(min_log_alpha) == -1) {
      mpfr_neg(alpha, alpha, MPFR_RNDN);
    }

    /* Compute theta. */
    mpfr_mul(theta, alpha, scale_factor_theta, MPFR_RNDN);
  
    /* Compute and store the probability. */
    mpfr_init2(norm_vector[i], PRECISION);

    switch (target) {
      case LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_D:
        linear_probability_d(norm_vector[i], theta, parameters);
        break;

      case LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_R:
        linear_probability_r(norm_vector[i], theta, parameters);
        break;
      
      default:
        critical("linear_distribution_slice_compute(): Unknown target.");
    }
  }

  /* Apply Simpson's method to each 3 element sub-vector, with main points X at  
   * the edges, to sum up the probabilities in the norm vector.
   * 
   * This is easy as we have stored the points in an interleaved fashion
   * 
   *      X A X : X
   * 
   * For each sub-vector X A X, we apply weights 1 4 1, sum up the results, and 
   * divide by 2 * 1 + 4 = 6. */

  slice->total_probability = 0;
  slice->total_error = 0; /* Not used but explicitly set to zero. */

  /* Bootstrap the exponentiation in alpha. */
  mpfr_set_ui_2exp(max_alpha, 1, abs_i(min_log_alpha), MPFR_RNDN);

  for (uint32_t i = 0; i < (2 * dimension); i += 2) {
    mpfr_set(min_alpha, max_alpha, MPFR_RNDN);
    mpfr_mul(max_alpha, min_alpha, pow_2step, MPFR_RNDN);
    mpfr_sub(alpha, max_alpha, min_alpha, MPFR_RNDN);
    
    /* Compute the norm using Simpson's rule. */
    mpfr_set_ui(avg_norm, 4, MPFR_RNDN);
    mpfr_mul(avg_norm, avg_norm, norm_vector[i + 1], MPFR_RNDN);
    mpfr_add(avg_norm, avg_norm, norm_vector[i    ], MPFR_RNDN);
    mpfr_add(avg_norm, avg_norm, norm_vector[i + 2], MPFR_RNDN);
    mpfr_div_ui(avg_norm, avg_norm, 6, MPFR_RNDN);

    /* Account for the number of alpha-values in this region. */
    mpfr_mul(avg_norm, avg_norm, alpha, MPFR_RNDN);
    
    if (LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_D == target) {
      /* There are 2^(l + kappa) pairs (j, k) that yield each admissible 
       * argument alpha. Only the arguments alpha that are multiples of 2^kappa  
       * are admissible. Hence we have a multiplicity of 2^(l+kappa), and a 
       * density of 2^kappa, so we have to multiply by 2^l to compensate below.
       * 
       * Recall that in the above context, kappa = kappa_d is the largest 
       * integer such that 2^kappa_d divides the logarithm d. */
      mpfr_mul(avg_norm, avg_norm, pow_2l, MPFR_RNDN);
    }
    
    /* Note that when performing order-finding, there are 2^(m+l-kappa) integers 
     * j that yield each admissible argument alpha. Only the arguments alpha 
     * that are multiples of 2^kappa are admissible. Hence we have a 
     * multiplicity of 2^kappa, and a density of 2^kappa, and these two factors
     * cancel. Hence, we need not multiply or divide to compensate.
     * 
     * Recall that in the order-finding context, kappa = kappa_r is the largest 
     * integer such that 2^kappa_r divides the order r. */

    /* Note that we do not account for the case where kappa is artifically 
     * large when constructing the histograms. If kappa is close to m in size, 
     * some bins will become unavailable, requiring us to compute the histogram 
     * in a slightly different way. This is expalined in the papers. 
     * 
     * Artifically large kappa is a special case. It does not occur in practice 
     * in cryptographic applications, as d be presumed to be random in such 
     * applications, and as the order r may in general be assumed to be prime.
     * Also d may be randomized if necessary. */

    const long double norm = mpfr_get_ld(avg_norm, MPFR_RNDN);
    slice->norm_vector[i / 2] = norm;
    slice->total_probability += norm;
  }

  slice->min_log_alpha = min_log_alpha;

  /* Set the method flag. */
  slice->flags &= ~(SLICE_FLAGS_MASK_METHOD);
  slice->flags |= SLICE_FLAGS_METHOD_SIMPSON;

  /* Clean up memory. */
  for (uint32_t i = 0; i <= (2 * dimension); i++) {
    mpfr_clear(norm_vector[i]);
  }
  free(norm_vector);

  mpfr_clear(scale_factor_theta);
  mpfr_clear(pow_2l);
  mpfr_clear(pow_2step);

  mpfr_clear(avg_norm);
  mpfr_clear(alpha);
  mpfr_clear(theta);
  mpfr_clear(min_alpha);
  mpfr_clear(max_alpha);
}