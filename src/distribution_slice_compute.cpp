/*!
 * \file    distribution_slice_compute.cpp
 * \ingroup two_dimensional_distribution_slice
 *
 * \brief   The definition of functions for computing slices in two-dimensional
 *          probability distributions.
 */

#include "distribution_slice.h"

#include "probability.h"
#include "parameters.h"
#include "errors.h"
#include "math.h"
#include "common.h"

#include <mpfr.h>
#include <gmp.h>

#include <math.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Note that the implementation of this function is fairly explicit. Performance 
 * and memory is sacrificed for explicitness and ease of verification.
 * 
 * For instance, we first compute norms (and errors when applicable) and store
 * them in a matrix to explicitly iterate over sub-matrices when performing 
 * Simpson's method. It would be enough to keep the current two endpoints and 
 * the midpoint in temporary registers, but since we are operating in two 
 * dimensions such a change would make the code harder to follow.
 * 
 *If we really want to optimize this, we could incrementally compute the 
 * probability, by using that we can keep theta_r fixed whilst varying theta_d
 * and so forth. We could also decrease the precision further in some operations
 * and we could integrate this function better with the Richardson extrapolation
 * in distribution_slice_compute_richardson(), as there is currently an overlap
 * causing the norm to be re-computed for a large fraction of theta-pairs. */
void distribution_slice_compute(
  Distribution_Slice * const slice,
  const Parameters * const parameters,
  const Distribution_Slice_Compute_Method method,
  const int32_t min_log_alpha_d,
  const int32_t min_log_alpha_r)
{
  /* Basic constants. */
  const uint32_t dimension = slice->dimension;

  const double step = (double) 1 / (double)(dimension);

  /* Allocate memory for norms computed before interpolating. */
  mpfr_t ** norm_matrix =
    (mpfr_t **)malloc((2 * dimension + 1) * sizeof(mpfr_t *));
  if (NULL == norm_matrix) {
    critical("distribution_slice_compute(): Failed to allocate memory.");
  }

  for (uint32_t i = 0; i <= (2 * dimension); i++) {
    norm_matrix[i] = (mpfr_t *)malloc((2 * dimension + 1) * sizeof(mpfr_t));

    if (NULL == norm_matrix[i]) {
      critical("distribution_slice_compute(): Failed to allocate memory.");
    }
  }

  /* Allocate memory for errors computed before interpolating. */
  mpfr_t ** error_matrix = NULL;

  bool bounded_error = TRUE;

  if (DISTRIBUTION_SLICE_COMPUTE_METHOD_QUICK != method) {
    error_matrix = (mpfr_t **)malloc((2 * dimension + 1) * sizeof(mpfr_t *));

    if (NULL == error_matrix) {
      critical("distribution_slice_compute(): Failed to allocate memory.");
    }

    for (uint32_t i = 0; i <= (2 * dimension); i++) {
      error_matrix[i] = (mpfr_t *)malloc((2 * dimension + 1) * sizeof(mpfr_t));

      if (NULL == error_matrix[i]) {
        critical("distribution_slice_compute(): Failed to allocate memory.");
      }
    }
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

  mpfr_t pow_2m;
  mpfr_init2(pow_2m, PRECISION);
  mpfr_set_ui_2exp(pow_2m, 1, (mpfr_exp_t)(parameters->m), MPFR_RNDN);

  mpfr_t pow_2step;
  mpfr_init2(pow_2step, PRECISION);
  mpfr_set_d(pow_2step, step, MPFR_RNDN);
  mpfr_exp2(pow_2step, pow_2step, MPFR_RNDN);

  /* Temporary variables. */
  mpfr_t avg_norm;
  mpfr_init2(avg_norm, PRECISION);

  mpfr_t avg_error;
  mpfr_init2(avg_error, PRECISION);

  mpfr_t alpha_d;
  mpfr_init2(alpha_d, PRECISION);

  mpfr_t alpha_r;
  mpfr_init2(alpha_r, PRECISION);

  mpfr_t theta_d;
  mpfr_init2(theta_d, PRECISION);

  mpfr_t theta_r;
  mpfr_init2(theta_r, PRECISION);

  mpfr_t min_alpha_d;
  mpfr_init2(min_alpha_d, PRECISION);

  mpfr_t max_alpha_d;
  mpfr_init2(max_alpha_d, PRECISION);

  mpfr_t min_alpha_r;
  mpfr_init2(min_alpha_r, PRECISION);

  mpfr_t max_alpha_r;
  mpfr_init2(max_alpha_r, PRECISION);

  /* Setup the sigma parameter. */
  uint32_t sigma;
  
  if (DISTRIBUTION_SLICE_COMPUTE_METHOD_HEURISTIC_SIGMA == method) {
    const uint32_t tau = 11;

    /* Note that we could extract tau from the parameters above. However, we  
     * have also imposed a hard limit on t not growing past 10 when generating 
     * the distribution, and this limit is not accounted for when reading from 
     * the parameters. Note also that we set tau to the maximum possible value 
     * above. It may be possible to do even better by varying tau. */

    sigma = round(((float)parameters->l + tau + 4 - 1.6515f) / 2.0f);
  }

  /* Evaluate the norm (and error when applicable, depending on the method) in 
   * the (dimension + 1)^2 main the points given by
   * 
   *      (alpha_d, alpha_r) = (2^(min_log_alpha_d + i * step),
   *                            2^(min_log_alpha_r + j * step)),
   * 
   * where 0 <= i, j <= dimension and step = 1 / dimension. Also evaluate in 
   * the average points inbetween the main points, given by
   * 
   * ((2^(min_log_alpha_d + i * step) + 2^(min_log_alpha_d + (i+1) * step))/2)
   *  (2^(min_log_alpha_r + j * step) + 2^(min_log_alpha_r + (j+1) * step))/2))
   * 
   * and store the results in an interleaved fashion in the norm (and error) 
   * matrix. For X the main points, and A the averages inbetween, store as:
   * 
   *      X A X : X
   *      A A A : A
   *      X A X : X
   *      A A A : A
   *      : : : : : 
   *      X A X : X
   *  
   * Note: To interleave we let i, j run through 0 <= i, j <= 2 * dimension 
   * below, and use for the main points that
   * 
   *      (alpha_d, alpha_r) = (2^(min_log_alpha_d + (i / 2) * step),
   *                            2^(min_log_alpha_r + (j / 2) * step)). */

  /* Bootstrap the exponentiation in alpha_d. */
  mpfr_set_ui_2exp(max_alpha_d, 1, abs_i(min_log_alpha_d), MPFR_RNDN);

  for (uint32_t i = 0; i <= (2 * dimension); i++) {
    /* Compute the absolute value of alpha_r in a main or average point. */
    if (0 == (i % 2)) {
      /* Compute a main point. */
      mpfr_set(min_alpha_d, max_alpha_d, MPFR_RNDN);
      mpfr_round(alpha_d, min_alpha_d);
    } else {
      /* Compute an average point. */
      mpfr_mul(max_alpha_d, min_alpha_d, pow_2step, MPFR_RNDN);

      mpfr_add(alpha_d, min_alpha_d, max_alpha_d, MPFR_RNDN);
      mpfr_div_ui(alpha_d, alpha_d, 2, MPFR_RNDN);
      mpfr_round(alpha_d, alpha_d);
    }

    /* Adjust the sign. */
    if (sgn_d(min_log_alpha_d) == -1) {
      mpfr_neg(alpha_d, alpha_d, MPFR_RNDN);
    }

    /* Compute theta_d. */
    mpfr_mul(theta_d, alpha_d, scale_factor_theta, MPFR_RNDN);

    /* Bootstrap the exponentiation in alpha_r. */
    mpfr_set_ui_2exp(max_alpha_r, 1, abs_i(min_log_alpha_r), MPFR_RNDN);

    for (uint32_t j = 0; j <= (2 * dimension); j++) {
      /* Compute the absolute value of alpha_r in a main or average point. */
      if (0 == (j % 2)) {
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

      /* Compute theta_r. */
      mpfr_mul(theta_r, alpha_r, scale_factor_theta, MPFR_RNDN);

      /* Compute and store the probability (and error when applicable). */
      mpfr_init2(norm_matrix[i][j], PRECISION);
      
      switch (method) {
        case DISTRIBUTION_SLICE_COMPUTE_METHOD_HEURISTIC_SIGMA:
          mpfr_init2(error_matrix[i][j], PRECISION);

          bounded_error &= probability_approx(
                              norm_matrix[i][j],
                              error_matrix[i][j],
                              sigma, /* selected above and held constant. */
                              theta_d,
                              theta_r,
                              parameters);
          
          break;

        case DISTRIBUTION_SLICE_COMPUTE_METHOD_OPTIMAL_LOCAL_SIGMA:
          mpfr_init2(error_matrix[i][j], PRECISION);
          
          if ((0 == i) && (0 == j)) {
            /* For the first sample in the slice we find the optimal sigma. */
            bounded_error &= probability_approx_optimal_sigma(
                                norm_matrix[i][j],
                                error_matrix[i][j],
                                sigma, /* selected by the function */
                                theta_d,
                                theta_r,
                                parameters);
          } else {
            /* For all samples in the slice but the first, we adaptively adjust 
            * sigma find the local optimum, by increasing or decreasing it. */
            bounded_error &= probability_approx_adjust_sigma(
                                norm_matrix[i][j],
                                error_matrix[i][j],
                                sigma, /* adjusted by the function */
                                theta_d,
                                theta_r,
                                parameters);
          }

          break;

        case DISTRIBUTION_SLICE_COMPUTE_METHOD_QUICK:
          probability_approx_quick(
            norm_matrix[i][j],
            theta_d,
            theta_r,
            parameters);
          
          break;
        
        default:
          critical("distribution_slice_compute(): "
            "Unknown method specified for computing the slice.");
      }
    }
  }

  /* Apply Simpson's method to each 9 x 9 sub-matrix, with main points X at the 
   * edges, to sum up the probabilities in the norm matrix, and when applicable 
   * also the error bounds in the error matrix.
   * 
   * This is easy as we have stored the points in an interleaved fashion
   * 
   *      X A X : X
   *      A A A : A
   *      X A X : X
   *      A A A : A
   *      : : : : : 
   *      X A X : X.
   * 
   * For each 3 x 3 sub-matrix, 
   * 
   *      X A X                           1  4 1
   *      A A A      we apply weights     4 16 4
   *      A A A                           1  4 1,
   * 
   * sum up the results, and divide by 4 * 1 + 4 * 4 + 16 = 36. */

  slice->total_probability = 0;
  slice->total_error = 0;
  
  /* Bootstrap the exponentiation in alpha_d. */
  mpfr_set_ui_2exp(max_alpha_d, 1, abs_i(min_log_alpha_d), MPFR_RNDN);

  for (uint32_t i = 0; i < (2 * dimension); i += 2) {
    mpfr_set(min_alpha_d, max_alpha_d, MPFR_RNDN);
    mpfr_mul(max_alpha_d, min_alpha_d, pow_2step, MPFR_RNDN);

    mpfr_sub(alpha_d, max_alpha_d, min_alpha_d, MPFR_RNDN);

    /* Bootstrap the exponentiation in alpha_r. */
    mpfr_set_ui_2exp(max_alpha_r, 1, abs_i(min_log_alpha_r), MPFR_RNDN);

    for (uint32_t j = 0; j < (2 * dimension); j += 2) {
      mpfr_set(min_alpha_r, max_alpha_r, MPFR_RNDN);
      mpfr_mul(max_alpha_r, min_alpha_r, pow_2step, MPFR_RNDN);

      mpfr_sub(alpha_r, max_alpha_r, min_alpha_r, MPFR_RNDN);

      /* Compute the norm in each 3 x 3 sub-matrix using Simpson's rule. */
      mpfr_set_ui(avg_norm, 4, MPFR_RNDN);
      mpfr_mul(avg_norm, avg_norm, norm_matrix[i + 1][j + 1], MPFR_RNDN);

      mpfr_add(avg_norm, avg_norm, norm_matrix[i    ][j + 1], MPFR_RNDN);
      mpfr_add(avg_norm, avg_norm, norm_matrix[i + 2][j + 1], MPFR_RNDN);
      mpfr_add(avg_norm, avg_norm, norm_matrix[i + 1][j    ], MPFR_RNDN);
      mpfr_add(avg_norm, avg_norm, norm_matrix[i + 1][j + 2], MPFR_RNDN);
      mpfr_mul_ui(avg_norm, avg_norm, 4, MPFR_RNDN);

      mpfr_add(avg_norm, avg_norm, norm_matrix[i    ][j    ], MPFR_RNDN);
      mpfr_add(avg_norm, avg_norm, norm_matrix[i + 2][j    ], MPFR_RNDN);
      mpfr_add(avg_norm, avg_norm, norm_matrix[i    ][j + 2], MPFR_RNDN);
      mpfr_add(avg_norm, avg_norm, norm_matrix[i + 2][j + 2], MPFR_RNDN);

      mpfr_div_ui(avg_norm, avg_norm, 36, MPFR_RNDN);

      if (NULL != error_matrix) {
        /* Compute the error in each 3 x 3 sub-matrix using Simpson's rule. */
        mpfr_set_ui(avg_error, 4, MPFR_RNDN);
        mpfr_mul(avg_error, avg_error, error_matrix[i + 1][j + 1], MPFR_RNDN);

        mpfr_add(avg_error, avg_error, error_matrix[i    ][j + 1], MPFR_RNDN);
        mpfr_add(avg_error, avg_error, error_matrix[i + 2][j + 1], MPFR_RNDN);
        mpfr_add(avg_error, avg_error, error_matrix[i + 1][j    ], MPFR_RNDN);
        mpfr_add(avg_error, avg_error, error_matrix[i + 1][j + 2], MPFR_RNDN);
        mpfr_mul_ui(avg_error, avg_error, 4, MPFR_RNDN);

        mpfr_add(avg_error, avg_error, error_matrix[i    ][j    ], MPFR_RNDN);
        mpfr_add(avg_error, avg_error, error_matrix[i + 2][j    ], MPFR_RNDN);
        mpfr_add(avg_error, avg_error, error_matrix[i    ][j + 2], MPFR_RNDN);
        mpfr_add(avg_error, avg_error, error_matrix[i + 2][j + 2], MPFR_RNDN);

        mpfr_div_ui(avg_error, avg_error, 36, MPFR_RNDN);
      }

      /* Account for the number of alpha-values in this square. */      
      mpfr_mul(avg_norm, avg_norm, alpha_d, MPFR_RNDN);
      mpfr_mul(avg_norm, avg_norm, alpha_r, MPFR_RNDN);
      
      /* Account for the density. */
      mpfr_div(avg_norm, avg_norm, pow_2m, MPFR_RNDN);

      if (NULL != error_matrix) {
        mpfr_mul(avg_error, avg_error, alpha_d, MPFR_RNDN);
        mpfr_mul(avg_error, avg_error, alpha_r, MPFR_RNDN);

        /* Account for the density. */
        mpfr_div(avg_error, avg_error, pow_2m, MPFR_RNDN);
      }

      const long double norm = mpfr_get_ld(avg_norm, MPFR_RNDN);
      slice->norm_matrix[(i / 2) + dimension * (j / 2)] = norm;
      slice->total_probability += norm;

      if (NULL != error_matrix) {
        slice->total_error += mpfr_get_ld(avg_error, MPFR_RNDN);
      }
    }
  }

  slice->min_log_alpha_d = min_log_alpha_d;
  slice->min_log_alpha_r = min_log_alpha_r;

  slice->flags &= ~(SLICE_FLAGS_MASK_METHOD);
  slice->flags |= SLICE_FLAGS_METHOD_SIMPSON;

  if (!bounded_error) {
    slice->flags |= SLICE_FLAGS_ERROR_BOUND_WARNING;
  }

  /* Clean up memory. */
  for (uint32_t i = 0; i <= (2 * dimension); i++) {
    for (uint32_t j = 0; j <= (2 * dimension); j++) {
      mpfr_clear(norm_matrix[i][j]);
    }
    free(norm_matrix[i]);
  }
  free(norm_matrix);

  if (NULL != error_matrix) {
    for (uint32_t i = 0; i <= (2 * dimension); i++) {
      for (uint32_t j = 0; j <= (2 * dimension); j++) {
        mpfr_clear(error_matrix[i][j]);
      }
      free(error_matrix[i]);
    }
    free(error_matrix);
  }

  mpfr_clear(scale_factor_theta);
  mpfr_clear(pow_2m);
  mpfr_clear(pow_2step);

  mpfr_clear(avg_norm);
  mpfr_clear(avg_error);
  mpfr_clear(alpha_d);
  mpfr_clear(alpha_r);
  mpfr_clear(theta_d);
  mpfr_clear(theta_r);
  mpfr_clear(min_alpha_d);
  mpfr_clear(max_alpha_d);
  mpfr_clear(min_alpha_r);
  mpfr_clear(max_alpha_r);
}
