/*!
 * \file    distribution_slice_compute_richardson.cpp
 * \ingroup two_dimensional_distribution_slice
 *
 * \brief   The definition of functions for computing slices in two-dimensional
 *          probability distributions with Richardson extrapolation.
 */

#include "distribution_slice.h"

#include "common.h"
#include "math.h"
#include "parameters.h"

#include <stdint.h>

void distribution_slice_compute_richardson(
  Distribution_Slice * const slice,
  const Parameters * const parameters,
  const Distribution_Slice_Compute_Method method,
  const int32_t min_log_alpha_d,
  const int32_t min_log_alpha_r)
{
  /* Extract the dimension. */
  const uint32_t dimension = slice->dimension;

  /* Allocate a slice with twice the resolution. */
  Distribution_Slice double_slice;
  distribution_slice_init(&double_slice, 2 * dimension);

  /* Compute the slice and the slice with twice the resolution. */
  distribution_slice_compute(
    slice,
    parameters,
    method,
    min_log_alpha_d,
    min_log_alpha_r);

  distribution_slice_compute(
    &double_slice,
    parameters,
    method,
    min_log_alpha_d,
    min_log_alpha_r);

  /* Perform Richardson extrapolation. */
  slice->total_probability = 0;

  for (uint32_t i = 0; i < dimension; i++) {
    for (uint32_t j = 0; j < dimension; j++) {
      const long double probability = slice->norm_matrix[dimension * j + i];

      const long double double_probability =
        double_slice.norm_matrix[(2 * dimension) * (2 * j) + (2 * i)] +
        double_slice.norm_matrix[(2 * dimension) * (2 * j) + (2 * i + 1)] +
        double_slice.norm_matrix[(2 * dimension) * (2 * j + 1) + (2 * i)] +
        double_slice.norm_matrix[(2 * dimension) * (2 * j + 1) + (2 * i + 1)];

      slice->norm_matrix[dimension * j + i] =
        2 * double_probability - probability;

      slice->total_probability += slice->norm_matrix[dimension * j + i];
    }
  }

  slice->total_error = 2 * double_slice.total_error - slice->total_error;

  /* Update the method flag. */
  slice->flags |= SLICE_FLAGS_METHOD_RICHARDSON;

  /* Clear memory. */
  distribution_slice_clear(&double_slice);
}
