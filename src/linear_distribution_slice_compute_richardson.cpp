/*!
 * \file    linear_distribution_slice_compute_richardson.cpp
 * \ingroup linear_distribution_slice
 *
 * \brief   The definition of functions for computing slices in linear
 *          probability distributions with Richardson extrapolation.
 */

#include "linear_distribution_slice.h"

#include "common.h"
#include "math.h"
#include "parameters.h"

#include <stdint.h>

void linear_distribution_slice_compute_richardson(
  Linear_Distribution_Slice * const slice,
  const Parameters * const parameters,
  const Linear_Distribution_Slice_Compute_Target target,
  const int32_t min_log_alpha)
{
  /* Extract the dimension. */
  const uint32_t dimension = slice->dimension;

  /* Allocate a slice with twice the resolution. */
  Linear_Distribution_Slice double_slice;
  linear_distribution_slice_init(&double_slice, 2 * dimension);

  /* Compute the slice and the slice with twice the resolution. */
  linear_distribution_slice_compute(
    slice,
    parameters,
    target,
    min_log_alpha);

  linear_distribution_slice_compute(
    &double_slice,
    parameters,
    target,
    min_log_alpha);

  /* Perform Richardson extrapolation. */
  slice->total_probability = 0;

  for (uint32_t i = 0; i < dimension; i++) {
    const long double probability = slice->norm_vector[i];

    const long double double_probability =
      double_slice.norm_vector[2 * i] +
      double_slice.norm_vector[2 * i + 1];

    slice->norm_vector[i] = 2 * double_probability - probability;

    slice->total_probability += slice->norm_vector[i];
  }

  slice->total_error = 2 * double_slice.total_error - slice->total_error;

  /* Update the method flag. */
  slice->flags |= SLICE_FLAGS_METHOD_RICHARDSON;

  /* Clear memory. */
  linear_distribution_slice_clear(&double_slice);
}
