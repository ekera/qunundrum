/*!
 * \file    diagonal_distribution_slice_compute_richardson.cpp
 * \ingroup diagonal_distribution_slice
 *
 * \brief   The definition of functions for computing slices in diagonal
 *          probability distributions with Richardson extrapolation.
 */

#include "diagonal_distribution_slice.h"

#include "common.h"
#include "diagonal_parameters.h"
#include "math.h"

#include <stdint.h>

void diagonal_distribution_slice_compute_richardson(
  Diagonal_Distribution_Slice * const slice,
  const Diagonal_Parameters * const parameters,
  const int32_t min_log_alpha_r)
{
  /* Extract the dimension. */
  const uint32_t dimension = slice->dimension;

  /* Allocate a slice with twice the resolution. */
  Diagonal_Distribution_Slice double_slice;
  diagonal_distribution_slice_init(&double_slice, 2 * dimension);

  /* Compute the slice and the slice with twice the resolution. */
  diagonal_distribution_slice_compute(
    slice,
    parameters,
    min_log_alpha_r);

  diagonal_distribution_slice_compute(
    &double_slice,
    parameters,
    min_log_alpha_r);

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
  diagonal_distribution_slice_clear(&double_slice);
}
