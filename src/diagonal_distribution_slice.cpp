/*!
 * \file    diagonal_distribution_slice.cpp
 * \ingroup diagonal_distribution_slice
 *
 * \brief   The definition of functions for manipulating slices in diagonal
 *          probability distributions.
 */

#include "diagonal_distribution_slice.h"

#include "diagonal_parameters.h"
#include "errors.h"
#include "common.h"
#include "sample.h"
#include "math.h"

#include <mpfr.h>
#include <gmp.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

Diagonal_Distribution_Slice * diagonal_distribution_slice_alloc()
{
  Diagonal_Distribution_Slice * slice =
    (Diagonal_Distribution_Slice *)malloc(sizeof(Diagonal_Distribution_Slice));
  if (NULL == slice) {
    critical("diagonal_distribution_slice_alloc(): Failed to allocate memory.");
  }

  memset(slice, 0, sizeof(Diagonal_Distribution_Slice));

  return slice;
}

void diagonal_distribution_slice_dealloc(
  Diagonal_Distribution_Slice ** const slice)
{
  free(*slice);
  *slice = NULL;
}

void diagonal_distribution_slice_init(
  Diagonal_Distribution_Slice * const slice,
  const uint32_t dimension)
{
  memset(slice, 0, sizeof(Diagonal_Distribution_Slice));

  slice->dimension = dimension;

  slice->norm_vector =
    (long double *)malloc(dimension * sizeof(long double));
  if (NULL == (slice->norm_vector)) {
    critical("distribution_slice_init(): Failed to allocate memory.");
  }

  memset(slice->norm_vector, 0, dimension * sizeof(long double));
}

void diagonal_distribution_slice_init_copy(
  Diagonal_Distribution_Slice * const dst,
  const Diagonal_Distribution_Slice * const src)
{
  diagonal_distribution_slice_init(dst, src->dimension);

  diagonal_distribution_slice_copy(dst, src);
}

void diagonal_distribution_slice_copy(
  Diagonal_Distribution_Slice * const dst_slice,
  const Diagonal_Distribution_Slice * const src_slice)
{
  const uint32_t dimension = src_slice->dimension;

  if (dimension != dst_slice->dimension) {
    diagonal_distribution_slice_clear(dst_slice);
    diagonal_distribution_slice_init(dst_slice, dimension);
  }

  dst_slice->min_log_alpha_r = src_slice->min_log_alpha_r;

  dst_slice->total_probability = src_slice->total_probability;
  dst_slice->total_error = src_slice->total_error;

  dst_slice->flags = src_slice->flags;

  memcpy(
    dst_slice->norm_vector,
    src_slice->norm_vector,
    dimension * sizeof(long double));
}

void diagonal_distribution_slice_copy_scale(
  Diagonal_Distribution_Slice * const dst_slice,
  const Diagonal_Distribution_Slice * const src_slice)
{
  if ((src_slice->dimension % dst_slice->dimension) != 0) {
    critical("diagonal_distribution_slice_copy_scale(): "
      "Incompatible dimensions.");
  }

  const uint32_t scale = src_slice->dimension / dst_slice->dimension;

  dst_slice->flags  = src_slice->flags;
  dst_slice->flags |= SLICE_FLAGS_SCALED;

  dst_slice->min_log_alpha_r = src_slice->min_log_alpha_r;

  dst_slice->total_probability = src_slice->total_probability;
  dst_slice->total_error = src_slice->total_error;

  for (uint32_t i = 0; i < dst_slice->dimension; i++) {
    dst_slice->norm_vector[i] = 0;

    for (uint32_t j = 0; j < scale; j++) {
      dst_slice->norm_vector[i] += src_slice->norm_vector[i * scale + j];
    }
  }
}

void diagonal_distribution_slice_clear(
  Diagonal_Distribution_Slice * const slice)
{
  free(slice->norm_vector);
  slice->norm_vector = NULL;

  memset(slice, 0, sizeof(Diagonal_Distribution_Slice));
}

void diagonal_distribution_slice_coordinates(
  const Diagonal_Distribution_Slice * const slice,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r)
{
  const double d_sgn = (double)sgn_i(slice->min_log_alpha_r);

  if (NULL != min_log_alpha_r) {
    (*min_log_alpha_r) = d_sgn * (abs_d((double)slice->min_log_alpha_r));
  }

  if (NULL != max_log_alpha_r) {
    (*max_log_alpha_r) = d_sgn * (abs_d((double)slice->min_log_alpha_r) + 1);
  }
}

void diagonal_distribution_slice_region_coordinates(
  const Diagonal_Distribution_Slice * const slice,
  const uint32_t j,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r)
{
  const double step = (double)1 / (double)(slice->dimension);

  const double d_sgn = (double)sgn_i(slice->min_log_alpha_r);

  if (NULL != min_log_alpha_r) {
    (*min_log_alpha_r) =
      abs_d((double)slice->min_log_alpha_r) + (double)j * step;
    (*min_log_alpha_r) *= d_sgn;
  }

  if (NULL != max_log_alpha_r) {
    (*max_log_alpha_r) =
      abs_d((double)slice->min_log_alpha_r) + (double)(j + 1) * step;
    (*max_log_alpha_r) *= d_sgn;
  }
}

void diagonal_distribution_slice_sample_region(
  const Diagonal_Distribution_Slice * const slice,
  Random_State * const random_state,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r)
{
  /* Select a pivot uniformly at random. */
  long double pivot = random_generate_pivot_inclusive(random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("diagonal_distribution_slice_sample_region(): "
    "Debug: Sampled pivot: %Lf\n", pivot);
  #endif

  pivot *= slice->total_probability;

  #ifdef DEBUG_TRACE_SAMPLING
  printf("diagonal_distribution_slice_sample_region(): "
    "Debug: Scaled pivot: %Lf\n", pivot);
  #endif

  /* Sample a region from the slice. */
  for (uint32_t i = 0; i < slice->dimension; i++) {
    pivot -= slice->norm_vector[i];

    #ifdef DEBUG_TRACE_SAMPLING
    printf("diagonal_distribution_slice_sample_region(): "
      "Debug: Decremented pivot to %Lf for region: %u\n", pivot, i);
    #endif

    if (pivot <= 0) {
      /* Select this region. */
      #ifdef DEBUG_TRACE_SAMPLING
      printf("diagonal_distribution_slice_sample_region(): "
        "Debug: Selected region: %u\n", i);
      #endif

      /* Compute the region coordinates. */
      diagonal_distribution_slice_region_coordinates(
        slice,
        i,
        min_log_alpha_r,
        max_log_alpha_r);

      #ifdef DEBUG_TRACE_SAMPLING
      printf("diagonal_distribution_slice_sample_region(): "
        "Debug: The region coordinates are: %f %f\n",
          *min_log_alpha_r, *max_log_alpha_r);
      #endif

      return;
    }
  }

  critical("diagonal_distribution_slice_sample_region(): "
    "Failed to sample a region from the slice.");
}
