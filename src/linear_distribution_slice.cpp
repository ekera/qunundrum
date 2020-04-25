/*!
 * \file    linear_distribution_slice.cpp
 * \ingroup linear_distribution_slice
 *
 * \brief   The definition of functions for manipulating slices in linear
 *          probability distributions.
 */

#include "linear_distribution_slice.h"

#include "sample.h"
#include "parameters.h"
#include "errors.h"
#include "common.h"
#include "sample.h"
#include "math.h"

#include <mpfr.h>
#include <gmp.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

Linear_Distribution_Slice * linear_distribution_slice_alloc()
{
  Linear_Distribution_Slice * slice =
    (Linear_Distribution_Slice *)malloc(sizeof(Linear_Distribution_Slice));
  if (NULL == slice) {
    critical("linear_distribution_slice_alloc(): Failed to allocate memory.");
  }

  memset(slice, 0, sizeof(Linear_Distribution_Slice));

  return slice;
}

void linear_distribution_slice_dealloc(
  Linear_Distribution_Slice ** const slice)
{
  free(*slice);
  *slice = NULL;
}

void linear_distribution_slice_init(
  Linear_Distribution_Slice * const slice,
  const uint32_t dimension)
{
  memset(slice, 0, sizeof(Linear_Distribution_Slice));

  slice->dimension = dimension;

  slice->norm_vector =
    (long double *)malloc(dimension * sizeof(long double));
  if (NULL == (slice->norm_vector)) {
    critical("distribution_slice_init(): Failed to allocate memory.");
  }

  memset(slice->norm_vector, 0, dimension * sizeof(long double));
}

void linear_distribution_slice_init_copy(
  Linear_Distribution_Slice * const dst,
  const Linear_Distribution_Slice * const src)
{
  linear_distribution_slice_init(dst, src->dimension);

  linear_distribution_slice_copy(dst, src);
}

void linear_distribution_slice_copy(
  Linear_Distribution_Slice * const dst_slice,
  const Linear_Distribution_Slice * const src_slice)
{
  const uint32_t dimension = src_slice->dimension;

  if (dimension != dst_slice->dimension) {
    linear_distribution_slice_clear(dst_slice);
    linear_distribution_slice_init(dst_slice, dimension);
  }

  dst_slice->min_log_alpha = src_slice->min_log_alpha;

  dst_slice->total_probability = src_slice->total_probability;
  dst_slice->total_error = src_slice->total_error;

  dst_slice->flags = src_slice->flags;

  memcpy(
    dst_slice->norm_vector,
    src_slice->norm_vector,
    dimension * sizeof(long double));
}

void linear_distribution_slice_copy_scale(
  Linear_Distribution_Slice * const dst_slice,
  const Linear_Distribution_Slice * const src_slice)
{
  if ((src_slice->dimension % dst_slice->dimension) != 0) {
    critical("distribution_slice_scale(): Incompatible dimensions.");
  }

  const uint32_t scale = src_slice->dimension / dst_slice->dimension;

  dst_slice->flags  = src_slice->flags;
  dst_slice->flags |= SLICE_FLAGS_SCALED;

  dst_slice->min_log_alpha = src_slice->min_log_alpha;

  dst_slice->total_probability = src_slice->total_probability;
  dst_slice->total_error = src_slice->total_error;

  for (uint32_t i = 0; i < dst_slice->dimension; i++) {
    dst_slice->norm_vector[i] = 0;

    for (uint32_t j = 0; j < scale; j++) {
      dst_slice->norm_vector[i] += src_slice->norm_vector[i * scale + j];
    }
  }
}

void linear_distribution_slice_clear(
  Linear_Distribution_Slice * const slice)
{
  free(slice->norm_vector);
  slice->norm_vector = NULL;

  memset(slice, 0, sizeof(Linear_Distribution_Slice));
}

void linear_distribution_slice_coordinates(
  const Linear_Distribution_Slice * const slice,
  double * const min_log_alpha,
  double * const max_log_alpha)
{
  const double d_sgn = (double)sgn_i(slice->min_log_alpha);

  if (NULL != min_log_alpha) {
    (*min_log_alpha) = d_sgn * (abs_d((double)slice->min_log_alpha));
  }

  if (NULL != max_log_alpha) {
    (*max_log_alpha) = d_sgn * (abs_d((double)slice->min_log_alpha) + 1);
  }
}

void linear_distribution_slice_region_coordinates(
  const Linear_Distribution_Slice * const slice,
  const uint32_t j,
  double * const min_log_alpha,
  double * const max_log_alpha)
{
  const double step = (double)1 / (double)(slice->dimension);

  const double d_sgn = (double)sgn_i(slice->min_log_alpha);

  if (NULL != min_log_alpha) {
    (*min_log_alpha) = abs_d((double)slice->min_log_alpha) + (double)j * step;
    (*min_log_alpha) *= d_sgn;
  }

  if (NULL != max_log_alpha) {
    (*max_log_alpha) = 
      abs_d((double)slice->min_log_alpha) + (double)(j + 1) * step;
    (*max_log_alpha) *= d_sgn;
  }
}

void linear_distribution_slice_sample_region(
  const Linear_Distribution_Slice * const slice,
  Random_State * const random_state,
  double * const min_log_alpha,
  double * const max_log_alpha)
{
  /* Select a pivot uniformly at random. */
  long double pivot = random_generate_pivot_inclusive(random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("linear_distribution_slice_sample_region(): "
    "Debug: Sampled pivot: %Lf\n", pivot);
  #endif
  
  pivot *= slice->total_probability;

  #ifdef DEBUG_TRACE_SAMPLING
  printf("linear_distribution_slice_sample_region(): "
    "Debug: Scaled pivot: %Lf\n", pivot);
  #endif

  /* Sample a region from the slice. */
  for (uint32_t i = 0; i < slice->dimension; i++) {
    pivot -= slice->norm_vector[i];

    #ifdef DEBUG_TRACE_SAMPLING
    printf("linear_distribution_slice_sample_region(): "
      "Debug: Decremented pivot to %Lf for region: %u\n", pivot, i);
    #endif

    if (pivot <= 0) {
      /* Select this region. */
      #ifdef DEBUG_TRACE_SAMPLING
      printf("linear_distribution_slice_sample_region(): "
        "Debug: Selected region: %u\n", i);
      #endif

      /* Compute the region coordinates. */
      linear_distribution_slice_region_coordinates(
        slice,
        i,
        min_log_alpha,
        max_log_alpha);

      #ifdef DEBUG_TRACE_SAMPLING
      printf("linear_distribution_slice_sample_region(): "
        "Debug: The region coordinates are: %f %f\n",
          *min_log_alpha, *max_log_alpha);
      #endif

      return;
    }
  }

  critical("linear_distribution_slice_sample_region(): "
    "Failed to sample a region from the slice.");
}
