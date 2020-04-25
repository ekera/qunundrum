/*!
 * \file    distribution_slice.cpp
 * \ingroup two_dimensional_distribution_slice
 *
 * \brief   The definition of functions for manipulating slices in
 *          two-dimensional probability distributions.
 */

#include "distribution_slice.h"

#include "parameters.h"
#include "errors.h"
#include "math.h"
#include "common.h"

#include <mpfr.h>
#include <gmp.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

Distribution_Slice * distribution_slice_alloc()
{
  Distribution_Slice * slice =
    (Distribution_Slice *)malloc(sizeof(Distribution_Slice));
  if (NULL == slice) {
    critical("distribution_slice_alloc(): Failed to allocate memory.");
  }

  memset(slice, 0, sizeof(Distribution_Slice));

  return slice;
}

void distribution_slice_dealloc(
  Distribution_Slice ** const slice)
{
  free(*slice);
  *slice = NULL;
}

void distribution_slice_init(
  Distribution_Slice * const slice,
  const uint32_t dimension)
{
  memset(slice, 0, sizeof(Distribution_Slice));

  slice->dimension = dimension;

  slice->norm_matrix =
    (long double *)malloc((dimension * dimension) * sizeof(long double));
  if (NULL == (slice->norm_matrix)) {
    critical("distribution_slice_init(): Failed to allocate memory.");
  }

  memset(slice->norm_matrix, 0, (dimension * dimension) * sizeof(long double));
}

void distribution_slice_clear(
  Distribution_Slice * const slice)
{
  free(slice->norm_matrix);
  slice->norm_matrix = NULL;

  memset(slice, 0, sizeof(Distribution_Slice));
}

void distribution_slice_copy(
  Distribution_Slice * const dst_slice,
  const Distribution_Slice * const src_slice)
{
  const uint32_t dimension = src_slice->dimension;

  if (dimension != dst_slice->dimension) {
    distribution_slice_clear(dst_slice);
    distribution_slice_init(dst_slice, dimension);
  }

  dst_slice->min_log_alpha_d = src_slice->min_log_alpha_d;
  dst_slice->min_log_alpha_r = src_slice->min_log_alpha_r;

  dst_slice->total_probability = src_slice->total_probability;
  dst_slice->total_error = src_slice->total_error;

  dst_slice->flags = src_slice->flags;

  memcpy(
    dst_slice->norm_matrix,
    src_slice->norm_matrix,
    dimension * dimension * sizeof(long double));
}

void distribution_slice_init_copy(
  Distribution_Slice * const slice,
  const Distribution_Slice * const src)
{
  distribution_slice_init(slice, src->dimension);

  distribution_slice_copy(slice, src);
}

void distribution_slice_coordinates(
  const Distribution_Slice * const slice,
  double * const min_log_alpha_d,
  double * const max_log_alpha_d,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r)
{
  const double d_sgn = (double)sgn_i(slice->min_log_alpha_d);

  if (NULL != min_log_alpha_d) {
    (*min_log_alpha_d) = d_sgn * (abs_d((double)slice->min_log_alpha_d));
  }

  if (NULL != max_log_alpha_d) {
    (*max_log_alpha_d) = d_sgn * (abs_d((double)slice->min_log_alpha_d) + 1);
  }

  const double r_sgn = (double)sgn_i(slice->min_log_alpha_r);

  if (NULL != min_log_alpha_r) {
    (*min_log_alpha_r) = r_sgn * (abs_d((double)slice->min_log_alpha_r));
  }

  if (NULL != max_log_alpha_r) {
    (*max_log_alpha_r) = r_sgn * (abs_d((double)slice->min_log_alpha_r) + 1);
  }
}

void distribution_slice_region_coordinates(
  const Distribution_Slice * const slice,
  const uint32_t j,
  double * const min_log_alpha_d,
  double * const max_log_alpha_d,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r)
{
  const double step = ((double)1) / (double)(slice->dimension);

  const double d_sgn = (double)sgn_i(slice->min_log_alpha_d);
  const double d_abs = abs_d((double)slice->min_log_alpha_d);

  const uint32_t d_j = j % (slice->dimension);

  if (NULL != min_log_alpha_d) {
    (*min_log_alpha_d)  = d_sgn * (d_abs + (double)d_j * step);
  }

  if (NULL != max_log_alpha_d) {
    (*max_log_alpha_d)  = d_sgn * (d_abs + (double)(d_j + 1) * step);
  }

  const double r_sgn = (double)sgn_i(slice->min_log_alpha_r);
  const double r_abs = abs_d((double)slice->min_log_alpha_r);

  const uint32_t r_j = j / (slice->dimension);

  if (NULL != min_log_alpha_r) {
    (*min_log_alpha_r) = r_sgn * (r_abs + (double)r_j * step);
  }

  if (NULL != max_log_alpha_r) {
    (*max_log_alpha_r) = r_sgn * (r_abs + (double)(r_j + 1) * step);
  }
}

void distribution_slice_sample_region(
  const Distribution_Slice * const slice,
  Random_State * const random_state,
  double * const min_log_alpha_d,
  double * const max_log_alpha_d,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r)
{
  /* Select a pivot uniformly at random. */
  long double pivot = random_generate_pivot_inclusive(random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("distribution_slice_sample_region(): "
    "Debug: Sampled pivot: %Lf\n", pivot);
  #endif
  
  pivot *= slice->total_probability;

  #ifdef DEBUG_TRACE_SAMPLING
  printf("distribution_slice_sample_region(): "
    "Debug: Scaled pivot: %Lf\n", pivot);
  #endif

  /* Sample a region from the slice. */
  for (uint32_t i = 0; i < (slice->dimension * slice->dimension); i++) {
    pivot -= slice->norm_matrix[i];

    #ifdef DEBUG_TRACE_SAMPLING
    printf("distribution_slice_sample_region(): "
      "Debug: Decremented pivot to %Lf for region: %u\n", pivot, i);
    #endif

    if (pivot <= 0) {
      /* Select this region. */
      #ifdef DEBUG_TRACE_SAMPLING
      printf("distribution_slice_sample_region(): "
        "Debug: Selected region: %u\n", i);
      #endif

      /* Compute the region coordinates. */
      distribution_slice_region_coordinates(
        slice,
        i,
        min_log_alpha_d,
        max_log_alpha_d,
        min_log_alpha_r,
        max_log_alpha_r);

      #ifdef DEBUG_TRACE_SAMPLING
      printf("distribution_slice_sample_region(): "
        "Debug: The region coordinates are: (%f %f, %f %f)\n",
          *min_log_alpha_d, *max_log_alpha_d,
            *min_log_alpha_r, *max_log_alpha_r);
      #endif

      return;
    }
  }

  critical("distribution_slice_sample_region(): "
    "Failed to sample a region from the slice.");
}

void distribution_slice_copy_scale(
  Distribution_Slice * const dst_slice,
  const Distribution_Slice * const src_slice)
{
  if ((src_slice->dimension % dst_slice->dimension) != 0) {
    critical("distribution_slice_copy_scale(): Incompatible dimensions.");
  }

  const uint32_t scale = src_slice->dimension / dst_slice->dimension;

  dst_slice->flags  = src_slice->flags;
  dst_slice->flags |= SLICE_FLAGS_SCALED;

  dst_slice->min_log_alpha_d = src_slice->min_log_alpha_d;
  dst_slice->min_log_alpha_r = src_slice->min_log_alpha_r;

  dst_slice->total_probability = src_slice->total_probability;
  dst_slice->total_error = src_slice->total_error;

  for (uint32_t i1 = 0; i1 < dst_slice->dimension; i1++) {
    for (uint32_t j1 = 0; j1 < dst_slice->dimension; j1++) {
      dst_slice->norm_matrix[i1 + dst_slice->dimension * j1] = 0;

      for (uint32_t i2 = 0; i2 < scale; i2++) {
        for (uint32_t j2 = 0; j2 < scale; j2++) {
          const uint32_t i = i1 * scale + i2;
          const uint32_t j = j1 * scale + j2;

          dst_slice->norm_matrix[i1 + dst_slice->dimension * j1] +=
            src_slice->norm_matrix[i + src_slice->dimension * j];
        }
      }
    }
  }
}
