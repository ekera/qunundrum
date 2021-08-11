/*!
 * \file    distribution_info.cpp
 * \ingroup two_dimensional_distribution
 *
 * \brief   The definitions of functions for exporting information about
 *          two-dimensional probability distributions.
 */

#include "distribution.h"

#include "common.h"
#include "errors.h"
#include "math.h"

#include <gmp.h>

#include <stdint.h>
#include <stdio.h>
#include <string.h>

/*!
* \brief   The maximum supported slice dimension when exporting information
*          on two-dimensional probability distributions.
 */
#define MAX_SLICE_DIMENSION                     16384

void distribution_export_info(
  const Distribution * const distribution,
  FILE * const file,
  const bool export_slices)
{
  fprintf(file, "Precision: %u bits\n\n", distribution->precision);

  const uint32_t m = distribution->parameters.m;
  const uint32_t s = distribution->parameters.s;
  const uint32_t l = distribution->parameters.l;

  fprintf(file, "Parameters:\n");
  fprintf(file, " m: %u\n", m);
  fprintf(file, " s: %u\n", s);
  fprintf(file, " l: %u\n", l);

  gmp_fprintf(file, " r: %Zd\n", distribution->parameters.r);
  gmp_fprintf(file, " d: %Zd\n\n", distribution->parameters.d);

  /* Print information on the slices. */
  long double pivot = 0.1f;
  while (distribution->slices[0]->total_probability < pivot) {
    pivot /= (long double)10;
  }

  uint32_t begin_i = 0;
  uint32_t bound_error_count = 0;

  long double region_probability = 0;

  fprintf(file, "Slices:\n");

  uint32_t dimensions[MAX_SLICE_DIMENSION + 1];
  memset(dimensions, 0, sizeof(dimensions));

  for (uint32_t i = 0; i < distribution->count; ) {
    uint32_t count = 0;

    for (; i < distribution->count; i++) {
      if (distribution->slices[i]->total_probability < pivot) {
        break;
      }

      count++;
      region_probability += distribution->slices[i]->total_probability;

      if (0 != (distribution->slices[i]->flags &
          SLICE_FLAGS_ERROR_BOUND_WARNING))
      {
        bound_error_count++;
      }

      if (distribution->slices[i]->dimension > MAX_SLICE_DIMENSION) {
        critical("distribution_export_info(): The dimension (%u) of slice %u "
          "exceeds MAX_SLICE_DIMENSION (%u).", 
            distribution->slices[i]->dimension, i, MAX_SLICE_DIMENSION);
      }

      dimensions[distribution->slices[i]->dimension]++;
    }

    fprintf(file, " Probability <= %Lg: %u slice(s) (%u - %u) (%.8Lg) (%u)\n",
      pivot,
      count,
      begin_i,
      i - 1,
      region_probability,
      bound_error_count);

    for (uint32_t dimension = 0; dimension <= MAX_SLICE_DIMENSION; dimension++)
    {
      if (dimensions[dimension] != 0) {
        fprintf(file, "  Dimension %u: %u slice(s)\n",
          dimension,
          dimensions[dimension]);
      }
    }

    begin_i = i;
    pivot /= (long double)10;

    region_probability = 0;
    bound_error_count = 0;

    memset(dimensions, 0, sizeof(dimensions));
  }

  if (export_slices) {
    fprintf(file, "\n");

    long double filtered_probability = 0;
    long double filtered_error = 0;
    uint32_t filtered_count = 0;

    for (uint32_t i = 0; i < distribution->count; i++) {
      const uint32_t dimension = distribution->slices[i]->dimension;
      const long double probability =
        distribution->slices[i]->total_probability;
      const long double error = distribution->slices[i]->total_error;

      const int32_t alpha_d = distribution->slices[i]->min_log_alpha_d;
      const int32_t alpha_r = distribution->slices[i]->min_log_alpha_r;

      fprintf(file, "Slice: %u (alpha_d: %d, alpha_r: %d, "
        "dimension: %u, probability: %Lg, error: %Lg)",
          i, alpha_d, alpha_r, dimension, probability, error);

      if ((probability > 0) && ((error / probability) > 0.01f)) {
        fprintf(file, " **");
      } else {
        filtered_probability += probability;
        filtered_error += error;
        filtered_count++;
      }

      fprintf(file, "\n");
    }

    fprintf(file, "\nFiltered number of slices: %u / %u\n",
      filtered_count, distribution->count);
    fprintf(file, "Filtered probability: %.24Lg\n", filtered_probability);
    fprintf(file, "Filtered error: %.24Lg\n", filtered_error);
  }

  fprintf(file, "\nTotal count: %u\n", distribution->count);
  fprintf(file, "Total probability: %.24Lg\n", distribution->total_probability);
  fprintf(file, "Total error: %.24Lg\n", distribution->total_error);
}
