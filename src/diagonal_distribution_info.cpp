/*!
 * \file    diagonal_distribution_info.cpp
 * \ingroup diagonal_distribution
 *
 * \brief   The definitions of functions for exporting information about
 *          linear probability distributions.
 */

#include "diagonal_distribution.h"

#include "errors.h"
#include "math.h"

#include <gmp.h>

#include <stdint.h>
#include <stdio.h>
#include <string.h>

/*!
 * \brief   The maximum supported slice dimension when exporting information
 *          on linear probability distributions.
 */
#define MAX_SLICE_DIMENSION                     16384

void diagonal_distribution_export_info(
  const Diagonal_Distribution * const distribution,
  FILE * const file,
  const bool export_slices,
  const bool export_regions)
{
  fprintf(file, "Precision: %u bits\n\n", distribution->precision);

  fprintf(file, "Parameters:\n");
  fprintf(file, " m:         %u\n", distribution->parameters.m);
  fprintf(file, " sigma:     %u\n", distribution->parameters.sigma);
  fprintf(file, " s:         %u\n", distribution->parameters.s);
  fprintf(file, " l:         %u\n\n", distribution->parameters.l);

  fprintf(file, " eta_bound: %u\n\n", distribution->parameters.eta_bound);

  gmp_fprintf(file, " r:         %Zd\n", distribution->parameters.r);
  gmp_fprintf(file, " d:         %Zd\n\n", distribution->parameters.d);

  fprintf(file, " Flags:     %.8x\n\n", distribution->flags);

  /* Print information on the slices. */
  long double pivot = 0.1f;
  while (distribution->slices[0]->total_probability < pivot) {
    pivot /= (long double)10;
  }

  uint32_t begin_i = 0;

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

      if (distribution->slices[i]->dimension > MAX_SLICE_DIMENSION) {
        critical("diagonal_distribution_export_info(): The dimension (%u) of "
          "slice %u exceeds MAX_SLICE_DIMENSION (%u).",
            distribution->slices[i]->dimension, i, MAX_SLICE_DIMENSION);
      }

      dimensions[distribution->slices[i]->dimension]++;
    }

    fprintf(file, " Probability <= %Lg: %u slice(s) (%u - %u) (%.8Lg)\n",
      pivot,
      count,
      begin_i,
      i - 1,
      region_probability);

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

    memset(dimensions, 0, sizeof(dimensions));
  }

  if (export_slices) {
    fprintf(file, "\n");

    for (uint32_t i = 0; i < distribution->count; i++) {
      const uint32_t dimension = distribution->slices[i]->dimension;
      const long double probability =
        distribution->slices[i]->total_probability;

      const int32_t alpha_r = distribution->slices[i]->min_log_alpha_r;
      const int32_t eta = distribution->slices[i]->eta;

      fprintf(file,
        "Slice: %u (alpha_r: %d, eta: %d, dimension: %u, probability: %Lg)\n",
          i, alpha_r, eta, dimension, probability);
    }
  }

  if (export_regions) {
    fprintf(file, "\n");

    for (uint32_t i = 0; i < distribution->count; i++) {
      const uint32_t dimension = distribution->slices[i]->dimension;
      const long double probability =
        distribution->slices[i]->total_probability;

      const int32_t alpha_r = distribution->slices[i]->min_log_alpha_r;
      const int32_t eta = distribution->slices[i]->eta;

      fprintf(file,
        "Slice: %u (alpha_r: %d, eta: %d, dimension: %u, probability: %Lg)\n",
          i, alpha_r, eta, dimension, probability);

      for (uint32_t j = 0; j < dimension; j++) {
        double region_alpha_r = (double)(abs_i(alpha_r));
        region_alpha_r += ((double)j) / (double)dimension;
        region_alpha_r *= ((double)sgn_i(alpha_r));

        fprintf(file, " Region: %u (alpha_r: %f, probability: %Lg)\n",
          j, region_alpha_r, distribution->slices[i]->norm_vector[j]);
      }
    }
  }

  fprintf(file, "\nTotal count: %u\n", distribution->count);
  fprintf(file, "Total probability: %.24Lg\n", distribution->total_probability);
}
