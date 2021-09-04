/*!
 * \file    plot_distribution.cpp
 * \ingroup plot_distribution
 *
 * \brief   The definition of functions for plotting distributions.
 */

#include "plot_distribution.h"

#include "common.h"
#include "distribution.h"
#include "distribution_slice.h"
#include "math.h"
#include "plot_distribution_axis.h"
#include "plot_distribution_common.h"

#include <stdint.h>
#include <stdio.h>

void plot_detailed_distribution(
  Distribution * const distribution,
  FILE * const file)
{
  const uint32_t m = distribution->parameters.m;

  if (distribution->count > 0) {
    /* Sort the distribution slices. */
    distribution_sort_slices(distribution);

    /* The maximum slice probability is the probability of the first slice. */
    const long double max_slice_probability =
      distribution->slices[0]->total_probability;

    /* Determine the maximum region probability. */
    long double max_region_probability = 0;

    for (uint32_t i = 0; i < distribution->count; i++) {
      /* Extract the slice. */
      const Distribution_Slice * slice = distribution->slices[i];

      /* If the total probability of the slice is smaller than the maximum
       * region probability already observed then we may stop. */
      if (slice->total_probability < max_region_probability) {
        break; /* No need to continue; all regions are less probable. */
      }

      /* Test all regions in the slice. */
      const uint32_t dimension = slice->dimension;

      for (uint32_t j = 0; j < (dimension * dimension); j++) {
        if (slice->norm_matrix[j] > max_region_probability) {
          max_region_probability = slice->norm_matrix[j];
        }
      }
    }

    /* Process the slices. */
    for (uint32_t i = 0; i < distribution->count; i++) {
      const Distribution_Slice * const slice = distribution->slices[i];

      const double relative_slice_probability =
        slice->total_probability / max_slice_probability;

      if (0 == (uint32_t)(100 * relative_slice_probability)) {
        continue; /* Do not draw white slices. */
      }

      /* Iterate over all regions in the slice. */
      const uint32_t dimension = slice->dimension;

      for (uint32_t j = 0; j < (dimension * dimension); j++) {
        const double relative_region_probability =
          slice->norm_matrix[j] / max_region_probability;

        const uint32_t color = (uint32_t)(100 * relative_region_probability);

        if (0 == color) {
          continue; /* Do not draw white regions. */
        }

        /* Draw the slice. */
        double alpha_d_min;
        double alpha_d_max;
        double alpha_r_min;
        double alpha_r_max;

        distribution_slice_region_coordinates(
          slice,
          j,
          &alpha_d_min,
          &alpha_d_max,
          &alpha_r_min,
          &alpha_r_max);

        fprintf(
          file,
          "\\fill[blue!%u] (%f, %f) rectangle (%f, %f); %% (%d, %d)\n",
          color,
          plot_distribution_coordinate(alpha_d_min, m),
          plot_distribution_coordinate(alpha_r_min, m),
          plot_distribution_coordinate(alpha_d_max, m),
          plot_distribution_coordinate(alpha_r_max, m),
          slice->min_log_alpha_d,
          slice->min_log_alpha_r);
      }

      #ifdef DRAW_SLICE_BOUNDARIES
      /* Draw the slice boundary for debug purposes. */
      double alpha_d_min;
      double alpha_d_max;
      double alpha_r_min;
      double alpha_r_max;

      distribution_slice_coordinates(
        slice,
        &alpha_d_min,
        &alpha_d_max,
        &alpha_r_min,
        &alpha_r_max);

      fprintf(
        file,
        "\\draw[red, thin] (%f, %f) rectangle (%f, %f); %% (%d, %d)\n",
        plot_distribution_coordinate(alpha_d_min, m),
        plot_distribution_coordinate(alpha_r_min, m),
        plot_distribution_coordinate(alpha_d_max, m),
        plot_distribution_coordinate(alpha_r_max, m),
        slice->min_log_alpha_d,
        slice->min_log_alpha_r);

      fprintf(file, "\n");
      #endif /* DRAW_SLICE_BOUNDARIES */
    }
  }

  /* Draw axes. */
  plot_distribution_axes(m, file);
}
