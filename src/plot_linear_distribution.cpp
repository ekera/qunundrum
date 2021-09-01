/*!
 * \file    plot_linear_distribution.cpp
 * \ingroup plot_distribution
 *
 * \brief   The definition of functions for plotting linear distributions.
 */

#include "plot_distribution.h"

#include "common.h"
#include "errors.h"
#include "linear_distribution.h"
#include "linear_distribution_slice.h"
#include "math.h"
#include "plot_distribution_axis.h"
#include "plot_distribution_common.h"

#include <stdint.h>
#include <stdio.h>

void plot_linear_distribution_horizontal(
  Linear_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const bool absolute)
{
  plot_linear_distribution_histogram_horizontal(
    distribution,
    offset_y,
    file,
    absolute);

  plot_linear_distribution_detailed_horizontal(
    distribution,
    offset_y - (PLOT_DISTRIBUTION_LINEAR_MAX_SIZE + 1) /
      PLOT_DISTRIBUTION_SCALE,
    file,
    absolute);

  /* Label the two axes with a single label. */
  const uint32_t m = distribution->parameters.m;
  const uint32_t flags = distribution->flags;

  if ((flags & LINEAR_DISTRIBUTION_FLAG_R) != 0) {
    if (absolute) {
      fprintf(file,  "\\draw (%f, %f) "
        "node[above, rotate=-90] "
          "{\\tiny $\\log_2(|\\alpha_r|)$};\n",
            plot_distribution_coordinate(m + PLOT_DISTRIBUTION_MAX_OFFSET_M, m),
            offset_y -
              (PLOT_DISTRIBUTION_LINEAR_MAX_SIZE + 1) /
                PLOT_DISTRIBUTION_SCALE / 2.0f);
    } else {
      fprintf(file,  "\\draw (%f, %f) "
        "node[above, rotate=-90] "
          "{\\tiny $\\text{sgn}(\\alpha_r) \\log_2(|\\alpha_r|)$};\n",
            plot_distribution_coordinate(m + PLOT_DISTRIBUTION_MAX_OFFSET_M, m),
            offset_y -
              (PLOT_DISTRIBUTION_LINEAR_MAX_SIZE + 1) /
                PLOT_DISTRIBUTION_SCALE / 2.0f);
    }
  } else {
    if (absolute) {
      fprintf(file,  "\\draw (%f, %f) "
        "node[above, rotate=-90] "
          "{\\tiny $\\log_2(|\\alpha_d|)$};\n",
            plot_distribution_coordinate(m + PLOT_DISTRIBUTION_MAX_OFFSET_M, m),
            offset_y -
              (PLOT_DISTRIBUTION_LINEAR_MAX_SIZE + 1) /
                PLOT_DISTRIBUTION_SCALE / 2.0f);
    } else {
      fprintf(file,  "\\draw (%f, %f) "
        "node[above, rotate=-90] "
          "{\\tiny $\\text{sgn}(\\alpha_d) \\log_2(|\\alpha_d|)$};\n",
            plot_distribution_coordinate(m + PLOT_DISTRIBUTION_MAX_OFFSET_M, m),
            offset_y -
              (PLOT_DISTRIBUTION_LINEAR_MAX_SIZE + 1) /
                PLOT_DISTRIBUTION_SCALE / 2.0f);
    }
  }
}

void plot_linear_distribution_detailed_horizontal(
    Linear_Distribution * const distribution,
    const double offset_y,
    FILE * const file,
    const bool absolute)
{
  const uint32_t m = distribution->parameters.m;

  if (distribution->count > 0) {
    /* Sort the distribution. */
    linear_distribution_sort_slices(distribution);

    /* Compute the maximum regional probability. */
    double max_region_probability = 0;

    for (uint32_t i = 0; i < distribution->count; i++) {
      const Linear_Distribution_Slice * const slice = distribution->slices[i];

      const uint32_t dimension = slice->dimension;

      for (uint32_t j = 0; j < dimension; j++) {
        if (slice->norm_vector[j] > max_region_probability) {
          max_region_probability = slice->norm_vector[j];
        }
      }
    }

    /* Plot the detailed linear distribution. */
    for (uint32_t i = 0; i < distribution->count; i++) {
      const Linear_Distribution_Slice * const slice = distribution->slices[i];

      const uint32_t dimension = slice->dimension;

      for (uint32_t j = 0; j < dimension; j++) {
        const double region_probability = slice->norm_vector[j];

        double y =
          PLOT_DISTRIBUTION_LINEAR_MAX_AMPLITUDE *
            (double)(region_probability / max_region_probability) /
              PLOT_DISTRIBUTION_SCALE;

        if (y < 0.0001) {
          continue;
        }

        double alpha_min;
        double alpha_max;

        linear_distribution_slice_region_coordinates(
          slice,
          j,
          &alpha_min,
          &alpha_max);

        if ((alpha_min < 0) && (alpha_max < 0)) {
          if (absolute) {
            continue;
          }

          if (alpha_min < -((double)m) - PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M)
          {
            continue;
          }

          if (alpha_max >= -((double)m) + PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M)
          {
            continue;
          }
        } else if ((alpha_min >= 0) && (alpha_max >= 0)) {
          if (alpha_min < ((double)m) - PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
            continue;
          }

          if (alpha_max >= ((double)m) + PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
            continue;
          }
        } else {
          critical("plot_linear_distribution_detailed_horizontal(): "
            "Unexpected region close to the origin.");
        }

        fprintf(
          file,
          "\\fill[blue] (%f, %f) rectangle (%f, %f); \n",
          plot_distribution_coordinate(alpha_min, m),
          offset_y,
          plot_distribution_coordinate(alpha_max, m),
          offset_y + y);
      }
    }
  }

  /* Draw axis. */
  plot_distribution_horizontal_axis(offset_y, m, file, absolute);

  fprintf(
    file,
    "\\draw[thin, -stealth] (%f, %f) -- (%f, %f);\n",
    0.0f,
    offset_y - PLOT_DISTRIBUTION_TICK_SIZE_MINOR,
    0.0f,
    offset_y +
      ((double)PLOT_DISTRIBUTION_LINEAR_MAX_AMPLITUDE + (double)1.0f) /
        (double)PLOT_DISTRIBUTION_SCALE);
}

void plot_linear_distribution_histogram_horizontal(
  Linear_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  bool absolute)
{
  const uint32_t m = distribution->parameters.m;

  if (distribution->count > 0) {
    /* Sort the distribution. */
    linear_distribution_sort_slices(distribution);

    /* Compute the maximum slice probability. */
    const double max_slice_probability =
      distribution->slices[0]->total_probability;

    /* Compute the maximum regional probability. */
    double max_region_probability = 0;

    for (uint32_t i = 0; i < distribution->count; i++) {
      const Linear_Distribution_Slice * const slice = distribution->slices[i];

      const uint32_t dimension = slice->dimension;

      for (uint32_t j = 0; j < dimension; j++) {
        if (slice->norm_vector[j] > max_region_probability) {
          max_region_probability = slice->norm_vector[j];
        }
      }
    }

    /* Plot the distribution histogram. */
    for (uint32_t i = 0; i < distribution->count; i++) {
      const Linear_Distribution_Slice * const slice = distribution->slices[i];

      double y =
        PLOT_DISTRIBUTION_LINEAR_MAX_AMPLITUDE *
          (slice->total_probability / max_slice_probability) /
            PLOT_DISTRIBUTION_SCALE;

      if (y < 0.0001) {
        continue;
      }

      double alpha_min;
      double alpha_max;

      linear_distribution_slice_coordinates(
        slice,
        &alpha_min,
        &alpha_max);

      if ((alpha_min < 0) && (alpha_max < 0)) {
        if (absolute) {
          continue;
        }

        if (alpha_min < -((double)m) - PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
          continue;
        }

        if (alpha_max >= -((double)m) + PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
          continue;
        }
      } else if ((alpha_min >= 0) && (alpha_max >= 0)) {
        if (alpha_min < ((double)m) - PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
          continue;
        }

        if (alpha_max >= ((double)m) + PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
          continue;
        }
      } else {
        critical("plot_linear_distribution_histogram_horizontal(): "
          "Unexpected region close to the origin.");
      }

      const double min_x =
        plot_distribution_coordinate(alpha_min, m) +
          sgn_d(alpha_min) * 0.2f / PLOT_DISTRIBUTION_SCALE;

      const double max_x =
        plot_distribution_coordinate(alpha_max, m) -
          sgn_d(alpha_min) * 0.2f / PLOT_DISTRIBUTION_SCALE;

      fprintf(
        file,
        "\\fill[blue] (%f, %f) rectangle (%f, %f); %% %f \n",
        min_x,
        offset_y,
        max_x,
        offset_y + y,
        alpha_min);

      double probability = (double)(slice->total_probability);
      if (absolute) {
        probability *= 2;
      }

      if (probability >= 0.001) {
        fprintf(
          file,
          "\\draw (%f, %f) "
            "node[anchor=west, left, rotate=-90, scale=0.75, black] "
            "{\\tiny %.3f};\n",
          (min_x + max_x) / 2,
          offset_y + y,
          probability);
      }
    }
  }

  /* Draw axis. */
  plot_distribution_horizontal_axis(offset_y, m, file, absolute);

  fprintf(
    file,
    "\\draw[thin, -stealth] (%f, %f) -- (%f, %f);\n",
    0.0f,
    offset_y - PLOT_DISTRIBUTION_TICK_SIZE_MINOR,
    0.0f,
    offset_y +
      ((double)PLOT_DISTRIBUTION_LINEAR_MAX_AMPLITUDE + (double)1.0f) /
        (double)PLOT_DISTRIBUTION_SCALE);
}

void plot_linear_distribution_vertical(
  Linear_Distribution * const distribution,
  const double offset_x,
  FILE * const file)
{
  plot_linear_distribution_histogram_vertical(
    distribution,
    offset_x,
    file);

  plot_linear_distribution_detailed_vertical(
    distribution,
    offset_x - (PLOT_DISTRIBUTION_LINEAR_MAX_SIZE + 1) /
      PLOT_DISTRIBUTION_SCALE,
    file);

  /* Label the two axes with a single label. */
  const uint32_t m = distribution->parameters.m;
  const uint32_t flags = distribution->flags;

  if ((flags & LINEAR_DISTRIBUTION_FLAG_R) != 0) {
    fprintf(file, "\\draw (%f, %f) node[above]"
      "{\\tiny $\\text{sgn}(\\alpha_r) \\log_2(|\\alpha_r|)$};\n",
        offset_x -
          (PLOT_DISTRIBUTION_LINEAR_MAX_SIZE + 1) /
            PLOT_DISTRIBUTION_SCALE / 2.0f,
        plot_distribution_coordinate(m + PLOT_DISTRIBUTION_MAX_OFFSET_M, m));
  } else {
    fprintf(file, "\\draw (%f, %f) node[above]"
      "{\\tiny $\\text{sgn}(\\alpha_d) \\log_2(|\\alpha_d|)$};\n",
        offset_x-
          (PLOT_DISTRIBUTION_LINEAR_MAX_SIZE + 1) /
            PLOT_DISTRIBUTION_SCALE / 2.0f,
        plot_distribution_coordinate(m + PLOT_DISTRIBUTION_MAX_OFFSET_M, m));
  }
}

void plot_linear_distribution_detailed_vertical(
  Linear_Distribution * const distribution,
  const double offset_x,
  FILE * const file)
{
  const uint32_t m = distribution->parameters.m;

  if (distribution->count > 0) {
    /* Sort the distribution slices. */
    linear_distribution_sort_slices(distribution);

    /* Compute the maximum regional probability. */
    double max_region_probability = 0;

    for (uint32_t i = 0; i < distribution->count; i++) {
      const Linear_Distribution_Slice * const slice = distribution->slices[i];

      const uint32_t dimension = slice->dimension;

      for (uint32_t j = 0; j < dimension; j++) {
        if (slice->norm_vector[j] > max_region_probability) {
          max_region_probability = slice->norm_vector[j];
        }
      }
    }

    /* Plot the detailed linear distribution. */
    for (uint32_t i = 0; i < distribution->count; i++) {
      const Linear_Distribution_Slice * const slice = distribution->slices[i];

      const uint32_t dimension = slice->dimension;

      for (uint32_t j = 0; j < dimension; j++) {
        const double region_probability = slice->norm_vector[j];

        const double x =
          PLOT_DISTRIBUTION_LINEAR_MAX_AMPLITUDE *
            (double)(region_probability / max_region_probability) /
              PLOT_DISTRIBUTION_SCALE;

        double alpha_min;
        double alpha_max;

        linear_distribution_slice_region_coordinates(
          slice,
          j,
          &alpha_min,
          &alpha_max);

        if ((alpha_min < 0) && (alpha_max < 0)) {
          if (alpha_min < -((double)m) - PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
            continue;
          }

          if (alpha_max >= -((double)m) + PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M)
          {
            continue;
          }
        } else if ((alpha_min >= 0) && (alpha_max >= 0)) {
          if (alpha_min < ((double)m) - PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
            continue;
          }

          if (alpha_max >= ((double)m) + PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M)
          {
            continue;
          }
        } else {
          critical("plot_linear_distribution_detailed_vertical(): "
            "Unexpected region close to the origin.");
        }

        fprintf(
          file,
          "\\fill[blue] (%f, %f) rectangle (%f, %f); \n",
          offset_x,
          plot_distribution_coordinate(alpha_min, m),
          offset_x + x,
          plot_distribution_coordinate(alpha_max, m));
      }
    }
  }

  /* Draw axis. */
  plot_distribution_vertical_axis(offset_x, m, file);

  fprintf(
    file,
    "\\draw[thin, -stealth] (%f, %f) -- (%f, %f);\n",
    offset_x - PLOT_DISTRIBUTION_TICK_SIZE_MINOR,
    0.0f,
    offset_x +
      ((double)PLOT_DISTRIBUTION_LINEAR_MAX_AMPLITUDE + (double)1.0f) /
        (double)PLOT_DISTRIBUTION_SCALE,
    0.0f);
}

void plot_linear_distribution_histogram_vertical(
  Linear_Distribution * const distribution,
  const double offset_x,
  FILE * const file)
{
  const uint32_t m = distribution->parameters.m;

  if (distribution->count > 0) {
    /* Sort the distribution slices. */
    linear_distribution_sort_slices(distribution);

    /* Compute the maximum slice probability. */
    const double max_slice_probability =
      distribution->slices[0]->total_probability;

    /* Compute the maximum regional probability. */
    double max_region_probability = 0;

    for (uint32_t i = 0; i < distribution->count; i++) {
      const Linear_Distribution_Slice * const slice = distribution->slices[i];

      const uint32_t dimension = slice->dimension;

      for (uint32_t j = 0; j < dimension; j++) {
        if (slice->norm_vector[j] > max_region_probability) {
          max_region_probability = slice->norm_vector[j];
        }
      }
    }

    /* Plot the distribution histogram */
    for (uint32_t i = 0; i < distribution->count; i++) {
      const Linear_Distribution_Slice * const slice = distribution->slices[i];

      double alpha_min;
      double alpha_max;

      linear_distribution_slice_coordinates(
        slice,
        &alpha_min,
        &alpha_max);

      if ((alpha_min < 0) && (alpha_max < 0)) {
        if (alpha_min < -((double)m) - PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
          continue;
        }

        if (alpha_max >= -((double)m) + PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
          continue;
        }
      } else if ((alpha_min >= 0) && (alpha_max >= 0)) {
        if (alpha_min < ((double)m) - PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
          continue;
        }

        if (alpha_max >= ((double)m) + PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M) {
          continue;
        }
      } else {
        critical("plot_linear_distribution_histogram_vertical(): "
          "Unexpected region close to the origin.");
      }

      const double x =
        PLOT_DISTRIBUTION_LINEAR_MAX_AMPLITUDE *
          (slice->total_probability / max_slice_probability) /
            PLOT_DISTRIBUTION_SCALE;

      const double min_y =
        plot_distribution_coordinate(alpha_min, m) +
          sgn_d(alpha_min) * 0.2f / PLOT_DISTRIBUTION_SCALE;

      const double max_y =
        plot_distribution_coordinate(alpha_max, m) -
          sgn_d(alpha_min) * 0.2f / PLOT_DISTRIBUTION_SCALE;

      fprintf(
        file,
        "\\fill[blue] (%f, %f) rectangle (%f, %f); %% %f \n",
        offset_x,
        min_y,
        offset_x + x,
        max_y,
        alpha_min);

      if (slice->total_probability >= 0.001) {
        fprintf(
          file,
          "\\draw (%f, %f) node[anchor=east, right, scale=0.75, black] "
            "{\\tiny %.3f}; \n",
          offset_x + x,
          (min_y + max_y) / 2,
          (double)(slice->total_probability));
      }
    }
  }

  /* Draw axis. */
  plot_distribution_vertical_axis(offset_x, m, file);

  fprintf(
    file,
    "\\draw[thin, -stealth] (%f, %f) -- (%f, %f);\n",
    offset_x - PLOT_DISTRIBUTION_TICK_SIZE_MINOR,
    0.0f,
    offset_x +
      ((double)PLOT_DISTRIBUTION_LINEAR_MAX_AMPLITUDE + (double)1.0f) /
        (double)PLOT_DISTRIBUTION_SCALE,
    0.0f);
}
