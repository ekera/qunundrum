/*!
 * \file    plot_distribution_axis.cpp
 * \ingroup plot_distribution
 *
 * \brief   The definition of functions for plotting distribution axes.
 */

#include "plot_distribution_axis.h"

#include "plot_distribution_common.h"

#include "common.h"

#include "math.h"

#include <stdio.h>
#include <stdint.h>

void plot_distribution_vertical_axis(
  const double x,
  const uint32_t m,
  FILE * const file)
{
  /* Draw the axis. */
  fprintf(file, "\\draw[-stealth, thin] (%f, %f)"
    " -- (%f, -0.20) -- (%f, -0.175) -- (%f, -0.125) -- (%f, -0.10)"
    " -- (%f, 0.10) -- (%f, 0.125) -- (%f, 0.175) -- (%f, 0.20)"
    " -- (%f, %f);\n",
      x,
      plot_distribution_coordinate(-((double)m) -
        PLOT_DISTRIBUTION_MAX_OFFSET_M, m),
      x, x + 0.025, x - 0.025, x,
      x, x - 0.025, x + 0.025, x,
      x,
      plot_distribution_coordinate(m + PLOT_DISTRIBUTION_MAX_OFFSET_M, m));

  /* Draw tick marks on the axis. */
  const int32_t min_tick = m - PLOT_DISTRIBUTION_MAX_OFFSET_M + 2;
  const int32_t max_tick = m + PLOT_DISTRIBUTION_MAX_OFFSET_M - 2;

  for (int32_t tick = min_tick; tick <= max_tick; tick++) {
    const double TICK_SIZE =
      (((tick - (int32_t)m) % 5) == 0) ?
        PLOT_DISTRIBUTION_TICK_SIZE_MAJOR :
        PLOT_DISTRIBUTION_TICK_SIZE_MINOR;

    fprintf(file, "\\draw[thin] (%f, %f) -- (%f, %f);\n",
      x,
      plot_distribution_coordinate(tick, m),
      x - TICK_SIZE,
      plot_distribution_coordinate(tick, m));

    fprintf(file, "\\draw[thin] (%f, %f) -- (%f, %f);\n",
      x,
      -plot_distribution_coordinate(tick, m),
      x - TICK_SIZE,
      -plot_distribution_coordinate(tick, m));
  }

  /* Label a subset of the tick marks. */
  fprintf(file, "\\draw (%f, %f) node[rotate=-90] {\\tiny $m - 5$};\n",
    x - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET,
    plot_distribution_coordinate(m - 5, m));
  fprintf(file, "\\draw (%f, %f) node[rotate=-90] {\\tiny $m$};\n",
    x - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET,
    plot_distribution_coordinate(m, m));
  fprintf(file, "\\draw (%f, %f) node[rotate=-90] {\\tiny $m + 5$};\n",
    x - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET,
    plot_distribution_coordinate(m + 5, m));

  fprintf(file, "\\draw (%f, %f) node[rotate=-90] {\\tiny $-m - 5$};\n",
    x - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET,
    plot_distribution_coordinate(-(double)(m + 5), m));
  fprintf(file, "\\draw (%f, %f) node[rotate=-90] {\\tiny $-m$};\n",
    x - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET,
    plot_distribution_coordinate(-(double)m, m));
  fprintf(file, "\\draw (%f, %f) node[rotate=-90] {\\tiny $-m + 5$};\n",
    x - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET,
    plot_distribution_coordinate(-(double)(m - 5), m));
}

void plot_distribution_horizontal_axis(
  const double y,
  const uint32_t m,
  FILE * const file,
  const bool absolute)
{
  /* Draw the axis. */
  if (absolute) {
    fprintf(file, "\\draw[-stealth, thin] (-0.10, %f)"
      " -- (0.10, %f) -- (0.125, %f) -- (0.175, %f) -- (0.20, %f)"
      " -- (%f, %f);\n",
        y,
        y, y - 0.025, y + 0.025, y,
        plot_distribution_coordinate(m + PLOT_DISTRIBUTION_MAX_OFFSET_M, m),
        y);
  } else {
    fprintf(file, "\\draw[-stealth, thin] (%f, %f)"
      " -- (-0.20, %f) -- (-0.175, %f) -- (-0.125, %f) -- (-0.10, %f)"
      " -- (0.10, %f) -- (0.125, %f) -- (0.175, %f) -- (0.20, %f)"
      " -- (%f, %f);\n",
        plot_distribution_coordinate(-((double)m) -
          PLOT_DISTRIBUTION_MAX_OFFSET_M, m),
        y,
        y, y + 0.025, y - 0.025, y,
        y, y - 0.025, y + 0.025, y,
        plot_distribution_coordinate(m + PLOT_DISTRIBUTION_MAX_OFFSET_M, m),
        y);
  }

  /* Draw tick marks on the axis. */
  const int32_t min_tick = m - PLOT_DISTRIBUTION_MAX_OFFSET_M + 2;
  const int32_t max_tick = m + PLOT_DISTRIBUTION_MAX_OFFSET_M - 2;

  for (int32_t tick = min_tick; tick <= max_tick; tick++) {
    const double TICK_SIZE =
      (((tick - (int32_t)m) % 5) == 0) ?
        PLOT_DISTRIBUTION_TICK_SIZE_MAJOR :
        PLOT_DISTRIBUTION_TICK_SIZE_MINOR;

    fprintf(file, "\\draw[thin] (%f, %f) -- (%f, %f);\n",
      plot_distribution_coordinate(tick, m),
      y,
      plot_distribution_coordinate(tick, m),
      y - TICK_SIZE);

    if (!absolute) {
      fprintf(file, "\\draw[thin] (%f, %f) -- (%f, %f);\n",
        -plot_distribution_coordinate(tick, m),
        y,
        -plot_distribution_coordinate(tick, m),
        y - TICK_SIZE);
    }
  }

  /* Label a subset of the tick marks. */
  fprintf(file, "\\draw (%f, %f) node {\\tiny $m - 5$};\n",
    plot_distribution_coordinate(m - 5, m),
    y - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET);
  fprintf(file, "\\draw (%f, %f) node {\\tiny $m$};\n",
    plot_distribution_coordinate(m, m),
    y - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET);
  fprintf(file, "\\draw (%f, %f) node {\\tiny $m + 5$};\n",
    plot_distribution_coordinate(m + 5, m),
    y - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET);

  if (!absolute) {
    fprintf(file, "\\draw (%f, %f) node {\\tiny $-m - 5$};\n",
      plot_distribution_coordinate(-(double)(m + 5), m),
      y - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET);
    fprintf(file, "\\draw (%f, %f) node {\\tiny $-m$};\n",
      plot_distribution_coordinate(-(double)m, m),
      y - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET);
    fprintf(file, "\\draw (%f, %f) node {\\tiny $-m + 5$};\n",
      plot_distribution_coordinate(-(double)(m - 5), m),
      y - PLOT_DISTRIBUTION_TICK_LABEL_OFFSET);
  }
}

void plot_distribution_axes(
  const uint32_t m,
  FILE * const file)
{
  plot_distribution_horizontal_axis(0, m, file); /* y = 0 */

  fprintf(file, "\\draw (%f, 0) node[above, rotate=-90] " /* y = 0 */
    "{\\tiny $\\text{sgn}(\\alpha_d) \\log_2(|\\alpha_d|)$};\n",
      plot_distribution_coordinate(m + PLOT_DISTRIBUTION_MAX_OFFSET_M, m));

  plot_distribution_vertical_axis(0, m, file); /* x = 0 */

  fprintf(file, "\\draw (0, %f) node[above]" /* x = 0 */
    "{\\tiny $\\text{sgn}(\\alpha_r) \\log_2(|\\alpha_r|)$};\n",
      plot_distribution_coordinate(m + PLOT_DISTRIBUTION_MAX_OFFSET_M, m));
}
