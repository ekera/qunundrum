/*!
 * \file    plot_distribution_axis.h
 * \ingroup plot_distribution
 *
 * \brief   The declaration of functions for plotting distribution axes.
 */

#ifndef PLOT_DISTRIBUTION_AXIS_H
#define PLOT_DISTRIBUTION_AXIS_H

#include <stdio.h>
#include <stdint.h>

#include "common.h"

/*!
 * \brief Plots a vertical axis at x.
 *
 * \param[in] x                 The x coordinate.
 * \param[in] m                 The parameter m.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the axes plotted.
 */
void plot_distribution_vertical_axis(
  const double x,
  const uint32_t m,
  FILE * const file);

/*!
 * \brief Plots a horizontal axis at y.
 *
 * \param[in] y                 The y coordinate.
 * \param[in] m                 The parameter m.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the axes plotted.
 * \param[in] absolute          Assume the distribution to be symmetric and
 *                              plot only the positive side of the axis.
 */
void plot_distribution_horizontal_axis(
  const double y,
  const uint32_t m,
  FILE * const file,
  const bool absolute = FALSE);

/*!
 * \brief Plots horizontal and vertical axes in a crosshair with the origin
 *        in (x, y) = (0, 0).
 *
 * This function draws an alpha_d label for the horizontal axis and a alpha_r
 * label for the vertical axis.
 *
 * \param[in] m                 The parameter m.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the axes plotted.
 */
void plot_distribution_axes(
  const uint32_t m,
  FILE * const file);

#endif /* PLOT_DISTRIBUTION_AXIS_H */
