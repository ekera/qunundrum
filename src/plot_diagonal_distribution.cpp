/*!
 * \file    plot_diagonal_distribution.cpp
 * \ingroup plot_distribution
 *
 * \brief   The definition of functions for plotting diagonal distributions.
 */

#include "plot_distribution.h"

#include "plot_distribution_axis.h"
#include "plot_distribution_common.h"

#include "distribution.h"
#include "distribution_slice.h"
#include "diagonal_distribution.h"
#include "diagonal_distribution_slice.h"
#include "parameters.h"
#include "common.h"
#include "errors.h"

#include "math.h"

#include <stdio.h>
#include <stdint.h>

/*!
 * \name Diagonal distributions
 * \{
 */

void plot_diagonal_distribution_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const int32_t offset_alpha_d,
  const bool absolute)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_extract_diagonal(
    &linear_distribution, distribution, offset_alpha_d);

  plot_linear_distribution_horizontal(
    &linear_distribution,
    offset_y,
    file,
    absolute);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

void plot_diagonal_distribution_detailed_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const int32_t offset_alpha_d,
  const bool absolute)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_extract_diagonal(
    &linear_distribution, distribution, offset_alpha_d);

  plot_linear_distribution_detailed_horizontal(
    &linear_distribution,
    offset_y,
    file,
    absolute);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

void plot_diagonal_distribution_histogram_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const int32_t offset_alpha_d,
  const bool absolute)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_extract_diagonal(
    &linear_distribution, distribution, offset_alpha_d);

  plot_linear_distribution_histogram_horizontal(
    &linear_distribution,
    offset_y,
    file,
    absolute);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

void plot_diagonal_distribution_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const int32_t offset_alpha_d)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_extract_diagonal(
    &linear_distribution, distribution, offset_alpha_d);

  plot_linear_distribution_vertical(
    &linear_distribution,
    offset_y,
    file);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

void plot_diagonal_distribution_detailed_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const int32_t offset_alpha_d)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_extract_diagonal(
    &linear_distribution, distribution, offset_alpha_d);

  plot_linear_distribution_detailed_vertical(
    &linear_distribution,
    offset_y,
    file);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

void plot_diagonal_distribution_histogram_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const int32_t offset_alpha_d)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_extract_diagonal(
    &linear_distribution, distribution, offset_alpha_d);

  plot_linear_distribution_histogram_vertical(
    &linear_distribution,
    offset_y,
    file);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

/*!
 * \}
 */

/*!
 * \name Collapsed diagonal distributions
 * \{
 */

void plot_collapsed_diagonal_distribution_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const bool absolute)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_collapse_diagonal(
    &linear_distribution, distribution);

  plot_linear_distribution_horizontal(
    &linear_distribution,
    offset_y,
    file,
    absolute);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

void plot_collapsed_diagonal_distribution_detailed_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const bool absolute)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_collapse_diagonal(
    &linear_distribution, distribution);

  plot_linear_distribution_detailed_horizontal(
    &linear_distribution,
    offset_y,
    file,
    absolute);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

void plot_collapsed_diagonal_distribution_histogram_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const bool absolute)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_collapse_diagonal(
    &linear_distribution, distribution);

  plot_linear_distribution_histogram_horizontal(
    &linear_distribution,
    offset_y,
    file,
    absolute);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

void plot_collapsed_diagonal_distribution_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_collapse_diagonal(
    &linear_distribution, distribution);

  plot_linear_distribution_vertical(
    &linear_distribution,
    offset_y,
    file);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

void plot_collapsed_diagonal_distribution_detailed_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_collapse_diagonal(
    &linear_distribution, distribution);

  plot_linear_distribution_detailed_vertical(
    &linear_distribution,
    offset_y,
    file);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

void plot_collapsed_diagonal_distribution_histogram_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file)
{
  Linear_Distribution linear_distribution;

  linear_distribution_init_collapse_diagonal(
    &linear_distribution, distribution);

  plot_linear_distribution_histogram_vertical(
    &linear_distribution,
    offset_y,
    file);

  /* Clean up. */
  linear_distribution_clear(&linear_distribution);
}

/*!
 * \}
 */
