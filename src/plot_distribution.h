/*!
 * \file    plot_distribution.h
 * \ingroup plot_distribution
 *
 * \brief   The declaration of functions for plotting distributions.
 */

/*!
 * \defgroup plot_distribution Plotting distributions
 * \ingroup  distribution
 *
 * \brief    A module for functions for plotting probability distributions.
 */

#ifndef PLOT_DISTRIBUTION_H
#define PLOT_DISTRIBUTION_H

#include "distribution.h"
#include "linear_distribution.h"
#include "diagonal_distribution.h"

#include <stdio.h>
#include <stdint.h>

/*!
 * \name Two-dimensional distributions
 * \{
 */

/*!
 * \brief  Plots a two-dimensional distribution in the signed logarithmic
 *         (alpha_d, alpha_r)-plane.
 *
 * \remark This function sorts the slices in the distribution in place.
 *
 * \param[in, out] distribution The distribution to plot in the plane.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 */
void plot_detailed_distribution(
  Distribution * const distribution,
  FILE * const file);

/*!
 * \}
 */

/*!
 * \name Linear distributions
 * \{
 */

/*!
 * \brief  Plots a detailed linear distribution, and a linear histogram
 *         distribution, along the signed logarithmic alpha_d or alpha_r axis
 *         with a horizontal orientation.
 *
 * \remark This function sorts the slices in the distribution in place.
 *
 * \param[in, out] distribution The linear distribution to plot on the axis.
 * \param[in]      offset_y     The offset of the horizontal axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] absolute          Assume the distribution to be symmetric and
 *                              plot only the positive side of the axis.
 */
void plot_linear_distribution_horizontal(
  Linear_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const bool absolute = FALSE);

/*!
 * \brief  Plots a linear detailed distribution along the signed logarithmic
 *         alpha_d or alpha_r axis with a horizontal orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * plot_linear_distribution_horizontal(), that draws both a histogram plot and
 * a detailed plot, and places one label for both axes that result.
 *
 * \remark This function sorts the slices in the distribution in place.
 *
 * \param[in, out] distribution The linear distribution to plot on the axis.
 * \param[in]      offset_y     The offset of the horizontal axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] absolute          Assume the distribution to be symmetric and
 *                              plot only the positive side of the axis.
 */
void plot_linear_distribution_detailed_horizontal(
  Linear_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const bool absolute = FALSE);

/*!
 * \brief  Plots a linear histogram distribution along the signed logarithmic
 *         alpha_d or alpha_r axis with a horizontal orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * plot_linear_distribution_horizontal(), that draws both a histogram plot and
 * a detailed plot, and places one label for both axes that result.
 *
 * \remark This function sorts the slices in the distribution in place.
 *
 * \param[in, out] distribution The linear distribution to plot.
 * \param[in]      offset_y     The offset of the horizontal axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] absolute          Assume the distribution to be symmetric and
 *                              plot only the positive side of the axis.
 */
void plot_linear_distribution_histogram_horizontal(
  Linear_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const bool absolute = FALSE);

/*!
 * \brief  Plots a detailed linear distribution, and a linear histogram
 *         distribution, along the signed logarithmic alpha_d or alpha_r axis
 *         with a vertical orientation.
 *
 * \remark This function sorts the slices in the distribution in place.
 *
 * \param[in, out] distribution The linear distribution to plot.
 * \param[in]      offset_x     The offset of the vertical axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 */
void plot_linear_distribution_vertical(
  Linear_Distribution * const distribution,
  const double offset_x,
  FILE * const file);

/*!
 * \brief  Plots a linear detailed distribution along the signed logarithmic
 *         alpha_d or alpha_r axis with a vertical orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * plot_linear_distribution_vertical(), that draws both a histogram plot and
 * a detailed plot, and places one label for both axes that result.
 *
 * \remark This function sorts the slices in the distribution in place.
 *
 * \param[in, out] distribution The linear distribution to plot.
 * \param[in]      offset_x     The offset of the vertical axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 */
void plot_linear_distribution_detailed_vertical(
  Linear_Distribution * const distribution,
  const double offset_x,
  FILE * const file);

/*!
 * \brief  Plots a linear histogram distribution along the signed logarithmic
 *         alpha_d or alpha_r axis with a vertical orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * plot_linear_distribution_vertical(), that draws both a histogram plot and
 * a detailed plot, and places one label for both axes that result.
 *
 * \remark This function sorts the slices in the distribution in place.
 *
 * \param[in, out] distribution The linear distribution to plot.
 * \param[in]      offset_x     The offset of the vertical axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 */
void plot_linear_distribution_histogram_vertical(
  Linear_Distribution * const distribution,
  const double offset_x,
  FILE * const file);

/*!
 * \}
 */

/*!
 * \name Diagonal distributions
 * \{
 */

/*!
 * \brief  Plots a detailed diagonal distribution, and a diagonal histogram
 *         distribution, along the signed logarithmic alpha_r axis with a
 *         horizontal orientation.
 *
 * \param[in, out] distribution The diagonal distribution to plot on the axis.
 * \param[in]      offset_y     The offset of the horizontal axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] offset_alpha_d    The offset in alpha_d.
 * \param[in] absolute          Assume the distribution to be symmetric and
 *                              plot only the positive side of the axis.
 */
void plot_diagonal_distribution_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const int32_t offset_alpha_d,
  const bool absolute = FALSE);

/*!
 * \brief  Plots a diagonal detailed distribution along the signed logarithmic
 *         alpha_r axis with a horizontal orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * other functions that draw multiple plots and label all axes.
 *
 * \param[in, out] distribution The diagonal distribution to plot.
 * \param[in]      offset_y     The offset of the horizontal axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] offset_alpha_d    The offset in alpha_d.
 * \param[in] absolute          Assume the distribution to be symmetric and
 *                              plot only the positive side of the axis.
 */
void plot_diagonal_distribution_detailed_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const int32_t offset_alpha_d,
  const bool absolute = FALSE);

/*!
 * \brief  Plots a diagonal histogram distribution along the signed logarithmic
 *         alpha_r axis with a horizontal orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * other functions that draw multiple plots and label all axes.
 *
 * \param[in, out] distribution The diagonal distribution to plot.
 * \param[in]      offset_y     The offset of the horizontal axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] offset_alpha_d    The offset in alpha_d.
 * \param[in] absolute          Assume the distribution to be symmetric and
 *                              plot only the positive side of the axis.
 */
void plot_diagonal_distribution_histogram_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const int32_t offset_alpha_d,
  const bool absolute = FALSE);

/*!
 * \brief  Plots a detailed diagonal distribution, and a diagonal histogram
 *         distribution, along the signed logarithmic alpha_r axis with a
 *         vertical orientation.
 *
 * \param[in, out] distribution The diagonal distribution to plot on the axis.
 * \param[in]      offset_x     The offset of the vertical axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] offset_alpha_d    The offset in alpha_d.
 */
void plot_diagonal_distribution_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_x,
  FILE * const file,
  const int32_t offset_alpha_d);

/*!
 * \brief  Plots a diagonal detailed distribution along the signed logarithmic
 *         alpha_r axis with a vertical orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * other functions that draw multiple plots and label all axes.
 *
 * \param[in, out] distribution The diagonal distribution to plot.
 * \param[in]      offset_x     The offset of the vertical axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] offset_alpha_d    The offset in alpha_d.
 */
void plot_diagonal_distribution_detailed_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_x,
  FILE * const file,
  const int32_t offset_alpha_d);

/*!
 * \brief  Plots a diagonal histogram distribution along the signed logarithmic
 *         alpha_r axis with a vertical orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * other functions that draw multiple plots and label all axes.
 *
 * \param[in, out] distribution The diagonal distribution to plot.
 * \param[in]      offset_x     The offset of the vertical axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] offset_alpha_d    The offset in alpha_d.
 */
void plot_diagonal_distribution_histogram_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_x,
  FILE * const file,
  const int32_t offset_alpha_d);

/*!
 * \}
 */

/*!
 * \name Collapsed diagonal distributions
 * \{
 */

/*!
 * \brief  Plots a detailed collapsed diagonal distribution, and a collapsed
 *         diagonal histogram distribution, along the signed logarithmic
 *         alpha_r axis with a horizontal orientation.
 *
 * \param[in, out] distribution The diagonal distribution to plot on the axis.
 * \param[in]      offset_y     The offset of the horizontal axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] absolute          Assume the distribution to be symmetric and
 *                              plot only the positive side of the axis.
 */
void plot_collapsed_diagonal_distribution_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const bool absolute = FALSE);

/*!
 * \brief  Plots a collapsed diagonal detailed distribution along the signed
 *         logarithmic alpha_r axis with a horizontal orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * other functions that draw multiple plots and label all axes.
 *
 * \param[in, out] distribution The diagonal distribution to plot.
 * \param[in]      offset_y     The offset of the horizontal axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] absolute          Assume the distribution to be symmetric and
 *                              plot only the positive side of the axis.
 */
void plot_collapsed_diagonal_distribution_detailed_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const bool absolute = FALSE);

/*!
 * \brief  Plots a collapsed diagonal histogram distribution along the signed
 *         logarithmic alpha_r axis with a horizontal orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * other functions that draw multiple plots and label all axes.
 *
 * \param[in, out] distribution The diagonal distribution to plot.
 * \param[in]      offset_y     The offset of the horizontal axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] absolute          Assume the distribution to be symmetric and
 *                              plot only the positive side of the axis.
 */
void plot_collapsed_diagonal_distribution_histogram_horizontal(
  Diagonal_Distribution * const distribution,
  const double offset_y,
  FILE * const file,
  const bool absolute = FALSE);

/*!
 * \brief  Plots a detailed collapsed diagonal distribution, and a collapsed
 *         diagonal histogram distribution, along the signed logarithmic
 *         alpha_r axis with a vertical orientation.
 *
 * \param[in, out] distribution The diagonal distribution to plot on the axis.
 * \param[in]      offset_x     The offset of the vertical axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 * \param[in] offset_alpha_d    The offset in alpha_d.
 */
void plot_collapsed_diagonal_distribution_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_x,
  FILE * const file,
  const int32_t offset_alpha_d);

/*!
 * \brief  Plots a collapsed diagonal detailed distribution along the signed
 *         logarithmic alpha_r axis with a vertical orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * other functions that draw multiple plots and label all axes.
 *
 * \param[in, out] distribution The diagonal distribution to plot.
 * \param[in]      offset_x     The offset of the vertical axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 */
void plot_collapsed_diagonal_distribution_detailed_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_x,
  FILE * const file);

/*!
 * \brief  Plots a collapsed diagonal histogram distribution along the signed
 *         logarithmic alpha_r axis with a vertical orientation.
 *
 * This function does not label the axis as it is designed to be called by
 * other functions that draw multiple plots and label all axes.
 *
 * \param[in, out] distribution The diagonal distribution to plot.
 * \param[in]      offset_x     The offset of the vertical axis.
 * \param[in, out] file         The file to which to write Latex source for
 *                              the figure plotted.
 */
void plot_collapsed_diagonal_distribution_histogram_vertical(
  Diagonal_Distribution * const distribution,
  const double offset_x,
  FILE * const file);

/*!
 * \}
 */

#endif /* PLOT_DISTRIBUTION_H */
