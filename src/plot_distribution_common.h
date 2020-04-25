/*!
 * \file    plot_distribution_common.h
 * \ingroup plot_distribution
 *
 * \brief   The definition of constants, and of an inline functions for mapping
 *          logarithmic arguments to coordinates.
 */

#ifndef PLOT_DISTRIBUTION_COMMON_H
#define PLOT_DISTRIBUTION_COMMON_H

#include "math.h"

#include <stdint.h>

/*!
 * \name Geometry
 * \{
 */

/*!
 * \brief The maximum offset from m covered on the signed logarithmic argument
 *        line, or in each quadrant in the two dimensional plot.
 *
 * This define is used to draw axes and to separate plots.
 *
 * See also #PLOT_DISTRIBUTION_LINEAR_MAX_SIZE.
 */
#define PLOT_DISTRIBUTION_MAX_OFFSET_M          11.0f

/*!
 * \brief The maximum offset from m covered on the signed logarithmic argument
 *        line, or in each quadrant in the two dimensional plot, for which
 *        there are slots and tick marks drawn.
 *
 * See also #PLOT_DISTRIBUTION_MAX_OFFSET_M.
 */
#define PLOT_DISTRIBUTION_MAX_OFFSET_SLOTS_M     9.0f

/*!
 * \brief The maximum amplitude in a linear plot.
 *
 * This is the maximum height or width of the blue curve or histogram plotted
 * in a horizontal or vertical plot, respectively.
 */
#define PLOT_DISTRIBUTION_LINEAR_MAX_AMPLITUDE   3.0f

/*!
 * \brief The maximum size of the linear plot.
 *
 * This is the maximum height of a horizontal plot, or the maximum width of a
 * vertical plot. This define is used to control e.g. spacing between plots.
 *
 * See also #PLOT_DISTRIBUTION_MAX_OFFSET_M.
 */
#define PLOT_DISTRIBUTION_LINEAR_MAX_SIZE        8.0f

/*!
 * \brief The scale factor used to scale the figure for the plot.
 */
#define PLOT_DISTRIBUTION_SCALE                  5.0f

/*!
 * \}
 */

/*!
 * \name Tick marks
 * \{
 */

/*!
 * \brief The size of a minor tick mark.
 */
#define PLOT_DISTRIBUTION_TICK_SIZE_MINOR       0.05f

/*!
 * \brief The weight of a major tick mark.
 */
#define PLOT_DISTRIBUTION_TICK_SIZE_MAJOR       0.10f

/*!
 * \brief The offset between the tick mark and the associated label.
 */
#define PLOT_DISTRIBUTION_TICK_LABEL_OFFSET     0.25f

/*!
 * \}
 */

/*!
 * \brief Given x = sgn(alpha) log_2(abs(alpha)) this function computes the
 *        coordinate of alpha on the signed logarithmic alpha axis.
 *
 * \param[in] x   The variable x = sgn(alpha) log_2(abs(alpha)).
 * \param[in] m   The parameter m.
 *
 * \return The coordinate of alpha on the signed logarithmic alpha axis.
 */
static inline double plot_distribution_coordinate(
  const double x,
  const uint32_t m)
{
  return (x - sgn_d(x) *
    ((double)m - PLOT_DISTRIBUTION_MAX_OFFSET_M)) / PLOT_DISTRIBUTION_SCALE;
}

#endif /* PLOT_DISTRIBUTION_COMMON_H */
