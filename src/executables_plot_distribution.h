/*!
 * \file    executables_plot_distribution.h
 * \ingroup plot_executable
 *
 * \brief   The declaration of constants shared by the executables for plotting
 *          probability distributions.
 */

/*!
 * \defgroup plot_executable Executables for plotting distributions
 * \ingroup  executable
 *
 * \brief    A module for executables for plotting probability distributions.
 */

#ifndef EXECUTABLES_PLOT_DISTRIBUTION_H
#define EXECUTABLES_PLOT_DISTRIBUTION_H

/*!
 * \brief   The dimension to which to scale two-dimensional probability
 *          distributions prior to plotting said distribution.
 *
 * If the distribution has smaller dimension, the distribution is not scaled.
 */
#define SCALED_DIMENSION 16

/*!
 * \brief   The dimension to which to scale linear probability distributions
 *          prior to plotting said distribution.
 *
 * If the distribution has smaller dimension, the distribution is not scaled.
 */
#define SCALED_DIMENSION_LINEAR 256

#endif /* EXECUTABLES_PLOT_DISTRIBUTION_H */
