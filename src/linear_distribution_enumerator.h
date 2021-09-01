/*!
 * \file    linear_distribution_enumerator.h
 * \ingroup linear_distribution_enumerator
 *
 * \brief   The declaration of functions for enumerating the slices of
 *          linear probability distributions.
 */

/*!
 * \defgroup linear_distribution_enumerator Linear distribution enumerator
 * \ingroup  linear_distribution
 *
 * \brief    A module for a enumerating over the slices that are part of a
 *           probability distribution.
 *
 * The enumerator is used to determine the order in which to compute the
 * slices that make up the probability distribution. The enumerator sorts
 * the slices in an order such that the slices that are likely to contain large
 * portions of the probability mass are processed first. This enables the
 * computation to be aborted once a certain fraction of the probability mass
 * has been captured.
 */

#ifndef LINEAR_DISTRIBUTION_ENUMERATOR_H
#define LINEAR_DISTRIBUTION_ENUMERATOR_H

#include "parameters.h"

#include <stdint.h>

/*!
 * \brief   A data structure used to identify the coordinates of a slice.
 */
typedef struct {
  /*!
   * \brief   The alpha_d- or alpha_r-coordinate of the slice.
   */
  int32_t min_log_alpha;

  /*!
   * \brief   A distance measure used to sort the coordinates so that
   *          coordinates for slices that are more likely to capture a large
   *          fraction of the total probability mass are processed first.
   *
   * The distance is abs(abs(min_log_alpha) - m).
   *
   * This entry is not really needed, as it can be computed from the values of
   * the parameter m and Linear_Distribution_Coordinate::min_log_alpha. The
   * entry is included to make the distance available to qsort(), as m cannot be
   * passed to qsort(), and as qsort_r() that does take an additional parameter
   * is non-portable in that it has different prototypes on different
   * platforms(!)
   *
   * Previously sorting was implemented natively using bubble sort, but this is
   * inefficient. In the future, we may consider adding a native efficient
   * sorting algorithm to get rid of this entry.
   */
  uint32_t distance;
} Linear_Distribution_Coordinate;

/*!
 * \brief   A data structure used as a handle for an enumeration operation.
 *
 * \ingroup linear_distribution_enumerator
 */
typedef struct {
  /*!
   * \brief   Pointers to the slice coordinates.
   */
  Linear_Distribution_Coordinate ** coordinates;

  /*!
   * \brief   The offset within the
   *          Linear_Distribution_Enumerator::coordinates vector.
   */
  uint32_t offset;

  /*!
   * \brief   The total number of entries in the
   *          Linear_Distribution_Enumerator::coordinates vector.
   */
  uint32_t count;
} Linear_Distribution_Enumerator;


/*!
 * \name Initialization
 * \{
 */

/*!
 * \brief   Initializes the enumerator.
 *
 * \param[in, out] enumerator The enumerator to be initialized.
 * \param[in] parameters      The distribution parameters.
 * \param[in] mirrored        A flag indicating if the distribution is to be
 *                            mirrored.
 */
void linear_distribution_enumerator_init(
  Linear_Distribution_Enumerator * const enumerator,
  const Parameters * const parameters,
  const bool mirrored);

/*!
 * \brief   Clears the enumerator.
 *
 * \param[in, out] enumerator The enumerator to be cleared.
 */
void linear_distribution_enumerator_clear(
  Linear_Distribution_Enumerator * const enumerator);

/*!
 * \}
 */

/*!
 * \name Enumerating
 * \{
 */

/*!
 * \brief   Fetches the coordinate of the next slice from the enumerator.
 *
 * \param[out] min_log_alpha  A pointer to an integer to set to the signed
 *                            logarithmic alpha_d- or alpha_r-coordinate.
 * \param[in] enumerator      An initialized enumerator.
 *
 * \return True if the coordinate of the next slice was returned, False
 *         otherwise in which case all slices have been enumerated.
 */
bool linear_distribution_enumerator_next(
  int32_t * const min_log_alpha,
  Linear_Distribution_Enumerator * const enumerator);

/*!
 * \brief   Counts the number of slices to be enumerated in total.
 *
 * \param[in] enumerator  An initialized enumerator.
 */
uint32_t linear_distribution_enumerator_count(
  const Linear_Distribution_Enumerator * const enumerator);

/*!
 * \}
 */

#endif /* LINEAR_DISTRIBUTION_ENUMERATOR_H */
