/*!
 * \file    distribution_enumerator.h
 * \ingroup two_dimensional_distribution_enumerator
 *
 * \brief   The declaration of functions for enumerating the slices of
 *          two-dimensional probability distributions, and the definition of
 *          data structures used in such enumerations.
 */

/*!
 * \defgroup  two_dimensional_distribution_enumerator \
 *            Two-dimensional distribution enumerator
 * \ingroup   two_dimensional_distribution
 *
 * \brief   A module for a enumerating over the slices that are part of a
 *          probability distribution.
 *
 * The enumerator is used to determine the order in which to compute the
 * slices that make up the probability distribution. The enumerator sorts
 * the slices in an order such that the slices that are likely to contain large
 * portions of the probability mass are processed first. This enables the
 * computation to be aborted once a certain fraction of the probability mass
 * has been captured.
 */

#ifndef DISTRIBUTION_ENUMERATOR_H
#define DISTRIBUTION_ENUMERATOR_H

#include "parameters.h"

#include <stdint.h>

/*!
 * \brief   A data structure used to identify the coordinates of a slice.
 */
typedef struct {
  /*!
   * \brief   The alpha_d-coordinate of the slice.
   */
  int32_t min_log_alpha_d;

  /*!
   * \brief   The alpha_r-coordinate of the slice.
   */
  int32_t min_log_alpha_r;

  /*!
   * \brief   A distance measure used to sort the coordinates so that
   *          coordinates for slices that are more likely to capture a large
   *          fraction of the total probability mass are processed first.
   *
   * The distance is sum_{x in {d, r}} abs(abs(min_log_alpha_x) - m)^2 with the
   * caveat that the distance is forced to zero for diagonal element.
   *
   * This entry is not really needed, as it can be computed from the values of
   * the parameter m, Linear_Distribution_Coordinate::min_log_alpha_d and
   * Linear_Distribution_Coordinate::min_log_alpha_r. The entry is included to
   * make the distance available to qsort(), as m cannot be passed to qsort(),
   * and as qsort_r() that does take an additional parameter is non-portable in
   * that it has different prototypes on different platforms(!)
   *
   * Previously sorting was implemented natively using bubble sort, but this is
   * inefficient. In the future, we may consider adding a native efficient
   * sorting algorithm to get rid of this entry.
   */
  uint32_t distance;
} Distribution_Coordinate;

/*!
 * \brief   A data structure used as a handle for an enumeration operation.
 *
 * \ingroup two_dimensional_distribution_enumerator
 */
typedef struct {
  /*!
   * \brief   Pointers to the slice coordinates.
   */
  Distribution_Coordinate ** coordinates;

  /*!
   * \brief   A flag indicating whether only on side of symmetric distributions
   *          are computed and then mirrored, or whether both sides are
   *          explicitly computed. In general, mirroring is used.
   */
  bool mirrored;

  /*!
   * \brief   The offset within the Distribution_Enumerator::coordinates vector.
   */
  uint32_t offset;

  /*!
   * \brief   The total number of entries in the
   *          Distribution_Enumerator::coordinates vector.
   */
  uint32_t count;
} Distribution_Enumerator;


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
void distribution_enumerator_init(
  Distribution_Enumerator * const enumerator,
  const Parameters * const parameters,
  const bool mirrored);

/*!
 * \brief   Clears the enumerator.
 *
 * \param[in, out] enumerator The enumerator to be cleared.
 */
void distribution_enumerator_clear(
  Distribution_Enumerator * const enumerator);

/*!
 * \}
 */

/*!
 * \name Enumerating
 * \{
 */

/*!
 * \brief   Fetches the coordinates of the next slice from the enumerator.
 *
 * \param[out] min_log_alpha_d  A pointer to an integer to set to the signed
 *                              logarithmic alpha_d-coordinate.
 * \param[out] min_log_alpha_r  A pointer to an integer to set to the signed
 *                              logarithmic alpha_r-coordinate.
 * \param[in, out] enumerator   An initialized enumerator.
 *
 * \return #TRUE if the coordinate of the next slice was returned, #FALSE
 *         otherwise in which case all slices have been enumerated.
 */
bool distribution_enumerator_next(
  int32_t * const min_log_alpha_d,
  int32_t * const min_log_alpha_r,
  Distribution_Enumerator * const enumerator);

/*!
 * \brief   Counts the number of slices to be enumerated in total.
 *
 * \param[in] enumerator  An initialized enumerator.
 */
uint32_t distribution_enumerator_count(
  const Distribution_Enumerator * const enumerator);

/*!
 * \}
 */

#endif /* DISTRIBUTION_ENUMERATOR_H */
