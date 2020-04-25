/*!
 * \file    diagonal_distribution_enumerator.cpp
 * \ingroup diagonal_distribution_enumerator
 *
 * \brief   The definition of functions for enumerating the slices of
 *          diagonal probability distributions.
 */

#include "diagonal_distribution_enumerator.h"

#include "diagonal_distribution.h"
#include "parameters.h"
#include "errors.h"
#include "math.h"

#include <stdint.h>
#include <stdlib.h>

/*!
 * \brief Compares two diagonal distribution slice coordinates with respect to 
 *        their distance from (abs(alpha_r), abs(offset_alpha_d)) = (m, 0), for 
 *        the purpose of enabling coordinates to be sorted with qsort().
 * 
 * A weight of three is use for the offset in alpha_d to reflect the fact that 
 * the probability drops off much more steeply in the offset from alpha than it
 * does in alpha_r. This heuristic seemingly produces reasonable results.
 * 
 * \param[in] a   A pointer to a pointer to the first coordinate.
 * \param[in] b   A pointer to a pointer to the second coordinate.
 * 
 * \return Returns -1 if the the distance for the first coordinate is greater
 *         than that of the second coordinate, -1 if the inverse is true and 
 *         0 if the distances are equal.
 */
static int diagonal_distribution_enumerator_sort_coordinates_cmp(
  const void * const a,
  const void * const b)
{
  const Diagonal_Distribution_Coordinate * const coordinate_a = 
    *((const Diagonal_Distribution_Coordinate * const * const)a);

  const Diagonal_Distribution_Coordinate * const coordinate_b = 
    *((const Diagonal_Distribution_Coordinate * const * const)b);
  
  if ((coordinate_a->distance) > (coordinate_b->distance)) {
    return 1;
  } else if ((coordinate_a->distance) < (coordinate_b->distance)) {
    return -1;
  } else {
    return 0;
  }
}

void diagonal_distribution_enumerator_init(
  Diagonal_Distribution_Enumerator * const enumerator,
  const Parameters * const parameters,
  const uint32_t bound_delta,
  const bool mirrored)
{
  if (bound_delta > DIAGONAL_DISTRIBUTION_MAX_DELTA) {
    critical("diagonal_distribution_enumerator_init(): "
      "The bound on delta must be less than or equal to %u.", 
        DIAGONAL_DISTRIBUTION_MAX_DELTA);
  }

  /* Setup the region in alpha. */
  int32_t min_alpha;

  if (parameters->t > parameters->m) {
    min_alpha = 0;
  } else {
    min_alpha = parameters->m - parameters->t;
  }

  int32_t max_alpha;

  if (parameters->t >= parameters->l) {
    /* The maximum alpha is 2^(m+l-1) - 1 this maximum is 2^(m+l-2) as we add
     * up to one to the exponent during the search. */
    max_alpha = parameters->m + parameters->l - 2;
  } else {
    max_alpha = parameters->m + parameters->t - 1;
  }

  /* Count the number of slices required. */
  uint32_t count_alpha = max_alpha - min_alpha + 1;

  if (!mirrored) {
    count_alpha *= 2;
  }
  
  const uint32_t count_offset = 2 * bound_delta + 1;

  const uint32_t count = count_alpha * count_offset;

  enumerator->count = count;
  enumerator->offset = 0;

  /* Allocate memory for the list of slice coordinates. */
  enumerator->coordinates =
    (Diagonal_Distribution_Coordinate **)malloc(
        count * sizeof(Diagonal_Distribution_Coordinate *));
  if (NULL == enumerator->coordinates) {
    critical("diagonal_distribution_enumerator_init(): "
      "Failed to allocate memory.");
  }

  uint32_t j = 0;

  for (uint32_t offset = 0; offset < count_offset; offset++) {
    for (uint32_t i = 0; i < count_alpha; i++, j++) {
      /* Allocate memory for a slice coordinate. */
      enumerator->coordinates[j] =
        (Diagonal_Distribution_Coordinate *)malloc(
          sizeof(Diagonal_Distribution_Coordinate));
      if (NULL == enumerator->coordinates[j]) {
        critical("diagonal_distribution_enumerator_init(): "
          "Failed to allocate memory.");
      }

      /* Compute and store the coordinates. */
      if (mirrored) {
        enumerator->coordinates[j]->min_log_alpha_r = min_alpha + i;
      } else {
        enumerator->coordinates[j]->min_log_alpha_r = min_alpha + (i >> 1);

        /* Use the least significant bit in the index to select the sign. */
        if ((i & 1) == 1) {
          enumerator->coordinates[j]->min_log_alpha_r *= -1;
        }
      }

      if (0 == offset) {
        enumerator->coordinates[j]->offset_alpha_d = 0;
      } else {
        enumerator->coordinates[j]->offset_alpha_d = ((offset - 1) >> 1) + 1;

        /* Use the least significant bit in the offset to select the sign. */
        if (((offset - 1) & 1) == 1) {
          enumerator->coordinates[j]->offset_alpha_d *= -1;
        }
      }

      /* Compute and store the distance. For more information on the distance, 
       * see documentation for Diagonal_Distribution_Coordinate::distance. */
      enumerator->coordinates[j]->distance =
        abs_i((int32_t)abs_i(enumerator->coordinates[j]->min_log_alpha_r) - 
          (int32_t)parameters->m) + 
            3 * abs_i(enumerator->coordinates[j]->offset_alpha_d);
    }
  }

  /* Sort the coordinates. */
  qsort(
    enumerator->coordinates, 
    enumerator->count,
    sizeof(Diagonal_Distribution_Coordinate *),
    diagonal_distribution_enumerator_sort_coordinates_cmp);
}

void diagonal_distribution_enumerator_clear(
  Diagonal_Distribution_Enumerator * const enumerator)
{
  free(enumerator->coordinates);
}

bool diagonal_distribution_enumerator_next(
  int32_t * const min_log_alpha_r,
  int32_t * const offset_alpha_d,
  Diagonal_Distribution_Enumerator * const enumerator)
{
  if (enumerator->offset >= enumerator->count) {
    return FALSE;
  }

  (*min_log_alpha_r) =
    enumerator->coordinates[enumerator->offset]->min_log_alpha_r;
  (*offset_alpha_d) =
    enumerator->coordinates[enumerator->offset]->offset_alpha_d;

  enumerator->offset++;

  return TRUE;
}

uint32_t diagonal_distribution_enumerator_count(
  const Diagonal_Distribution_Enumerator * const enumerator)
{
  return enumerator->count;
}
