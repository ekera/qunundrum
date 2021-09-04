/*!
 * \file    diagonal_distribution_enumerator.cpp
 * \ingroup diagonal_distribution_enumerator
 *
 * \brief   The definition of functions for enumerating the slices of
 *          diagonal probability distributions.
 */

#include "diagonal_distribution_enumerator.h"

#include "common.h"
#include "diagonal_parameters.h"
#include "errors.h"
#include "math.h"

#include <stdint.h>
#include <stdlib.h>

/*!
 * \brief Compares two diagonal distribution slice coordinates with respect to
 *        their distance, for the purpose of enabling coordinates to be sorted
 *        with qsort().
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
    *((const Diagonal_Distribution_Coordinate * const *)a);

  const Diagonal_Distribution_Coordinate * const coordinate_b =
    *((const Diagonal_Distribution_Coordinate * const *)b);

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
  const Diagonal_Parameters * const parameters,
  const bool mirrored)
{
  /* Setup the region in alpha. */
  int32_t min_alpha;

  if (parameters->t > parameters->m) {
    min_alpha = 0;
  } else {
    min_alpha = parameters->m - parameters->t;
  }

  int32_t max_alpha;

  if (parameters->t >= parameters->sigma) {
    /* The maximum alpha is 2^(m + sigma - 1) - 1. */
    max_alpha = parameters->m + parameters->sigma - 2;
  } else {
    max_alpha = parameters->m + parameters->t - 1;
  }

  /* Count the number of slices required. */
  uint32_t count = max_alpha - min_alpha + 1;

  if (!mirrored) {
    count *= 2;
  }

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

  for (uint32_t i = 0; i < count; i++) {
    /* Allocate memory for the slice coordinate. */
    enumerator->coordinates[i] =
      (Diagonal_Distribution_Coordinate *)malloc(
        sizeof(Diagonal_Distribution_Coordinate));
    if (NULL == enumerator->coordinates[i]) {
      critical("diagonal_distribution_enumerator_init(): "
        "Failed to allocate memory.");
    }

    /* Compute and store the coordinates. */
    if (mirrored) {
      enumerator->coordinates[i]->min_log_alpha_r = min_alpha + i;
    } else {
      enumerator->coordinates[i]->min_log_alpha_r = min_alpha + (i >> 1);

      /* Use the least significant bit in the index to select the sign. */
      if ((i & 1) == 1) {
        enumerator->coordinates[i]->min_log_alpha_r *= -1;
      }
    }

    /* Compute and store the distance. For more information on the distance,
     * see documentation for Diagonal_Distribution_Coordinate::distance. */
    enumerator->coordinates[i]->distance =
      abs_i((int32_t)abs_i(enumerator->coordinates[i]->min_log_alpha_r) -
        (int32_t)parameters->m);
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
  Diagonal_Distribution_Enumerator * const enumerator)
{
  if (enumerator->offset >= enumerator->count) {
    return FALSE;
  }

  (*min_log_alpha_r) =
    enumerator->coordinates[enumerator->offset]->min_log_alpha_r;

  enumerator->offset++;

  return TRUE;
}

uint32_t diagonal_distribution_enumerator_count(
  const Diagonal_Distribution_Enumerator * const enumerator)
{
  return enumerator->count;
}
