/*!
 * \file    distribution_enumerator.cpp
 * \ingroup two_dimensional_distribution_enumerator
 *
 * \brief   The definition of functions for enumerating the slices of
 *          two-dimensional probability distributions.
 */

#include "distribution_enumerator.h"

#include "parameters.h"
#include "errors.h"
#include "math.h"
#include "common.h"

#include <stdint.h>
#include <stdlib.h>

/*!
 * \brief Compares two distribution slice coordinates with respect to their 
 *        distance from (abs(alpha_d), abs(alpha_r)) = (m, m), for the 
 *        purpose of enabling coordinates to be sorted with qsort().
 * 
 * \param[in] a   A pointer to a pointer to the first coordinate.
 * \param[in] b   A pointer to a pointer to the second coordinate.
 * 
 * \return Returns -1 if the the distance for the first coordinate is greater
 *         than that of the second coordinate, -1 if the inverse is true and 
 *         0 if the distances are equal.
 */
static int distribution_enumerator_sort_coordinates_cmp(
  const void * const a,
  const void * const b)
{
  const Distribution_Coordinate * const coordinate_a = 
    *((const Distribution_Coordinate * const * const)a);

  const Distribution_Coordinate * const coordinate_b = 
    *((const Distribution_Coordinate * const * const)b);
  
  if ((coordinate_a->distance) > (coordinate_b->distance)) {
    return 1;
  } else if ((coordinate_a->distance) < (coordinate_b->distance)) {
    return -1;
  } else {
    return 0;
  }
}

void distribution_enumerator_init(
  Distribution_Enumerator * const enumerator,
  const Parameters * const parameters,
  const bool mirrored)
{
  /* Count the number of slices required. */
  uint32_t count;

  if (mirrored) {
    count = 2;
  } else {
    count = 4;
  }

  count *= (parameters->max_alpha_d - parameters->min_alpha_d + 1);
  count *= (parameters->max_alpha_r - parameters->min_alpha_r + 1);

  enumerator->mirrored = mirrored;
  enumerator->count = 0;
  enumerator->offset = 0;

  /* Allocate memory for the list of slice coordinates. */
  enumerator->coordinates =
    (Distribution_Coordinate **)malloc(
        count * sizeof(Distribution_Coordinate *));
  if (NULL == enumerator->coordinates) {
    critical("distribution_enumerator_init(): Failed to allocate memory.");
  }

  /* Initialize all pointers to NULL. */
  for (uint32_t i = 0; i < count; i++) {
    enumerator->coordinates[i] = NULL;
  }

  for (uint32_t i = 0; i < count; i++) {
    /* Compute the coordinates. */
    uint32_t size_d = parameters->max_alpha_d + 1 - parameters->min_alpha_d;

    uint32_t min_alpha_d;
    uint32_t min_alpha_r;

    if (mirrored) {
      min_alpha_d = (int32_t)(parameters->min_alpha_d + ((i >> 1) % size_d));
      min_alpha_r = (int32_t)(parameters->min_alpha_r + ((i >> 1) / size_d));
    } else {
      min_alpha_d = (int32_t)(parameters->min_alpha_d + ((i >> 2) % size_d));
      min_alpha_r = (int32_t)(parameters->min_alpha_r + ((i >> 2) / size_d));
    }

    /* Sanity checks. */
    if (min_alpha_d >= (parameters->m + parameters->l - 1)) {
      critical("distribution_enumerator_init(): "
        "The alpha_d coordinate is out of range.");
    }

    if (min_alpha_r >= (parameters->m + parameters->l - 1)) {
      critical("distribution_enumerator_init(): "
        "The alpha_r coordinate is out of range.");
    }

    /* The coordinate is in range. Allocate memory for the coordinate. */
    Distribution_Coordinate * const coordinate =
      (Distribution_Coordinate *)malloc(sizeof(Distribution_Coordinate));
    if (NULL == coordinate) {
      critical("distribution_enumerator_init(): Failed to allocate memory.");
    }

    /* Store the coordinates. */
    coordinate->min_log_alpha_d = (int32_t)min_alpha_d;
    coordinate->min_log_alpha_r = (int32_t)min_alpha_r;

    /* Use the two lower bits in the index to select the sign. */
    if ((i & 1) == 1) {
      coordinate->min_log_alpha_d *= -1;
    }

    if ((!mirrored) && ((i & 2) == 2)) {
      coordinate->min_log_alpha_r *= -1;
    }

    /* Compute and store the distance. For more information on the distance, 
     * see documentation for Distribution_Coordinate::distance. */
    uint32_t tmp;

    tmp = abs_i((int32_t)abs_i(coordinate->min_log_alpha_d) - 
      (int32_t)parameters->m);
    coordinate->distance  = tmp * tmp;
   
    tmp = abs_i((int32_t)abs_i(coordinate->min_log_alpha_r) - 
      (int32_t)parameters->m);
    coordinate->distance += tmp * tmp;

    /* Zeroize the norms of elements on the diagonal to process these first. */
    if (abs_i(coordinate->min_log_alpha_d - coordinate->min_log_alpha_r) <= 1) {
        if (sgn_i(coordinate->min_log_alpha_d) == 
          sgn_i(coordinate->min_log_alpha_r))
        {
          if (abs_i(coordinate->min_log_alpha_d) > parameters->m) {
            coordinate->distance = 0;
          }
        }
    }

    /* Insert the coordinate into the enumerator and increment the count. */
    enumerator->coordinates[enumerator->count++] = coordinate;
  }

  /* Sort the coordinates. */
  qsort(
    enumerator->coordinates, 
    enumerator->count,
    sizeof(Distribution_Coordinate *),
    distribution_enumerator_sort_coordinates_cmp);
}

void distribution_enumerator_clear(
  Distribution_Enumerator * const enumerator)
{
  for (uint32_t i = 0; i < enumerator->count; i++) {
    free(enumerator->coordinates[i]);
    enumerator->coordinates[i] = NULL;
  }

  free(enumerator->coordinates);
}

bool distribution_enumerator_next(
  int32_t * const min_log_alpha_d,
  int32_t * const min_log_alpha_r,
  Distribution_Enumerator * const enumerator)
{
  if (enumerator->offset >= enumerator->count) {
    return FALSE;
  }

  (*min_log_alpha_d) = 
    enumerator->coordinates[enumerator->offset]->min_log_alpha_d;
  (*min_log_alpha_r) =
    enumerator->coordinates[enumerator->offset]->min_log_alpha_r;

  enumerator->offset++;

  return TRUE;
}

uint32_t distribution_enumerator_count(
  const Distribution_Enumerator * const enumerator)
{
  return enumerator->count;
}
