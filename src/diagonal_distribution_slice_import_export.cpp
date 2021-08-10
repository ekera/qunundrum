/*!
 * \file    diagonal_distribution_slice_import_export.cpp
 * \ingroup diagonal_distribution_slice
 *
 * \brief   The definition of functions for importing and exporting slices in
 *          diagonal probability distributions.
 */

#include "diagonal_distribution_slice.h"

#include "errors.h"
#include "math.h"
#include "common.h"

#include <mpfr.h>
#include <gmp.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void diagonal_distribution_slice_import_common(
  Diagonal_Distribution_Slice * const slice,
  FILE * const file)
{
  if (1 != fscanf(file, "%d\n", &(slice->min_log_alpha_r))) {
    critical("diagonal_distribution_slice_import_common(): "
      "Failed to import min_log_alpha_r.");
  }

  if (1 != fscanf(file, "%x\n", &(slice->flags))) {
    critical("diagonal_distribution_slice_import_common(): "
      "Failed to import flags.");
  }

  slice->total_probability = 0;

  for (uint32_t i = 0; i < slice->dimension; i++) {
    if (1 != fscanf(file, "%Lg\n", &(slice->norm_vector[i]))) {
      critical("diagonal_distribution_slice_import_common(): "
        "Failed to import an element in the norm vector.");
    }

    slice->total_probability += slice->norm_vector[i];
  }

  if (1 != fscanf(file, "%Lg\n", &(slice->total_error))) {
    critical("diagonal_distribution_slice_import_common(): "
      "Failed to import the total error.");
  }
}

void diagonal_distribution_slice_import(
  Diagonal_Distribution_Slice * const slice,
  FILE * const file)
{
  uint32_t dimension;

  if (1 != fscanf(file, "%u\n", &dimension)) {
    critical("diagonal_distribution_slice_import(): "
      "Failed to import the dimension.");
  }

  if (dimension != slice->dimension) {
    diagonal_distribution_slice_clear(slice);
    diagonal_distribution_slice_init(slice, dimension);
  }

  diagonal_distribution_slice_import_common(slice, file);
}

void diagonal_distribution_slice_init_import(
  Diagonal_Distribution_Slice * const slice,
  FILE * const file)
{
  uint32_t dimension;

  if (1 != fscanf(file, "%u\n", &dimension)) {
    critical("diagonal_distribution_slice_init_import(): "
      "Failed to import the dimension.");
  }

  diagonal_distribution_slice_init(slice, dimension);

  diagonal_distribution_slice_import_common(slice, file);
}

void diagonal_distribution_slice_export(
  const Diagonal_Distribution_Slice * const slice,
  FILE * const file)
{
  fprintf(file, "%u\n", slice->dimension);
  fprintf(file, "%d\n", slice->min_log_alpha_r);
  fprintf(file, "%.8x\n", slice->flags);

  for (uint32_t i = 0; i < slice->dimension; i++) {
    fprintf(file, "%.24Lg\n", slice->norm_vector[i]);
  }

  fprintf(file, "%.24Lg\n", slice->total_error);
}
