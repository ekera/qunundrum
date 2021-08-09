/*!
 * \file    linear_distribution_slice_import_export.cpp
 * \ingroup linear_distribution_slice
 *
 * \brief   The definition of functions for importing and exporting slices in
 *          linear probability distributions.
 */

#include "linear_distribution_slice.h"

#include "errors.h"
#include "math.h"
#include "common.h"

#include <mpfr.h>
#include <gmp.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void linear_distribution_slice_import_common(
  Linear_Distribution_Slice * const slice,
  FILE * const file)
{
  if (1 != fscanf(file, "%d\n", &(slice->min_log_alpha))) {
    critical("linear_distribution_slice_import(): "
      "Failed to import min_log_alpha.");
  }

  if (1 != fscanf(file, "%x\n", &(slice->flags))) {
    critical("linear_distribution_slice_import(): "
      "Failed to import flags.");
  }

  slice->total_probability = 0;

  for (uint32_t i = 0; i < slice->dimension; i++) {
    if (1 != fscanf(file, "%Lg\n", &(slice->norm_vector[i]))) {
      critical("linear_distribution_slice_import(): "
        "Failed to import vector element.");
    }

    slice->total_probability += slice->norm_vector[i];
  }

  if (1 != fscanf(file, "%Lg\n", &(slice->total_error))) {
    critical("linear_distribution_slice_import(): "
      "Failed to import total error.");
  }
}

void linear_distribution_slice_import(
  Linear_Distribution_Slice * const slice,
  FILE * const file)
{
  uint32_t dimension;

  if (1 != fscanf(file, "%u\n", &dimension)) {
    critical("linear_distribution_slice_import(): "
      "Failed to import dimension.");
  }

  if (dimension != slice->dimension) {
    linear_distribution_slice_clear(slice);
    linear_distribution_slice_init(slice, dimension);
  }

  linear_distribution_slice_import_common(slice, file);
}

void linear_distribution_slice_init_import(
  Linear_Distribution_Slice * const slice,
  FILE * const file)
{
  uint32_t dimension;

  if (1 != fscanf(file, "%u\n", &dimension)) {
    critical("linear_distribution_slice_import(): "
      "Failed to import dimension.");
  }

  linear_distribution_slice_init(slice, dimension);

  linear_distribution_slice_import_common(slice, file);
}

void linear_distribution_slice_export(
  const Linear_Distribution_Slice * const slice,
  FILE * const file)
{
  fprintf(file, "%u\n", slice->dimension);
  fprintf(file, "%d\n", slice->min_log_alpha);
  fprintf(file, "%.8x\n", slice->flags);

  for (uint32_t i = 0; i < slice->dimension; i++) {
    fprintf(file, "%.24Lg\n", slice->norm_vector[i]);
  }

  fprintf(file, "%.24Lg\n", slice->total_error);
}
