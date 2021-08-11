/*!
 * \file    distribution_slice_import_export.cpp
 * \ingroup two_dimensional_distribution_slice
 *
 * \brief   The definition of functions for importing and exporting slices in
 *          two-dimensional probability distributions.
 */

#include "distribution_slice.h"

#include "common.h"
#include "errors.h"
#include "math.h"

#include <stdint.h>
#include <stdio.h>

static void distribution_slice_import_common(
  Distribution_Slice * const slice,
  FILE * const file)
{
  if (1 != fscanf(file, "%d\n", &(slice->min_log_alpha_d))) {
    critical("distribution_slice_import_common(): "
      "Failed to import min_log_alpha_d.");
  }

  if (1 != fscanf(file, "%d\n", &(slice->min_log_alpha_r))) {
    critical("distribution_slice_import_common(): "
      "Failed to import min_log_alpha_r.");
  }

  if (1 != fscanf(file, "%x\n", &(slice->flags))) {
    critical("distribution_slice_import_common(): "
      "Failed to import flags.");
  }

  slice->total_probability = 0;

  for (uint32_t i = 0; i < (slice->dimension * slice->dimension); i++) {
    if (1 != fscanf(file, "%Lg\n", &(slice->norm_matrix[i]))) {
      critical("distribution_slice_import_common(): "
        "Failed to import an element in the norm matrix.");
    }

    slice->total_probability += slice->norm_matrix[i];
  }

  if (1 != fscanf(file, "%Lg\n", &(slice->total_error))) {
    critical("distribution_slice_import_common(): "
      "Failed to import the total error.");
  }
}

void distribution_slice_import(
  Distribution_Slice * const slice,
  FILE * const file)
{
  uint32_t dimension;

  if (1 != fscanf(file, "%u\n", &dimension)) {
    critical("distribution_slice_import(): "
      "Failed to import the dimension.");
  }

  if (dimension != slice->dimension) {
    distribution_slice_clear(slice);
    distribution_slice_init(slice, dimension);
  }

  distribution_slice_import_common(slice, file);
}

void distribution_slice_init_import(
  Distribution_Slice * const slice,
  FILE * const file)
{
  uint32_t dimension;

  if (1 != fscanf(file, "%u\n", &dimension)) {
    critical("distribution_slice_init_import(): "
      "Failed to import the dimension.");
  }

  distribution_slice_init(slice, dimension);

  distribution_slice_import_common(slice, file);
}

void distribution_slice_export(
  const Distribution_Slice * const slice,
  FILE * const file)
{
  fprintf(file, "%u\n", slice->dimension);
  fprintf(file, "%d\n", slice->min_log_alpha_d);
  fprintf(file, "%d\n", slice->min_log_alpha_r);
  fprintf(file, "%.8x\n", slice->flags);

  for (uint32_t i = 0; i < (slice->dimension * slice->dimension); i++) {
    fprintf(file, "%.24Lg\n", slice->norm_matrix[i]);
  }

  fprintf(file, "%.24Lg\n", slice->total_error);
}
