/*!
 * \file    linear_distribution.cpp
 * \ingroup linear_distribution
 *
 * \brief   The definition of functions for manipulating linear probability
 *          distributions.
 */

#include "linear_distribution.h"

#include "common.h"
#include "distribution.h"
#include "diagonal_distribution.h"
#include "errors.h"
#include "linear_distribution_slice.h"
#include "math.h"
#include "parameters.h"
#include "random.h"
#include "sample.h"

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Linear_Distribution * linear_distribution_alloc()
{
  Linear_Distribution * distribution =
    (Linear_Distribution *)malloc(sizeof(Linear_Distribution));
  if (NULL == distribution) {
    critical("linear_distribution_alloc(): Failed to allocate memory.");
  }

  return distribution;
}

void linear_distribution_dealloc(
  Linear_Distribution ** const distribution)
{
  free(*distribution);
  (*distribution) = NULL;
}

void linear_distribution_init(
  Linear_Distribution * const distribution,
  const Parameters * const parameters,
  const uint32_t flags,
  const uint32_t capacity)
{
  /* Zeroize the distribution. */
  memset(distribution, 0, sizeof(Linear_Distribution));

  /* Allocate memory for the list of slices. */
  distribution->slices =
    (Linear_Distribution_Slice **)malloc(
      sizeof(Linear_Distribution_Slice *) * capacity);
  if (NULL == distribution->slices) {
    critical("linear_distribution_init(): Failed to allocate memory.");
  }

  for (uint32_t i = 0; i < capacity; i++) {
    distribution->slices[i] = NULL;
  }

  /* Copy the parameters to the distribution. */
  parameters_init(&(distribution->parameters));
  parameters_copy(&(distribution->parameters), parameters);

  /* Store other parameters in the distribution. */
  distribution->precision = PRECISION;
  distribution->capacity = capacity;
  distribution->count = 0;
  distribution->total_probability = 0;
  distribution->total_error = 0;
  distribution->flags = flags;
}

void linear_distribution_init_copy_scale(
  Linear_Distribution * const dst_distribution,
  const Linear_Distribution * const src_distribution,
  const uint32_t dimension)
{
  linear_distribution_init(
    dst_distribution,
    &(src_distribution->parameters),
    src_distribution->flags,
    src_distribution->capacity);

  for (uint32_t i = 0; i < src_distribution->count; i++) {
    Linear_Distribution_Slice * const slice = linear_distribution_slice_alloc();

    linear_distribution_slice_init(slice, dimension);
    linear_distribution_slice_copy_scale(slice, src_distribution->slices[i]);
    linear_distribution_insert_slice(dst_distribution, slice);
  }
}

void linear_distribution_init_import(
  Linear_Distribution * const distribution,
  FILE * const file)
{
  /* Import the precision. */
  uint32_t precision;

  if (1 != fscanf(file, "%u\n", &precision)) {
    critical("linear_distribution_init_import(): "
      "Failed to import the precision.");
  }

  /* Import the flags. */
  uint32_t flags;

  if (1 != fscanf(file, "%x\n", &flags)) {
    critical("linear_distribution_init_import(): Failed to import flags.");
  }

  /* Import the parameters and use them to initialize the distribution. */
  Parameters parameters;

  parameters_init(&parameters);
  parameters_import(&parameters, file);

  /* Read the slice count. */
  uint32_t count;

  if (1 != fscanf(file, "%u\n", &count)) {
    critical("linear_distribution_init_import(): "
      "Failed to read the slice count.");
  }

  /* Initialize the distribution with the minimum required capacity. */
  linear_distribution_init(distribution, &parameters, flags, count);

  /* Import the slices. */
  for (uint32_t i = 0; i < count; i++) {
    Linear_Distribution_Slice * slice = linear_distribution_slice_alloc();

    linear_distribution_slice_init_import(slice, file);

    /* The slice pointer is copied the distribution. The slice is deallocated
     * when the distribution is deallocated. */
    linear_distribution_insert_slice(distribution, slice);
  }

  /* Clear memory. */
  parameters_clear(&parameters);
}

void linear_distribution_init_collapse_d(
  Linear_Distribution * const dst,
  const Distribution * const src)
{
  /* Initialize the distribution. */
  uint32_t flags;

  flags  = LINEAR_DISTRIBUTION_FLAG_D;
  flags |= LINEAR_DISTRIBUTION_FLAG_COLLAPSED;

  linear_distribution_init(dst, &(src->parameters), flags, src->count);

  /* If there are no slices in the distribution, we need not do anything. */
  if (0 == src->count) {
    return;
  }

  /* Find the maximum dimension of slices in the distribution. */
  uint32_t max_dimension = src->slices[0]->dimension;

  for (uint32_t i = 1; i < src->count; i++) {
    if (src->slices[i]->dimension > max_dimension) {
      max_dimension = src->slices[i]->dimension;
    }
  }

  /* Assert that all slice dimensions divides the maximum dimension. */
  for (uint32_t i = 0; i < src->count; i++) {
    if ((max_dimension % src->slices[i]->dimension) != 0) {
      critical("linear_distribution_init_collapse_d(): All slices in the "
        "source distribution must be of dimension that divides the maximum "
        "slice dimension.");
    }
  }

  /* Iterate over all slices in the source distribution. */
  for (uint32_t i = 0; i < src->count; i++) {
    /* Check if there is a slice in dst for this alpha_d-range. */
    uint32_t j;

    for (j = 0; j < dst->count; j++) {
      if (src->slices[i]->min_log_alpha_d == dst->slices[j]->min_log_alpha) {
        break;
      }
    }

    Linear_Distribution_Slice * slice = NULL;

    if (j == dst->count) {
      /* No slice was found. Create and insert a slice. */
      slice = linear_distribution_slice_alloc();
      linear_distribution_slice_init(slice, max_dimension);

      slice->min_log_alpha = src->slices[i]->min_log_alpha_d;

      linear_distribution_insert_slice(dst, slice);
    } else {
      /* An existing slice was found. */
      slice = dst->slices[j];
    }

    /* Add total probability and error to the slice. */
    slice->total_probability += src->slices[i]->total_probability;
    slice->total_error += src->slices[i]->total_error;

    /* Add the norm matrix to the slice. */
    const uint32_t dimension = src->slices[i]->dimension;
    const uint32_t divisor = max_dimension / dimension;

    for (uint32_t x = 0; x < dimension; x++) {
      for (uint32_t y = 0; y < dimension; y++) {
        const long double probability =
          src->slices[i]->norm_matrix[x + y * dimension];

        for (uint32_t j = 0; j < divisor; j++) {
          slice->norm_vector[divisor * x + j] +=
            probability / (long double)divisor;
        }

      }
    }
  }

  dst->total_probability = src->total_probability;
  dst->total_error = src->total_error;
}

void linear_distribution_init_collapse_r(
  Linear_Distribution * const dst,
  const Distribution * const src)
{
  /* Initialize the distribution. */
  uint32_t flags;

  flags  = LINEAR_DISTRIBUTION_FLAG_R;
  flags |= LINEAR_DISTRIBUTION_FLAG_COLLAPSED;

  linear_distribution_init(dst, &(src->parameters), flags, src->count);

  /* If there are no slices in the distribution, we need not do anything. */
  if (0 == src->count) {
    return;
  }

  /* Find the maximum dimension of slices in the distribution. */
  uint32_t max_dimension = src->slices[0]->dimension;

  for (uint32_t i = 1; i < src->count; i++) {
    if (src->slices[i]->dimension > max_dimension) {
      max_dimension = src->slices[i]->dimension;
    }
  }

  /* Assert that all slice dimensions divides the maximum dimension. */
  for (uint32_t i = 0; i < src->count; i++) {
    if ((max_dimension % src->slices[i]->dimension) != 0) {
      critical("linear_distribution_init_collapse_r(): All slices in the "
        "source distribution must be of dimension that divides the maximum "
        "slice dimension.");
    }
  }

  /* Iterate over all slices in the source distribution. */
  for (uint32_t i = 0; i < src->count; i++) {
    /* Check if there is a slice in dst for this alpha_r-range. */
    uint32_t j;

    for (j = 0; j < dst->count; j++) {
      if (src->slices[i]->min_log_alpha_r == dst->slices[j]->min_log_alpha) {
        break;
      }
    }

    Linear_Distribution_Slice * slice = NULL;

    if (j == dst->count) {
      /* No slice was found. Create and insert a slice. */
      slice = linear_distribution_slice_alloc();
      linear_distribution_slice_init(slice, max_dimension);

      slice->min_log_alpha = src->slices[i]->min_log_alpha_r;

      linear_distribution_insert_slice(dst, slice);
    } else {
      /* An existing slice was found. */
      slice = dst->slices[j];
    }

    /* Add total probability and error to the slice. */
    slice->total_probability += src->slices[i]->total_probability;
    slice->total_error += src->slices[i]->total_error;

    /* Add the norm matrix to the slice. */
    const uint32_t dimension = src->slices[i]->dimension;
    const uint32_t divisor = max_dimension / dimension;

    for (uint32_t x = 0; x < dimension; x++) {
      for (uint32_t y = 0; y < dimension; y++) {
        const long double probability =
          src->slices[i]->norm_matrix[x + y * dimension];

        for (uint32_t j = 0; j < divisor; j++) {
          slice->norm_vector[divisor * y + j] +=
            probability / (long double)divisor;
        }

      }
    }
  }

  dst->total_probability = src->total_probability;
  dst->total_error = src->total_error;
}

void linear_distribution_init_collapse_diagonal(
  Linear_Distribution * const dst,
  const Diagonal_Distribution * const src,
  const uint32_t eta_bound)
{
  /* Initialize the distribution. */
  uint32_t flags;

  flags  = LINEAR_DISTRIBUTION_FLAG_R;
  flags |= LINEAR_DISTRIBUTION_FLAG_COLLAPSED;

  Parameters parameters;
  parameters_init(&parameters);
  parameters.m = src->parameters.m;
  parameters.l = src->parameters.sigma;
  parameters.s = 0;

  mpz_set(parameters.d, src->parameters.d);
  mpz_set(parameters.r, src->parameters.r);

  parameters.min_alpha_d = 0;
  parameters.max_alpha_d = 0;

  parameters.min_alpha_r = src->parameters.min_alpha_r;
  parameters.max_alpha_r = src->parameters.max_alpha_r;

  parameters.t = src->parameters.t;

  linear_distribution_init(dst, &parameters, flags, src->count);

  /* Clear memory. */
  parameters_clear(&parameters);

  /* If there are no slices in the distribution, we need not do anything. */
  if (0 == src->count) {
    return;
  }

  /* Assert that all slice in the source dimension have the same dimension. */
  uint32_t dimension = src->slices[0]->dimension;

  for (uint32_t i = 1; i < src->count; i++) {
    if (src->slices[i]->dimension != dimension) {
      critical("linear_distribution_init_collapse_diagonal(): All slices in "
        "the source distribution must be of the same dimension.");
    }
  }

  /* Iterate over all slices in the source distribution. */
  for (uint32_t i = 0; i < src->count; i++) {
    /* Check if there is a slice in dst for this alpha_r-range. */
    uint32_t j;

    if (abs_i(src->slices[i]->eta) > eta_bound) {
      continue;
    }

    for (j = 0; j < dst->count; j++) {
      if (src->slices[i]->min_log_alpha_r == dst->slices[j]->min_log_alpha) {
        break;
      }
    }

    Linear_Distribution_Slice * slice = NULL;

    if (j == dst->count) {
      /* No slice was found. Create and insert a slice. */
      slice = linear_distribution_slice_alloc();
      linear_distribution_slice_init(slice, dimension);

      slice->min_log_alpha = src->slices[i]->min_log_alpha_r;

      linear_distribution_insert_slice(dst, slice);
    } else {
      /* An existing slice was found. */
      slice = dst->slices[j];
    }

    /* Add total probability and error to the slice. */
    slice->total_probability += src->slices[i]->total_probability;
    slice->total_error += src->slices[i]->total_error;

    /* Add the norm vector to the slice. */
    for (uint32_t x = 0; x < dimension; x++) {
      slice->norm_vector[x] += src->slices[i]->norm_vector[x];
    }
  }

  dst->total_probability = src->total_probability;
  dst->total_error = src->total_error;
}

void linear_distribution_clear(
  Linear_Distribution * const distribution)
{
  for (uint32_t i = 0; i < distribution->count; i++) {
    linear_distribution_slice_clear(distribution->slices[i]);
    linear_distribution_slice_dealloc(&(distribution->slices[i]));
  }

  free(distribution->slices);
  distribution->slices = NULL;

  parameters_clear(&(distribution->parameters));

  memset(distribution, 0, sizeof(Linear_Distribution));
}

void linear_distribution_insert_slice(
  Linear_Distribution * const distribution,
  Linear_Distribution_Slice * const slice)
{
  /* Ensure that there is room for the slice. */
  if (distribution->count >= distribution->capacity) {
    critical("linear_distribution_insert_slice(): "
      "Inserting the slice would exceed the distribution capacity.");
  }

  /* Insert the slice. */
  distribution->slices[distribution->count++] = slice;

  /* Add the total error and probability in the slice to the overall total
   * probability and error in the distribution. */
  distribution->total_probability += slice->total_probability;
  distribution->total_error += slice->total_error;
}

/*!
 * \brief Compares two linear slices with respect to the total probability
 *        they capture, for the purpose of enabling slices to be sorted in
 *        descending order with qsort().
 *
 * \param[in] a   A pointer to a pointer to the first slice.
 * \param[in] b   A pointer to a pointer to the second slice.
 *
 * \return Returns -1 if the total probability captured by the first slice is
 *         greater than that captured by the second slice, 1 if the inverse is
 *         true, and 0 if the probabilities are equal.
 */
static int linear_distribution_sort_slices_cmp(
  const void * const a,
  const void * const b)
{
  const Linear_Distribution_Slice * const slice_a =
    *((const Linear_Distribution_Slice * const *)a);

  const Linear_Distribution_Slice * const slice_b =
    *((const Linear_Distribution_Slice * const *)b);

  if ((slice_a->total_probability) > (slice_b->total_probability)) {
    return -1;
  } else if ((slice_a->total_probability) < (slice_b->total_probability)) {
    return 1;
  } else {
    return 0;
  }
}

void linear_distribution_sort_slices(
  Linear_Distribution * const distribution)
{
  if (0 == distribution->count) {
    return;
  }

  qsort(
    distribution->slices,
    distribution->count,
    sizeof(Linear_Distribution_Slice *),
    linear_distribution_sort_slices_cmp);
}

void linear_distribution_export(
  const Linear_Distribution * const distribution,
  FILE * const file)
{
  /* Export the precision. */
  fprintf(file, "%u\n", distribution->precision);

  /* Export the flags. */
  fprintf(file, "%.8x\n", distribution->flags);

  /* Export the parameters. */
  parameters_export(&(distribution->parameters), file);

  /* Export the number of slices. */
  uint32_t count = distribution->count;
  fprintf(file, "%u\n", count);

  /* Export the slices. */
  for (uint32_t i = 0; i < count; i++) {
    linear_distribution_slice_export(distribution->slices[i], file);
  }
}

const Linear_Distribution_Slice * linear_distribution_sample_slice(
  const Linear_Distribution * const distribution,
  Random_State * const random_state)
{
  /* Select a pivot uniformly at random. */
  long double pivot = random_generate_pivot_inclusive(random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("linear_distribution_sample_slice(): "
    "Debug: Sampled pivot: %Lf\n", pivot);
  #endif

  /* Normalize if the total probability exceeds one. Of course, this should
   * never occur in practice, unless the dimension that controls the resolution
   * of the distribution is set very low. This check handles such cases. */
  if (distribution->total_probability > 1) {
    pivot *= distribution->total_probability;

    #ifdef DEBUG_TRACE_SAMPLING
    printf("linear_distribution_sample_slice(): "
      "Debug: Warning: Scaled pivot: %Lf\n", pivot);
    #endif
  }

  /* Sample a slice. */
  for (uint32_t i = 0 ; i < distribution->count; i++) {
    const Linear_Distribution_Slice * const slice = distribution->slices[i];

    pivot -= slice->total_probability;

    #ifdef DEBUG_TRACE_SAMPLING
    printf("linear_distribution_sample_slice(): "
      "Debug: Decremented pivot to %Lf for slice: %u\n", pivot, i);
    #endif

    if (pivot <= 0) {
      #ifdef DEBUG_TRACE_SAMPLING
      printf("linear_distribution_sample_slice(): "
        "Debug: Selected slice: %u (%d)\n", i, slice->min_log_alpha);
      #endif

      return slice;
    }
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("linear_distribution_sample_slice(): Debug: Out of bounds.\n");
  #endif

  /* Signal out of bounds. */
  return NULL;
}

bool linear_distribution_sample_region(
  const Linear_Distribution * const distribution,
  Random_State * const random_state,
  double * const min_log_alpha,
  double * const max_log_alpha)
{
  const Linear_Distribution_Slice * const slice =
    linear_distribution_sample_slice(distribution, random_state);

  if (NULL == slice) {
    #ifdef DEBUG_TRACE_SAMPLING
    printf("linear_distribution_sample_region(): "
      "Debug: Sampled slice via linear_distribution_sample_slice(): "
        "Out of bounds.\n");
    #endif

    (*min_log_alpha) = 0;
    (*max_log_alpha) = 0;

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("linear_distribution_sample_region(): "
    "Debug: Sampled slice via linear_distribution_sample_slice(): %d\n",
      slice->min_log_alpha);
  #endif

  linear_distribution_slice_sample_region(
    slice,
    random_state,
    min_log_alpha,
    max_log_alpha);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("linear_distribution_sample_region(): "
    "Debug: Sampled region via linear_distribution_slice_sample_region(): "
      "%f %f\n", *min_log_alpha, *max_log_alpha);
  #endif

  /* Signal success. */
  return TRUE;
}

bool linear_distribution_sample_approximate_alpha(
  const Linear_Distribution * const distribution,
  Random_State * const random_state,
  mpfr_t alpha)
{
  /* Sample a region from the distribution. */
  double min_log_alpha;
  double max_log_alpha;

  bool result;

  result = linear_distribution_sample_region(
    distribution,
    random_state,
    &min_log_alpha,
    &max_log_alpha);
  if (FALSE == result) {
    #ifdef DEBUG_TRACE_SAMPLING
    printf("linear_distribution_sample_approximate_alpha(): "
      "Debug: Sampled region via linear_distribution_sample_region(): "
        "Out of bounds.\n");
    #endif

    mpfr_set_inf(alpha, 0); /* positive infinity */

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("linear_distribution_sample_approximate_alpha(): "
    "Debug: Sampled region via linear_distribution_sample_region(): "
      "%f %f\n", min_log_alpha, max_log_alpha);
  #endif

  sample_approximate_alpha_from_region(
    alpha,
    min_log_alpha,
    max_log_alpha,
    random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("linear_distribution_sample_approximate_alpha(): "
    "Debug: Sampled alpha via sample_approximate_alpha_from_region(): "
      "%Zd\n", alpha);
  #endif

  /* Signal success. */
  return TRUE;
}

bool linear_distribution_sample_alpha(
  const Linear_Distribution * const distribution,
  Random_State * const random_state,
  mpz_t alpha)
{
  /* Sample a region from the distribution. */
  double min_log_alpha;
  double max_log_alpha;

  bool result;

  result = linear_distribution_sample_region(
    distribution,
    random_state,
    &min_log_alpha,
    &max_log_alpha);
  if (FALSE == result) {
    #ifdef DEBUG_TRACE_SAMPLING
    printf("linear_distribution_sample_alpha(): "
      "Debug: Sampled region via linear_distribution_sample_region(): "
        "Out of bounds.\n");
    #endif

    mpz_set_ui(alpha, 0); /* no notion of infinity exists in GMP */

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("linear_distribution_sample_alpha(): "
    "Debug: Sampled region via linear_distribution_sample_region(): "
      "%f %f\n", min_log_alpha, max_log_alpha);
  #endif

  uint32_t kappa_d_or_r;

  if ((distribution->flags & LINEAR_DISTRIBUTION_FLAG_R) != 0) {
    kappa_d_or_r = kappa(distribution->parameters.r);
  } else {
    kappa_d_or_r = kappa(distribution->parameters.d);
  }

  sample_alpha_from_region(
    alpha,
    min_log_alpha,
    max_log_alpha,
    kappa_d_or_r,
    random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("linear_distribution_sample_alpha(): "
    "Debug: Kappa: %u\n", kappa_d_or_r);
  gmp_printf("linear_distribution_sample_alpha(): "
    "Debug: Sampled alpha via sample_alpha_from_region(): %Zd\n", alpha);
  #endif

  /* Signal success. */
  return TRUE;
}
