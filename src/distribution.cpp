/*!
 * \file    distribution.cpp
 * \ingroup two_dimensional_distribution
 *
 * \brief   The definition of functions for manipulating two-dimensional
 *          probability distributions.
 */

#include "distribution.h"

#include "common.h"
#include "distribution_slice.h"
#include "errors.h"
#include "lattice_sample.h"
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

Distribution * distribution_alloc()
{
  Distribution * distribution =
    (Distribution *)malloc(sizeof(Distribution));
  if (NULL == distribution) {
    critical("distribution_alloc(): Failed to allocate memory.");
  }

  return distribution;
}

void distribution_dealloc(
  Distribution ** const distribution)
{
  free(*distribution);
  (*distribution) = NULL;
}

void distribution_init(
  Distribution * const distribution,
  const Parameters * const parameters,
  const uint32_t capacity)
{
  /* Zeroize the distribution. */
  memset((void *)distribution, 0, sizeof(Distribution));

  /* Allocate memory for the list of slices. */
  distribution->slices =
    (Distribution_Slice **)malloc(sizeof(Distribution_Slice *) * capacity);
  if (NULL == distribution->slices) {
    critical("distribution_init(): Failed to allocate memory.");
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

  /* Compute the alpha lattice to prepare for future sampling requests. */
  lattice_alpha_init(
    &(distribution->lattice_alpha),
    &(distribution->parameters));
}

void distribution_init_copy_scale(
  Distribution * const dst_distribution,
  const Distribution * const src_distribution,
  const uint32_t dimension)
{
  distribution_init(
    dst_distribution,
    &(src_distribution->parameters),
    src_distribution->capacity);

  for (uint32_t i = 0; i < src_distribution->count; i++) {
    Distribution_Slice * const slice = distribution_slice_alloc();

    distribution_slice_init(slice, dimension);
    distribution_slice_copy_scale(slice, src_distribution->slices[i]);
    distribution_insert_slice(dst_distribution, slice);
  }
}

void distribution_init_import(
  Distribution * const distribution,
  FILE * const file)
{
  /* Import the precision. */
  if (1 != fscanf(file, "%u\n", &(distribution->precision))) {
    critical("distribution_init_import(): Failed to import the precision.");
  }

  /* Import the parameters and use them to initialize the distribution. */
  Parameters parameters;

  parameters_init(&parameters);
  parameters_import(&parameters, file);

  /* Read the slice count. */
  uint32_t count;

  if (1 != fscanf(file, "%u\n", &count)) {
    critical("distribution_init_import(): Failed to read the slice count.");
  }

  /* Initialize the distribution with the minimum required capacity. */
  distribution_init(distribution, &parameters, count);

  /* Import the slices. */
  for (uint32_t i = 0; i < count; i++) {
    Distribution_Slice * slice = distribution_slice_alloc();

    distribution_slice_init_import(slice, file);

    /* The slice pointer is copied the distribution. The slice is deallocated
     * when the distribution is deallocated. */
    distribution_insert_slice(distribution, slice);
  }

  /* Clear memory. */
  parameters_clear(&parameters);
}

void distribution_clear(
  Distribution * const distribution)
{
  for (uint32_t i = 0; i < distribution->count; i++) {
    distribution_slice_clear(distribution->slices[i]);
    distribution_slice_dealloc(&(distribution->slices[i]));
  }

  free(distribution->slices);
  distribution->slices = NULL;

  parameters_clear(&(distribution->parameters));

  lattice_alpha_clear(&(distribution->lattice_alpha));

  memset((void *)distribution, 0, sizeof(Distribution));
}

void distribution_insert_slice(
  Distribution * const distribution,
  Distribution_Slice * const slice)
{
  /* Ensure that there is room for the slice. */
  if (distribution->count >= distribution->capacity) {
    critical("distribution_insert_slice(): "
      "Inserting the slice would exceed the distribution capacity.");
  }

  /* Insert the slice. */
  distribution->slices[distribution->count++] = slice;

  /* Add the total error and probability in the slice to the overall total
   * probability and error in the distribution. */
  distribution->total_probability += slice->total_probability;
  distribution->total_error += slice->total_error;
}

void distribution_remove_slice(
  Distribution * const distribution,
  const uint32_t index)
{
  if (index >= distribution->count) {
    critical("distribution_remove_slice(): Slice index out of range.");
  }

  distribution_slice_clear(distribution->slices[index]);
  distribution_slice_dealloc(&(distribution->slices[index]));

  for (uint32_t i = index + 1; i < distribution->count; i++) {
    distribution->slices[i - 1] = distribution->slices[i];
  }

  distribution->count--;
}

bool distribution_replace_slice(
  Distribution * const distribution,
  Distribution_Slice * const slice)
{
  for (uint32_t i = 0; i < distribution->count; i++) {
    if ((distribution->slices[i]->min_log_alpha_d == slice->min_log_alpha_d) &&
        (distribution->slices[i]->min_log_alpha_r == slice->min_log_alpha_r))
    {
      /* Update the total probability and total error in the distribution. */
      distribution->total_probability -=
        distribution->slices[i]->total_probability;
      distribution->total_error -=
        distribution->slices[i]->total_error;

      distribution->total_probability += slice->total_probability;
      distribution->total_error += slice->total_error;

      /* Clear and deallocate the existing slice that is being replaced. */
      distribution_slice_clear(distribution->slices[i]);
      distribution_slice_dealloc(&(distribution->slices[i]));

      /* Insert the new slice in place of the existing slice. */
      distribution->slices[i] = slice;

      return TRUE;
    }
  }

  return FALSE;
}

/*!
 * \brief Compares two slices with respect to the total probability they
 *        capture, for the purpose of enabling slices to be sorted in descending
 *        order with qsort().
 *
 * \param[in] a   A pointer to a pointer to the first slice.
 * \param[in] b   A pointer to a pointer to the second slice.
 *
 * \return Returns -1 if the total probability captured by the first slice is
 *         greater than that captured by the second slice, 1 if the inverse is
 *         true, and 0 if the probabilities are equal.
 */
static int distribution_sort_slices_cmp(
  const void * const a,
  const void * const b)
{
  const Distribution_Slice * const slice_a =
    *((const Distribution_Slice * const *)a);

  const Distribution_Slice * const slice_b =
    *((const Distribution_Slice * const *)b);

  if ((slice_a->total_probability) > (slice_b->total_probability)) {
    return -1;
  } else if ((slice_a->total_probability) < (slice_b->total_probability)) {
    return 1;
  } else {
    return 0;
  }
}

void distribution_sort_slices(
  Distribution * const distribution)
{
  if (0 == distribution->count) {
    return;
  }

  qsort(
    distribution->slices,
    distribution->count,
    sizeof(Distribution_Slice *),
    distribution_sort_slices_cmp);
}

void distribution_filter_slices(
  Distribution * const distribution)
{
  /* Filter the slices in the distribution. */
  for (uint32_t i = 0; i < distribution->count; ) {
    const double total_probability = distribution->slices[i]->total_probability;
    const double total_error = distribution->slices[i]->total_error;

    if ((total_probability == 0) || (total_error / total_probability > 0.01f)) {
      /* Remove this slice. */
      distribution_remove_slice(distribution, i);
      continue;
    }

    i++; /* Process the next index. */
  }
}

bool distribution_is_filtered(
  const Distribution * const distribution)
{
  for (uint32_t i = 0; i < distribution->count; i++) {
    const double total_probability = distribution->slices[i]->total_probability;
    const double total_error = distribution->slices[i]->total_error;

    if ((total_probability == 0) || (total_error / total_probability > 0.01f)) {
      return FALSE;
    }
  }

  return TRUE;
}

void distribution_export(
  const Distribution * const distribution,
  FILE * const file)
{
  /* Export the precision. */
  fprintf(file, "%u\n", distribution->precision);

  /* Export the parameters. */
  parameters_export(&(distribution->parameters), file);

  /* Export the number of slices. */
  uint32_t count = distribution->count;
  fprintf(file, "%u\n", count);

  /* Export the slices. */
  for (uint32_t i = 0; i < count; i++) {
    distribution_slice_export(distribution->slices[i], file);
  }
}

void distribution_export_clear_dealloc(
  Distribution * const distribution,
  FILE * const file)
{
  /* Export the precision. */
  fprintf(file, "%u\n", distribution->precision);

  /* Export the parameters. */
  parameters_export(&(distribution->parameters), file);

  /* Export the number of slices. */
  uint32_t count = distribution->count;
  fprintf(file, "%u\n", count);

  /* Export the slices. */
  for (uint32_t i = 0; i < count; i++) {
    distribution_slice_export(distribution->slices[i], file);

    /* Clear and deallocate the slice. */
    distribution_slice_clear(distribution->slices[i]);
    distribution_slice_dealloc(&(distribution->slices[i]));
  }

  /* Deallocate the list of slices. */
  free(distribution->slices);
  distribution->slices = NULL;

  /* Clear the parameters. */
  parameters_clear(&(distribution->parameters));

  /* Zeroize the distribution. */
  memset((void *)distribution, 0, sizeof(Distribution));
}

const Distribution_Slice * distribution_sample_slice(
  const Distribution * const distribution,
  Random_State * const random_state)
{
  /* Select a pivot uniformly at random. */
  long double pivot = random_generate_pivot_inclusive(random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("distribution_sample_slice(): "
    "Debug: Sampled pivot: %Lf\n", pivot);
  #endif

  /* Normalize if the total probability exceeds one. Of course, this should
   * never occur in practice, unless the dimension that controls the resolution
   * of the distribution is set very low. This check handles such cases. */
  if (distribution->total_probability > 1) {
    pivot *= distribution->total_probability;

    #ifdef DEBUG_TRACE_SAMPLING
    printf("distribution_sample_slice(): "
      "Debug: Warning: Scaled pivot: %Lf\n", pivot);
    #endif
  }

  /* Sample a slice. */
  for (uint32_t i = 0 ; i < distribution->count; i++) {
    const Distribution_Slice * const slice = distribution->slices[i];

    pivot -= slice->total_probability;

    #ifdef DEBUG_TRACE_SAMPLING
    printf("distribution_sample_slice(): "
      "Debug: Decremented pivot to %Lf for slice: %u\n", pivot, i);
    #endif

    if (pivot <= 0) {
      #ifdef DEBUG_TRACE_SAMPLING
      printf("distribution_sample_slice(): "
        "Debug: Selected slice: %u (%d, %d)\n", i,
          slice->min_log_alpha_d, slice->min_log_alpha_r);
      #endif

      return slice;
    }
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("distribution_sample_slice(): Debug: Out of bounds.\n");
  #endif

  /* Signal out of bounds. */
  return NULL;
}

bool distribution_sample_region(
  const Distribution * const distribution,
  Random_State * const random_state,
  double * const min_log_alpha_d,
  double * const max_log_alpha_d,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r)
{
  /* Sample a slice. */
  const Distribution_Slice * const slice =
    distribution_sample_slice(distribution, random_state);

  if (NULL == slice) {
    #ifdef DEBUG_TRACE_SAMPLING
    printf("distribution_sample_region(): "
      "Debug: Sampled slice via distribution_sample_slice(): Out of bounds.\n");
    #endif

    (*min_log_alpha_d) = 0;
    (*max_log_alpha_d) = 0;
    (*min_log_alpha_r) = 0;
    (*max_log_alpha_r) = 0;

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("distribution_sample_region(): "
    "Debug: Sampled slice via distribution_sample_slice(): (%d, %d)\n",
      slice->min_log_alpha_d, slice->min_log_alpha_r);
  #endif

  /* Sample a region from the slice. */
  distribution_slice_sample_region(
    slice,
    random_state,
    min_log_alpha_d,
    max_log_alpha_d,
    min_log_alpha_r,
    max_log_alpha_r);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("distribution_sample_region(): "
    "Debug: Sampled region via distribution_slice_sample_region(): "
      "(%f %f, %f %f)\n", *min_log_alpha_d, *max_log_alpha_d,
        *min_log_alpha_r, *max_log_alpha_r);
  #endif

  /* Signal success. */
  return TRUE;
}

bool distribution_sample_approximate_alpha_d_r(
  const Distribution * const distribution,
  Random_State * const random_state,
  mpfr_t alpha_d,
  mpfr_t alpha_r)
{
  /* Sample a region from the distribution. */
  double min_log_alpha_d;
  double max_log_alpha_d;
  double min_log_alpha_r;
  double max_log_alpha_r;

  bool result;

  result = distribution_sample_region(
    distribution,
    random_state,
    &min_log_alpha_d,
    &max_log_alpha_d,
    &min_log_alpha_r,
    &max_log_alpha_r);
  if (FALSE == result) {
    #ifdef DEBUG_TRACE_SAMPLING
    printf("distribution_sample_approximate_alpha_d_r(): "
      "Debug: Sampled region via distribution_sample_region(): "
        "Out of bounds.\n");
    #endif

    mpfr_set_inf(alpha_d, 0); /* positive infinity */
    mpfr_set_inf(alpha_r, 0); /* positive infinity */

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("distribution_sample_approximate_alpha_d_r(): "
    "Debug: Sampled region via distribution_sample_region(): "
      "(%f %f, %f %f)\n", min_log_alpha_d, max_log_alpha_d,
        min_log_alpha_r, max_log_alpha_r);
  #endif

  /* Sample an argument pair from the region. */
  sample_approximate_alpha_from_region(
    alpha_d,
    min_log_alpha_d,
    max_log_alpha_d,
    random_state);

  sample_approximate_alpha_from_region(
    alpha_r,
    min_log_alpha_r,
    max_log_alpha_r,
    random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("distribution_sample_approximate_alpha_d_r(): "
    "Debug: Sampled (alpha_d, alpha_r) via "
      "sample_approximate_alpha_from_region():\n"
        "   alpha_d: %Zd\n   alpha_r: %Zd\n", alpha_d, alpha_r);
  #endif

  /* Signal success. */
  return TRUE;
}

bool distribution_sample_alpha_d_r(
  const Distribution * const distribution,
  Random_State * const random_state,
  mpz_t alpha_d,
  mpz_t alpha_r)
{
  /* Sample a region from the distribution. */
  double min_log_alpha_d;
  double max_log_alpha_d;
  double min_log_alpha_r;
  double max_log_alpha_r;

  bool result;

  result = distribution_sample_region(
    distribution,
    random_state,
    &min_log_alpha_d,
    &max_log_alpha_d,
    &min_log_alpha_r,
    &max_log_alpha_r);
  if (FALSE == result) {
    #ifdef DEBUG_TRACE_SAMPLING
    printf("distribution_sample_alpha_d_r(): "
      "Debug: Sampled region via distribution_sample_region(): "
        "Out of bounds.\n");
    #endif

    mpz_set_ui(alpha_d, 0); /* no notion of infinity exists in GMP */
    mpz_set_ui(alpha_r, 0); /* no notion of infinity exists in GMP */

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("distribution_sample_alpha_d_r(): "
    "Debug: Sampled region via distribution_sample_region(): "
      "(%f %f, %f %f)\n", min_log_alpha_d, max_log_alpha_d,
        min_log_alpha_r, max_log_alpha_r);
  #endif

  /* Sample an argument pair from the region. */
  const uint32_t kappa_d = kappa(distribution->parameters.d);

  sample_alpha_from_region(
    alpha_d,
    min_log_alpha_d,
    max_log_alpha_d,
    kappa_d,
    random_state);

  const uint32_t kappa_r = kappa(distribution->parameters.r);

  sample_alpha_from_region(
    alpha_r,
    min_log_alpha_r,
    max_log_alpha_r,
    kappa_r,
    random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("distribution_sample_alpha_d_r(): "
    "Debug: Kappa: (%u, %u)\n", kappa_d, kappa_r);
  gmp_printf("distribution_sample_alpha_d_r(): "
    "Debug: Sampled non-admissible alpha via sample_alpha_from_region():\n"
      "   alpha_d: %Zd\n   alpha_r: %Zd\n", alpha_d, alpha_r);
  #endif

  /* Map the argument pair to the closest admissible argument pair. */
  lattice_alpha_map(
    alpha_d,
    alpha_r,
    &(distribution->lattice_alpha),
    &(distribution->parameters));

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("distribution_sample_alpha_d_r(): "
    "Debug: Mapped to admissible alpha via lattice_alpha_map():\n"
      "   alpha_d: %Zd\n   alpha_r: %Zd\n", alpha_d, alpha_r);
  #endif

  /* Signal success. */
  return TRUE;
}

bool distribution_sample_pair_j_k(
  const Distribution * const distribution,
  Random_State * const random_state,
  mpz_t j,
  mpz_t k)
{
  /* Declare variables. */
  mpz_t alpha_d;
  mpz_init(alpha_d);

  mpz_t alpha_r;
  mpz_init(alpha_r);

  /* Sample an argument pair (alpha_d, alpha_r). */
  bool result;

  result = distribution_sample_alpha_d_r(
      distribution,
      random_state,
      alpha_d,
      alpha_r);
  if (FALSE == result) {
    #ifdef DEBUG_TRACE_SAMPLING
    gmp_printf("distribution_sample_pair_j_k(): "
      "Debug: Sampled (alpha_d, alpha_r) via "
        "distribution_sample_alpha_d_r():\n Out of bounds.\n");
    #endif

    mpz_clear(alpha_d);
    mpz_clear(alpha_r);

    mpz_set_ui(j, 0);
    mpz_set_ui(k, 0);

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("distribution_sample_pair_j_k(): "
    "Debug: Sampled (alpha_d, alpha_r) via distribution_sample_alpha_d_r():\n"
      "   alpha_d: %Zd\n   alpha_r: %Zd\n", alpha_d, alpha_r);
  #endif

  /* Sample an integer pair (j, k) from (alpha_d, alpha_r). */
  sample_j_k_from_alpha_d_r(
    j,
    k,
    alpha_d,
    alpha_r,
    &(distribution->parameters),
    random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("distribution_sample_pair_j_k(): "
    "Debug: Sampled (j, k) via sample_j_k_from_alpha_d_r():\n"
      "   j: %Zd\n   k: %Zd\n", j, k);
  #endif

  /* Clear memory. */
  mpz_clear(alpha_d);
  mpz_clear(alpha_r);

  /* Signal success. */
  return TRUE;
}
