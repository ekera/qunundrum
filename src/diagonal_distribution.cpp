/*!
 * \file    diagonal_distribution.cpp
 * \ingroup diagonal_distribution
 *
 * \brief   The definition of functions for manipulating diagonal probability
 *          distributions.
 */

#include "diagonal_distribution.h"

#include "distribution.h"
#include "diagonal_distribution_slice.h"
#include "diagonal_distribution.h"
#include "parameters.h"
#include "probability.h"
#include "random.h"
#include "errors.h"
#include "sample.h"
#include "math.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

Diagonal_Distribution * diagonal_distribution_alloc()
{
  Diagonal_Distribution * distribution =
    (Diagonal_Distribution *)malloc(sizeof(Diagonal_Distribution));
  if (NULL == distribution) {
    critical("diagonal_distribution_alloc(): Failed to allocate memory.");
  }

  return distribution;
}

void diagonal_distribution_dealloc(
  Diagonal_Distribution ** const distribution)
{
  free(*distribution);
  (*distribution) = NULL;
}

void diagonal_distribution_init(
  Diagonal_Distribution * const distribution,
  const Parameters * const parameters,
  const uint32_t flags,
  const uint32_t capacity)
{
  /* Zeroize the distribution. */
  memset(distribution, 0, sizeof(Diagonal_Distribution));

  /* Allocate memory for the list of slices. */
  distribution->slices =
    (Diagonal_Distribution_Slice **)malloc(
      sizeof(Diagonal_Distribution_Slice *) * capacity);
  if (NULL == distribution->slices) {
    critical("diagonal_distribution_init(): Failed to allocate memory.");
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

void diagonal_distribution_init_copy_scale(
  Diagonal_Distribution * const dst_distribution,
  const Diagonal_Distribution * const src_distribution,
  const uint32_t dimension)
{
  diagonal_distribution_init(
    dst_distribution,
    &(src_distribution->parameters),
    src_distribution->flags,
    src_distribution->capacity);

  for (uint32_t i = 0; i < src_distribution->count; i++) {
    Diagonal_Distribution_Slice * const slice = 
      diagonal_distribution_slice_alloc();

    diagonal_distribution_slice_init(slice, dimension);
    diagonal_distribution_slice_copy_scale(slice, src_distribution->slices[i]);
    diagonal_distribution_insert_slice(dst_distribution, slice);
  }
}

void diagonal_distribution_init_import(
  Diagonal_Distribution * const distribution,
  FILE * const file)
{
  /* Import the precision. */
  uint32_t precision;

  if (1 != fscanf(file, "%u\n", &precision)) {
    critical("diagonal_distribution_init_import(): "
      "Failed to import the distribution precision.");
  }

  /* Import the flags. */
  uint32_t flags;

  if (1 != fscanf(file, "%x\n", &flags)) {
    critical("diagonal_distribution_init_import(): Failed to import flags.");
  }

  /* Import the parameters and use them to initialize the distribution. */
  Parameters parameters;

  parameters_init(&parameters);
  parameters_import(&parameters, file);

  /* Read the slice count. */
  uint32_t count;

  if (1 != fscanf(file, "%u\n", &count)) {
    critical("diagonal_distribution_init_import(): "
      "Failed to read the distribution slice count.");
  }

  /* Initialize the distribution with the minimum required capacity. */
  diagonal_distribution_init(distribution, &parameters, flags, count);

  /* Import the slices. */
  for (uint32_t i = 0; i < count; i++) {
    Diagonal_Distribution_Slice * slice = diagonal_distribution_slice_alloc();

    diagonal_distribution_slice_init_import(slice, file);

    /* The slice pointer is copied the distribution. The slice is deallocated
     * when the distribution is deallocated. */
    diagonal_distribution_insert_slice(distribution, slice);
  }

  /* Clear memory. */
  parameters_clear(&parameters);
}

void diagonal_distribution_clear(
  Diagonal_Distribution * const distribution)
{
  for (uint32_t i = 0; i < distribution->count; i++) {
    diagonal_distribution_slice_clear(distribution->slices[i]);
    diagonal_distribution_slice_dealloc(&(distribution->slices[i]));
  }

  free(distribution->slices);
  distribution->slices = NULL;

  parameters_clear(&(distribution->parameters));

  memset(distribution, 0, sizeof(Diagonal_Distribution));
}

void diagonal_distribution_insert_slice(
  Diagonal_Distribution * const distribution,
  Diagonal_Distribution_Slice * const slice)
{
  /* Ensure that there is room for the slice. */
  if (distribution->count >= distribution->capacity) {
    critical("diagonal_distribution_insert_slice(): "
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
 * \brief Compares two diagonal slices with respect to the total probability  
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
static int diagonal_distribution_sort_slices_cmp(
  const void * const a,
  const void * const b)
{
  const Diagonal_Distribution_Slice * const slice_a = 
    *((const Diagonal_Distribution_Slice * const *)a);

  const Diagonal_Distribution_Slice * const slice_b = 
    *((const Diagonal_Distribution_Slice * const *)b);

  if ((slice_a->total_probability) > (slice_b->total_probability)) {
    return -1;
  } else if ((slice_a->total_probability) < (slice_b->total_probability)) {
    return 1;
  } else {
    return 0;
  }
}

void diagonal_distribution_sort_slices(
  Diagonal_Distribution * const distribution)
{
  if (0 == distribution->count) {
    return;
  }

  qsort(
    distribution->slices, 
    distribution->count, 
    sizeof(Diagonal_Distribution_Slice *), 
    diagonal_distribution_sort_slices_cmp);
}

void diagonal_distribution_export(
  const Diagonal_Distribution * const distribution,
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
    diagonal_distribution_slice_export(distribution->slices[i], file);
  }
}

const Diagonal_Distribution_Slice * diagonal_distribution_sample_slice(
  const Diagonal_Distribution * const distribution,
  Random_State * const random_state)
{
  /* Select a pivot uniformly at random. */
  long double pivot = random_generate_pivot_inclusive(random_state);
  
  #ifdef DEBUG_TRACE_SAMPLING
  printf("diagonal_distribution_sample_slice(): "
    "Debug: Sampled pivot: %Lf\n", pivot);
  #endif

  /* Normalize if the total probability exceeds one. Of course, this should 
   * never occur in practice, unless the dimension that controls the resolution
   * of the distribution is set very low. This check handles such cases. */
  if (distribution->total_probability > 1) {
    pivot *= distribution->total_probability;
  
    #ifdef DEBUG_TRACE_SAMPLING
    printf("diagonal_distribution_sample_slice(): "
      "Debug: Warning: Scaled pivot: %Lf\n", pivot);
    #endif
  }

  /* Sample a slice. */
  for (uint32_t i = 0 ; i < distribution->count; i++) {
    const Diagonal_Distribution_Slice * const slice = distribution->slices[i];

    pivot -= slice->total_probability;

    #ifdef DEBUG_TRACE_SAMPLING
    printf("diagonal_distribution_sample_slice(): "
      "Debug: Decremented pivot to %Lf for slice: %u\n", pivot, i);
    #endif
  
    if (pivot <= 0) {
      #ifdef DEBUG_TRACE_SAMPLING
      printf("diagonal_distribution_sample_slice(): "
        "Debug: Selected slice: %u (%d %d)\n", i, 
          slice->min_log_alpha_r, slice->offset_alpha_d);
      #endif

      return slice;
    }
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("diagonal_distribution_sample_slice(): Debug: Out of bounds.\n");
  #endif

  /* Signal out of bounds. */
  return NULL;
}

bool diagonal_distribution_sample_region(
  const Diagonal_Distribution * const distribution,
  Random_State * const random_state,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r,
  int32_t * const offset_alpha_d)
{
  const Diagonal_Distribution_Slice * const slice =
    diagonal_distribution_sample_slice(distribution, random_state);

  if (NULL == slice) {
    #ifdef DEBUG_TRACE_SAMPLING
    printf("diagonal_distribution_sample_region(): "
      "Debug: Sampled slice via diagonal_distribution_sample_slice(): "
        "Out of bounds.\n");
    #endif

    (*min_log_alpha_r) = 0;
    (*max_log_alpha_r) = 0;
    (*offset_alpha_d) = 0;

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("diagonal_distribution_sample_region(): "
    "Debug: Sampled slice via diagonal_distribution_sample_slice(): %d %d\n",
      slice->min_log_alpha_r, slice->offset_alpha_d);
  #endif

  diagonal_distribution_slice_sample_region(
    slice,
    random_state,
    min_log_alpha_r,
    max_log_alpha_r,
    offset_alpha_d);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("diagonal_distribution_sample_region(): "
    "Debug: Sampled region via diagonal_distribution_slice_sample_region(): "
      "%f %f %d\n", *min_log_alpha_r, *max_log_alpha_r, *offset_alpha_d);
  #endif

  /* Signal success. */
  return TRUE;
}

bool diagonal_distribution_sample_alpha_d_r(
  const Diagonal_Distribution * const distribution,
  Random_State * const random_state,
  mpz_t alpha_d,
  mpz_t alpha_r)
{
  /* Sample a region from the distribution. */
  double min_log_alpha_r;
  double max_log_alpha_r;
  int32_t offset_alpha_d;

  bool result;

  result = diagonal_distribution_sample_region(
    distribution,
    random_state,
    &min_log_alpha_r,
    &max_log_alpha_r,
    &offset_alpha_d);
  if (FALSE == result) {
    #ifdef DEBUG_TRACE_SAMPLING
    printf("diagonal_distribution_sample_alpha_d_r(): "
      "Debug: Sampled region via diagonal_distribution_sample_region(): "
        "Out of bounds.\n");
    #endif

    mpz_set_ui(alpha_r, 0); /* no notion of infinity exists in GMP */
    mpz_set_ui(alpha_d, 0);

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("diagonal_distribution_sample_alpha_d_r(): "
    "Debug: Sampled region via diagonal_distribution_sample_region(): "
      "%f %f %d\n", min_log_alpha_r, max_log_alpha_r, offset_alpha_d);
  #endif

  uint32_t kappa_r;

  kappa_r = kappa(distribution->parameters.r);

  sample_alpha_from_region(
    alpha_r,
    min_log_alpha_r,
    max_log_alpha_r,
    kappa_r,
    random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("diagonal_distribution_sample_alpha_d_r(): "
    "Debug: Kappa: %u\n", kappa_r);
  gmp_printf("diagonal_distribution_sample_alpha_d_r(): "
    "Debug: Sampled alpha_r via sample_alpha_from_region(): %Zd\n", alpha_r);
  #endif

  /* Compute alpha_d given alpha_r and the offset in alpha_d. */
  const uint32_t precision =
    2 * (distribution->parameters.m + distribution->parameters.l);

  mpfr_t f_alpha_d;
  mpfr_init2(f_alpha_d, precision);

  mpfr_set_z(f_alpha_d, distribution->parameters.d, MPFR_RNDN);
  mpfr_div_z(f_alpha_d, f_alpha_d, distribution->parameters.r, MPFR_RNDN);
  mpfr_mul_z(f_alpha_d, f_alpha_d, alpha_r, MPFR_RNDN);
  mpfr_round(f_alpha_d, f_alpha_d);
  mpfr_get_z(alpha_d, f_alpha_d, MPFR_RNDN);

  mpfr_clear(f_alpha_d);

  if (offset_alpha_d < 0) {
    mpz_sub_ui(alpha_d, alpha_d, (uint32_t)(-offset_alpha_d));
  } else {
    mpz_add_ui(alpha_d, alpha_d, offset_alpha_d);
  }

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("diagonal_distribution_sample_alpha_d_r(): "
    "Debug: Computed alpha_d from alpha_r (with offset %d): %Zd\n",
      offset_alpha_d, alpha_d);
  #endif

  /* Signal success. */
  return TRUE;
}

bool diagonal_distribution_sample_pair_j_k(
  const Diagonal_Distribution * const distribution,
  Random_State * const random_state,
  mpz_t j,
  mpz_t k)
{
  mpz_t alpha_d;
  mpz_init(alpha_d);

  mpz_t alpha_r;
  mpz_init(alpha_r);

  bool result;

  result = diagonal_distribution_sample_alpha_d_r(
              distribution,
              random_state,
              alpha_d,
              alpha_r);
  if (FALSE == result) {
    #ifdef DEBUG_TRACE_SAMPLING
    gmp_printf("diagonal_distribution_sample_pair_j_k(): "
      "Debug: Sampled (alpha_d, alpha_r) via "
        "diagonal_distribution_sample_alpha_d_r():\n Out of bounds.\n");
    #endif

    mpz_clear(alpha_d);
    mpz_clear(alpha_r);

    mpz_set_ui(j, 0);
    mpz_set_ui(k, 0);

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("diagonal_distribution_sample_pair_j_k(): "
    "Debug: Sampled (alpha_d, alpha_r) via "
      "diagonal_distribution_sample_alpha_d_r():\n"
        "   alpha_d: %Zd\n   alpha_r: %Zd\n", alpha_d, alpha_r);
  #endif

  /* Sample an integer pair (j, k) from (alpha_d, alpha_r). */
  sample_j_k_from_diagonal_alpha_d_r(
    j,
    k,
    alpha_d,
    alpha_r,
    &(distribution->parameters),
    random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("diagonal_distribution_sample_pair_j_k(): "
    "Debug: Sampled (j, k) via sample_j_k_from_diagonal_alpha_d_r():\n"
      "   j: %Zd\n   k: %Zd\n", j, k);
  #endif

  /* Clear memory. */
  mpz_clear(alpha_d);
  mpz_clear(alpha_r);

  /* Signal success. */
  return TRUE;
}
