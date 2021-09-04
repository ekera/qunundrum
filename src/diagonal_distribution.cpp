/*!
 * \file    diagonal_distribution.cpp
 * \ingroup diagonal_distribution
 *
 * \brief   The definition of functions for manipulating diagonal probability
 *          distributions.
 */

#include "diagonal_distribution.h"

#include "common.h"
#include "diagonal_distribution_slice.h"
#include "diagonal_parameters.h"
#include "diagonal_probability.h"
#include "errors.h"
#include "math.h"
#include "random.h"
#include "sample.h"

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*!
 * \brief   The bound on the absolute value of Delta when creating the histogram
 *          in Delta used to sample the probability distribution.
 */
#define BOUND_DELTA   1000

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
  const Diagonal_Parameters * const parameters,
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
  diagonal_parameters_init(&(distribution->parameters));
  diagonal_parameters_copy(&(distribution->parameters), parameters);

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
  Diagonal_Parameters parameters;

  diagonal_parameters_init(&parameters);
  diagonal_parameters_import(&parameters, file);

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
  diagonal_parameters_clear(&parameters);
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

  diagonal_parameters_clear(&(distribution->parameters));

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
  diagonal_parameters_export(&(distribution->parameters), file);

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
        "Debug: Selected slice: %u (%d)\n", i, slice->min_log_alpha_r);
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
  double * const max_log_alpha_r)
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

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("diagonal_distribution_sample_region(): "
    "Debug: Sampled slice via diagonal_distribution_sample_slice(): %d\n",
      slice->min_log_alpha_r);
  #endif

  diagonal_distribution_slice_sample_region(
    slice,
    random_state,
    min_log_alpha_r,
    max_log_alpha_r);

  #ifdef DEBUG_TRACE_SAMPLING
  printf("diagonal_distribution_sample_region(): "
    "Debug: Sampled region via diagonal_distribution_slice_sample_region(): "
      "%f %fn", *min_log_alpha_r, *max_log_alpha_r);
  #endif

  /* Signal success. */
  return TRUE;
}

bool diagonal_distribution_sample_alpha_r(
  const Diagonal_Distribution * const distribution,
  Random_State * const random_state,
  mpz_t alpha_r)
{
  /* Sample a region from the distribution. */
  double min_log_alpha_r;
  double max_log_alpha_r;

  bool result;

  result = diagonal_distribution_sample_region(
    distribution,
    random_state,
    &min_log_alpha_r,
    &max_log_alpha_r);
  if (FALSE == result) {
    #ifdef DEBUG_TRACE_SAMPLING
    printf("diagonal_distribution_sample_alpha_r(): "
      "Debug: Sampled region via diagonal_distribution_sample_region(): "
        "Out of bounds.\n");
    #endif

    mpz_set_ui(alpha_r, 0); /* no notion of infinity exists in GMP */

    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  printf("diagonal_distribution_sample_alpha_r(): "
    "Debug: Sampled region via diagonal_distribution_sample_region(): "
      "%f %f\n", min_log_alpha_r, max_log_alpha_r);
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
  printf("diagonal_distribution_sample_alpha_r(): "
    "Debug: Kappa: %u\n", kappa_r);
  gmp_printf("diagonal_distribution_sample_alpha_r(): "
    "Debug: Sampled alpha_r via sample_alpha_from_region(): %Zd\n", alpha_r);
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
  const uint32_t precision =
    2 * (max_ui(distribution->parameters.m +
      distribution->parameters.sigma, PRECISION));

  mpz_t alpha_r;
  mpz_init(alpha_r);

  bool result;

  result = diagonal_distribution_sample_alpha_r(
              distribution,
              random_state,
              alpha_r);
  if (FALSE == result) {
    #ifdef DEBUG_TRACE_SAMPLING
    gmp_printf("diagonal_distribution_sample_pair_j_k(): "
      "Debug: Sampled alpha_r via "
        "diagonal_distribution_sample_alpha_r():\n Out of bounds.\n");
    #endif

    mpz_clear(alpha_r);

    mpz_set_ui(j, 0);
    mpz_set_ui(k, 0);

    /* Signal sampling error in f. */
    return FALSE;
  }

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("diagonal_distribution_sample_pair_j_k(): "
    "Debug: Sampled alpha_r via "
      "diagonal_distribution_sample_alpha_r():\n"
        "   alpha_r: %Zd\n", alpha_r);
  #endif

  /* Sample an integer j given alpha_r. */
  sample_j_from_diagonal_alpha_r(
    j,
    alpha_r,
    &(distribution->parameters),
    random_state);

  #ifdef DEBUG_TRACE_SAMPLING
  gmp_printf("diagonal_distribution_sample_pair_j_k(): "
    "Debug: Sampled (j, k) via sample_j_k_from_diagonal_alpha_d_r():\n"
      "   j: %Zd\n   k: %Zd\n", j, k);
  #endif

  /* Initialize additional variables. */
  mpz_t k0;
  mpz_init(k0);

  mpz_t tmp_z;
  mpz_init(tmp_z);

  mpfr_t tmp;
  mpfr_init2(tmp, precision);

  mpfr_t tmp2;
  mpfr_init2(tmp2, precision);

  /* Compute: ((d / r) alpha_r - d * j) / 2^(m + sigma - l) */
  mpfr_set_z(tmp, alpha_r, MPFR_RNDN); /* tmp = alpha_r */
  mpfr_mul_z(tmp, tmp, distribution->parameters.d, MPFR_RNDN);
    /* tmp = alpha_r * d */
  mpfr_div_z(tmp, tmp, distribution->parameters.r, MPFR_RNDN);
    /* tmp = alpha_r * d / r */

  mpz_mul(tmp_z, distribution->parameters.d, j); /* tmp_z = d * j */
  mpfr_sub_z(tmp, tmp, tmp_z, MPFR_RNDN); /* tmp = alpha_r * d / r - d * j */

  mpz_set_ui(tmp_z, 0); /* tmp_z = 0 */
  mpz_setbit(tmp_z,
    distribution->parameters.m +
    distribution->parameters.sigma -
    distribution->parameters.l); /* tmp_z = 2^(m + sigma - l) */
  mpfr_div_z(tmp, tmp, tmp_z, MPFR_RNDN);
    /* tmp = (alpha_r * d / r - d * j) / 2^(m + sigma - l) */

  mpfr_round(tmp, tmp);
  mpfr_get_z(k0, tmp, MPFR_RNDN); /* k0 = round(tmp) */

  /* Sample pivot. */
  long double pivot = random_generate_pivot_inclusive(random_state);

  mpz_t alpha_d;
  mpz_init(alpha_d);

  for (uint32_t delta = 0; delta < BOUND_DELTA; delta++) {
    for (int32_t sign = 1; sign >= -1; sign -= 2) {
      /* Compute alpha_d. */
      if (1 == sign) {
        mpz_add_ui(k, k0, delta);
      } else {
        if (0 == delta) {
          continue;
        }

        mpz_sub_ui(k, k0, delta);
      }

      mpz_set_ui(tmp_z, 0); /* tmp_z = 0 */
      mpz_setbit(tmp_z, distribution->parameters.l); /* tmp_z = 2^l */
      mpz_mod(k, k, tmp_z);

      mpz_set_ui(tmp_z, 0); /* tmp_z = 0 */
      mpz_setbit(tmp_z,
        distribution->parameters.m +
          distribution->parameters.sigma -
            distribution->parameters.l); /* tmp_z = 2^(m + sigma - l) */
      mpz_mul(tmp_z, tmp_z, k); /* tmp_z = 2^(m + sigma - l) k */

      mpz_mul(alpha_d, distribution->parameters.d, j); /* alpha_d = dj */
      mpz_add(alpha_d, alpha_d, tmp_z);
        /* alpha_d = dj + 2^(m + sigma - l) k */

      mpz_set_ui(tmp_z, 0); /* tmp_z = 0 */
      mpz_setbit(tmp_z,
        distribution->parameters.m + distribution->parameters.sigma);
          /* tmp_z = 2^(m + sigma) */
      mod_reduce(alpha_d, tmp_z);
        /* alpha_d = {dj + 2^(m + sigma - l) k}_{2^(m + sigma)} */

      mpfr_set_z(tmp, distribution->parameters.d, MPFR_RNDN); /* tmp = d */
      mpfr_mul_z(tmp, tmp, alpha_r, MPFR_RNDN); /* tmp = d * alpha_r */
      mpfr_div_z(tmp, tmp, distribution->parameters.r, MPFR_RNDN);
        /* tmp = (d / r) * alpha_r */

      mpfr_set_z(tmp2, alpha_d, MPFR_RNDN); /* tmp2 = alpha_d */
      mpfr_sub(tmp, tmp2, tmp, MPFR_RNDN); /* tmp = alpha_d - (d / r) alpha_r */

      mpfr_const_pi(tmp2, MPFR_RNDN); /* tmp2 = pi */
      mpfr_mul_ui(tmp2, tmp2, 2, MPFR_RNDN); /* tmp2 = 2 pi */
      mpfr_div_z(tmp2, tmp2, tmp_z, MPFR_RNDN);
        /* tmp2 = 2 pi / 2^(m + sigma) */
      mpfr_mul(tmp, tmp2, tmp, MPFR_RNDN);
        /* tmp = (2 pi / 2^(m + sigma)) (alpha_d - (d / r) alpha_r) */

      diagonal_probability_approx_h(tmp2, tmp, &(distribution->parameters));

      pivot -= mpfr_get_ld(tmp2, MPFR_RNDN);

      if (pivot <= 0) {
        break;
      }
    }

    if (pivot <= 0) {
      break;
    }
  }

  /* Clear memory. */
  mpz_clear(alpha_r);
  mpz_clear(alpha_d);

  mpz_clear(k0);
  mpz_clear(tmp_z);

  mpfr_clear(tmp);
  mpfr_clear(tmp2);

  if (pivot > 0) {
    mpz_set_ui(j, 0);
    mpz_set_ui(k, 0);

    /* Signal sampling error in h. */
    return FALSE;
  }

  /* Signal success. */
  return TRUE;
}
