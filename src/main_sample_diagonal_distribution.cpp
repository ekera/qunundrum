/*!
 * \file    main_sample_diagonal_distribution.cpp
 * \ingroup sample_diagonal_distribution_exe
 *
 * \brief   The definition of the main entry point to the
 *          sample_diagonal_distribution executable, and of associated
 *          functions.
 */

/*!
 * \defgroup sample_diagonal_distribution_exe \
 *           The sample_diagonal_distribution executable
 * \ingroup  sample_executable
 *
 * \brief    A module for the sample_diagonal_distribution executable.
 */

#include "common.h"
#include "diagonal_distribution.h"
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

#include <unistd.h>

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  printf("Synopsis: sample_diagonal_distribution \\\n"
    "   [ -delta-bound <delta-bound> ] [ -eta-bound <eta-bound> ] \\\n"
    "      <distribution> [ <n> ]\n");

  fprintf(file, "\n");
  fprintf(file, "Delta bound: -- defaults to %u\n", BOUND_DELTA);
  fprintf(file,
    " -delta-bound  explicitly set the delta bound to <delta-bound>\n");

  fprintf(file, "\n");
  fprintf(file, "Eta bound: -- "
    "defaults to the eta bound for the distribution\n");
  fprintf(file,
    " -eta-bound    explicitly set the eta bound to <eta-bound>\n");

  fprintf(file, "\n");
  fprintf(file, "Number of samples n: -- defaults to 1000\n");
  fprintf(file, " <n>  The number of samples to draw from <distribution>.\n");
}

/*!
 * \brief The main entry point to the sample_diagonal_distribution executable.
 *
 * This is a conveniency executable provide for increased accessibility.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \return Zero upon successful execution, non-zero otherwise.
 */
int main(int argc, char ** argv) {
  mpfr_set_default_prec(PRECISION);

  if (argc < 2) {
    print_synopsis(stdout);
    return 0;
  }

  /* Parse the command line arguments. */
  int i = 1;

  uint32_t eta_bound = 0;
  uint32_t delta_bound = BOUND_DELTA;

  bool delta_bound_specified = FALSE;
  bool eta_bound_specified = FALSE;

  for (i = 1; i < argc; i++) {
    /* Parse the delta bound. */
    if (0 == strcmp(argv[i], "-delta-bound")) {
      /* Check that a delta bound has not already been specified. */
      if (FALSE != delta_bound_specified) {
        fprintf(stderr, "Error: The delta bound cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr, "Error: Expected <delta-bound> to follow after "
          "-delta-bound.\n");
        return FALSE;
      }

      const int x = atoi(argv[i + 1]);
      if (x < 0) {
        fprintf(stderr, "Error: The <delta-bound> passed to -delta-bound must "
          "be non-negative.\n");
        return FALSE;
      }

      /* Store the delta bound. */
      delta_bound = (uint32_t)x;
      delta_bound_specified = TRUE;

      i++;

      continue;
    }

    /* Parse the eta bound. */
    if (0 == strcmp(argv[i], "-eta-bound")) {
      /* Check that an eta bound has not already been specified. */
      if (FALSE != eta_bound_specified) {
        fprintf(stderr, "Error: The eta bound cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr, "Error: Expected <eta-bound> to follow after "
          "-eta-bound.\n");
        return FALSE;
      }

      const int x = atoi(argv[i + 1]);
      if ((x < 0) || (x > 256)) {
        fprintf(stderr, "Error: The <eta-bound> passed to -eta-bound must be "
          "on the interval [0, 256].\n");
        return FALSE;
      }

      /* Store the eta bound. */
      eta_bound = (uint32_t)x;
      eta_bound_specified = TRUE;

      i++;

      continue;
    }

    /* Stop parsing. */
    break;
  }

  if ((i + 1 != argc) && (i + 2 != argc)) {
    fprintf(stderr, "Error: Incorrect command line arguments.\n");
    return -1;
  }

  const char * const path = argv[i];

  if (0 != access(path, F_OK)) {
    fprintf(stderr, "Error: The distribution \"%s\" does not exist.\n", path);
    return -1;
  }

  if (0 != access(path, R_OK)) {
    fprintf(stderr, "Error: The distribution \"%s\" is not readable.\n", path);
    return -1;
  }

  /* Parse n. */
  uint32_t n = 1000;

  if (i + 2 == argc) {
    const int tmp = atoi(argv[i + 1]);

    if (tmp < 1) {
      fprintf(stderr, "Error: Failed to parse <n>.\n");
      return -1;
    } else {
      n = (uint32_t)tmp;
    }
  }

  /* Setup a random state. */
  Random_State random_state;
  random_init(&random_state);

  /* Load the distribution. */
  Diagonal_Distribution distribution;

  fprintf(stderr, "Importing the distribution from \"%s\"...\n", path);

  FILE * file = fopen(path, "rb");
  if (NULL == file) {
    fprintf(stderr, "Error: Failed to open \"%s\".", path);
    return -1;
  }

  diagonal_distribution_init_import(&distribution, file);
  fclose(file);

  if (TRUE == eta_bound_specified) {
    if (eta_bound > distribution.parameters.eta_bound) {
      critical("main(): The eta bound specified via -eta-bound is greater than "
        "the eta bound used to generate the distribution.");
    }
  }

  /* Sample from the distribution. */
  mpz_t alpha_d;
  mpz_init(alpha_d);

  mpz_t alpha_r;
  mpz_init(alpha_r);

  mpz_t j;
  mpz_init(j);

  mpz_t k;
  mpz_init(k);

  mpz_t tmp;
  mpz_init(tmp);

  for (uint32_t iteration = 0; iteration < n; iteration++) {
    printf("sample: %u / %u\n", iteration + 1, n);

    int32_t eta;

    bool result;

    result = diagonal_distribution_sample_pair_j_k(
      &distribution,
      &random_state,
      delta_bound,
      j,
      k,
      &eta);

    if (TRUE != result) {
      printf("*** out of bounds\n\n");
      continue;
    }

    if (TRUE == eta_bound_specified) {
      if (abs_i(eta) > eta_bound) {
        printf("*** out of bounds\n\n");
        continue;
      }
    }

    /* Compute alpha_r. */
    mpz_mul(alpha_r, distribution.parameters.r, j);
    mpz_set_ui(tmp, 0);
    mpz_setbit(tmp, distribution.parameters.m + distribution.parameters.sigma);
    mod_reduce(alpha_r, tmp);

    /* Compute alpha_d. */
    mpz_mul(alpha_d, distribution.parameters.d, j);

    mpz_set_ui(tmp, 0);
    mpz_setbit(tmp, distribution.parameters.m +
                    distribution.parameters.sigma -
                    distribution.parameters.l);
    mpz_mul(tmp, tmp, k);

    mpz_add(alpha_d, alpha_d, tmp);

    mpz_set_ui(tmp, 0);
    mpz_setbit(tmp, distribution.parameters.m + distribution.parameters.sigma);
    mod_reduce(alpha_d, tmp);

    /* Perform the printout. */
    gmp_printf("alpha_d: %Zd\n", alpha_d);
    gmp_printf("alpha_r: %Zd\n", alpha_r);

    gmp_printf("j: %Zd\n", j);
    gmp_printf("k: %Zd\n", k);
    gmp_printf("eta: %d\n", eta);

    printf("\n");
    fflush(stdout);
  }

  /* Close the random state. */
  random_close(&random_state);

  /* Clear memory. */
  mpz_clear(alpha_d);
  mpz_clear(alpha_r);
  mpz_clear(j);
  mpz_clear(k);

  mpz_clear(tmp);

  diagonal_distribution_clear(&distribution);

  return 0;
}
