/*!
 * \file    main_sample_distribution.cpp
 * \ingroup sample_distribution_exe
 *
 * \brief   The definition of the main entry point to the sample_distribution
 *          executable, and of associated functions.
 */

/*!
 * \defgroup sample_distribution_exe \
 *           The sample_distribution executable
 * \ingroup  sample_executable
 *
 * \brief    A module for the sample_distribution executable.
 */

#include "distribution.h"
#include "sample.h"
#include "random.h"
#include "common.h"

#include "errors.h"

#include <mpfr.h>
#include <gmp.h>

#include <stdio.h>
#include <stdint.h>
#include <stdint.h>

#include <unistd.h>

#include <sys/stat.h>
#include <sys/types.h>

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  printf("Synopsis: sample_distribution <distribution> [ <n> ]\n");
  fprintf(file, "\n");
  fprintf(file, "Number of samples n: -- defaults to 1000\n");
  fprintf(file, " <n>  The number of samples to draw from <distribution>.\n");
}

/*!
 * \brief The main entry point to the sample_distribution executable.
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

  if (1 == argc) {
    print_synopsis(stdout);
    return 0;
  }

  /* Parse the command line arguments. */
  if ((2 != argc) && (3 != argc)) {
    fprintf(stderr, "Error: Incorrect command line arguments.\n");
    return -1;
  }

  if (0 != access(argv[1], F_OK)) {
    fprintf(stderr, "Error: The distribution \"%s\" does not exist.\n",
      argv[1]);
    return -1;
  }

  if (0 != access(argv[1], R_OK)) {
    fprintf(stderr, "Error: The distribution \"%s\" is not readable.\n",
      argv[1]);
    return -1;
  }

  /* Parse n. */
  uint32_t n = 1000;

  if (3 == argc) {
    int tmp = atoi(argv[2]);

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
  Distribution distribution;

  fprintf(stderr, "Importing the distribution from \"%s\"...\n", argv[1]);

  FILE * file = fopen(argv[1], "rb");
  if (NULL == file) {
    fprintf(stderr, "Error: Failed to open \"%s\".", argv[1]);
    return -1;
  }

  distribution_init_import(&distribution, file);
  fclose(file);

  /* Sample from the distribution. */
  mpz_t alpha_d;
  mpz_init(alpha_d);

  mpz_t alpha_r;
  mpz_init(alpha_r);

  mpz_t j;
  mpz_init(j);

  mpz_t k;
  mpz_init(k);

  for (uint32_t i = 0; i < n; i++) {
    printf("sample: %u / %u\n", i + 1, n);

    bool result;

    result = distribution_sample_alpha_d_r(
      &distribution,
      &random_state,
      alpha_d,
      alpha_r);

    if (TRUE != result) {
      printf("*** out of bounds\n\n");
      continue;
    }

    gmp_printf("alpha_d: %Zd\n", alpha_d);
    gmp_printf("alpha_r: %Zd\n", alpha_r);

    sample_j_k_from_alpha_d_r(
      j,
      k,
      alpha_d,
      alpha_r,
      &(distribution.parameters),
      &random_state);

    gmp_printf("j: %Zd\n", j);
    gmp_printf("k: %Zd\n", k);

    printf("\n");
    fflush(stdout);
  }

  /* Clean up. */
  mpz_clear(alpha_d);
  mpz_clear(alpha_r);
  mpz_clear(j);
  mpz_clear(k);

  distribution_clear(&distribution);

  random_close(&random_state);

  return 0;
}
