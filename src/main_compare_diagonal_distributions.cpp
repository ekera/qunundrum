/*!
 * \file    main_compare_diagonal_distributions.cpp
 * \ingroup compare_diagonal_distributions_exe
 *
 * \brief   The definition of the main entry point to the
 *          compare_diagonal_distributions executable, and of associated
 *          functions.
 */

/*!
 * \defgroup compare_diagonal_distributions_exe \
 *           The compare_diagonal_distributions executable
 * \ingroup  manipulate_executable
 *
 * \brief    A module for the compare_diagonal_distributions executable.
 */

#include "diagonal_distribution.h"
#include "parameters.h"
#include "errors.h"
#include "math.h"

#include <mpfr.h>

#include <stdio.h>
#include <stdint.h>

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  fprintf(file, "Synopsis: compare_diagonal_distributions "
    "<distribution1> <distribution2>\n");
}

/*!
 * \brief The main entry point to the compare_diagonal_distributions executable.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \return Zero upon successful execution, non-zero otherwise.
 */
int main(int argc, char ** argv) {
  if (argc <= 1) {
    print_synopsis(stdout);
    return 0;
  }

  if (3 != argc) {
    fprintf(stderr, "Error: Incorrect command line arguments.\n");
    return -1;
  }

  FILE * file;

  Diagonal_Distribution distribution_a;

  file = fopen(argv[1], "rb");
  if (NULL == file) {
    fprintf(stderr, "Error: Failed to open \"%s\" for reading.\n", argv[1]);
    return -1;
  }

  printf("Loading distribution \"%s\"...\n", argv[1]);
  diagonal_distribution_init_import(&distribution_a, file);

  fclose(file);
  file = NULL;


  Diagonal_Distribution distribution_b;

  file = fopen(argv[2], "rb");
  if (NULL == file) {
    fprintf(stderr, "Error: Failed to open \"%s\" for reading.\n", argv[2]);
    return -1;
  }

  printf("Loading distribution \"%s\"...\n", argv[2]);
  diagonal_distribution_init_import(&distribution_b, file);

  fclose(file);
  file = NULL;

  printf("\n");

  long double total_probability_a = 0;
  long double total_probability_b = 0;

  long double total_abs_delta = 0;
  long double max_abs_delta = 0;

  for (int32_t t = -30; t <= 30; t++) {
    int32_t t_a = distribution_a.parameters.m + t;
    int32_t t_b = distribution_b.parameters.m + t;

    long double total_slice_probability_a = 0;
    bool found_slice_a = FALSE;

    for (uint32_t i = 0; i < distribution_a.count; i++) {
      if (distribution_a.slices[i]->min_log_alpha_r == t_a) {
        total_slice_probability_a += 
          distribution_a.slices[i]->total_probability;
        found_slice_a = TRUE;
      }
    }

    long double total_slice_probability_b = 0;
    bool found_slice_b = FALSE;

    for (uint32_t i = 0; i < distribution_b.count; i++) {
      if (distribution_b.slices[i]->min_log_alpha_r == t_b) {
        total_slice_probability_b += 
          distribution_b.slices[i]->total_probability;
        found_slice_b = TRUE;
      }
    }

    if ((FALSE == found_slice_a) && (FALSE == found_slice_b)) {
      continue; /* Skip this slice. */
    }

    if ((FALSE == found_slice_a) || (FALSE == found_slice_b)) {
      fprintf(stderr,
        "Warning: Failed to find slice %d in both distributions. This may be "
          "due to a difference in the number of padding bits.\n", t);
      continue;
    }

    const long double abs_delta =
      abs_ld(total_slice_probability_a - total_slice_probability_b);

    printf("%d %d: %Lg %Lg -- %Lg\n",
      t_a,
      t_b,
      2 * total_slice_probability_a,
      2 * total_slice_probability_b,
      2 * abs_delta);

    total_probability_a += total_slice_probability_a;
    total_probability_b += total_slice_probability_b;

    total_abs_delta += abs_delta;
    if (abs_delta > max_abs_delta) {
      max_abs_delta = abs_delta;
    }
  }

  printf("\n%Lf %Lf -- %Lg (max: %Lg)\n",
    2 * total_probability_a, /* double probabilities due to symmetry */
    2 * total_probability_b,
    2 * total_abs_delta,
    2 * max_abs_delta);

  diagonal_distribution_clear(&distribution_a);
  diagonal_distribution_clear(&distribution_b);

  return 0;
}
