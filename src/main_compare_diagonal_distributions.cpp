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

#include "common.h"
#include "diagonal_distribution.h"
#include "errors.h"
#include "math.h"

#include <stdint.h>
#include <stdio.h>

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

  const uint32_t t_a = distribution_a.parameters.t;
  const uint32_t t_b = distribution_b.parameters.t;

  if ((t_a > 256) || (t_b > 256)) {
    critical("main(): The range covered by the distribution is too large.");
  }

  const int32_t t_max = (int32_t)max_ui(t_a, t_b);

  const uint32_t eta_bound_a = distribution_a.parameters.eta_bound;
  const uint32_t eta_bound_b = distribution_b.parameters.eta_bound;

  const uint32_t eta_bound_min = min_ui(eta_bound_a, eta_bound_b);

  for (uint32_t eta_abs = 0; eta_abs <= eta_bound_min; eta_abs++) {
    for (int32_t eta_sgn = 1; eta_sgn >= -1; eta_sgn -= 2) {

      if ((0 == eta_abs) && (-1 == eta_sgn)) {
        continue;
      }

      const int32_t eta = eta_sgn * ((int32_t)eta_abs);

      printf("-- slices for eta: %d\n", eta);

      for (int32_t t = -t_max; t <= t_max; t++) {
        int32_t t_a = ((int32_t)distribution_a.parameters.m) + t;
        int32_t t_b = ((int32_t)distribution_b.parameters.m) + t;

        long double total_slice_probability_a = 0;
        bool found_slice_a = FALSE;

        for (uint32_t i = 0; i < distribution_a.count; i++) {
          if ((distribution_a.slices[i]->min_log_alpha_r == t_a) &&
              (distribution_a.slices[i]->eta == eta))
          {
            total_slice_probability_a +=
              distribution_a.slices[i]->total_probability;
            found_slice_a = TRUE;
          }
        }

        long double total_slice_probability_b = 0;
        bool found_slice_b = FALSE;

        for (uint32_t i = 0; i < distribution_b.count; i++) {
          if ((distribution_b.slices[i]->min_log_alpha_r == t_b) &&
              (distribution_b.slices[i]->eta == eta))
          {
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
            "Warning: Failed to find slice %d in one distribution. This may "
              "be due to a difference in the number of padding bits.\n", t);
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

      printf("\n");
    }
  }

  if (eta_bound_a > eta_bound_min) {
    fprintf(stderr, "Warning: Slices for eta in [%d, %d] have been ignored "
      "in the above comparison. These slices are only present in the first "
        "distribution.\n", eta_bound_min + 1, eta_bound_a);
  } else if (eta_bound_b > eta_bound_min) {
    fprintf(stderr, "Warning: Slices for eta in [%d, %d] have been ignored "
      "in the above comparison. These slices are only present in the second "
        "distribution.\n", eta_bound_min + 1, eta_bound_b);
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
