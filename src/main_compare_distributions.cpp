/*!
 * \file    main_compare_distributions.cpp
 * \ingroup compare_distributions_exe
 *
 * \brief   The definition of the main entry point to the
 *          compare_distributions executable, and of associated
 *          functions.
 */

/*!
 * \defgroup compare_distributions_exe \
 *           The compare_distributions executable
 * \ingroup  manipulate_executable
 *
 * \brief    A module for the compare_distributions executable.
 */

#include "linear_distribution.h"
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
  fprintf(file, "Synopsis: compare_distributions "
    "<distribution1> <distribution2>\n");
}

/*!
 * \brief The main entry point to the compare_distributions executable.
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

  Distribution distribution_a;

  file = fopen(argv[1], "rb");
  if (NULL == file) {
    fprintf(stderr, "Error: Failed to open \"%s\" for reading.\n", argv[1]);
    return -1;
  }

  printf("Loading distribution \"%s\"...\n", argv[1]);
  distribution_init_import(&distribution_a, file);

  fclose(file);
  file = NULL;

  Distribution distribution_b;

  file = fopen(argv[2], "rb");
  if (NULL == file) {
    fprintf(stderr, "Error: Failed to open \"%s\" for reading.\n", argv[2]);
    return -1;
  }

  printf("Loading distribution \"%s\"...\n", argv[2]);
  distribution_init_import(&distribution_b, file);

  fclose(file);
  file = NULL;

  printf("\n");

  long double total_probability_a = 0;
  long double total_probability_b = 0;

  long double total_abs_delta = 0;
  long double max_abs_delta = 0;

  for (int32_t td = -30; td <= 30; td++) {
    for (int32_t sgnd = 0; sgnd <= 1; sgnd++) {
      for (int32_t tr = -30; tr <= 30; tr++) {
        for (int32_t sgnr = 0; sgnr <= 1; sgnr++) {
          Distribution_Slice * slice_a = NULL;
          Distribution_Slice * slice_b = NULL;

          int32_t td_a = distribution_a.parameters.m + td;
          int32_t td_b = distribution_b.parameters.m + td;

          if (0 != sgnd) {
            td_a = -td_a;
            td_b = -td_b;
          }

          int32_t tr_a = distribution_a.parameters.m + tr;
          int32_t tr_b = distribution_b.parameters.m + tr;

          if (0 != sgnr) {
            tr_a = -tr_a;
            tr_b = -tr_b;
          }

          for (uint32_t i = 0; i < distribution_a.count; i++) {
            if ((distribution_a.slices[i]->min_log_alpha_d == td_a) &&
                (distribution_a.slices[i]->min_log_alpha_r == tr_a))
            {
              slice_a = distribution_a.slices[i];
            }
          }

          for (uint32_t i = 0; i < distribution_b.count; i++) {
            if ((distribution_b.slices[i]->min_log_alpha_d == td_b) &&
                (distribution_b.slices[i]->min_log_alpha_r == tr_b))
            {
              slice_b = distribution_b.slices[i];
            }
          }

          if ((NULL == slice_a) && (NULL == slice_b)) {
            continue; /* Skip this slice. */
          }

          if ((NULL == slice_a) || (NULL == slice_b)) {
            fprintf(stderr, "Warning: "
              "Failed to find slice (%d, %d) in both distributions.\n", td, tr);
            continue;
          }

          const long double abs_delta =
            abs_ld(slice_a->total_probability - slice_b->total_probability);

          printf("(%d, %d) (%d, %d): %Lg %Lg -- %Lg\n",
            slice_a->min_log_alpha_d,
            slice_a->min_log_alpha_r,
            slice_b->min_log_alpha_d,
            slice_b->min_log_alpha_r,
            slice_a->total_probability,
            slice_b->total_probability,
            abs_delta);

          total_probability_a += slice_a->total_probability;
          total_probability_b += slice_b->total_probability;

          total_abs_delta += abs_delta;
          if (abs_delta > max_abs_delta) {
            max_abs_delta = abs_delta;
          }
        }
      }
    }
  }

  printf("\n%Lf %Lf -- %Lg (max: %Lg)\n",
    total_probability_a,
    total_probability_b,
    total_abs_delta,
    max_abs_delta);

  distribution_clear(&distribution_a);
  distribution_clear(&distribution_b);

  return 0;
}
