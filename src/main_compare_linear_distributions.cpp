/*!
 * \file    main_compare_linear_distributions.cpp
 * \ingroup compare_linear_distributions_exe
 *
 * \brief   The definition of the main entry point to the
 *          compare_linear_distributions executable, and of associated
 *          functions.
 */

/*!
 * \defgroup compare_linear_distributions_exe \
 *           The compare_linear_distributions executable
 * \ingroup  manipulate_executable
 *
 * \brief    A module for the compare_linear_distributions executable.
 */

#include "linear_distribution.h"
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
  fprintf(file, "Synopsis: compare_linear_distributions "
    "<distribution1> <distribution2>\n");
}

/*!
 * \brief The main entry point to the compare_linear_distributions executable.
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

  Linear_Distribution distribution_a;

  file = fopen(argv[1], "rb");
  if (NULL == file) {
    fprintf(stderr, "Error: Failed to open \"%s\" for reading.\n", argv[1]);
    return -1;
  }

  printf("Loading distribution \"%s\"...\n", argv[1]);
  linear_distribution_init_import(&distribution_a, file);

  fclose(file);
  file = NULL;

  Linear_Distribution distribution_b;

  file = fopen(argv[2], "rb");
  if (NULL == file) {
    fprintf(stderr, "Error: Failed to open \"%s\" for reading.\n", argv[2]);
    return -1;
  }

  printf("Loading distribution \"%s\"...\n", argv[2]);
  linear_distribution_init_import(&distribution_b, file);

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

  for (int32_t t = -t_max; t <= t_max; t++) {
    Linear_Distribution_Slice * slice_a = NULL;
    Linear_Distribution_Slice * slice_b = NULL;

    int32_t t_a = distribution_a.parameters.m + t;
    int32_t t_b = distribution_b.parameters.m + t;

    for (uint32_t i = 0; i < distribution_a.count; i++) {
      if (distribution_a.slices[i]->min_log_alpha == t_a) {
        slice_a = distribution_a.slices[i];
      }
    }

    for (uint32_t i = 0; i < distribution_b.count; i++) {
      if (distribution_b.slices[i]->min_log_alpha == t_b) {
        slice_b = distribution_b.slices[i];
      }
    }

    if ((NULL == slice_a) && (NULL == slice_b)) {
      continue; /* Skip this slice. */
    }

    if ((NULL == slice_a) || (NULL == slice_b)) {
      fprintf(stderr,
        "Warning: Failed to find slice %d in one distribution.\n", t);
      continue;
    }

    const long double abs_delta =
      abs_ld(slice_a->total_probability - slice_b->total_probability);

    printf("%d %d: %Lg %Lg -- %Lg\n",
      slice_a->min_log_alpha,
      slice_b->min_log_alpha,
      2 * slice_a->total_probability,
      2 * slice_b->total_probability,
      2 * abs_delta);

    total_probability_a += slice_a->total_probability;
    total_probability_b += slice_b->total_probability;

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

  linear_distribution_clear(&distribution_a);
  linear_distribution_clear(&distribution_b);

  return 0;
}
