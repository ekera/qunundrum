/*!
 * \file    main_plot_diagonal_distribution.cpp
 * \ingroup plot_diagonal_distribution_exe
 *
 * \brief   The definition of the main entry point to the
 *          plot_diagonal_distribution executable, and of associated functions.
 */

/*!
 * \defgroup plot_diagonal_distribution_exe \
 *           The plot_diagonal_distribution executable
 * \ingroup  plot_executable
 *
 * \brief    A module for the plot_diagonal_distribution executable.
 */

#include "executables.h"
#include "executables_plot_distribution.h"

#include "diagonal_distribution.h"
#include "diagonal_distribution_slice.h"
#include "diagonal_distribution_enumerator.h"
#include "random.h"
#include "common.h"
#include "string_utilities.h"

#include "plot_distribution.h"
#include "plot_distribution_common.h"
#include "plot_distribution_axis.h"

#include "errors.h"
#include "math.h"

#include <mpfr.h>

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include <unistd.h>

#include <sys/stat.h>
#include <sys/types.h>

static void plot_diagonal_distribution(
  Diagonal_Distribution * const distribution,
  const bool absolute)
{
  /* Assert that the distribution is non-empty. */
  if (0 == distribution->count) {
    critical("process_diagonal_distribution(): The distribution is empty.");
  }

  /* Assert that all slices are of the same dimension. */
  for (uint32_t i = 0; i < distribution->count; i++) {
    const uint32_t dimension = distribution->slices[i]->dimension;

    if ((0 == dimension) || ((dimension % SCALED_DIMENSION_LINEAR) != 0)) {
      critical("process_diagonal_distribution(): The distribution contains "
        "slices of dimension not divisible by %u.", SCALED_DIMENSION_LINEAR);
    }
  }

  /* Extract m and s. */
  const uint32_t m = distribution->parameters.m;
  const uint32_t l = distribution->parameters.l;

  /* Check that the output directory exists and is accessible. */
  if (0 != access(PLOTS_DIRECTORY, F_OK)) {
    if (0 != mkdir(PLOTS_DIRECTORY, DIRECTORY_PERMISSIONS)) {
      critical("plot_diagonal_distribution(): The output directory \"%s\" does "
        "not exist and could not be created.\n", PLOTS_DIRECTORY);
    }
  }

  /* Setup the path to the output file. */
  char path[MAX_SIZE_PATH_BUFFER];
  safe_snprintf(path, MAX_SIZE_PATH_BUFFER,
    "%s/plot-diagonal-distribution-m-%u-l-%u.tex",
    PLOTS_DIRECTORY, m, l);

  /* Check that the output file does not already exist. */
  if (0 == access(path, F_OK)) {
    critical("plot_diagonal_distribution(): "
      "The plot \"%s\" already exists.", path);
  }

  /* Open the output file. */
  FILE * file = fopen(path, "wb");
  if (NULL == file) {
    critical("plot_diagonal_distribution(): Failed to open \"%s\".", path);
  }

  printf("Writing the distribution to \"%s\"...\n", path);

  /* Write the file header. */
  fprintf(file, "%% m = %u\n", m);
  fprintf(file, "%% l = %u\n", l);
  gmp_fprintf(file, "%% d = %Zd\n", distribution->parameters.d);
  gmp_fprintf(file, "%% r = %Zd\n", distribution->parameters.r);
  fprintf(file, "%% dimension = %u\n", SCALED_DIMENSION_LINEAR);
  fprintf(file, "%% flags = %.8x\n\n", distribution->flags);

  fprintf(file, "\\documentclass[crop, tikz]{standalone}\n");
  fprintf(file, "\\usepackage{pgf}\n");
  fprintf(file, "\\usepackage{tikz}\n");
  fprintf(file, "\\usepackage{amsmath}\n");
  fprintf(file, "\\usepackage{amssymb}\n\n");

  fprintf(file, "\\begin{document}\n");
  fprintf(file, "\\begin{tikzpicture}\n\n");

  /* Plot the linear distribution along the horizontal axis. */
  Diagonal_Distribution scaled_distribution;

  diagonal_distribution_init_copy_scale(
    &scaled_distribution,
    distribution,
    SCALED_DIMENSION_LINEAR);

  plot_collapsed_diagonal_distribution_horizontal(
    &scaled_distribution,
    0, /* offset_y */
    file,
    absolute);

  /* Write the file footer. */
  fprintf(file, "\n");
  fprintf(file, "\\end{tikzpicture}\n");
  fprintf(file, "\\end{document}\n");

  /* Close the file. */
  fclose(file);
  file = NULL;

  /* Clear memory. */
  diagonal_distribution_clear(&scaled_distribution);
}

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  fprintf(file, "Synopsis: plot_diagonal_distribution "
    "[ -sgn | -abs ] <distribution> { <distribution> }\n");
  fprintf(file, "\n");
  fprintf(file, "Absolute or signed: -- defaults to -sgn\n");
  fprintf(file, " -abs  "
    "Assume symmetry around zero. Plot only the positive side of the axis.\n");
  fprintf(file, " -sgn  "
    "Independently plot both the negative and positive sides of the axis.\n");
}

/*!
 * \brief The main entry point to the plot_diagonal_distribution executable.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \return Zero upon successful execution, non-zero otherwise.
 */
int main(int argc, char ** argv) {
  mpfr_set_default_prec(PRECISION);

  /* Parse command line arguments. */
  if (argc < 2) {
    print_synopsis((argc == 1) ? stdout : stderr);
    return (argc == 1) ? 0 : -1;
  }

  bool absolute = FALSE;

  int i = 1;

  if (strcmp(argv[i], "-sgn") == 0) {
    absolute = FALSE;
    i++;
  } else if (strcmp(argv[i], "-abs") == 0) {
    absolute = TRUE;
    i++;
  }

  for ( ; i < argc; i++) {
    /* Load the distribution. */
    Diagonal_Distribution distribution;

    printf("Importing the distribution from \"%s\"...\n", argv[i]);
    fflush(stdout);

    FILE * file = fopen(argv[i], "rb");
    if (NULL == file) {
      critical("main(): Failed to open \"%s\".", argv[i]);
    }

    diagonal_distribution_init_import(&distribution, file);
    fclose(file);

    /* Process the distribution. */
    plot_diagonal_distribution(&distribution, absolute);

    /* Clear memory. */
    diagonal_distribution_clear(&distribution);
  }

  printf("Done.\n");
  fflush(stdout);

  return 0;
}
