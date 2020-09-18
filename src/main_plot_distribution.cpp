/*!
 * \file    main_plot_distribution.cpp
 * \ingroup plot_distribution_exe
 *
 * \brief   The definition of the main entry point to the plot_distribution
 *          executable, and of associated functions.
 */

/*!
 * \defgroup plot_distribution_exe \
 *           The plot_distribution executable
 * \ingroup  plot_executable
 *
 * \brief    A module for the plot_distribution executable.
 */

#include "executables.h"
#include "executables_plot_distribution.h"

#include "distribution.h"
#include "linear_distribution.h"
#include "linear_distribution_slice.h"
#include "linear_distribution_enumerator.h"
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

#include <unistd.h>

#include <sys/stat.h>
#include <sys/types.h>

/*!
 * \brief   Plots a two-dimensional probability distribution to file.
 *
 * \param[in] distribution    The distribution to plot to file.
 */
static void plot_distribution(
  const Distribution * const distribution)
{
  /* Assert that the distribution is non-empty. */
  if (0 == distribution->count) {
    critical("plot_distribution(): The distribution is empty.");
    return;
  }

  /* Assert that all slice dimensions are multiples of SCALED_DIMENSION. */
  for (uint32_t i = 0; i < distribution->count; i++) {
    const uint32_t dimension = distribution->slices[i]->dimension;

    if ((0 == dimension) || ((dimension % SCALED_DIMENSION) != 0)) {
      critical("plot_distribution(): The distribution contains slices of "
        "dimension not divisible by %u.", SCALED_DIMENSION);
    }
  }

  /* Extract m and s. */
  const uint32_t m = distribution->parameters.m;
  const uint32_t s = distribution->parameters.s;
  const uint32_t l = distribution->parameters.l;

  /* Check that the output directory exists and is accessible. */
  if (0 != access(PLOTS_DIRECTORY, F_OK)) {
    if (0 != mkdir(PLOTS_DIRECTORY, DIRECTORY_PERMISSIONS)) {
      critical("plot_distribution(): The output directory \"%s\" does not "
        "exist and could not be created.\n", PLOTS_DIRECTORY);
    }
  }

  /* Setup the path to the output file. */
  char path[MAX_SIZE_PATH_BUFFER];
  safe_snprintf(path, MAX_SIZE_PATH_BUFFER,
    "%s/plot-distribution-m-%u-%c-%u.tex",
    PLOTS_DIRECTORY,
    m,
    (0 != s) ? 's' : 'l',
    (0 != s) ?  s :   l);

  /* Check that the output file does not already exist. */
  if (0 == access(path, F_OK)) {
    critical("plot_distribution(): The plot \"%s\" already exists.", path);
  }

  /* Open the output file. */
  FILE * file = fopen(path, "wb");
  if (NULL == file) {
    critical("plot_distribution(): Failed to open \"%s\" for writing.", path);
  }

  printf("Writing the distribution to \"%s\"...\n", path);

  fprintf(file, "%% Use e.g. pdflatex, xelatex or lualatex to compile this "
    "file to a PDF.\n");
  fprintf(file, "%% Using lualatex may be preferable as it allocates memory "
    "dynamically.\n");
  fprintf(file, "%% When using lualatex, you may need to uncomment the below "
    "line:\n");
  fprintf(file, "%% \\RequirePackage{luatex85}\n\n");

  /* Write the file header. */
  fprintf(file, "%% m = %u\n", m);
  fprintf(file, "%% s = %u\n", s);
  fprintf(file, "%% l = %u\n", l);
  gmp_fprintf(file, "%% d = %Zd\n", distribution->parameters.d);
  gmp_fprintf(file, "%% r = %Zd\n", distribution->parameters.r);
  fprintf(file, "%% dimension = %u\n\n", SCALED_DIMENSION);

  fprintf(file, "\\documentclass[crop, tikz]{standalone}\n");
  fprintf(file, "\\usepackage{pgf}\n");
  fprintf(file, "\\usepackage{tikz}\n");
  fprintf(file, "\\usepackage{amsmath}\n");
  fprintf(file, "\\usepackage{amssymb}\n\n");

  fprintf(file, "\\begin{document}\n");
  fprintf(file, "\\begin{tikzpicture}\n\n");

  /* Plot the distribution. */
  Distribution scaled_distribution;

  distribution_init_copy_scale(
    &scaled_distribution,
    distribution,
    SCALED_DIMENSION);

  plot_detailed_distribution(&scaled_distribution, file);

  /* Plot the collapsed distribution in alpha_d. */
  Linear_Distribution linear_distribution_d;

  linear_distribution_init_collapse_d(
    &linear_distribution_d,
    distribution);

  Linear_Distribution linear_distribution_d_scaled;

  linear_distribution_init_copy_scale(
    &linear_distribution_d_scaled,
    &linear_distribution_d,
    SCALED_DIMENSION_LINEAR);

  const double offset_y =
    plot_distribution_coordinate(-(double)(
      m + PLOT_DISTRIBUTION_MAX_OFFSET_M +
        PLOT_DISTRIBUTION_LINEAR_MAX_SIZE), m);

  plot_linear_distribution_horizontal(
    &linear_distribution_d_scaled,
    offset_y,
    file);

  /* Plot the collapsed distribution in alpha_r. */
  Linear_Distribution linear_distribution_r;

  linear_distribution_init_collapse_r(
    &linear_distribution_r,
    distribution);

  Linear_Distribution linear_distribution_r_scaled;

  linear_distribution_init_copy_scale(
    &linear_distribution_r_scaled,
    &linear_distribution_r,
    SCALED_DIMENSION_LINEAR);

  const double offset_x =
    plot_distribution_coordinate(-(double)(
      m + PLOT_DISTRIBUTION_MAX_OFFSET_M +
        PLOT_DISTRIBUTION_LINEAR_MAX_SIZE), m);

  plot_linear_distribution_vertical(
    &linear_distribution_r_scaled,
    offset_x,
    file);

  /* Write the file footer. */
  fprintf(file, "\n");
  fprintf(file, "\\end{tikzpicture}\n");
  fprintf(file, "\\end{document}\n");

  /* Close the file. */
  fclose(file);
  file = NULL;

  /* Clear memory. */
  distribution_clear(&scaled_distribution);

  linear_distribution_clear(&linear_distribution_d);
  linear_distribution_clear(&linear_distribution_d_scaled);

  linear_distribution_clear(&linear_distribution_r);
  linear_distribution_clear(&linear_distribution_r_scaled);
}

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  fprintf(file, "Synopsis: plot_distribution "
    "<distribution> { <distribution> }\n");
}

/*!
 * \brief The main entry point to the plot_distribution executable.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \return Zero upon successful execution, non-zero otherwise.
 */
int main(int argc, char ** argv) {
  mpfr_set_default_prec(PRECISION);

  if (argc < 2) {
    print_synopsis((argc == 1) ? stdout : stderr);
    return (argc == 1) ? 0 : -1;
  }

  for (int i = 1; i < argc; i++) {
    /* Load the distribution. */
    Distribution distribution;

    printf("Importing the distribution from \"%s\"...\n", argv[i]);
    fflush(stdout);

    FILE * file = fopen(argv[i], "rb");
    if (NULL == file) {
      critical("main(): Failed to open \"%s\".", argv[i]);
    }

    distribution_init_import(&distribution, file);

    /* Plot the distribution. */
    plot_distribution(&distribution);

    /* Close the file. */
    fclose(file);
    file = NULL;

    /* Clear memory. */
    distribution_clear(&distribution);
  }

  printf("Done.\n");
  fflush(stdout);

  return 0;
}
