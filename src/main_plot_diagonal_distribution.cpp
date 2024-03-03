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

#include "common.h"
#include "diagonal_distribution.h"
#include "errors.h"
#include "math.h"
#include "plot_distribution.h"
#include "plot_distribution_axis.h"
#include "plot_distribution_common.h"
#include "string_utilities.h"

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <unistd.h>
#include <sys/stat.h>

static void plot_diagonal_distribution(
  Diagonal_Distribution * const distribution,
  const bool absolute,
  const uint32_t eta_bound = UINT32_MAX)
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

  /* Extract m, sigma, s and l. */
  const uint32_t m = distribution->parameters.m;
  const uint32_t sigma = distribution->parameters.sigma;
  const uint32_t s = distribution->parameters.s;
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
    "%s/plot-diagonal-distribution-m-%u-sigma-%u-%c-%u.tex",
    PLOTS_DIRECTORY,
    m,
    sigma,
    (0 != s) ? 's' : 'l',
    (0 != s) ?  s  :  l);

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

  printf("Writing the plot to \"%s\"...\n", path);

  /* Write the file header. */
  fprintf(file, "%% m = %u\n", m);
  fprintf(file, "%% sigma = %u\n", sigma);
  fprintf(file, "%% s = %u\n", s);
  fprintf(file, "%% l = %u\n", l);
  gmp_fprintf(file, "%% d = %Zd\n", distribution->parameters.d);
  gmp_fprintf(file, "%% r = %Zd\n", distribution->parameters.r);
  fprintf(file, "%% dimension = %u\n", SCALED_DIMENSION_LINEAR);
  fprintf(file, "%% eta-bound = %u\n", distribution->parameters.eta_bound);
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

  plot_diagonal_distribution_horizontal_collapsed(
    &scaled_distribution,
    0, /* offset_y */
    file,
    absolute,
    eta_bound);

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
  fprintf(file, "Synopsis: plot_diagonal_distribution \\\n"
    "   [ -sgn | -abs ] [ -eta-bound <eta-bound> ] \\\n"
      "      <distribution> { <distribution> }\n");

  fprintf(file, "\n");
  fprintf(file, "Absolute or signed: -- defaults to -sgn\n");
  fprintf(file, " -abs  "
    "Assume symmetry around zero. Plot only the positive side of the axis.\n");
  fprintf(file, " -sgn  "
    "Independently plot both the negative and positive sides of the axis.\n");

  fprintf(file, "\n");
  fprintf(file, "Eta bound: -- "
    "defaults to the eta bound for the distribution\n");
  fprintf(file,
    " -eta-bound    explicitly set the eta bound to <eta-bound>\n");
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
  uint32_t eta_bound = UINT32_MAX;

  bool absolute_specified = FALSE;
  bool eta_bound_specified = FALSE;

  int i = 1;

  for (i = 1; i < argc; i++) {
    /* Parse signed and absolute flags. */
    if ((strcmp(argv[i], "-sgn") == 0) || (strcmp(argv[i], "-abs") == 0)) {
      if (TRUE == absolute_specified) {
        fprintf(stderr, "Error: The -sgn and -abs flags cannot be twice "
          "specified, or both be specified.\n");
        return -1;
      }

      if (strcmp(argv[i], "-sgn") == 0) {
        absolute = FALSE;
      } else if (strcmp(argv[i], "-abs") == 0) {
        absolute = TRUE;
      }

      absolute_specified = TRUE;

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

    break;
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

    /* Plot the distribution. */
    if (FALSE == eta_bound_specified) {
      eta_bound = distribution.parameters.eta_bound;
    } else {
      if (eta_bound > distribution.parameters.eta_bound) {
        critical("main(): The eta bound specified via -eta-bound is greater "
          "than the eta bound used to generate the distribution.");
      }
    }

    plot_diagonal_distribution(&distribution, absolute, eta_bound);

    /* Clear memory. */
    diagonal_distribution_clear(&distribution);
  }

  printf("Done.\n");
  fflush(stdout);

  return 0;
}
