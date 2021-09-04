/*!
 * \file    main_info_diagonal_distribution.cpp
 * \ingroup info_diagonal_distribution_exe
 *
 * \brief   The definition of the main entry point to the
 *          info_diagonal_distribution executable, and of associated functions.
 */

/*!
 * \defgroup info_diagonal_distribution_exe \
 *           The info_diagonal_distribution executable
 * \ingroup  manipulate_executable
 *
 * \brief    A module for the info_diagonal_distribution executable.
 */

#include "common.h"
#include "diagonal_distribution.h"
#include "errors.h"
#include "math.h"

#include <mpfr.h>

#include <stdio.h>

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  fprintf(file, "Synopsis: info_diagonal_distribution <distribution>\n");
}

/*!
 * \brief The main entry point to the info_diagonal_distribution executable.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \return Zero upon successful execution, non-zero otherwise.
 */
int main(int argc, char ** argv) {
  /* Set precision. */
  mpfr_set_default_prec(PRECISION);

  /* Verify synopsis. */
  if (2 != argc) {
    print_synopsis((argc == 1) ? stdout : stderr);
    return (argc == 1) ? 0 : -1;
  }

  /* Load the distribution. */
  Diagonal_Distribution distribution;

  printf("Importing the distribution from \"%s\"...\n", argv[1]);

  FILE * file = fopen(argv[1], "rb");
  if (NULL == file) {
    critical("main(): Unable to open \"%s\".", argv[1]);
  }

  diagonal_distribution_init_import(&distribution, file);
  fclose(file);

  printf("Done importing the distribution.\n\n");
  fflush(stdout);

  /* Export distribution information. */
  diagonal_distribution_export_info(&distribution, stdout, TRUE, FALSE);

  /* Clear memory. */
  diagonal_distribution_clear(&distribution);

  return 0;
}
