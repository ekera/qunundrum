/*!
 * \file    main_filter_distribution.cpp
 * \ingroup filter_distribution_exe
 *
 * \brief   The definition of the main entry point to the filter_distribution
 *          executable, and of associated functions.
 */

/*!
 * \defgroup filter_distribution_exe \
 *           The filter_distribution executable
 * \ingroup  filter_executable
 *
 * \brief    A module for the filter_distribution executable.
 */

#include "executables.h"
#include "executables_filter_distribution.h"

#include "common.h"
#include "distribution.h"
#include "errors.h"
#include "string_utilities.h"

#include <mpfr.h>

#include <stdio.h>

#include <unistd.h>
#include <sys/stat.h>

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  fprintf(file, "Synopsis: filter_distribution "
    "<distribution> { <distribution> }\n");
}

/*!
 * \brief The main entry point to the filter_distribution executable.
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
  
  /* Check that the output directory exists and is accessible. */
  if (0 != access(DISTRIBUTIONS_DIRECTORY, F_OK)) {
    if (0 != mkdir(DISTRIBUTIONS_DIRECTORY, DIRECTORY_PERMISSIONS)) {
      fprintf(stderr, "Error: The output directory \"%s\" does not exist.\n",
        DISTRIBUTIONS_DIRECTORY);
      return FALSE;
    }
  }

  if (0 != access(DISTRIBUTIONS_DIRECTORY, R_OK | W_OK | X_OK)) {
    fprintf(stderr, "Error: The output directory \"%s\" is not readable, "
      "writeable and executable.\n", DISTRIBUTIONS_DIRECTORY);
    return FALSE;
  }

  /* Load the distributions. */
  Distribution distribution;

  for (int i = 1; i < argc; i++) {
    printf("Processing \"%s\"...\n", argv[i]);
  
    /* Import the distribution. */
    FILE * file = fopen(argv[i], "rb");
    if (NULL == file) {
      critical("main(): Unable to open \"%s\" for reading.", argv[i]);
    }
  
    distribution_init_import(&distribution, file);

    fclose(file);
    file = NULL;

    /* Filter the distribution. */    
    distribution_filter_slices(&distribution);
    distribution_sort_slices(&distribution);

    /* Construct the new path. */
    const char * const filename = truncate_path(argv[i]);

    char path[MAX_SIZE_PATH_BUFFER];
    safe_snprintf(path,
                  MAX_SIZE_PATH_BUFFER,
                  "%s/filtered-%s",
                  DISTRIBUTIONS_DIRECTORY,
                  filename);
    
    /* Export the filtered distribution. */
    file = fopen(path, "wb+");
    if (NULL == file) {
      critical("main(): Unable to open \"%s\" for writing.", path);
    }

    distribution_export(&distribution, file);

    fclose(file);
    file = NULL;

    /* Clear the distribution. */
    distribution_clear(&distribution);
  }

  printf("Done.\n");

  /* Signal success. */
  return 0;
}
