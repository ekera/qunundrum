/*!
 * \file    main_generate_distribution.cpp
 * \ingroup generate_distribution_exe
 *
 * \brief   The definition of the main entry point to the generate_distribution
 *          executable, and of associated functions.
 */

/*!
 * \defgroup generate_distribution_exe \
 *           The generate_distribution executable
 * \ingroup  generate_executable
 *
 * \brief    A module for the generate_distribution executable.
 */

#include "executables.h"
#include "executables_generate_distribution.h"

#include "distribution.h"
#include "distribution_slice.h"
#include "distribution_enumerator.h"
#include "linear_distribution.h"
#include "parameters.h"
#include "parameters_selection.h"
#include "thread_pool.h"
#include "probability.h"
#include "math.h"
#include "errors.h"
#include "random.h"
#include "string_utilities.h"

#include <gmp.h>
#include <mpfr.h>
#include <mpi.h>

#include <pthread.h>

#include <memory.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <unistd.h>

#include <sys/stat.h>
#include <sys/types.h>

/*!
 * \brief   The minimum probability mass in a slice for it to be kept.
 *
 * If the total probability mass contained in a distribution slice is found to
 * be lower than the minimum slice probability, the slice is dropped.
 *
 * This allows the size of the distribution to be reduced.
 */
#define MIN_SLICE_PROBABILITY                 1e-16

/*!
 * \brief   The maximum dimension of slices not causing scaling to occur.
 *
 * If the slice dimension exceeds the maximum dimension, it will be scaled to
 * the maximum dimension before being inserted into the distribution.
 *
 * This allows the size of the distribution to be reduced.
 */
#define MAX_SLICE_DIMENSION                   256

/*!
 * \brief  An invalid dimension used to indicate that no fixed dimension or
 *         heuristic has been specified when parsing the command line arguments.
 */
#define DIMENSION_UNDEFINED                   0xffffffff

/*!
 * \brief  An invalid dimension used to indicate that the dimension should be
 *         heuristically selected, as opposed to being fixed.
 */
#define DIMENSION_HEURISTIC                   0

/*!
 * \brief   A data structure representing argument entries in the form of
 *          parsed \<m\> \<s\> or \<m\> \<l\> tuples from the command line
 *          arguments.
 *
 * \ingroup generate_distribution_exe
 */
typedef struct {
  /*!
   * \brief   The parameter m.
   */
  uint32_t m;

  /*!
   * \brief   The parameter s.
   */
  uint32_t s;

  /*!
   * \brief   The parameter l.
   */
  uint32_t l;

  /*!
   * \brief   The logarithm d.
   */
  mpz_t d;

  /*!
   * \brief   The order r.
   */
  mpz_t r;

  /*!
   * \brief   The path to which to export this distribution.
   */
  char path[MAX_SIZE_PATH_BUFFER];

  /*!
   * \brief   The base name from which to additional paths may be constructed.
   */
  char name[MAX_SIZE_PATH_BUFFER];
} Generate_Distribution_Arguments_Entry;

/*!
 * \brief   A data structure representing parsed command line arguments.
 *
 * \ingroup generate_distribution_exe
 */
typedef struct {
  /*!
   * \brief   The selection method for the logarithm d and order r.
   */
  Selection_Method selection_method;

  /*!
   * \brief   The tradeoff method.
   */
  Tradeoff_Method tradeoff_method;

  /*!
   * \brief   The method to use to select sigma in the probability estimate.
   *
   * Has no effect for #PROBABILITY_ESTIMATE_QUICK.
   */
  Sigma_Selection_Method sigma_method;

  /*!
   * \brief   The probability estimate.
   */
  Probability_Estimate probability_estimate;

  /*!
   * \brief   The explicitly specified fixed dimension.
   *
   * If the dimension is not explicitly specified, this entry is set to zero.
   */
  uint32_t dimension;

  /*!
   * \brief   The explicitly specified value of the logarithm d.
   *
   * If d is not explicitly specified, this entry is set to zero.
   */
  mpz_t explicit_d;

  /*!
   * \brief   The explicitly specified value of the order r.
   *
   * If r is not explicitly specified, this entry is set to zero.
   */
  mpz_t explicit_r;

  /*!
   * \brief   The number of \<m\> \<s\> or \<m\> \<l\> tuples passed to the
   *          executable.
   *
   * This number also corresponds to the number of distributions to generate.
   */
  uint32_t count;

  /*!
   * \brief   A vector of entries in the form of parsed tuples.
   *
   * Each entry corresponds to a distribution to generate.
   */
  Generate_Distribution_Arguments_Entry * entries;
} Generate_Distribution_Arguments;

/*!
 * \brief   A data structure representing an export job.
 *
 * \ingroup generate_distribution_exe
 */
typedef struct {
  /*!
   * \brief   A pointer to the linear distribution to export.
   */
  Distribution * distribution;

  /*!
   * \brief   The path to which to export the distribution.
   */
  char path[MAX_SIZE_PATH_BUFFER];

  /*!
   * \brief   The base name from which to additional paths may be constructed.
   */
  char name[MAX_SIZE_PATH_BUFFER];
} Distribution_Export_Job;


/*!
 * \brief   Parses the command line arguments.
 *
 * \param[in, out] arguments   The argument data structure in which to store
 *                             the parsed command line arguments.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \remark   So as to allow parallelized executables to be shut down gracefully,
 *           this function does not call critical() to signal a critical error.
 *           Rather it prints an informative error message to stderr and returns
 *           False. When returning False, this function expects the caller to
 *           terminate the executable.
 *
 * \return   True if the command line arguments were successfully parsed, False
 *           otherwise. If False is returned, the data structure may be only
 *           partially initialized and memory only partially allocated.
 */
static bool arguments_init_parse_command_line(
  Generate_Distribution_Arguments * const arguments,
  const int argc,
  char ** argv)
{
  /* Initialize the arguments data structure. */
  arguments->selection_method = SELECTION_METHOD_UNDEFINED;
  arguments->tradeoff_method = TRADEOFF_METHOD_UNDEFINED;
  arguments->sigma_method = SIGMA_METHOD_UNDEFINED;
  arguments->probability_estimate = PROBABILITY_ESTIMATE_UNDEFINED;
  arguments->dimension = DIMENSION_UNDEFINED;
  arguments->entries = NULL;
  arguments->count = 0;

  mpz_init(arguments->explicit_d);
  mpz_set_ui(arguments->explicit_d, 0);
  mpz_init(arguments->explicit_r);
  mpz_set_ui(arguments->explicit_r, 0);

  /* Iterate over the command line arguments. */
  int i;

  for (i = 1; i < argc; i++) {
    /* Parse the selection method. */
    Selection_Method selection_method = SELECTION_METHOD_UNDEFINED;

    if (0 == strcmp(argv[i], "-det")) {
      selection_method = SELECTION_METHOD_DETERMINISTIC;
    } else if (0 == strcmp(argv[i], "-rnd")) {
      selection_method = SELECTION_METHOD_RANDOM;
    } else if (0 == strcmp(argv[i], "-exp")) {
      selection_method = SELECTION_METHOD_EXPLICIT;
    }

    if (SELECTION_METHOD_UNDEFINED != selection_method) {
      /* Check that a selection method has not already been specified. */
      if (SELECTION_METHOD_UNDEFINED != (arguments->selection_method)) {
        fprintf(stderr, "Error: The selection method cannot be twice "
          "specified.\n");
        return FALSE;
      }

      /* Store the selection method. */
      arguments->selection_method = selection_method;

      if (SELECTION_METHOD_EXPLICIT == selection_method) {
        if ((i + 2) >= argc) {
          fprintf(stderr, "Error: Expected <d> <r> to follow after -exp.\n");
          return FALSE;
        }

        if (0 != mpz_set_str(arguments->explicit_d, argv[i + 1], 10)) {
          fprintf(stderr, "Error: Failed to parse <d> after -exp as an "
            "integer.\n");
          return FALSE;
        }

        if (0 != mpz_set_str(arguments->explicit_r, argv[i + 2], 10)) {
          fprintf(stderr, "Error: Failed to parse <r> after -exp as an "
            "integer.\n");
          return FALSE;
        }

        /* Check that d > 0. */
        if (mpz_cmp_ui(arguments->explicit_d, 0) <= 0) {
          fprintf(stderr, "Error: The value of <d> after -exp must be "
            "positive.\n");
          return FALSE;
        }

        /* Check that d < r. Then 0 < d < r by the above check. */
        if (mpz_cmp(arguments->explicit_d, arguments->explicit_r) >= 0) {
          fprintf(stderr, "Error: The value of <d> after -exp must be strictly "
            "less than <r>.\n");
          return FALSE;
        }

        i += 2;
      }

      continue;
    }

    /* Parse the tradeoff method. */
    Tradeoff_Method tradeoff_method = TRADEOFF_METHOD_UNDEFINED;

    if (0 == strcmp(argv[i], "-s")) {
      tradeoff_method = TRADEOFF_METHOD_FACTOR;
    } else if (0 == strcmp(argv[i], "-l")) {
      tradeoff_method = TRADEOFF_METHOD_EXPLICIT;
    }

    if (TRADEOFF_METHOD_UNDEFINED != tradeoff_method) {
      /* Check that a tradeoff method has not already been specified. */
      if (TRADEOFF_METHOD_UNDEFINED != (arguments->tradeoff_method)) {
        fprintf(stderr, "Error: The tradeoff method cannot be twice "
          "specified.\n");
        return FALSE;
      }

      /* Store the tradeoff method. */
      arguments->tradeoff_method = tradeoff_method;

      continue;
    }

    /* Parse the sigma selection method. */
    Sigma_Selection_Method sigma_method = SIGMA_METHOD_UNDEFINED;

    if (0 == strcmp(argv[i], "-sigma-heuristic")) {
      sigma_method = SIGMA_METHOD_HEURISTIC;
    } else if (0 == strcmp(argv[i], "-sigma-optimal")) {
      sigma_method = SIGMA_METHOD_OPTIMAL;
    }

    if (SIGMA_METHOD_UNDEFINED != sigma_method) {
      /* Check that a sigma selection method has not already been specified. */
      if (SIGMA_METHOD_UNDEFINED != (arguments->sigma_method)) {
        fprintf(stderr, "Error: The sigma selection method cannot be twice "
          "specified.\n");
        return FALSE;
      }

      /* Store the sigma method. */
      arguments->sigma_method = sigma_method;

      continue;
    }

    /* Parse the probability estimate. */
    Probability_Estimate probability_estimate = PROBABILITY_ESTIMATE_UNDEFINED;

    if (0 == strcmp(argv[i], "-approx-with-error-bound")) {
      probability_estimate = PROBABILITY_ESTIMATE_BOUNDED_ERROR;
    } else if (0 == strcmp(argv[i], "-approx-quick")) {
      probability_estimate = PROBABILITY_ESTIMATE_QUICK;
    }

    if (PROBABILITY_ESTIMATE_UNDEFINED != probability_estimate) {
      /* Check that a probability estimate has not already been specified. */
      if (PROBABILITY_ESTIMATE_UNDEFINED != (arguments->probability_estimate)) {
        fprintf(stderr, "Error: The probability estimate method cannot be "
          "twice specified.\n");
        return FALSE;
      }

      /* Store the probability estimate. */
      arguments->probability_estimate = probability_estimate;

      continue;
    }

    /* Parse the fixed dimension. */
    if (0 == strcmp(argv[i], "-dim")) {
      /* Check that a dimension has not already been specified. */
      if (DIMENSION_UNDEFINED != arguments->dimension) {
        fprintf(stderr, "Error: The dimension cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr, "Error: Expected <dimension> to follow after -dim.\n");
        return FALSE;
      }

      const int x = atoi(argv[i+1]);
      if ((x < 16) || (x > 8192) || (!is_pow2((uint32_t)x))) {
        fprintf(stderr, "Error: The <dimension> passed to -dim must be a power "
          "of two on the interval [16, 8192].\n");
        return FALSE;
      }

      /* Store the dimension. */
      arguments->dimension = (uint32_t)x;

      i++;

      continue;
    } else if (0 == strcmp(argv[i], "-heuristic-dim")) {
      /* Check that a dimension has not already been specified. */
      if (DIMENSION_UNDEFINED != arguments->dimension) {
        fprintf(stderr, "Error: The dimension cannot be twice specified.\n");
        return FALSE;
      }

      /* Store the dimension. */
      arguments->dimension = DIMENSION_HEURISTIC;

      continue;
    }

    /* Stop parsing. */
    break;
  }

  /* Set default parameters if arguments where not explicitly specified. */
  if (DIMENSION_UNDEFINED == arguments->dimension) {
    arguments->dimension = DIMENSION_HEURISTIC;
  }

  if (SELECTION_METHOD_UNDEFINED == arguments->selection_method) {
    arguments->selection_method = SELECTION_METHOD_DETERMINISTIC;
  }

  if (TRADEOFF_METHOD_UNDEFINED == arguments->tradeoff_method) {
    arguments->tradeoff_method = TRADEOFF_METHOD_FACTOR;
  }

  if (PROBABILITY_ESTIMATE_UNDEFINED == arguments->probability_estimate) {
    arguments->probability_estimate = PROBABILITY_ESTIMATE_BOUNDED_ERROR;
  }

  if (PROBABILITY_ESTIMATE_BOUNDED_ERROR == arguments->probability_estimate) {
    if (SIGMA_METHOD_UNDEFINED == arguments->sigma_method) {
      arguments->sigma_method = SIGMA_METHOD_HEURISTIC;
    }
  } else {
    if (SIGMA_METHOD_UNDEFINED != arguments->sigma_method) {
      fprintf(stderr, "Error: The flags -sigma-heuristic and "
        "-sigma-optimal are only compatible with "
          "-probability-with-error-bound.\n");
    }
  }

  /* Parse tuples { <m> <s> } or { <m> <l> }. */
  if (((argc - i) <= 0) || (0 != ((argc - i) % 2))) {
    fprintf(stderr, "Error: Incorrect command line arguments; expected tuples "
      "but found an odd number of arguments.\n");
    return FALSE;
  }

  arguments->count = (uint32_t)((argc - i) / 2);

  arguments->entries =
    (Generate_Distribution_Arguments_Entry *)malloc(
      (arguments->count) * sizeof(Generate_Distribution_Arguments_Entry));
  if (NULL == arguments->entries) {
    fprintf(stderr, "Error: Failed to allocate memory for argument entries.\n");
    return FALSE;
  }

  for (uint32_t j = 0; j < arguments->count; j++) {
    mpz_init(arguments->entries[j].d);
    mpz_set_ui(arguments->entries[j].d, 0);

    mpz_init(arguments->entries[j].r);
    mpz_set_ui(arguments->entries[j].r, 0);
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

  /* Iterate over the tuples. */
  for (uint32_t j = 0; j < arguments->count; j++, i += 2) {
    const int m = atoi(argv[i]);
    if (m <= 1) {
      fprintf(stderr, "Error: Failed to parse <m>.\n");
      return FALSE;
    }

    arguments->entries[j].m = (uint32_t)m;

    if (TRADEOFF_METHOD_FACTOR == arguments->tradeoff_method) {
      const int s = atoi(argv[i + 1]);
      if (s <= 0) {
        fprintf(stderr, "Error: Failed to parse <s>.\n");
        return FALSE;
      }

      arguments->entries[j].s = (uint32_t)s;
      arguments->entries[j].l = (uint32_t)ceil((double)m / (double)s);
    } else {
      const int l = atoi(argv[i + 1]);
      if (l <= 0) {
        fprintf(stderr, "Error: Failed to parse <l>.\n");
        return FALSE;
      }

      arguments->entries[j].s = 0;
      arguments->entries[j].l = (uint32_t)l;
    }

    char suffix[MAX_SIZE_PATH_BUFFER];
    if (TRADEOFF_METHOD_FACTOR == arguments->tradeoff_method) {
      safe_snprintf(suffix, MAX_SIZE_PATH_BUFFER, "m-%u-s-%u",
        arguments->entries[j].m, arguments->entries[j].s);
    } else if (TRADEOFF_METHOD_EXPLICIT == arguments->tradeoff_method) {
      safe_snprintf(suffix, MAX_SIZE_PATH_BUFFER, "m-%u-l-%u",
        arguments->entries[j].m, arguments->entries[j].l);
    } else {
      critical("arguments_init_parse_command_line(): "
        "Undefined tradeoff method.");
    }

    char sigma_prefix[MAX_SIZE_PATH_BUFFER];
    if (PROBABILITY_ESTIMATE_BOUNDED_ERROR == arguments->probability_estimate) {
      if (SIGMA_METHOD_OPTIMAL == arguments->sigma_method) {
        safe_snprintf(sigma_prefix, MAX_SIZE_PATH_BUFFER, "sigma-optimal");
      } else if (SIGMA_METHOD_HEURISTIC == arguments->sigma_method) {
        safe_snprintf(sigma_prefix, MAX_SIZE_PATH_BUFFER, "sigma-heuristic");
      } else {
        critical("arguments_init_parse_command_line(): "
          "Undefined sigma selection method.");
      }
    } else if (PROBABILITY_ESTIMATE_QUICK == arguments->probability_estimate) {
      safe_snprintf(sigma_prefix, MAX_SIZE_PATH_BUFFER, "quick");
    } else {
      critical("arguments_init_parse_command_line(): "
        "Undefined probability estimate.");
    }

    char dim_prefix[MAX_SIZE_PATH_BUFFER];
    if (DIMENSION_HEURISTIC == arguments->dimension) {
      safe_snprintf(dim_prefix, MAX_SIZE_PATH_BUFFER, "dim-heuristic");
    } else if (DIMENSION_UNDEFINED != arguments->dimension){
      safe_snprintf(dim_prefix, MAX_SIZE_PATH_BUFFER, "dim-%u",
        arguments->dimension);
    } else {
      critical("arguments_init_parse_command_line(): Undefined dimension.");
    }

    /* Set the value. */
    if (SELECTION_METHOD_DETERMINISTIC == arguments->selection_method) {
      if (m > 8192) {
        fprintf(stderr,
          "Error: When using the -det flag it is required that m <= 8192.\n");
        return FALSE;
      }

      parameters_selection_deterministic_d_r(
          arguments->entries[j].d,
          arguments->entries[j].r,
          (uint32_t)m);

      safe_snprintf(arguments->entries[j].name, MAX_SIZE_PATH_BUFFER,
        "distribution-det-%s-%s-%s", dim_prefix, sigma_prefix, suffix);
    } else if (SELECTION_METHOD_RANDOM == arguments->selection_method) {
      if (j > 0) { /* Check if we should copy the previous value */
        if ((arguments->entries[j-1].m) == (arguments->entries[j].m)) {
          mpz_set(arguments->entries[j].d, arguments->entries[j-1].d);
          mpz_set(arguments->entries[j].r, arguments->entries[j-1].r);
        }
      }

      /* If we did not copy the previous value, select the value at random. */
      if (0 == mpz_cmp_ui(arguments->entries[j].d, 0)) {
        parameters_selection_random_d_and_r(
            arguments->entries[j].d,
            arguments->entries[j].r,
            arguments->entries[j].m);
      }

      safe_snprintf(arguments->entries[j].name, MAX_SIZE_PATH_BUFFER,
        "distribution-rnd-%s-%s-%s", dim_prefix, sigma_prefix, suffix);
    } else if (SELECTION_METHOD_EXPLICIT == arguments->selection_method) {
      /* Check that the explicitly specified value of r is less than 2^m. */
      mpz_setbit(arguments->entries[j].r, arguments->entries[j].m);

      if (mpz_cmp(arguments->explicit_r, arguments->entries[j].r) >= 0) {
        fprintf(stderr, "Error: The explicitly specified <r> with -exp is "
          "greater than or equal to 2^m for at least one <m>.\n");
        return FALSE;
      }

      /* Check that the specified value of r is greater than 2^(m-1). */
      mpz_set_ui(arguments->entries[j].r, 0);
      mpz_setbit(arguments->entries[j].r, arguments->entries[j].m - 1);

      if (mpz_cmp(arguments->explicit_r, arguments->entries[j].r) <= 0)
      {
        fprintf(stderr, "Error: The explicitly specified <r> with -exp is "
          "less than or equal to 2^(m-1) for at least one <m>.\n");
        return FALSE;
      }

      /* Note: It has already been checked that 0 < d < r when parsing the
       *       command line arguments so we need not do it again. */

      mpz_set(arguments->entries[j].d, arguments->explicit_d);
      mpz_set(arguments->entries[j].r, arguments->explicit_r);

      safe_snprintf(arguments->entries[j].name, MAX_SIZE_PATH_BUFFER,
        "distribution-exp-%s-%s-%s", dim_prefix, sigma_prefix, suffix);
    } else {
      fprintf(stderr, "Error: Internal error: Unknown selection method.\n");
      return FALSE;
    }

    /* Construct the path to the output file. */
    safe_snprintf(
      arguments->entries[j].path,
      MAX_SIZE_PATH_BUFFER,
      "%s/%s.txt",
      DISTRIBUTIONS_DIRECTORY,
      arguments->entries[j].name);

    /* Check that the output file does not already exist. */
    if (0 == access(arguments->entries[j].path, F_OK)) {
      fprintf(stderr, "Error: The distribution \"%s\" already exists.\n",
        arguments->entries[j].path);
      return FALSE;
    }
  }

  /* Signal success. */
  return TRUE;
}

/*!
 * \brief   Clears an initialized command line arguments data structure.
 * 
 * \param[in, out] arguments   The argument data structure to clear.
 */
static void arguments_clear(
  Generate_Distribution_Arguments * const arguments)
{
  mpz_clear(arguments->explicit_d);
  mpz_clear(arguments->explicit_r);

  if (NULL != arguments->entries) {
    for (uint32_t i = 0; i < arguments->count; i++) {
      mpz_clear(arguments->entries[i].d);
      mpz_clear(arguments->entries[i].r);
    }
          
    free(arguments->entries);
    arguments->entries = NULL;
  }

  memset(arguments, 0, sizeof(Generate_Distribution_Arguments));
}

/*!
 * \brief   Exports human-readable information about a distribution.
 *
 * \param[in] distribution  The distribution for which to export information.
 * \param[in] name          The base name from which to construct the path.
 */
static void main_server_export_distribution_info(
  const Distribution * const distribution,
  const char * const name)
{
  FILE * file = NULL;

  char path[MAX_SIZE_PATH_BUFFER];

  safe_snprintf(path, MAX_SIZE_PATH_BUFFER,
    "%s/%s.info", DISTRIBUTIONS_DIRECTORY, name);

  printf("Exporting distribution information to \"%s\"...\n", path);
  file = fopen(path, "wb+");
  if (NULL == file) {
    critical("main_server_export_distribution_info(): "
      "Failed to open \"%s\" for writing.", path);
  }

  distribution_export_info(distribution, file);

  fclose(file);
  file = NULL;
}

/*!
 * \brief   Collapses a two-dimensional distribution to two margin
 *          distributions and exports these distributions.
 *
 * \param[in] distribution  The distribution to collapse and export.
 * \param[in] name          The base name from which to construct the path.
 */
static void main_server_export_collapsed_distributions(
  const Distribution * const distribution,
  const char * const name)
{
  /* Export the collapsed distribution with respect to alpha_d. */
  Linear_Distribution linear_distribution_d;
  linear_distribution_init_collapse_d(&linear_distribution_d, distribution);
  linear_distribution_sort_slices(&linear_distribution_d);

  FILE * file = NULL;

  char path[MAX_SIZE_PATH_BUFFER];

  safe_snprintf(path, MAX_SIZE_PATH_BUFFER, "%s/collapsed-d-%s.txt",
    DISTRIBUTIONS_DIRECTORY, name);

  printf("Exporting collapsed distribution to \"%s\"...\n", path);
  file = fopen(path, "wb+");
  if (NULL == file) {
    critical("main_server_export_collapsed_distributions(): "
      "Failed to open \"%s\" for writing.", path);
  }

  linear_distribution_export(&linear_distribution_d, file);

  fclose(file);
  file = NULL;

  /* Export the collapsed distribution with respect to alpha_r. */
  Linear_Distribution linear_distribution_r;
  linear_distribution_init_collapse_r(&linear_distribution_r, distribution);
  linear_distribution_sort_slices(&linear_distribution_r);

  safe_snprintf(path, MAX_SIZE_PATH_BUFFER, "%s/collapsed-r-%s.txt",
    DISTRIBUTIONS_DIRECTORY, name);

  printf("Exporting collapsed distribution to \"%s\"...\n", path);

  file = fopen(path, "wb+");
  if (NULL == file) {
    critical("main_server_export_collapsed_distributions(): "
      "Failed to open \"%s\" for writing.", path);
  }

  linear_distribution_export(&linear_distribution_r, file);

  fclose(file);
  file = NULL;

  /* Clear memory. */
  linear_distribution_clear(&linear_distribution_d);
  linear_distribution_clear(&linear_distribution_r);
}

/*!
 * \brief   Exports the computed probability distribution to file.
 *
 * This function is intended to run in a separate thread.
 *
 * \param[in, out] ptr    A pointer to a Distribution_Export_Job data structure
 *                        that describes the export job.
 *
 * \remark   The export job and the distribution it contains will be cleared
 *           and deallocated by this function.
 */
static void * main_server_export_distribution(void * ptr)
{
  /* Map the input to data structures. */
  Distribution_Export_Job * job = (Distribution_Export_Job *)ptr;

  Distribution * distribution = job->distribution;

  const char * path = job->path;
  const char * name = job->name;

  /* Sort the slices in the distribution. */
  printf("Sorting the slices in the distribution...\n");
  distribution_sort_slices(distribution);

  /* Export distribution information. */
  main_server_export_distribution_info(distribution, name);

  /* Collapse and export collapsed distributions. */
  main_server_export_collapsed_distributions(distribution, name);

  /* Export the two-dimensional distribution. */
  printf("Exporting the distribution to \"%s\"...\n", path);

  FILE * file = fopen(path, "wb+");
  if (NULL == file) {
    critical("main_server_export_distribution(): Failed to open \"%s\".", path);
  }
  
  if (distribution_is_filtered(distribution)) {
    distribution_export_clear_dealloc(distribution, file);
    printf("Finished exporting the distribution to \"%s\".\n", path);
  } else {
    distribution_export(distribution, file);
    printf("Finished exporting the distribution to \"%s\".\n", path);

    fclose(file);
    file = NULL;

    /* Filter the slices in the distribution. */
    printf("Filtering the slices in the distribution...\n");
    distribution_filter_slices(distribution);

    printf("Sorting the slices in the distribution...\n");
    distribution_sort_slices(distribution); /* Not strictly necessary. */
  
    /* Setup the name of the filtered distribution. */
    char filtered_name[MAX_SIZE_PATH_BUFFER];
    safe_snprintf(filtered_name, MAX_SIZE_PATH_BUFFER, "filtered-%s", name);

    /* Export distribution information. */
    main_server_export_distribution_info(distribution, filtered_name);

    /* Collapse and export collapsed filtered distributions. */
    main_server_export_collapsed_distributions(distribution, filtered_name);

    /* Setup the path to the filtered distribution. */
    char filtered_path[MAX_SIZE_PATH_BUFFER];
    safe_snprintf(filtered_path, MAX_SIZE_PATH_BUFFER, "%s/%s.txt",
      DISTRIBUTIONS_DIRECTORY, filtered_name);

    /* Export the filtered two-dimensional distribution. */
    file = fopen(filtered_path, "wb+");
    if (NULL == file) {
      critical("main_server_export_distribution(): Failed to open \"%s\".",
        filtered_path);
    }

    distribution_export_clear_dealloc(distribution, file);

    printf("Finished exporting the distribution to \"%s\".\n", filtered_path);
  }

  fclose(file);
  file = NULL;

  free(job);
  job = NULL;
  ptr = NULL;

  /* Finished. */
  return NULL;
}

/*!
 * \brief   The main function on the server node.
 *
 * This function is called once by main() for each distribution to generate.
 *
 * \param[in] arguments   The parsed command line arguments.
 * \param[in] entry       The \<m\> \<s\> or \<m\> \<l\> tuple to process.
 * \param[in] mpi_size    The number of nodes.
 * \param[in] pool        The thread pool from which to spawn threads for
 *                        export operations that run in the background.
 */
static void main_server(
  const Generate_Distribution_Arguments * const arguments,
  const Generate_Distribution_Arguments_Entry * const entry,
  int mpi_size,
  Thread_Pool * const pool)
{
  /* Setup the parameters. */
  const uint32_t t = 30;

  Parameters parameters;
  parameters_init(&parameters);

  if (TRADEOFF_METHOD_FACTOR == arguments->tradeoff_method) {
    parameters_explicit_m_s(
      &parameters,
      entry->d,
      entry->r,
      entry->m,
      entry->s,
      t);
  } else {
    parameters_explicit_m_l(
      &parameters,
      entry->d,
      entry->r,
      entry->m,
      entry->l,
      t);
  }

  /* Send broadcast of the parameters. */
  parameters_bcast_send(&parameters, MPI_RANK_ROOT);

  /* Setup the distribution enumerator. */
  const bool mirrored = TRUE;

  Distribution_Enumerator enumerator;
  distribution_enumerator_init(&enumerator, &parameters, mirrored);

  /* Get the required distribution capacity in slices. */
  uint32_t capacity = distribution_enumerator_count(&enumerator);

  if (mirrored) {
    capacity *= 2;
  }

  /* Setup the distribution. */
  Distribution * distribution = distribution_alloc();
  distribution_init(distribution, &parameters, capacity);

  /* Keep track of the number of nodes remaining. */
  int mpi_nodes_remaining = mpi_size;

  while (mpi_nodes_remaining > 1) {
    /* Listen for a notification. */
    MPI_Status status;

    uint32_t notification;

    if (MPI_SUCCESS != MPI_Recv(
        &notification,
        1, /* count */
        MPI_UNSIGNED,
        MPI_ANY_SOURCE, /* source */
        MPI_TAG_NOTIFY,
        MPI_COMM_WORLD,
        &status))
    {
      critical("main_server(): Failed to receive notification.");
    }

    if (MPI_NOTIFY_READY == notification) {
      int32_t min_log_alpha[2];

      bool found = distribution_enumerator_next(
                      &min_log_alpha[0],
                      &min_log_alpha[1],
                      &enumerator);
      if (found) {
        printf("Processing slice: %u / %u (%d, %d)\n",
          enumerator.offset,
          enumerator.count,
          min_log_alpha[0],
          min_log_alpha[1]);

        uint32_t job = MPI_JOB_PROCESS_SLICE;

        if (MPI_SUCCESS != MPI_Send(
            &job,
            1, /* count */
            MPI_UNSIGNED,
            status.MPI_SOURCE, /* destination */
            MPI_TAG_JOB,
            MPI_COMM_WORLD))
        {
          critical("main_server(): Failed to send job.");
        }

        if (MPI_SUCCESS != MPI_Send(
            min_log_alpha,
            2, /* count */
            MPI_INT,
            status.MPI_SOURCE, /* destination */
            MPI_TAG_SLICE_MIN_LOG_ALPHA,
            MPI_COMM_WORLD))
        {
          critical("main_server(): "
            "Failed to send min_log_alpha_d and min_log_alpha_r.");
        }
      } else {
        uint32_t job = MPI_JOB_STOP;

        if (MPI_SUCCESS != MPI_Send(
            &job,
            1, /* count */
            MPI_UNSIGNED,
            status.MPI_SOURCE, /* destination */
            MPI_TAG_JOB,
            MPI_COMM_WORLD))
        {
          critical("main_server(): Failed to send job.");
        }

        mpi_nodes_remaining--;

        printf("Stopping node %d with %u node(s) remaining...\n",
          status.MPI_SOURCE,
          mpi_nodes_remaining);
      }

    } else if (MPI_NOTIFY_SLICE_DONE == notification) {
      /* The node has finished processing a slice. Receive the slice. */
      Distribution_Slice * slice = distribution_slice_alloc();

      distribution_slice_init_recv(slice, status.MPI_SOURCE);

      distribution_insert_slice(distribution, slice);

      if (mirrored) {
        Distribution_Slice * mirrored_slice = distribution_slice_alloc();
        distribution_slice_init_copy(mirrored_slice, slice);

        mirrored_slice->min_log_alpha_d *= -1; /* Mirror alpha_d. */
        mirrored_slice->min_log_alpha_r *= -1; /* Mirror alpha_r. */

        mirrored_slice->flags |= SLICE_FLAGS_MIRRORED;

        distribution_insert_slice(distribution, mirrored_slice);
      }

      /* Print some feedback to the user. */
      printf("Received slice: (%d, %d)\n",
        slice->min_log_alpha_d,
        slice->min_log_alpha_r);

      printf("Slice total probability is: %.24LG\n", slice->total_probability);

      if (PROBABILITY_ESTIMATE_BOUNDED_ERROR ==
        arguments->probability_estimate)
      {
          printf("Slice total error is: %.24LG\n", slice->total_error);
      }

      printf("Slice dimension is: %u\n", slice->dimension);

      printf("Total probability is now: %.24LG\n",
        distribution->total_probability);

      if (PROBABILITY_ESTIMATE_BOUNDED_ERROR ==
        arguments->probability_estimate)
      {
        printf("Total error is now: %.24LG\n", distribution->total_error);
      }

      /* Note: The slices are inserted by reference and will be deallocated by
       *       the call to distribution_dealloc(). They must not be deallocated
       *       here as doing so would lead to memory corruption. */
    } else if (MPI_NOTIFY_SLICE_SKIP == notification) {
      /* Ignore. */
    } else {
      critical("main_server(): Unknown notification.");
    }
  }

  /* Clear memory. */
  distribution_enumerator_clear(&enumerator);

  parameters_clear(&parameters);

  /* Setup an export job. */
  Distribution_Export_Job * export_job =
    (Distribution_Export_Job *)malloc(sizeof(Distribution_Export_Job));
  if (NULL == export_job) {
    critical("main_server(): Failed to allocate memory.");
  }

  export_job->distribution = distribution;

  safe_strlcpy(export_job->path, entry->path, MAX_SIZE_PATH_BUFFER);
  safe_strlcpy(export_job->name, entry->name, MAX_SIZE_PATH_BUFFER);

  /* Spawn a thread to sort, export and deallocate the distribution. */
  thread_pool_spawn(pool, main_server_export_distribution, export_job);
}

/*!
 * \brief Computes a slice using a given probability estimate and sigma method.
 *
 * This function always compute the slice using Simpson's method with
 * Richardson extrapolation.
 *
 * \param[in, out] slice        The slice to compute.
 * \param[in, out] parameters   The parameters for which to compute the slice.
 *
 * \param[in, out] min_log_alpha_d  The signed logarithmic alpha_d-coordinate of
 *                                  the slice to compute.
 * \param[in, out] min_log_alpha_r  The signed logarithmic alpha_r-coordinate of
 *                                  the slice to compute.
 *
 * \param[in] probability_estimate  The probability estimate to use.
 * \param[in] sigma_method          The sigma method to use.
 */
static void main_client_compute_slice(
  Distribution_Slice * const slice,
  const Parameters * const parameters,
  const int32_t min_log_alpha_d,
  const int32_t min_log_alpha_r,
  const Probability_Estimate probability_estimate,
  const Sigma_Selection_Method sigma_method)
{
  Distribution_Slice_Compute_Method method = 
    DISTRIBUTION_SLICE_COMPUTE_METHOD_QUICK;

  if (PROBABILITY_ESTIMATE_BOUNDED_ERROR == probability_estimate) {
    if (SIGMA_METHOD_OPTIMAL == sigma_method) {
      method = DISTRIBUTION_SLICE_COMPUTE_METHOD_OPTIMAL_LOCAL_SIGMA;
    } else if (SIGMA_METHOD_HEURISTIC == sigma_method) {
      method = DISTRIBUTION_SLICE_COMPUTE_METHOD_HEURISTIC_SIGMA;
    } else {
        critical("main_client_compute_slice(): Internal error: "
          "Unknown sigma method.");
    }
  } else if (PROBABILITY_ESTIMATE_QUICK != probability_estimate) {
    critical("main_client_compute_slice(): Internal error: "
      "Unknown probability estimate.");
  }

  distribution_slice_compute_richardson(
    slice,
    parameters,
    method,
    min_log_alpha_d,
    min_log_alpha_r);
}

/*!
 * \brief   The main function on the client node.
 *
 * This function is called once by main() for each distribution to generate.
 */
static void main_client(
  const Probability_Estimate probability_estimate,
  const Sigma_Selection_Method sigma_method,
  const uint32_t dimension)
{
  Parameters parameters;
  parameters_init(&parameters);

  /* Receive broadcast of the distribution parameters. */
  parameters_bcast_recv(&parameters, MPI_RANK_ROOT);

  while (TRUE) {
    /* Notify the server that we are ready to receive a job. */
    uint32_t notification = MPI_NOTIFY_READY;

    if (MPI_SUCCESS != MPI_Send(
        &notification,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_TAG_NOTIFY,
        MPI_COMM_WORLD))
    {
      critical("main_client(): Failed to send notification.");
    }

    /* Receive the response. */
    uint32_t job;

    MPI_Status status;

    if (MPI_SUCCESS != MPI_Recv(
        &job,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_TAG_JOB,
        MPI_COMM_WORLD,
        &status))
    {
      critical("main_client(): Failed to receive job.");
    }

    if (MPI_JOB_STOP == job) {
      break;
    } else if (MPI_JOB_PROCESS_SLICE != job) {
      critical("main_client(): Unknown job (%u).", job);
    }

    /* Receive the slice region. */
    uint32_t min_log_alpha[2];

    if (MPI_SUCCESS != MPI_Recv(
        &min_log_alpha,
        2, /* count */
        MPI_INT,
        MPI_RANK_ROOT,
        MPI_TAG_SLICE_MIN_LOG_ALPHA,
        MPI_COMM_WORLD,
        &status))
    {
      critical("main_client(): "
        "Failed to receive min_log_alpha_d and min_log_alpha_r.");
    }

    const int32_t min_log_alpha_d = min_log_alpha[0];
    const int32_t min_log_alpha_r = min_log_alpha[1];

    /* Skip this slice if the region is out of bounds. */
    const int32_t max_alpha =
      max_i(abs_i(min_log_alpha_d), abs_i(min_log_alpha_r));

    const int32_t m = (int32_t)(parameters.m);

    if (max_alpha > m + 10) {
      /* Notify the server that we are skipping this slice. */
      notification = MPI_NOTIFY_SLICE_SKIP;

      if (MPI_SUCCESS != MPI_Send(
          &notification,
          1, /* count */
          MPI_UNSIGNED,
          MPI_RANK_ROOT,
          MPI_TAG_NOTIFY,
          MPI_COMM_WORLD))
      {
        critical("main_client(): Failed to send notification.");
      }

      continue;
    }

    /* Compute the slice. */
    Distribution_Slice slice;

    if (DIMENSION_HEURISTIC == dimension) {
      /* Establish the initial dimension requirements heuristically. */
      uint32_t required_dimension = 256;

      /* The tail along the diagonal requires high resolution. */
      if (abs_i(min_log_alpha_d - min_log_alpha_r) <= 1) {
        if (sgn_i(min_log_alpha_d) == sgn_i(min_log_alpha_r)) {
          if (max_alpha >= m + 3) {
            if (required_dimension < 1024) {
              required_dimension = 1024;
            }
          } else if (max_alpha >= m) {
            if (required_dimension < 512) {
              required_dimension = 512;
            }
          }
        }
      }

      /* Compute the slice for the selected dimension. Note that we divide the 
       * dimension by two below as it will be doubled in the Richardson
       * extrapolation when the slice is computed. */
      distribution_slice_init(&slice, required_dimension / 2);

      main_client_compute_slice(
        &slice,
        &parameters,
        min_log_alpha_d,
        min_log_alpha_r,
        probability_estimate,
        sigma_method);

      /* Check if the dimension should be incremented given the result. */

      /* To eliminate noise in the tail we need high resolution not only for
       * the tail itself but also the low probability area around it; otherwise
       * these areas will erroneously sum to a significant probability. The 
       * below inferred and fairly coarse heuristic fixes such problems. */
      
      bool dimension_updated = FALSE;
  
      if (slice.total_probability >= 1e-7) {
        if (max_alpha >= m) {
          if (required_dimension < 512) {
            required_dimension = 512;
            dimension_updated = TRUE;
          }
        }
      }

      if (slice.total_probability >= 1e-10) {
        if (max_alpha >= m + 10) {
          if (required_dimension < 1024) {
            required_dimension = 1024;
            dimension_updated = TRUE;
          }
        }
      }

      /* Recompute the slice if the required dimension was updated above. */
      if (dimension_updated) {
        distribution_slice_clear(&slice);
        distribution_slice_init(&slice, required_dimension / 2);

        main_client_compute_slice(
          &slice,
          &parameters,
          min_log_alpha_d,
          min_log_alpha_r,
          probability_estimate,
          sigma_method);
      }
    } else {
      /* Compute slice with fixed dimension. */
      distribution_slice_init(&slice, dimension / 2);

      main_client_compute_slice(
        &slice,
        &parameters,
        min_log_alpha_d,
        min_log_alpha_r,
        probability_estimate,
        sigma_method);
    }

    if ((0 == MIN_SLICE_PROBABILITY) ||
      (slice.total_probability >= MIN_SLICE_PROBABILITY))
    {
      /* Notify the server we are done computing the slice. */
      notification = MPI_NOTIFY_SLICE_DONE;

      if (MPI_SUCCESS != MPI_Send(
          &notification,
          1, /* count */
          MPI_UNSIGNED,
          MPI_RANK_ROOT,
          MPI_TAG_NOTIFY,
          MPI_COMM_WORLD))
      {
        critical("main_client(): Failed to send notification.");
      }

      /* Scale the slice if its dimension exceeds MAX_SLICE_DIMENSION. */
      if ((0 != MAX_SLICE_DIMENSION) &&
        (slice.dimension > MAX_SLICE_DIMENSION))
      {
        Distribution_Slice scaled_slice;
        distribution_slice_init(&scaled_slice, MAX_SLICE_DIMENSION);

        /* Scale the slice to dimension MAX_SLICE_DIMENSION. */
        distribution_slice_copy_scale(&scaled_slice, &slice);

        /* Send back the slice. */
        distribution_slice_send(&scaled_slice, MPI_RANK_ROOT);

        /* Clear memory. */
        distribution_slice_clear(&scaled_slice);
      } else {
        /* Send back the slice. */
        distribution_slice_send(&slice, MPI_RANK_ROOT);
      }
    } else {
      /* Notify the server that we are skipping this slice. */
      notification = MPI_NOTIFY_SLICE_SKIP;

      if (MPI_SUCCESS != MPI_Send(
          &notification,
          1, /* count */
          MPI_UNSIGNED,
          MPI_RANK_ROOT,
          MPI_TAG_NOTIFY,
          MPI_COMM_WORLD))
      {
        critical("main_client(): Failed to send notification.");
      }
    }

    /* Clear memory. */
    distribution_slice_clear(&slice);
  }

  /* Clear memory. */
  parameters_clear(&parameters);
}

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  fprintf(file, "Synopsis: mpirun generate_distribution \\\n"
          "   [ -det | -rnd | -exp <d> <r> ] "
            "[ -dim-heuristic | -dim <dimension> ] \\\n"
          "      [ -approx-quick | "
            "[ -sigma-heuristic | -sigma-optimal ] ] \\\n"
          "         ( [ -s ] <m> <s> { <m> <s> } | -l <m> <l> { <m> <l> } )\n");

  fprintf(file, "\n");
  fprintf(file, "Selection method for d and r: -- defaults to -det\n");
  fprintf(file,
    " -det  select d and r deterministically from Catalan's constant\n");
  fprintf(file,
    " -rnd  select r uniformly at random from (2^(m-1), 2^m)\n");
  fprintf(file,
    "       and d uniformly at random from [r/2, r)\n");
  fprintf(file,
      " -exp  explicitly set d and r to <d> and <r>\n");

  fprintf(file, "\n");
  fprintf(file,
     "Approximation method: -- defaults to approximation with error bound\n");
  fprintf(file,
    " -approx-quick             use the quick approximation for which no "
      "error\n"
    "                           bound is known (incompatible with -sigma-*)\n");
  fprintf(file,
    " -approx-with-error-bound  explicitly specify use of the error-bounded\n"
    "                           approximation involving sigma "
      "(see -sigma-*) \n");

    fprintf(file, "\n");
    fprintf(file,
       "Sigma selection method: -- defaults to heuristic selection\n");
    fprintf(file,
      " -sigma-optimal    adaptively search for the optimal sigma "
        "minimizing the\n"
      "                   error bound in the approximation\n");
    fprintf(file,
      " -sigma-heuristic  explicitly specify heuristic selection of sigma\n");

  fprintf(file, "\n");
  fprintf(file, "Dimension: -- heuristically selected by default\n");
  fprintf(file,
    " -dim <dimension>  explicitly set the slice dimension to <dimension>\n");
  fprintf(file,
    " -heuristic-dim    explicitly specifies heuristic selection\n");

  fprintf(file, "\n");
  fprintf(file,
    "Tuples <m> <s> or <m> <l>: -- one distribution is generated for each "
      "tuple\n");
  fprintf(file,
    " <m>   the length m in bits of r\n");
  fprintf(file,
    " <s>   the tradeoff factor s; used to set l = ceil(m / s)\n");
  fprintf(file,
    " <l>   the parameter l\n");

  fprintf(file,
    "\nNote: If -rnd is specified and m is kept constant for consecutive "
      "tuples\n");
  fprintf(file,
    "<m> <s> or <m> <l>, the same values of d and r will be re-used.\n");

  fprintf(file,
    "\nNote: This implementation is optimized for small to medium kappa. If ");
  fprintf(file,
    "\nkappa is very large, it may fail to capture the probability mass.\n");


  fprintf(file,
    "\nNote: Once computed, slices of dimension exceeding 256 are scaled to ");
  fprintf(file,
    "\ndimension 256 when inserted into the distribution to save space.\n");
}

/*!
 * \brief The main entry point to the generate_distribution executable.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \return Zero upon successful execution, non-zero otherwise.
 */
int main(int argc, char ** argv) {
  mpfr_set_default_prec(PRECISION);

  int provided;

  if (MPI_SUCCESS != MPI_Init_thread(
                        &argc,
                        &argv,
                        MPI_THREAD_FUNNELED, /* requested level of support */
                        &provided))
  {
    critical("main(): Failed to initialize MPI.");
  }

  int mpi_rank;
  int mpi_size;

  if (MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank)) {
    critical("main(): Failed to get the MPI world communicator rank.");
  }

  if (MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD, &mpi_size)) {
    critical("main(): Failed to get the MPI world communicator size.");
  }

  if (1 == argc) {
    if (MPI_RANK_ROOT == mpi_rank) {
      print_synopsis(stdout);
    }

    /* Exit gracefully. */
    MPI_Finalize();
    return 0;
  }

  if (mpi_size < 2) {
    if (MPI_RANK_ROOT == mpi_rank) {
      fprintf(stderr, "Incorrect environment: "
        "At least two MPI processors are required.\n");
    }

    /* Exit gracefully. */
    MPI_Finalize();
    return 0;
  }

  /* Parse the command line arguments and signal errors in a graceful way. */
  uint32_t result;

  Generate_Distribution_Arguments arguments;

  if (MPI_RANK_ROOT == mpi_rank) {
    result =
      (uint32_t)arguments_init_parse_command_line(&arguments, argc, argv);
  }

  if (MPI_SUCCESS != 
    MPI_Bcast(&result, 1, MPI_UNSIGNED, MPI_RANK_ROOT, MPI_COMM_WORLD))
  {
    critical("main(): "
      "Failed to broadcast the result of parsing the arguments.");
  };

  if (!result) {
    /* Exit gracefully if an error occurred. */
    MPI_Finalize();
    return 0;
  }

  /* Begin the generation process. */
  if (MPI_RANK_ROOT == mpi_rank) {
    printf("Number of MPI processors: %u\n", mpi_size);

    /* Report back information on the level of MPI thread support provided. */
    switch (provided) {
      case MPI_THREAD_SINGLE:
        printf("MPI level of thread support provided: "
          "MPI_THREAD_SINGLE\n");
        break;

      case MPI_THREAD_FUNNELED:
        printf("MPI level of thread support provided: "
          "MPI_THREAD_FUNNELED\n");
        break;

      case MPI_THREAD_SERIALIZED:
        printf("MPI level of thread support provided: "
          "MPI_THREAD_SERIALIZED\n");
        break;

      case MPI_THREAD_MULTIPLE:
        printf("MPI level of thread support provided: "
          "MPI_THREAD_MULTIPLE\n");
        break;

      default:
        printf("Level of thread support provided: "
          "Unknown (%d)\n", provided);
        break;
    }

    printf("\n");

    /* Broadcast the arguments count. */
    if (MPI_SUCCESS !=
      MPI_Bcast(
        &(arguments.count),
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): Failed to broadcast the arguments count.");
    }

    /* Broadcast the dimension. */
    if (MPI_SUCCESS !=
      MPI_Bcast(
        &(arguments.dimension),
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): Failed to broadcast the dimension.");
    }

    /* Broadcast the sigma method. */
    uint32_t sigma_method = arguments.sigma_method;

    if (MPI_SUCCESS !=
      MPI_Bcast(
        &sigma_method,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): Failed to broadcast the sigma method.");
    }

    /* Broadcast the probability estimate. */
    uint32_t probability_estimate = arguments.probability_estimate;

    if (MPI_SUCCESS !=
      MPI_Bcast(
        &probability_estimate,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): Failed to broadcast the probability estimate.");
    }

    /* Process each distribution. */
    Thread_Pool pool;
    thread_pool_init(&pool);

    for (uint32_t i = 0; i < arguments.count; i++) {
      main_server(&arguments, &(arguments.entries[i]), mpi_size, &pool);
    }

    printf("Waiting for all export operations to finish...\n");

    thread_pool_join(&pool);
    thread_pool_clear(&pool);

    arguments_clear(&arguments);
  } else {
    /* Broadcast the arguments count. */
    uint32_t count;

    if (MPI_SUCCESS !=
      MPI_Bcast(
        &count,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): Failed to broadcast the arguments count.");
    }

    /* Broadcast the dimension. */
    uint32_t dimension;

    if (MPI_SUCCESS !=
      MPI_Bcast(
        &dimension,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): Failed to broadcast the dimension.");
    }

    /* Broadcast the sigma method. */
    uint32_t tmp;

    if (MPI_SUCCESS !=
      MPI_Bcast(
        &tmp,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): Failed to broadcast the sigma method.");
    }

    Sigma_Selection_Method sigma_method = (Sigma_Selection_Method)tmp;

    /* Broadcast the probability estimate. */
    if (MPI_SUCCESS !=
      MPI_Bcast(
        &tmp,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): Failed to broadcast the probability estimate.");
    }

    Probability_Estimate probability_estimate = (Probability_Estimate)tmp;

    for (uint32_t i = 0; i < count; i++) {
      main_client(
        probability_estimate,
        sigma_method,
        dimension);
    }
  }

  if (MPI_SUCCESS != MPI_Finalize()) {
    critical("main(): Failed to finalize MPI.");
  }

  /* Signal success. */
  return 0;
}
