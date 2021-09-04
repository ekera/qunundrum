/*!
 * \file    main_generate_diagonal_distribution.cpp
 * \ingroup generate_diagonal_distribution_exe
 *
 * \brief   The definition of the main entry point to the
 *          generate_diagonal_distribution executable, and of associated
 *          functions.
 */

/*!
 * \defgroup generate_diagonal_distribution_exe \
 *           The generate_diagonal_distribution executable
 * \ingroup  generate_executable
 *
 * \brief    A module for the generate_diagonal_distribution executable.
 */

#include "executables.h"
#include "executables_generate_distribution.h"

#include "common.h"
#include "diagonal_distribution.h"
#include "diagonal_distribution_enumerator.h"
#include "diagonal_distribution_slice.h"
#include "diagonal_parameters.h"
#include "errors.h"
#include "math.h"
#include "parameters_selection.h"
#include "string_utilities.h"
#include "thread_pool.h"

#include <gmp.h>
#include <mpfr.h>

#include <mpi.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <sys/stat.h>

/*!
 * \brief   A data structure representing argument entries in the form of parsed
 *          \<m\> \<sigma\> \<s\> or \<m\> \<sigma\> \<l\> tuples from the
 *          command line arguments.
 *
 * \ingroup generate_diagonal_distribution_exe
 */
typedef struct {
  /*!
   * \brief   The parameter m.
   */
  uint32_t m;

  /*!
   * \brief   The parameter sigma.
   */
  uint32_t sigma;

  /*!
   * \brief   The tradeoff factor s.
   */
  uint32_t s;

  /*!
   * \brief   The parameter l.
   */
  uint32_t l;

  /*!
   * \brief   The value of d for the selection method.
   *
   * The selection method is common for all argument tuples, see
   * Generate_Diagonal_Distribution_Arguments::selection_method.
   */
  mpz_t d;

  /*!
   * \brief   The value of r for the selection method.
   *
   * The selection method is common for all argument tuples, see
   * Generate_Diagonal_Distribution_Arguments::selection_method.
   */
  mpz_t r;

  /*!
   * \brief   The path to which to export this distribution.
   */
  char path[MAX_SIZE_PATH_BUFFER];
} Generate_Diagonal_Distribution_Arguments_Entry;

/*!
 * \brief   A data structure representing parsed command line arguments.
 *
 * \ingroup generate_diagonal_distribution_exe
 */
typedef struct {
  /*!
   * \brief   The selection method for the logarithm d or order r.
   */
  Selection_Method selection_method;

  /*!
   * \brief   The tradeoff method.
   */
  Tradeoff_Method tradeoff_method;

  /*!
   * \brief   The explicitly specified value of the logarithm d.
   *
   * If the value is not explicitly specified, this entry is set to zero.
   */
  mpz_t explicit_d;

  /*!
   * \brief   The explicitly specified value of the order r.
   *
   * If the value is not explicitly specified, this entry is set to zero.
   */
  mpz_t explicit_r;

  /*!
   * \brief   The dimension of the slices in the distribution.
   */
  uint32_t dimension;

  /*!
   * \brief   The number of \<m\> \<sigma\> \<s\> or \<m\> \<sigma\> \<l\>
   *          tuples passed to the executable on the command line.
   *
   * This number also corresponds to the number of distributions to generate.
   */
  uint32_t count;

  /*!
   * \brief   A vector of entries in the form of parsed tuples.
   *
   * Each entry corresponds to a distribution to generate.
   */
  Generate_Diagonal_Distribution_Arguments_Entry * entries;
} Generate_Diagonal_Distribution_Arguments;

/*!
 * \brief   A data structure representing an export job.
 *
 * \ingroup generate_diagonal_distribution_exe
 */
typedef struct {
  /*!
   * \brief   A pointer to the  distribution to export.
   */
  Diagonal_Distribution * distribution;

  /*!
   * \brief   The path to which to export the distribution.
   */
  char path[MAX_SIZE_PATH_BUFFER];
} Diagonal_Distribution_Export_Job;

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
 *           #FALSE. When returning #FALSE, this function expects the caller to
 *           terminate the executable.
 *
 * \return   #TRUE if the command line arguments were successfully parsed,
 *           #FALSE otherwise. If #FALSE is returned, the data structure may be
 *           only partially initialized and memory only partially allocated.
 */
static bool arguments_init_parse_command_line(
  Generate_Diagonal_Distribution_Arguments * const arguments,
  const int argc,
  char ** argv)
{
  /* Initialize the arguments data structure. */
  arguments->selection_method = SELECTION_METHOD_UNDEFINED;
  arguments->tradeoff_method = TRADEOFF_METHOD_UNDEFINED;
  arguments->dimension = 0;
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

        if (mpz_cmp_ui(arguments->explicit_d, 0) <= 0) {
          fprintf(stderr, "Error: The value of <d> after -exp must be "
            "positive.\n");
          return FALSE;
        }

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

    /* Parse the dimension. */
    if (0 == strcmp(argv[i], "-dim")) {
      /* Check that a dimension has not already been specified. */
      if (0 != arguments->dimension) {
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
    }

    /* Stop parsing. */
    break;
  }

  /* Set default parameters if arguments where not explicitly specified. */
  if (0 == arguments->dimension) {
    arguments->dimension = 2048;
  }

  if (SELECTION_METHOD_UNDEFINED == arguments->selection_method) {
    arguments->selection_method = SELECTION_METHOD_DETERMINISTIC;
  }

  if (TRADEOFF_METHOD_UNDEFINED == arguments->tradeoff_method) {
    arguments->tradeoff_method = TRADEOFF_METHOD_FACTOR;
  }

  /* Parse tuples { <m> <sigma> <s> } or { <m> <sigma> <l> }. */
  if (((argc - i) <= 0) || (0 != ((argc - i) % 3))) {
    fprintf(stderr, "Error: Incorrect command line arguments; expected tuples "
      "but found an incorrect number of arguments.\n");
    return FALSE;
  }

  arguments->count = (uint32_t)((argc - i) / 3);

  arguments->entries =
    (Generate_Diagonal_Distribution_Arguments_Entry *)malloc(
      (arguments->count) * sizeof(
        Generate_Diagonal_Distribution_Arguments_Entry));
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
  for (uint32_t j = 0; j < arguments->count; j++, i += 3) {
    const int m = atoi(argv[i]);
    if (m <= 1) {
      fprintf(stderr, "Error: Failed to parse <m>.\n");
      return FALSE;
    }

    arguments->entries[j].m = (uint32_t)m;

    const int sigma = atoi(argv[i + 1]);
    if (sigma < 0) {
      fprintf(stderr, "Error: Failed to parse <sigma>.\n");
      return FALSE;
    }

    arguments->entries[j].sigma = (uint32_t)sigma;

    if (arguments->tradeoff_method == TRADEOFF_METHOD_FACTOR) {
      const int s = atoi(argv[i + 2]);
      if (s < 1) {
        fprintf(stderr, "Error: Failed to parse <s>.\n");
        return FALSE;
      }

      arguments->entries[j].s = (uint32_t)s;
      arguments->entries[j].l =
        (uint32_t)ceil((double)(m + sigma) / (double)s);
    } else {
      const int l = atoi(argv[i + 2]);
      if (l < 0) {
        fprintf(stderr, "Error: Failed to parse <l>.\n");
        return FALSE;
      }

      arguments->entries[j].s = 0;
      arguments->entries[j].l = (uint32_t)l;
    }

    char suffix[MAX_SIZE_PATH_BUFFER];
    if (arguments->tradeoff_method == TRADEOFF_METHOD_FACTOR) {
      safe_snprintf(suffix,
        MAX_SIZE_PATH_BUFFER,
        "dim-%u-m-%u-sigma-%u-s-%u.txt",
        arguments->dimension,
        arguments->entries[j].m,
        arguments->entries[j].sigma,
        arguments->entries[j].s);
    } else {
      safe_snprintf(suffix,
        MAX_SIZE_PATH_BUFFER,
        "dim-%u-m-%u-sigma-%u-l-%u.txt",
        arguments->dimension,
        arguments->entries[j].m,
        arguments->entries[j].sigma,
        arguments->entries[j].l);
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

      safe_snprintf(arguments->entries[j].path, MAX_SIZE_PATH_BUFFER,
        "%s/diagonal-distribution-det-%s", DISTRIBUTIONS_DIRECTORY, suffix);
    } else if (SELECTION_METHOD_RANDOM == arguments->selection_method) {
      if (j > 0) { /* Check if we should copy the previous value */
        if ((arguments->entries[j-1].m) == (arguments->entries[j].m)) {
          mpz_set(arguments->entries[j].d, arguments->entries[j-1].d);
          mpz_set(arguments->entries[j].r, arguments->entries[j-1].r);
        }
      }

      /* If we did not copy the previous value, select the value at random. */
      if (0 == mpz_cmp_ui(arguments->entries[j].r, 0)) {
        parameters_selection_random_d_and_r(
          arguments->entries[j].d,
          arguments->entries[j].r,
          arguments->entries[j].m);
      }

      safe_snprintf(arguments->entries[j].path, MAX_SIZE_PATH_BUFFER,
        "%s/diagonal-distribution-rnd-%s", DISTRIBUTIONS_DIRECTORY, suffix);
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

      mpz_set(arguments->entries[j].r, arguments->explicit_r);
      mpz_set(arguments->entries[j].d, arguments->explicit_d);

      safe_snprintf(arguments->entries[j].path, MAX_SIZE_PATH_BUFFER,
        "%s/diagonal-distribution-exp-%s", DISTRIBUTIONS_DIRECTORY, suffix);
    } else {
      fprintf(stderr, "Error: Internal error: Unknown selection method.\n");
      return FALSE;
    }

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
  Generate_Diagonal_Distribution_Arguments * const arguments)
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

  memset(arguments, 0, sizeof(Generate_Diagonal_Distribution_Arguments));
}

/*!
 * \brief   Exports the computed probability distribution to file.
 *
 * This function is intended to run in a separate thread.
 *
 * \param[in, out] ptr    A pointer to a Diagonal_Distribution_Export_Job data
 *                        structure that describes the export job.
 *
 * \remark   The export job and the distribution it contains will be cleared
 *           and deallocated by this function.
 */
static void * main_server_export_distribution(void * ptr)
{
  /* Map the input to data structures. */
  Diagonal_Distribution_Export_Job * job =
    (Diagonal_Distribution_Export_Job *)ptr;

  Diagonal_Distribution * distribution = job->distribution;
  const char * path = job->path;

  /* Sort the slices in the distribution. */
  printf("Sorting the slices in the distribution...\n");
  diagonal_distribution_sort_slices(distribution);

  /* Export, clear and deallocate the distribution. */
  printf("Exporting the distribution to \"%s\"...\n", path);

  FILE * file = fopen(path, "wb+");
  if (NULL == file) {
    critical("main_server_export_distribution(): Failed to open \"%s\".", path);
  }

  diagonal_distribution_export(distribution, file);
  fclose(file);

  diagonal_distribution_clear(distribution);
  diagonal_distribution_dealloc(&distribution);

  printf("Finished exporting the distribution to \"%s\".\n", path);

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
 * \param[in] entry       The \<m\> \<sigma\> \<s\> or \<m\> \<sigma\> \<l\>
 *                        tuple to process.
 * \param[in] mpi_size    The number of nodes.
 * \param[in] pool        The thread pool from which to spawn threads for
 *                        export operations that run in the background.
 */
static void main_server(
  const Generate_Diagonal_Distribution_Arguments * const arguments,
  const Generate_Diagonal_Distribution_Arguments_Entry * const entry,
  int mpi_size,
  Thread_Pool * const pool)
{
  /* Setup the parameters. */
  const uint32_t t = 30;

  Diagonal_Parameters parameters;
  diagonal_parameters_init(&parameters);

  if (TRADEOFF_METHOD_FACTOR == arguments->tradeoff_method) {
    diagonal_parameters_explicit_m_s(
      &parameters,
      entry->d,
      entry->r,
      entry->m,
      entry->sigma,
      entry->s,
      t);
  } else {
    diagonal_parameters_explicit_m_l(
      &parameters,
      entry->d,
      entry->r,
      entry->m,
      entry->sigma,
      entry->l,
      t);
  }

  /* Send broadcast of the parameters. */
  diagonal_parameters_bcast_send(&parameters, MPI_RANK_ROOT);

  /* Send broadcast of the dimension. */
  uint32_t dimension = arguments->dimension;

  if (MPI_SUCCESS !=
    MPI_Bcast(&dimension,
    1, /* count */
    MPI_UNSIGNED,
    MPI_RANK_ROOT,
    MPI_COMM_WORLD))
  {
    critical("main_server(): Failed to send broadcast of the dimension.");
  }

  /* Setup the distribution enumerator. */
  const bool mirrored = TRUE;

  Diagonal_Distribution_Enumerator enumerator;
  diagonal_distribution_enumerator_init(&enumerator, &parameters, mirrored);

  /* Get the required distribution capacity in slices. */
  uint32_t capacity = diagonal_distribution_enumerator_count(&enumerator);

  if (mirrored) {
    capacity *= 2;
  }

  /* Setup the distribution. */
  uint32_t flags = 0;

  if (SELECTION_METHOD_DETERMINISTIC == arguments->selection_method) {
    flags |= DIAGONAL_DISTRIBUTION_FLAG_DETERMINISTIC;
  } else if (SELECTION_METHOD_EXPLICIT == arguments->selection_method) {
    flags |= DIAGONAL_DISTRIBUTION_FLAG_EXPLICIT;
  } else if (SELECTION_METHOD_RANDOM == arguments->selection_method) {
    flags |= DIAGONAL_DISTRIBUTION_FLAG_RANDOM;
  } else {
    critical("main_server(): Unknown selection method.");
  }

  Diagonal_Distribution * distribution = diagonal_distribution_alloc();
  diagonal_distribution_init(distribution, &parameters, flags, capacity);

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
      int32_t min_log_alpha_r;

      bool found = diagonal_distribution_enumerator_next(
                      &min_log_alpha_r,
                      &enumerator);
      if (found) {
        printf("Processing slice: %u / %u\n",
          enumerator.offset,
          enumerator.count);

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
            &min_log_alpha_r,
            1, /* count */
            MPI_INT,
            status.MPI_SOURCE, /* destination */
            MPI_TAG_SLICE_MIN_LOG_ALPHA,
            MPI_COMM_WORLD))
        {
          critical("main_server(): Failed to send min_log_alpha_r.");
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
      Diagonal_Distribution_Slice * slice = diagonal_distribution_slice_alloc();

      diagonal_distribution_slice_init_recv(slice, status.MPI_SOURCE);

      diagonal_distribution_insert_slice(distribution, slice);

      if (mirrored) {
        Diagonal_Distribution_Slice * mirrored_slice =
          diagonal_distribution_slice_alloc();
        diagonal_distribution_slice_init_copy(mirrored_slice, slice);

        mirrored_slice->min_log_alpha_r *= -1; /* Mirror alpha_r. */

        mirrored_slice->flags |= SLICE_FLAGS_MIRRORED;

        diagonal_distribution_insert_slice(distribution, mirrored_slice);
      }

      /* Print some feedback to the user. */
      printf("Slice total probability is: %.24LG\n", slice->total_probability);
      printf("Slice dimension is: %u\n", slice->dimension);

      printf("Total probability is now: %.24LG\n",
        distribution->total_probability);

      /* Note: The slice is inserted by reference and will be deallocated by
       *       the call to diagonal_distribution_dealloc(). It must not be
       *       deallocated here as doing so would lead to memory corruption. */
    } else if (MPI_NOTIFY_SLICE_SKIP == notification) {
      /* Ignore. */
    } else {
      critical("main_server(): Unknown notification.");
    }
  }

  /* Clear memory. */
  diagonal_distribution_enumerator_clear(&enumerator);

  diagonal_parameters_clear(&parameters);

  /* Setup an export job. */
  Diagonal_Distribution_Export_Job * export_job =
    (Diagonal_Distribution_Export_Job *)
      malloc(sizeof(Diagonal_Distribution_Export_Job));
  if (NULL == export_job) {
    critical("main_server(): Failed to allocate memory.");
  }

  export_job->distribution = distribution;
  safe_strlcpy(export_job->path, entry->path, MAX_SIZE_PATH_BUFFER);

  /* Spawn a thread to sort, export and deallocate the distribution. */
  thread_pool_spawn(pool, main_server_export_distribution, export_job);
}

/*!
 * \brief   The main function on the client node.
 *
 * This function is called once by main() for each distribution to generate.
 */
static void main_client()
{
  Diagonal_Parameters parameters;
  diagonal_parameters_init(&parameters);

  /* Receive broadcast of the distribution parameters. */
  diagonal_parameters_bcast_recv(&parameters, MPI_RANK_ROOT);

  /* Receive broadcast of the dimension. */
  uint32_t dimension;

  if (MPI_SUCCESS !=
    MPI_Bcast(
      &dimension,
      1, /* count */
      MPI_UNSIGNED,
      MPI_RANK_ROOT,
      MPI_COMM_WORLD))
  {
    critical("main_client(): Failed to receive broadcast of the dimension.");
  }

  /* Process jobs. */
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
    int32_t min_log_alpha_r;

    if (MPI_SUCCESS != MPI_Recv(
        &min_log_alpha_r,
        1, /* count */
        MPI_INT,
        MPI_RANK_ROOT,
        MPI_TAG_SLICE_MIN_LOG_ALPHA,
        MPI_COMM_WORLD,
        &status))
    {
      critical("main_client(): Failed to receive min_log_alpha_r.");
    }

    /* Compute slice. */
    Diagonal_Distribution_Slice slice;
    diagonal_distribution_slice_init(&slice, dimension);

    diagonal_distribution_slice_compute_richardson(
        &slice,
        &parameters,
        min_log_alpha_r);

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

    /* Send back the slice. */
    diagonal_distribution_slice_send(&slice, MPI_RANK_ROOT);

    /* Clear memory. */
    diagonal_distribution_slice_clear(&slice);
  }

  /* Clear memory. */
  diagonal_parameters_clear(&parameters);
}

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  fprintf(file, "Synopsis: mpirun generate_diagonal_distribution \\\n"
          "   [ -dim <dimension> ] [ -det | -rnd | -exp <d> <r> ] \\\n"
            "      ( [ -s ] <m> <sigma> <s> { <m> <sigma> <s> } | \\\n"
            "          -l   <m> <sigma> <l> { <m> <sigma> <l> } )\n");

  fprintf(file, "\n");
  fprintf(file, "Selection method for d and r: -- defaults to -det\n");

  fprintf(file,
    " -det     select d and r deterministically from Catalan's "
      "constant\n");
  fprintf(file,
    " -rnd     select r uniformly at random from (2^(m-1), 2^m)\n");
  fprintf(file,
    "          and d uniformly at random from [r/2, r)\n");
  fprintf(file,
    " -exp     explicitly set d and r to <d> and <r>\n");

  fprintf(file, "\n");
  fprintf(file,
    "Tuples <m> <sigma> <l> or <m> <sigma> <s>:\n");
  fprintf(file,
    " <m>      the length m in bits of r\n");
  fprintf(file,
    " <sigma>  the padding length sigma \n");
  fprintf(file,
    " <s>      the tradeoff factor s; used to set l = ceil((m + sigma) / s)\n");
  fprintf(file,
    " <l>      the parameter l\n");

  fprintf(file, "\n");
  fprintf(file, "Dimension: -- defaults to 2048\n");
  fprintf(file,
    " -dim     explicitly set the slice dimension to <dimension>\n");

  fprintf(file,
    "\nNote: If -rnd is specified and m is kept constant for consecutive "
      "tuples\n");
  fprintf(file,
    "<m> <sigma> <s> or <l>, the same values of d and r will be re-used.\n");

  fprintf(file,
    "\nNote: This implementation is optimized for small to medium kappa_r. If");
  fprintf(file,
    "\nkappa_r is very large, it may fail to capture the probability mass.\n");
}

/*!
 * \brief   The main entry point to the generate_diagonal_distribution
 *          executable.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \return  Zero upon successful execution, non-zero otherwise.
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

  Generate_Diagonal_Distribution_Arguments arguments;

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
        printf("Level of thread support provided: "
          "MPI_THREAD_SINGLE\n");
        break;

      case MPI_THREAD_FUNNELED:
        printf("Level of thread support provided: "
          "MPI_THREAD_FUNNELED\n");
        break;

      case MPI_THREAD_SERIALIZED:
        printf("Level of thread support provided: "
          "MPI_THREAD_SERIALIZED\n");
        break;

      case MPI_THREAD_MULTIPLE:
        printf("Level of thread support provided: "
          "MPI_THREAD_MULTIPLE\n");
        break;

      default:
        printf("Level of thread support provided: "
          "Unknown (%d)\n", provided);
        break;
    }

    printf("\n");

    /* Broadcast the number of distributions. */
    if (MPI_SUCCESS !=
      MPI_Bcast(
        &(arguments.count),
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): "
        "Failed to broadcast the number of distributions.");
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
    uint32_t count;

    if (MPI_SUCCESS !=
      MPI_Bcast(
        &count,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): "
        "Failed to broadcast the number of distributions.");
    }

    for (uint32_t i = 0; i < count; i++) {
      main_client();
    }
  }

  if (MPI_SUCCESS != MPI_Finalize()) {
    critical("main(): Failed to finalize MPI.");
  }

  /* Signal success. */
  return 0;
}
